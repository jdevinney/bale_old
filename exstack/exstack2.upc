/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
// 
 *****************************************************************/ 
/*! \file exstack2.upc
 * \brief A library to do asynchronous buffered communications in parallel programs.
 */
#include "libgetput.h" // from exstack library
#include "exstack.h" // from exstack library

/* internal routines */
int64_t exstack2_flush_needed(exstack2_t * Xstk2);

/*! \brief Initialize a exstack2_t structure.
 * \param buf_cnt The number of packages in each of the send and the receive stacks (buffer).
 * \param pkg_size The number of bytes in an package.
 * \return 0 on success, non-0 on error.
 * \ingroup exstackgrp
 */
exstack2_t * exstack2_init(int64_t buf_cnt, size_t pkg_size)
{
  int64_t i;
  if(!buf_cnt) return(NULL);

  exstack2_t * XS2 = calloc(1, sizeof(exstack2_t));
  
  XS2->buf_cnt    = buf_cnt;   // in pkgs
  XS2->pkg_size   = pkg_size;       // in bytes
  
  // set up send and receive buffers
  int64_t buffer_bytes = pkg_size*buf_cnt;

  XS2->s_snd_buffer = lgp_all_alloc (buffer_bytes*(size_t)THREADS*(size_t)THREADS, sizeof(char));
  if(XS2->s_snd_buffer  == NULL) return(NULL);

  XS2->s_rcv_buffer = lgp_all_alloc (buffer_bytes*(size_t)THREADS*(size_t)THREADS, sizeof(char));
  if(XS2->s_rcv_buffer  == NULL) return(NULL);

  XS2->l_snd_buffer  = calloc(THREADS, sizeof(char*));  
  if(XS2->l_snd_buffer == NULL) return(NULL);

  XS2->l_rcv_buffer  = calloc(THREADS, sizeof(char*));
  if(XS2->l_rcv_buffer == NULL) return(NULL);
  
  XS2->l_snd_buffer[0] = lgp_local_part(char, XS2->s_snd_buffer);
  XS2->l_rcv_buffer[0] = lgp_local_part(char, XS2->s_rcv_buffer);
  for(int64_t th=1; th<THREADS; th++) {
    XS2->l_snd_buffer[th] = XS2->l_snd_buffer[th-1] + buffer_bytes;
    XS2->l_rcv_buffer[th] = XS2->l_rcv_buffer[th-1] + buffer_bytes;    
  }
    
  // the message queue is a circular buffer, indices into the queue are ever increasing 
  // int64_t that are masked down the length of the queue
  int64_t Q_len = 1;
  while( Q_len < 2L*(size_t)THREADS ) // find is smallest Power of 2 larger then 2*THREADS
    Q_len = (Q_len << 1);
  XS2->msg_Q_mask = Q_len-1;

  XS2->s_msg_queue = lgp_all_alloc(Q_len*(size_t)THREADS, sizeof(int64_t));
  if(XS2->s_msg_queue == NULL) return(NULL);
  XS2->l_msg_queue = lgp_local_part(volatile int64_t, XS2->s_msg_queue);
  for(i = 0; i < Q_len; i++)
    XS2->l_msg_queue[i] = -1L;
  
  XS2->s_num_msgs      = lgp_all_alloc(THREADS, sizeof(int64_t));
  if(XS2->s_num_msgs   == NULL) return(NULL);
  XS2->l_num_msgs = lgp_local_part(volatile int64_t, XS2->s_num_msgs);
  
  XS2->num_popped = 0L;
  XS2->l_num_msgs[0] = 0L;
  
  XS2->active_buffer_queue = calloc(THREADS, sizeof(int64_t));
  XS2->num_active_buffers = 0L;
  XS2->current_active_index = -1L;
  XS2->num_made_active = 0L;

  // allocate and initialize s_can_send 
  XS2->s_can_send = lgp_all_alloc((size_t)THREADS*(size_t)THREADS, sizeof(uint64_t));
  if(XS2->s_can_send == NULL) return(NULL);
  XS2->l_can_send = lgp_local_part(volatile int64_t,  XS2->s_can_send);
  for(i = 0; i < THREADS; i++)
    XS2->l_can_send[i] = 1L;

  // allocate and initialize push and pop ptrs
  XS2->push_cnt = calloc(THREADS, sizeof(int64_t));
  if(XS2->push_cnt == NULL) return(NULL);
  XS2->push_ptr = calloc(THREADS, sizeof(char *));
  if(XS2->push_ptr == NULL) return(NULL);
  for(i = 0; i < THREADS; i++){
    XS2->push_ptr[i] = XS2->l_snd_buffer[i];
    XS2->push_cnt[i] = 0L;
  }
  XS2->pop_cnt = calloc(THREADS, sizeof(int64_t));
  XS2->pop_ptr = calloc(THREADS, sizeof(char*));
  XS2->pop_pe = -1L;

  // linked list trick for flush
  int64_t t, s, swap;
  XS2->flush_perm = calloc(THREADS, sizeof(int64_t));
  if(XS2->flush_perm == NULL) return(NULL);
  for(t = 0; t < THREADS; t++){
     XS2->flush_perm[t] = t;
  }
  srand48( MYTHREAD*MYTHREAD*MYTHREAD*313 + 909091 );
  for(t=THREADS-1; t>1; t--){
    s = lrand48() % t; 
    swap = XS2->flush_perm[t];
    XS2->flush_perm[t] = XS2->flush_perm[s];
    XS2->flush_perm[s] = swap;
  }

  XS2->flush_order = calloc(THREADS+1, sizeof(int64_t)); 
  if(XS2->flush_order == NULL) return(NULL); 

  t = THREADS;
  for(i = 0; i < THREADS; i++){
    XS2->flush_order[t] = XS2->flush_perm[i];
    t = XS2->flush_perm[i];
  }
  XS2->flush_order[t] = THREADS;

  XS2->all_done = 0L;
  XS2->num_done_sending = 0;     // incremented by exstack2_pop, for each done pushing message

  lgp_barrier();
  
  
  return(XS2);
}

/*! \brief This function should be called regularly after pushes and pops to trigger the exstack2 endgame.
 *
 * One of exstack2_proceed's jobs is to tell it's caller when this exstack2  is finished.
 * An exstack2 is finished from this threads's perspective once...
 * a) all the PE's have nothing left to push onto send stacks
 * b) all the PE's send stacks are empty
 * c) this PE's receive stacks are empty
 *
 * \param Xstk2 The exstack2 struct.
 * \param done_pushing Set this to one to signal to the exstack2 that you will no
 * longer be calling push in this context. Once you set done_pushing to one, this
 * triggers the "exstack2 endgame". You cannot undo this action and so it is an error
 * to try to push more onto this exstack2 once you set this flag (unless you call exstack2_reset).
 * \return 0 if the Xstk2 is finished, non-zero otherwise
 * \ingroup exstackgrp
 */
int64_t exstack2_proceed(exstack2_t *Xstk2, int done_pushing)
{
  if( done_pushing ) {
    if( exstack2_flush_needed(Xstk2) ) {
      return(2L);
    } else {
      if( Xstk2->num_done_sending < THREADS ) 
        return(4L);
      else {
        if ( Xstk2->num_popped == Xstk2->l_num_msgs[0] ){
          Xstk2->all_done = 1L;
          return(0L);
        } else
          return(8L);
      }
    }
  } else {
    return(1L);
  }
}

/*!\brief Try to send all of your unsent buffers.
 * \param Xstk2 The exstack2 struct.
 * \return 0 if you have successfully flushed all buffers, 1 otherwise.
 */
int64_t exstack2_flush_needed(exstack2_t * Xstk2)
{
  int64_t i, t;
  int64_t ret = 1L;

  t = THREADS;
  i = Xstk2->flush_order[t];
  if(i == THREADS) 
    return(0L);
  while( i != THREADS ) {
    if( exstack2_send(Xstk2, i, 1) == 0L ) {  // could not send
      t = i;
      ret = 1L;
    } else {  // send worked, so take this out of the linked list
      Xstk2->flush_order[t] = Xstk2->flush_order[i];
    }
    i = Xstk2->flush_order[t];
  }
  
  return(ret);
}

/******************************************************************************/
/*! \brief This is how you push items onto the exstack2 struct.
 *
 * When you push an item onto an exstack2, it first attempts to 
 * push the item onto a send buffer. If there is no more room
 * on that buffer, the buffer will attempt to send itself to its 
 * destination. If the destination is not ready to receive this buffer
 * then this function returns 1 (otherwise 0).

 * \param Xstk2 The exstack2 struct.
 * \param pkg A pointer to memory of the package you want to push
 * \param pe The destination pe.
 * \return A LIE
 * 1 upon success or 0 if this item cannot be put on a stack at this time.
 * \ingroup exstackgrp
 */
#if 0
int64_t exstack2_push(exstack2_t * Xstk2, void *pkg, int64_t pe)
{
  if(Xstk2->push_cnt[pe] < Xstk2->buf_cnt){           // copy to the buffer
    memcpy(Xstk2->push_ptr[pe], (char*)pkg, Xstk2->pkg_size);
    Xstk2->push_cnt[pe]++;
    Xstk2->push_ptr[pe] += Xstk2->pkg_size;
  }
  if(Xstk2->push_cnt[pe] == Xstk2->buf_cnt){   // this push filled a buffer, try to send it
    if(exstack2_send(Xstk2, pe, 0) == 0){     // but we could maybe push to other stacks
      return(0L);
    } 
  }
  return(1L);
}
int64_t 
exstack2_push(exstack2_t * Xstk2, void *pkg, int64_t pe)
{
  if(Xstk2->push_cnt[pe] < Xstk2->buf_cnt){
    memcpy(Xstk2->push_ptr[pe], (char*)pkg, Xstk2->pkg_size);
    Xstk2->push_cnt[pe]++;
    Xstk2->push_ptr[pe] += Xstk2->pkg_size;
    return(1L);
  } else { 
    if( exstack2_send(Xstk2, pe, 0) ) {
      memcpy(Xstk2->push_ptr[pe], (char*)pkg, Xstk2->pkg_size);
      Xstk2->push_cnt[pe]++;
      Xstk2->push_ptr[pe] += Xstk2->pkg_size;
      return(1L);
    } else {
      return(0L);
    }
  }
}
#else
int64_t
exstack2_push(exstack2_t * Xstk2, void *pkg, int64_t pe)
{
  int64_t head = Xstk2->buf_cnt - Xstk2->push_cnt[pe]; 
  if( head ){
    memcpy(Xstk2->push_ptr[pe], (char*)pkg, Xstk2->pkg_size);
    Xstk2->push_cnt[pe]++;
    Xstk2->push_ptr[pe] += Xstk2->pkg_size;
  } else { 
    if( exstack2_send(Xstk2, pe, 0) ) {
      head = Xstk2->buf_cnt; 
      memcpy(Xstk2->push_ptr[pe], (char*)pkg, Xstk2->pkg_size);
      Xstk2->push_cnt[pe]++;
      Xstk2->push_ptr[pe] += Xstk2->pkg_size;
    } 
  }
  return(head);
} 
#endif


/******************************************************************************/
/*! \brief This subroutine sends a stack.
 * If it can send the stack it sets _can_send to 0 and claims a spot on the 
 * receiver's pop_from queue.
 * \param Xstk2 The exstack2 struct
 * \param pe The destination pe.
 * \param islast Signal that this is the last buffer that we will be
 *        sending to this pe in this context.
 * \return 0 if is not able to be sent at this time, * \return 1 if the buffer is sent successfully
 */
int64_t 
exstack2_send(exstack2_t * Xstk2, int64_t pe, int islast)
{
  uint64_t pos;

  if(Xstk2->l_can_send[pe] == 0L) // not safe to send ... return 1 
      return(0);
  if(Xstk2->push_cnt[pe] > 0L){      // flushing might call send with an empty buffer
    lgp_memput(Xstk2->s_rcv_buffer,
               (const void * restrict)Xstk2->l_snd_buffer[pe],
               Xstk2->push_cnt[pe]*Xstk2->pkg_size,
               ((size_t)THREADS*MYTHREAD*(Xstk2->pkg_size*Xstk2->buf_cnt) + pe)
               );
  }
  
  // mark this buffer so that it can't be sent before the receiver has emptied the corresponding recv buffer 
  // this may have to be done with atomics...
  //_amo_aadd((volatile int64_t *)&Xstk2->l_can_send[pe], -1L); 
  Xstk2->l_can_send[pe] = 0L;
  
  lgp_fence();
  
  // send a message to the receiver pe, that a buffer has been delivered
  //#pragma pgas defer_sync
  //pos = _amo_afadd((SHARED int64_t *)&Xstk2->s_num_msgs[pe], 1L);  //get local index on remote pe for this message  
  pos = lgp_fetch_and_inc((SHARED int64_t *)Xstk2->s_num_msgs, pe);
  pos = pos & Xstk2->msg_Q_mask;
  //Xstk2->s_msg_queue[pos*THREADS + pe] = msg_pack( Xstk2->push_cnt[pe], islast );
  lgp_put_int64((SHARED int64_t *)Xstk2->s_msg_queue, pos*THREADS + pe, msg_pack( Xstk2->push_cnt[pe], islast ));

  // reset the buffer to continue pushing to it
  Xstk2->push_cnt[pe] = 0L;
  Xstk2->push_ptr[pe] = Xstk2->l_snd_buffer[pe];

  lgp_fence();
  
  return(1L);
}

/*! \brief This is the way we pop items off of the receive stacks in a exstack2.
 * \param Xstk2 The exstack2 struct
 * \param pkg A pointer to a pointer to memory that holds the item from the buffer. Note, that this memory is not copied
 *  from the buffer, it is just pointed to by pkg.
 * \param from_pe if non-NULL returns the pe from which the item was popped
 * \return 1 if something was successfully popped, 0 otherwise.
 * Note there is a slight advantage to just assigning a pointer to pkg, rather than copying the package as
 *  we do here.  This seems safer and a lot less ugly.
 * \ingroup exstackgrp
 */
int64_t exstack2_pop(exstack2_t * Xstk2, void *pkg, int64_t *from_pe)
{
  int64_t msg_index, msg, pe, cnt;

  lgp_poll();

  int64_t s2l_num_msgs = Xstk2->l_num_msgs[0];

  if(Xstk2->all_done){
    return(0L);
  }
  
  if(Xstk2->pop_pe != -1L){ // We have a buffer queued up... 
    if(Xstk2->pop_cnt[Xstk2->pop_pe] > 0L){  // there is something to pop

      memcpy((char *)pkg, Xstk2->pop_ptr[Xstk2->pop_pe], Xstk2->pkg_size);
      if( from_pe != NULL )
        *from_pe = Xstk2->pop_pe;
      
      // move the pop pointer and count
      Xstk2->pop_ptr[Xstk2->pop_pe] += Xstk2->pkg_size;
      Xstk2->pop_cnt[Xstk2->pop_pe]--;

      return(1L);
    
    }else{  // the buffer we are working on is now empty 

      // tell whoever sent this buffer, it is safe to start send again 
#pragma pgas defer_sync
      //_amo_aadd((SHARED uint64_t *)&Xstk2->s_can_send[MYTHREAD*THREADS + Xstk2->curr_pop_pe], 1L);
      lgp_put_int64((SHARED int64_t *)Xstk2->s_can_send, (MYTHREAD*THREADS + Xstk2->pop_pe), 1L);
      //Xstk2->s_can_send[MYTHREAD*THREADS + Xstk2->pop_pe] = 1L;
      lgp_fence();  // Necessary?, Harmful
      
      Xstk2->num_popped++;
      Xstk2->num_active_buffers--;
      
      /* swap in the last active buffer from the active buffer queue to fill in the hole left by this empty buffer. */
      assert(Xstk2->current_active_index >= 0);
      Xstk2->active_buffer_queue[Xstk2->current_active_index] = Xstk2->active_buffer_queue[Xstk2->num_active_buffers];
      Xstk2->pop_pe = -1L;
    }
  }

  //   Use this opportunity to activate new messages
  while(Xstk2->num_made_active < s2l_num_msgs){
    msg_index = (Xstk2->num_made_active & Xstk2->msg_Q_mask);
    msg = Xstk2->l_msg_queue[msg_index];
    if(msg == -1L) break;                 // the message is still in flight
    pe = msg_pe(msg);
    cnt = msg_cnt(msg);
    if( msg_islast(msg) ) 
      Xstk2->num_done_sending++;
    if( cnt ){        // this is not an empty "I'm done sending" message
      assert(Xstk2->pop_cnt[pe] == 0L);
      assert(cnt > 0);
      Xstk2->pop_cnt[pe] = cnt;
      lgp_poll();
      Xstk2->pop_ptr[pe] = Xstk2->l_rcv_buffer[pe];
    }
    Xstk2->active_buffer_queue[Xstk2->num_active_buffers] = pe;
    Xstk2->num_active_buffers++;
    Xstk2->num_made_active++;
    Xstk2->l_msg_queue[msg_index] = -1L;
  }
  
  /***************************************************/
  /*  If there are any buffers waiting to be popped  */
  /*  pick one (randomly or the 'next' one) and      */
  /*  the call pop again.                            */
  /***************************************************/
  if(Xstk2->num_active_buffers){
    Xstk2->current_active_index = rand() % Xstk2->num_active_buffers;
    Xstk2->pop_pe = Xstk2->active_buffer_queue[Xstk2->current_active_index];
    return(exstack2_pop(Xstk2, pkg, from_pe));   // we will succeed this time
  }

  return(0L);
}

/******************************************************************************/
/*! \brief Unpop the previous pop in a exstack2.  It is only legal to call after a successful pop. 
 * \param Xstk2 The exstack2 struct
 */
void exstack2_unpop(exstack2_t * Xstk2)
{
  Xstk2->pop_ptr[Xstk2->pop_pe] -= Xstk2->pkg_size;
  Xstk2->pop_cnt[Xstk2->pop_pe]++;
  Xstk2->pop_pe = -1L; // pick a new buffer next time
}

/******************************************************************************/
/*! \brief This is an exstack2 pop that returns a pointer the package instead of copying the 
 *  to a supplied pointer the way pop does.
 * \param Xstk2 The exstack2 struct
 * \param from_pe if non-NULL returns the pe from which the item was popped
 * \return a pointer into the exstack2 buffers for the popped item, or NULL if there is nothing to pop
 * \ingroup exstackgrp
 */
void *exstack2_pull(exstack2_t * Xstk2, int64_t *from_pe) // sets pointer to pkg, not pkg
{
  void *ret;
  int64_t msg_index, msg, pe, cnt;

  lgp_poll();

  int64_t s2l_num_msgs = Xstk2->l_num_msgs[0];

  if(Xstk2->all_done){
    return(0L);
  }
  
  if(Xstk2->pop_pe != -1L){ // We have a buffer queued up... 
    if(Xstk2->pop_cnt[Xstk2->pop_pe] > 0L){  // there is something to pop

      //memcpy((char *)pkg, Xstk2->pop_ptr[Xstk2->pop_pe], Xstk2->pkg_size);
      ret =  Xstk2->pop_ptr[Xstk2->pop_pe];
      if( from_pe != NULL )
        *from_pe = Xstk2->pop_pe;
      
      // move the pop pointer and count
      Xstk2->pop_ptr[Xstk2->pop_pe] += Xstk2->pkg_size;
      Xstk2->pop_cnt[Xstk2->pop_pe]--;

      return(ret);
    
    }else{  // the buffer we are working on is now empty 

      // tell whoever sent this buffer, it is safe to start send again 
#pragma pgas defer_sync
      //_amo_aadd((SHARED uint64_t *)&Xstk2->s_can_send[MYTHREAD*THREADS + Xstk2->curr_pop_pe], 1L);
      //Xstk2->s_can_send[MYTHREAD*THREADS + Xstk2->pop_pe] = 1L;
      lgp_put_int64((SHARED int64_t *)Xstk2->s_can_send, MYTHREAD*THREADS + Xstk2->pop_pe, 1L);
      lgp_fence(); //Necessary?, Harmful?
      
      Xstk2->num_popped++;
      Xstk2->num_active_buffers--;
      
      /* swap in the last active buffer from the active buffer queue to fill in the hole left by this empty buffer. */
      assert(Xstk2->current_active_index >= 0);
      Xstk2->active_buffer_queue[Xstk2->current_active_index] = Xstk2->active_buffer_queue[Xstk2->num_active_buffers];
      Xstk2->pop_pe = -1L;
    }
  }
  
  //   Use this opportunity to activate new messages
  while(Xstk2->num_made_active < s2l_num_msgs){
    msg_index = (Xstk2->num_made_active & Xstk2->msg_Q_mask);
    msg = Xstk2->l_msg_queue[msg_index];
    if(msg == -1L) break;                 // the message is still in flight
    pe = msg_pe(msg);
    cnt = msg_cnt(msg);
    if( msg_islast(msg) ) 
      Xstk2->num_done_sending++;
    if( cnt ){        // this is not an empty "I'm done sending" message
      lgp_poll();
      assert(Xstk2->pop_cnt[pe] == 0L);
      assert(cnt > 0);
      Xstk2->pop_cnt[pe] = cnt;
      //Xstk2->pop_ptr[pe] = Xstk2->l_rcv_buffer[pe] + (cnt-1)*Xstk2->pkg_size;
      Xstk2->pop_ptr[pe] = Xstk2->l_rcv_buffer[pe];
    }
    Xstk2->active_buffer_queue[Xstk2->num_active_buffers] = pe;
    Xstk2->num_active_buffers++;
    Xstk2->num_made_active++;
    Xstk2->l_msg_queue[msg_index] = -1L;
  }
  
  /***************************************************/
  /*  If there are any buffers waiting to be popped  */
  /*  pick one (randomly or the 'next' one) and      */
  /*  the call pop again.                            */
  /***************************************************/
  if(Xstk2->num_active_buffers){
    Xstk2->current_active_index = rand() % Xstk2->num_active_buffers;
    Xstk2->pop_pe = Xstk2->active_buffer_queue[Xstk2->current_active_index];
    return(exstack2_pull(Xstk2, from_pe));   // we will succeed this time
  }

  return(NULL);
}

/*! \brief Unpull. Same as unpop. Provided to match the syntax for exstack2_pull 
 * \param Xstk2 The exstack2 struct
 * \ingroup exstackgrp
 */
void exstack2_unpull(exstack2_t * Xstk2)
{
  Xstk2->pop_ptr[Xstk2->pop_pe] -= Xstk2->pkg_size;
  Xstk2->pop_cnt[Xstk2->pop_pe]++;
  Xstk2->pop_pe = -1L; // pick a new buffer next time
}

/******************************************************************************/
/*! \brief Reset this exstack2 to be used again. This routine is collective.
 * \param Xstk2 The exstack2 struct
 * \ingroup exstackgrp
 */
void exstack2_reset(exstack2_t * Xstk2)
{
  int64_t i, t;

  lgp_barrier();
  Xstk2->all_done = 0L;
  Xstk2->num_done_sending = 0;

  assert( Xstk2->num_popped == Xstk2->l_num_msgs[0]);
  
  Xstk2->num_popped = 0L;
  Xstk2->l_num_msgs[0] = 0L;
  Xstk2->num_active_buffers = 0;
  Xstk2->current_active_index = -1L;
  Xstk2->num_made_active = 0;

  for(i = 0; i < THREADS; i++) {
    assert( Xstk2->l_can_send[i] == 1L );   // Maybe this shouldn't be an error
  }

  /* reset flush_order */
  t = THREADS;
  for(i = 0; i < THREADS; i++){
    Xstk2->flush_order[t] = Xstk2->flush_perm[i];
    t = Xstk2->flush_perm[i];
  }
  Xstk2->flush_order[t] = THREADS;

  lgp_barrier();
}

/******************************************************************************/
/*! \brief Clear up space allocated by the exstack2.
 * \ingroup exstackgrp
 */
void exstack2_clear(exstack2_t * Xstk2)
{
  lgp_all_free((SHARED char *)Xstk2->s_snd_buffer);
  lgp_all_free((SHARED char *)Xstk2->s_rcv_buffer);
  lgp_all_free((SHARED int64_t *)Xstk2->s_msg_queue);
  lgp_all_free((SHARED int64_t *)Xstk2->s_num_msgs);
  lgp_all_free((SHARED int64_t *)Xstk2->s_can_send);

  free(Xstk2->l_snd_buffer);
  free(Xstk2->l_rcv_buffer);
  free(Xstk2->push_cnt);
  free(Xstk2->pop_cnt);
  free(Xstk2->push_ptr);
  free(Xstk2->pop_ptr);
  free(Xstk2->active_buffer_queue);
  free(Xstk2->flush_order);
  free(Xstk2->flush_perm);
  lgp_barrier();
}

