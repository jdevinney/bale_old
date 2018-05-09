/******************************************************************
 * Copyright 2014, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 * POC: Bale <bale@super.org>
 * Please contact the POC before disseminating this code.
 *****************************************************************/ 
/*! \file exstack.upc
 * \brief A library to do bulk synchronous buffered communications in parallel programs.
 */

#include "exstack.h"
#include "libgetput.h"
//#include <upc_castable.h>

/*! \brief Initialize an exstack_t.
 * \param buf_cnt the number of package that will fit into a send or receive stack
 * \param pkg_size The number of bytes per package
 * \return pointer to the initialized exstack_t
 * \ingroup exstackgrp
 */
exstack_t * exstack_init( int64_t buf_cnt, size_t pkg_size)
{
  uint64_t th , buf_siz_alloc ;
  int i , j , k , tmp ;

  if(!buf_cnt) return(NULL);
  exstack_t *XStk = calloc(1, sizeof(exstack_t));

  srand(THREADS*MYTHREAD + 1);                    // create  a different random
  XStk->put_order = calloc(THREADS, sizeof(int));     // memput sequence for each thread
  if(XStk->put_order == NULL) 
    return(NULL);
  for(i=0; i<THREADS; i++) 
    XStk->put_order[i] = i;
  for(j=THREADS-1; j>1; j--) {
    k               = rand()%(j-1);
    tmp             = XStk->put_order[k];
    XStk->put_order[k] = XStk->put_order[j];
    XStk->put_order[j] = tmp;
  }

  /* each thread has a send and receive buffer for every other thread   */
  buf_siz_alloc = sizeof(uint64_t) + pkg_size*buf_cnt;
  XStk->snd_buf = lgp_all_alloc (buf_siz_alloc*THREADS*THREADS, sizeof(char));
  if(XStk->snd_buf   == NULL) return(NULL);
  XStk->rcv_buf = lgp_all_alloc (buf_siz_alloc*THREADS*THREADS, sizeof(char));
  if(XStk->rcv_buf   == NULL) return(NULL);

  XStk->l_snd_buf = calloc(THREADS,sizeof(SHARED char*));
  if(XStk->l_snd_buf == NULL) return(NULL);
  XStk->l_rcv_buf = calloc(THREADS,sizeof(SHARED char*));
  if(XStk->l_rcv_buf == NULL) return(NULL);
  XStk->fifo_ptr  = calloc(THREADS,sizeof(char*));
  if(XStk->fifo_ptr  == NULL) return(NULL);
  XStk->push_ptr  = calloc(THREADS,sizeof(char*));
  if(XStk->push_ptr  == NULL) return(NULL);

  XStk->l_snd_buf[0] = lgp_local_part( char, XStk->snd_buf);
  XStk->l_rcv_buf[0] = lgp_local_part( char, XStk->rcv_buf);
  for(th=1; th<THREADS; th++) {
    XStk->l_snd_buf[th] = XStk->l_snd_buf[th-1] + pkg_size*buf_cnt + sizeof(uint64_t);
    XStk->l_rcv_buf[th] = XStk->l_rcv_buf[th-1] + pkg_size*buf_cnt + sizeof(uint64_t);
  }

  for(th=0; th<THREADS; th++) {
    // counters live in the first 8 bytes of the buffers
    *((uint64_t*)XStk->l_snd_buf[th])  = 0L;
    *((uint64_t*)XStk->l_rcv_buf[th])  = 0L;
    // set the push and pop ptrs to the item in the bufers
    XStk->push_ptr[th]  = &XStk->l_snd_buf[th][sizeof(uint64_t)];
    XStk->fifo_ptr[th]  = &XStk->l_rcv_buf[th][sizeof(uint64_t)];
  }

  XStk->buf_cnt               = buf_cnt;
  XStk->pkg_size              = pkg_size;
  XStk->first_ne_rcv          = 0;

  XStk->wait_done = lgp_all_alloc(THREADS*THREADS,sizeof(SHARED int64_t));
  if(XStk->wait_done == NULL) return(NULL);             
  
  XStk->l_wait_done = lgp_local_part(int64_t, XStk->wait_done);
  for(th=0; th<THREADS; th++) 
    XStk->l_wait_done[th] = 0;
  XStk->notify_done = 0;

  lgp_barrier();
  return XStk;
}


/*! \brief proceeds until all threads have
    reported that they have hit the done condition.
 * \param Xstk a pointer to the exstack struct
 * \param im_done a flag that indicates that this thread has finished
 * \return 0 if all threads have signaled done,  1 otherwise.
 * \ingroup exstackgrp
 */
int64_t exstack_proceed(exstack_t *Xstk , int im_done) {
  int i;
  int64_t breakout;
  if( im_done && !Xstk->notify_done ) {
    for(i=0; i<THREADS; i++) 
      lgp_put_int64(Xstk->wait_done, MYTHREAD*THREADS + i, 1L);
      //Xstk->wait_done[MYTHREAD*THREADS + i] = 1L;      
    Xstk->notify_done = 1;
  }
  lgp_barrier();
  
  breakout = 0;
  
  if( Xstk->notify_done ){ 
    breakout = 1L; 
    for(i=0; i<THREADS; i++)
      breakout &= Xstk->l_wait_done[i];
  }
  
  if(breakout){ 
    return(0); 
  }else{
    return(1);
  } 
}

/*! \brief push a pkg from a thread onto a send stack
 * \param Xstk pointer to an exstack_t
 * \param th_num the id of the thread to which we are sending these package
 * \param push_item a pointer to the package being pushed
 * \return headroom the amount of room (in packages) that was on the stack when called
 * Note: hence returns 0 if push fails and non-zero otherwise
 * \ingroup exstackgrp
 */
int64_t exstack_push(exstack_t *Xstk, void *push_item, int64_t th_num)
{
  // first 64 bits of the buffer hold the current buffer count, regardless of the package size.
  uint64_t h ,*snd_buf_ct_ptr = (uint64_t*)Xstk->l_snd_buf[th_num];   // first address of the buffer (a 64-bit int)

  h = Xstk->buf_cnt - *snd_buf_ct_ptr;    // the headroom in packages for this buffer
  if( h ) {
    memcpy(Xstk->push_ptr[th_num] , (char*)push_item , Xstk->pkg_size);
    (*snd_buf_ct_ptr)++;
    Xstk->push_ptr[th_num] += Xstk->pkg_size;
  }
  return(h);
}


/*! \brief Uses memputs to essentially do any an all_to_all of 
     send stacks to receive stacks
   Note: this is a variable length all_to_all, since we on send whats 
     needed from each stack
 * \param Xstk pointer to an exstack_t
 * \ingroup exstackgrp
 */
void exstack_exchange(exstack_t *Xstk )
{
  uint64_t ran , th , snd_buf_ct;

  // copy the contents of all non-empty send buffers from my thread
  // to the receive buffers of all other threads
  for(ran=0; ran<THREADS; ran++) {
    th = (uint64_t)Xstk->put_order[ran];
    snd_buf_ct = *((uint64_t*)Xstk->l_snd_buf[th]);
    if(snd_buf_ct) {                                           // only send if buffer has something in it
      lgp_memput(Xstk->rcv_buf,
                  Xstk->l_snd_buf[th],
                  snd_buf_ct*Xstk->pkg_size + sizeof(int64_t),
                  THREADS*MYTHREAD*(Xstk->pkg_size*Xstk->buf_cnt + sizeof(int64_t)) + th);

    }
  }

  //upc_synci();
  lgp_barrier();

  // clear send buffers for my thread and 
  // reset crnt_min_headroom and first_ne_rcv
  for(th=0; th<THREADS; th++) {
    *((uint64_t*)Xstk->l_snd_buf[th]) = 0L;
    Xstk->push_ptr[th] = &Xstk->l_snd_buf[th][sizeof(int64_t)];
  }
  Xstk->first_ne_rcv = 0;

  // Note the barrier is after the memputs so that each thread can start as soon as its stacks are ready
  lgp_barrier();

  // set rcv buffer pointers after we know that all the memputs have completed
  for(th=0; th<THREADS; th++) {
    Xstk->fifo_ptr[th]  = Xstk->l_rcv_buf[th] + sizeof(int64_t);
  }
}

/*! \brief pops a package from specified thread
 * \param Xstk pointer to the exstack_t
 * \param pop_item pointer to object holding the package being pop'd
 * \param th_num thread id of the stack to be popped.
 * \return 1 on success, 0 if the stack was empty
 * Note: this uses the buffer as a fifo
 * \ingroup exstackgrp
 */
int64_t exstack_pop_thread(exstack_t *Xstk, void *pop_item, int64_t th_num)
{
  uint64_t i;

  if( (*((uint64_t*)Xstk->l_rcv_buf[th_num])) == 0) 
    return(0);

  memcpy((char*)pop_item , Xstk->fifo_ptr[th_num] , Xstk->pkg_size);
  Xstk->fifo_ptr[th_num] += Xstk->pkg_size;
  (*((uint64_t*)Xstk->l_rcv_buf[th_num]))--;
  return(1);
}

/*! \brief unpops an item after a exstack_pop_thread.
 * \param Xstk pointer to an exstack_t
 * \param th_num thread id of item that gets unpopped
 * Note: it is only safe to unpop from an exstack immediately after 
 * a successful pop.
 * \ingroup exstackgrp
 */
void exstack_unpop_thread(exstack_t *Xstk , int64_t th_num)
/**********************************************************************/
{
  Xstk->fifo_ptr[th_num] -= Xstk->pkg_size;
  (*((uint64_t*)Xstk->l_rcv_buf[th_num]))++;
}


/*! \brief first-in-first-out pop from any thread, ie the next available package
 * \param Xstk pointer to an exstack_t
 * \param pop_item pointer to object holding the package being popped
 * \param from_th  pointer to thread id of item that gets popped or NULL if you don't care
 * \return 1 on success, 0 if all stacks are empty
 * \ingroup exstackgrp
 */
int64_t exstack_pop(exstack_t *Xstk, void *pop_item ,  int64_t *from_th)
{
  int64_t th;

  for(th=Xstk->first_ne_rcv; th<THREADS; th++) {
    if(*((uint64_t*)Xstk->l_rcv_buf[th]) == 0L) continue;

    // non-empty rcv buffer found: pop work item:
    memcpy((char*)pop_item , Xstk->fifo_ptr[th] , Xstk->pkg_size);
    Xstk->fifo_ptr[th] += Xstk->pkg_size;
    (*((uint64_t*)Xstk->l_rcv_buf[th]))--;
    Xstk->first_ne_rcv = th;
    if( from_th != NULL )
       *from_th = th;
    return(1);
  }

  // all receive buffers empty:  get ready for next exstack_memcpy and  return failure
  Xstk->first_ne_rcv = 0;
  for(th=0; th<THREADS; th++) {
    Xstk->fifo_ptr[th] = &Xstk->l_rcv_buf[th][sizeof(uint64_t)];
  }
  return (0);
}

/*! \brief unpops an item after a exstack_pop
 * \param Xstk pointer to an exstack_t
 * Note, in order to unpop a pop_any 
 * use the from_thread filled in by the that pop command.
 * Note: it is only safe to unpop from an exstack immediately after 
 * a successful pop.
 * \ingroup exstackgrp
 */
void exstack_unpop(exstack_t *Xstk)
/**********************************************************************/
{
  // ONLY A GUESS AT WHAT MIGHT WORK 
  Xstk->fifo_ptr[Xstk->first_ne_rcv] -= Xstk->pkg_size;
  (*((uint64_t*)Xstk->l_rcv_buf[Xstk->first_ne_rcv]))++;
}

/*! \brief first-in-first-out pop from any thread, ie the next available package
 * \param Xstk pointer to an exstack_t
 * \param from_th  pointer to thread id of item that gets popped or NULL if you don't care
 * \return a ptr to the popped item with the exstack, or on failure
 * \ingroup exstackgrp
 */
void *exstack_pull(exstack_t *Xstk, int64_t *from_th)
{
  void *ret;
  int64_t th;

  for(th=Xstk->first_ne_rcv; th<THREADS; th++) {
    if(*((uint64_t*)Xstk->l_rcv_buf[th]) == 0L) continue;

    // non-empty rcv buffer found: pop work item:
    ret = Xstk->fifo_ptr[th];

    Xstk->fifo_ptr[th] += Xstk->pkg_size;
    (*((uint64_t*)Xstk->l_rcv_buf[th]))--;
    Xstk->first_ne_rcv = th;
    if( from_th != NULL )
       *from_th = th;

    return(ret);
  }

  // all receive buffers empty:  get ready for next exstack_memcpy and  return failure
  Xstk->first_ne_rcv = 0;
  for(th=0; th<THREADS; th++) {
    Xstk->fifo_ptr[th] = &Xstk->l_rcv_buf[th][sizeof(uint64_t)];
  }
  return (NULL);
}

/*! \brief unpulls an item after a exstack_pull
 * \param Xstk pointer to an exstack_t
 * Note: it is only safe to unpull from an exstack immediately after 
 * a successful pull.
 * \ingroup exstackgrp
 */
void exstack_unpull(exstack_t *Xstk)
/**********************************************************************/
{
  // ONLY A GUESS AT WHAT MIGHT WORK 
  Xstk->fifo_ptr[Xstk->first_ne_rcv] -= Xstk->pkg_size;
  (*((uint64_t*)Xstk->l_rcv_buf[Xstk->first_ne_rcv]))++;
}



/*! \brief returns the minimum headroom over all of my send buffers
 * \param Xstk is the pointer to the exstack
 * \return the min headroom (in number of packages) across all send buffers
 * \ingroup exstackgrp
 */
int64_t exstack_min_headroom(exstack_t *Xstk)
{
  uint64_t t, ret, *snd_buf_ct_ptr;
  // first 64 bits of the buffer hold the current buffer count, regardless of the package size.
  ret = 0;
  for(t=0; t<THREADS; t++) {
    snd_buf_ct_ptr = (uint64_t*)Xstk->l_snd_buf[t];   
    ret += (Xstk->buf_cnt - *snd_buf_ct_ptr); // the headroom in packages going to MYTHREAD
  }
  return(ret);      // min across all my send buffers
}

/*! \brief returns the headroom for a send buffer
 * \param Xstk is the pointer to the exstack
 * \param th_num is the desired thread number
 * \return the headroom (in number of packages) for the send buffer for th_num
 * \ingroup exstackgrp
 */
int64_t exstack_headroom(exstack_t *Xstk , int64_t th_num)
{
  int64_t *snd_buf_ct_ptr = (int64_t*)Xstk->l_snd_buf[th_num];
  return( Xstk->buf_cnt - *snd_buf_ct_ptr);
}


/*! \brief frees all memory allocated by exstack_init
 * \param Xstk is the pointer to the exstack stack
 * \ingroup exstackgrp
 */
void exstack_clear(exstack_t *Xstk )
{
  lgp_all_free(Xstk->wait_done);
  lgp_all_free(Xstk->snd_buf);
  lgp_all_free(Xstk->rcv_buf);
  free(Xstk->put_order);

  free(Xstk->l_snd_buf);
  free(Xstk->l_rcv_buf);
  free(Xstk->fifo_ptr );
  free(Xstk->push_ptr );
}

/*! \brief resets the exstack buffers and resets the end game 
 * \param Xstk is the pointer to the exstack stack
 * \ingroup exstackgrp
 */
void exstack_reset(exstack_t * Xstk){
  
  for(int th=0; th<THREADS; th++) {
    *((uint64_t*)Xstk->l_snd_buf[th])  = 0L;
    *((uint64_t*)Xstk->l_rcv_buf[th])  = 0L;
    Xstk->push_ptr[th]  = &Xstk->l_snd_buf[th][sizeof(SHARED uint64_t)];
    Xstk->fifo_ptr[th]  = &Xstk->l_rcv_buf[th][sizeof(SHARED uint64_t)];
  }
  Xstk->first_ne_rcv = 0;
  for(int th=0; th<THREADS; th++) 
    Xstk->l_wait_done[th] = 0;
  Xstk->notify_done = 0;
}

