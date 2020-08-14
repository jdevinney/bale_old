use Random;
use CyclicDist;
use BlockDist;
use PrivateDist;
use Barriers;
use Time;

//
// This is mostly a direct copy of topo.chpl with triangle counting loops
// replacing the toposort algorithm.


// How data is laid out
// nonzero array is BLOCK distributed
// ROWS are cyclic distributed among locales
// offset array is cyclic distributed
// rows are BLOCK distributed amongs tasks (by convention)


// QUESTIONS, thoughts, suggestions:
// 1. Review iterator story... in transpose and permute_matrix I could use an iterator that iterates over all nonzeros properly
// 2. Learn more about returning a tuple (toposort)
// 3. it would be nice if you could demote atomics so that they were no longer atomic and could be treated normally
// 4. Is there a way to do a global (+ reduce x) in the middle of a SPMD like loop?
// 5. I think it would be nice to have a better way to get task locale subdomains


class Sparsemat{
  // you could try to make this a record if performance suffers
  var nr: int;
  var nc: int;
  var nnz: int;
  var off_dom;
  var nnz_dom;  
  var offset : [off_dom] int; 
  var nonzero : [nnz_dom] int;
  var first_row_per_task;
  var num_rows_per_task;
  var row_stride_per_task;
  
  /******************************************************************/
  proc init(nr:int, nc:int){
    var th_per_locale = dataParTasksPerLocale;
    var ntasks = th_per_locale*numLocales;
    this.nr  = nr;
    this.nc = nc;
    this.off_dom = {0..0};
    this.nnz_dom = {0..0};
    this.first_row_per_task  = newCyclicArr(0..#ntasks, int);
    this.num_rows_per_task   = newCyclicArr(0..#ntasks, int);
    this.row_stride_per_task = newCyclicArr(0..#ntasks, int);
    coforall l in Locales do on l{
        var l_nr = my_num_things(nr, l.id, numLocales);
        coforall t in 0..#th_per_locale{
          var my_rank = l.id*th_per_locale + t;
          this.first_row_per_task[my_rank] = l.id + my_first_thing(l_nr, t, th_per_locale, 0)*numLocales;
          this.num_rows_per_task[my_rank] = my_num_things(l_nr, t, th_per_locale);
          this.row_stride_per_task[my_rank] = numLocales;
        }
      }
  }
  /******************************************************************/
  proc init(nr:int, nc:int, nnz_per_locale:[]int){
    var max_nnz_per_locale = (max reduce nnz_per_locale);   
    var th_per_locale = dataParTasksPerLocale;
    var ntasks = th_per_locale*numLocales;
    this.nr  = nr;
    this.nc = nc;
    this.nnz = (+ reduce nnz_per_locale);
    this.off_dom = newCyclicDom(0..#(nr + numLocales));
    this.nnz_dom = newBlockDom(0..#(max_nnz_per_locale*numLocales));
    this.first_row_per_task  = newCyclicArr(0..#ntasks, int);
    this.num_rows_per_task   = newCyclicArr(0..#ntasks, int);
    this.row_stride_per_task = newCyclicArr(0..#ntasks, int);
    coforall l in Locales do on l{
        var l_nr = my_num_things(nr, l.id, numLocales);
        var local_off_dom = off_dom.localSubdomain();
        coforall t in 0..#th_per_locale{
          var my_rank = l.id*th_per_locale + t;
          this.first_row_per_task[my_rank] = local_off_dom.low + my_first_thing(l_nr, t, th_per_locale, 0)*local_off_dom.stride;
          this.num_rows_per_task[my_rank] = my_num_things(l_nr, t, th_per_locale);
          this.row_stride_per_task[my_rank] = local_off_dom.stride;
        }
      }
  }
  /******************************************************************/
  proc getrowcounts(): []int{
    var rc = newCyclicArr(0..#this.nr, int);
    //forall row in 0..#this.nr{
    //rc[row] = this.row_end(row) - this.row_start(row) + 1;
    //}
    forall row in this.rows(){
      rc[row] = this.row_end(row) - this.row_start(row);
    }
    return(rc);
  }
  /******************************************************************/
  proc get_offsets_from_rowcounts(rowcounts){
    //nonzero is BLOCK distributed
    coforall l in Locales do on l{
        ref lrowcounts = rowcounts.localSlice[rowcounts.localSubdomain()];
        const local_offset_inds = this.offset.localSubdomain();
        ref loffset = this.offset.localSlice[local_offset_inds((l.id + numLocales)..)];
        loffset = (+ scan lrowcounts);
        this.offset.localSlice[local_offset_inds] += this.nonzero.localSubdomain().low; //add block offset
      }    
  }
  /******************************************************************/
  proc transpose(){
    var th_per_locale = dataParTasksPerLocale;
    
    // get column counts
    var D = newCyclicDom(0..#this.nc);
    var colcounts: [D] atomic int;
    var cc: [D] int;

    forall nz in this{
      colcounts[nz[2]].add(1);
    }

    /*
    coforall l in Locales do on l{
        coforall t in 0..#th_per_locale{
          var my_rank = l.id*th_per_locale + t;
          for row in this.my_task_rows(my_rank){
            for col in this.row_nz(row){
              colcounts[col].add(1);
            }
          }
        }     
    }
    */
    
    // copy colcounts into non-atomic array
    forall i in 0..#colcounts.size{
      cc[i] = colcounts[i].peek();
    }

    // allocate transpose matrix
    var lnnz: [0..#numLocales] int;
    coforall l in Locales do on l{
        lnnz[l.id] = (+ reduce cc.localSlice[cc.localSubdomain()]);
      }
    var At = new owned Sparsemat(this.nc, this.nr, lnnz);
    At.get_offsets_from_rowcounts(cc);

    // get an temp. offset array that is atomic
    D = newCyclicDom(0..#(At.nr + numLocales));
    var tmpoffset: [D] atomic int;
    forall i in 0..#At.offset.size{
      tmpoffset[i].write(At.offset[i]);
    }

    // populate the nonzero array of At
    forall row in this.rows(){
      for col in this.row_nz(row){
        var pos = tmpoffset[col].fetchAdd(1);
        At.nonzero[pos] = row;
      }    
    }
    /*
    coforall l in Locales do on l{
        coforall t in 0..#th_per_locale{
          var my_rank = l.id*th_per_locale + t;
          for row in this.my_task_rows(my_rank){
            for col in this.row_nz(row){
              var pos = tmpoffset[col].fetchAdd(1);
              At.nonzero[pos] = row;
            }
          }
        }
      }
    */
    return(At);
  }

  /******************************************************************/
  proc is_unit_upper_triangular(): bool{
    var error = 0;
    //forall row in 0..#this.nr with (+ reduce error){
    forall row in this.rows() with (+ reduce error){
      var pivot = false;
      for col in this.row_nz(row){
        if(col < row){error += 1;}
        if(col == row){pivot = true;}
      }
      if(!pivot){error += 1;}
    }
    if(error){return(false);}
    return(true);
  }
  
  /******************************************************************/
  proc row_start(row: int){
    return(this.offset[row]);
  }
  /******************************************************************/
  proc row_end(row: int){
    return(this.offset[row + numLocales]); // this depends on the rows being CYCLIC dist
  }
  /******************************************************************/
  /*                           ITERATORS                            */
  /******************************************************************/
  iter row_nz(row: int) ref{ //serial iterator over nonzeros in a row
    var rs = this.row_start(row);
    var re = this.row_end(row);
    if(re < rs){
      writeln("Hey!!!", rs, " ", re);
    }
    for i in (rs..(re - 1)){ // this depends on nonzero being BLOCK distributed
      ref nz = this.nonzero[i];
      yield(nz);
    }
  }

  /******************************************************************/
  iter my_task_rows(my_rank:int){
    for rowind in (0..#this.num_rows_per_task[my_rank]){
      yield(rowind*this.row_stride_per_task[my_rank] + this.first_row_per_task[my_rank]);
    }
  }
  /******************************************************************/

  iter these(){ //serial class iterator
    writeln("Using serial CLASS iterator!");
    for row in (0..#this.nr){
      for col in this.row_nz(row){
        yield((row,col));
      }
    }
  }
  
  /******************************************************************/
  iter rows(){ //serial row iterator
    writeln("Using serial row iterator!");
    for row in (0..#this.nr){
      yield(row);
    }
  }

  /******************************************************************/
  iter nzs() ref{ //serial nonzero (column indices) iterator
    writeln("Using serial nonzeros iterator");
    for r in this.rows(){
      for c in this.row_nz(r){        
        yield(c);
      }
    }
  }
  /******************************************************************/
  iter rows(param tag: iterKind) where tag == iterKind.standalone{ //parallel row iterator
    var th_per_locale = dataParTasksPerLocale;
    coforall l in Locales do on l{
        coforall t in (0..#th_per_locale){
          var my_rank = l.id*th_per_locale + t;
          for row in this.my_task_rows(my_rank){
            yield(row);
          }
        }
        //forall row in (l.id..(this.nr-1) by numLocales){
        //yield(row);
        //}
      }
  }
  
  /******************************************************************/
  iter these(param tag: iterKind) where tag == iterKind.standalone{ //parallel class iterator
    forall row in this.rows(){
      for col in this.row_nz(row){
        yield(row,col);
      }
    }
  }
  /*coforall l in Locales do on l{
    forall row in (l.id..(this.nr-1) by numLocales){
    for col in this.row_nz(row){
    yield(row, col);
    }
    }
    }  
  */
  /******************************************************************/
  iter nzs(param tag: iterKind) ref where tag == iterKind.standalone{ // parallel nonzero (column indices) iterator
    coforall l in Locales do on l{
        forall row in (l.id..(this.nr-1) by numLocales){
          for col in this.row_nz(row){
            yield(col);
          }
        }
      }
  }

  
  /******************************************************************/
  proc print(detail=1){
    if(detail == 1){
      writeln("MATRIX: nr = ", this.nr, " nc = ", this.nc, " nnz = ", this.nnz);
    }else if(detail > 1){
      writeln("****************************************************");
      writeln("MATRIX: nr = ", this.nr, " nc = ", this.nc, " nnz = ", this.nnz);
      writeln("first_row_per_task: ", this.first_row_per_task);
      writeln("num_rows_per_task: ", this.num_rows_per_task);
      writeln("row_stride_per_task: ", this.row_stride_per_task);
      writeln("OFFSET: ", this.offset);
      writeln("NONZERO: ", this.nonzero);
      
      if(detail == 3){
        for row in this.rows(){
          stdout.write("row ", row, ": ");
          var rs = this.offset[row];
          var re = this.offset[row];

          for col in this.row_nz(row){
            stdout.write(col, ",");
          }
          writeln();
        }
      }
      writeln("****************************************************");
    }
  }
}

// this function is agnostic of layout (BLOCK or CYCLIC)
proc my_num_things(nthings:int, my_rank:int, total_ranks:int){
  return((nthings + total_ranks - my_rank - 1)/total_ranks);
}

proc my_first_thing(nthings:int, my_rank:int, total_ranks:int, layout:int){
  if(layout==0){//BLOCK
    return(my_rank*(nthings/total_ranks) + min(my_rank, mod(nthings,total_ranks)));
  }else{//CYCLIC
    return(my_rank);
  }
}

// We are generating an Erdos-Renyi random graph using
// a method inspired by the paper "Efficient Generation
// of Large Random Networks" by Batagelj and Brandes (2005).
// We observe that since the wait time between nonzeros
// is geometric and since the geometric distribution is memoryless,
// we can just start each row from scratch rather than doing it
// the way the paper does it (this may lead to more random
// draws, but it makes it easier to parallelize!).
// type = 0: upper triangular
// type = 1; lower triangular
// type = 2; symmetric
// unit_diagonal = true -> means all ones on the diagonal

//proc generate_erdos_renyi_random_adj_matrix(n: int, p:real, seed=1208, type=0, unit_diagonal=true){

//}

proc generate_erdos_renyi_half(const n: int, const p:real, const seed=1208, const upper=true, const unit_diagonal=true){
  var th_per_locale = dataParTasksPerLocale;
  var rowcounts = newCyclicArr(0..#n, int);  
  rowcounts = 0;
  var lnnz = newCyclicArr(0..#numLocales, int);
  var DUMMY = new owned Sparsemat(n,n);
  coforall l in Locales do on l{
      var tot: [0..#dataParTasksPerLocale] int;
      ref lrc = rowcounts.localSlice[rowcounts.localSubdomain()];
      coforall t in 0..#dataParTasksPerLocale{
        var my_rank = l.id*dataParTasksPerLocale + t;
        var rndstream = new owned RandomStream(real, seed + my_rank);
        var lq = log(1 - p);
        var bound = n;
        for row in DUMMY.my_task_rows(my_rank){
          var cnt = if unit_diagonal then 1 else 0;
          var col = if upper then row+1 else 0;
          if(!upper){ bound = row-1; }
          while(1){
            var r = rndstream.getNext();
            col += floor(log(1 - r)/lq): int;
            if(col < bound){
              cnt += 1;
              col += 1;
            }else{
              break;
            }
          }
          lrc[row] = cnt;
          tot[t] += cnt;
        }
      }
      lnnz[l.id] = (+ reduce tot);
    }
  
  var A = new owned Sparsemat(n, n, lnnz);
  A.get_offsets_from_rowcounts(rowcounts);
  writeln(rowcounts);
  coforall l in Locales do on l{
      const local_nonzero_inds = A.nonzero.localSubdomain();
      ref lnonzero = A.nonzero.localSlice[local_nonzero_inds];
      coforall t in 0..#dataParTasksPerLocale{
        var my_rank = l.id*dataParTasksPerLocale + t;
        var rndstream = new owned RandomStream(real, seed + my_rank);
        var lq = log(1 - p);
        var bound = n;
        for row in A.my_task_rows(my_rank){
          var nzpos = A.offset[row];
          var col = if upper then row+1 else 0;
          if(!upper){ bound = row-1; }
          if(unit_diagonal){
            lnonzero[nzpos] = row; //diagonal element
            nzpos += local_nonzero_inds.stride;
          }
          while(1){
            var r = rndstream.getNext();
            col += floor(log(1 - r)/lq): int;
            if(col < bound){
              lnonzero[nzpos] = col;
              nzpos += local_nonzero_inds.stride;
              col += 1;
            }else{
              break;
            }
          }
        }
      }
    }

  return(A);
  
}

/***************************************************************/
proc rand_permp(N:int, seed:int) : []int{
  var D = newCyclicDom(0..#(2*N));
  var target: [D] atomic int;
  var perm = newCyclicArr(0..#N, int);
  forall p in target{
    p.write(-1);
  }
  if(0){
    coforall l in Locales do on l{
        var th_per_locale = dataParTasksPerLocale;      
        coforall t in (0..#th_per_locale){
          var my_rank = l.id*th_per_locale + t;
          var rndstream = new owned RandomStream(real, seed + my_rank);
          var i = my_rank;
          var npes = th_per_locale*numLocales;
          while( i < N ){
            var r = (rndstream.getNext() * 2*N): int;
            if(target[r].compareExchange(-1, i)){
              i += npes;
            }
          }
        }
      }

    // count how many darts landed on each Locale
    var cnt: [0..#(numLocales+1)] int;
    coforall l in Locales do on l{
        var tmpcnt = 0;
        forall p in target.localSlice[target.localSubdomain()] with (+ reduce tmpcnt){
          if(p.read() > -1){
            tmpcnt += 1;
          }
        }
        cnt[l.id + 1] = tmpcnt;
      }

    var offset = (+ scan cnt);

    // this loop is only parallel over locales
    coforall l in Locales do on l{
        var pos = 0;
        for p in target.localSlice[target.localSubdomain()]{
          if (p.read() > -1){
            perm[offset[l.id] + pos] = p.read();
            pos += 1;
          }
        }
      }
  }else{
    // new try
    var thpl = dataParTasksPerLocale;
    var b = new Barrier(numLocales*thpl);
    var offset = newBlockArr(0..numLocales*thpl, int);
    coforall l in Locales do on l{
        var th_per_locale = dataParTasksPerLocale;
        coforall t in (0..#th_per_locale){
          var npes = numLocales*th_per_locale;
          var my_rank = l.id*th_per_locale + t;
          var rndstream = new owned RandomStream(real, seed + my_rank);
          // throw the darts
          for i in my_rank..(N-1) by npes{
            do{ var r = (rndstream.getNext() * 2*N): int;
            }while(!target[r].compareExchange(-1, i));
          }
        
          b.barrier();

          //count the number of darts that landed in my task domain
          for i in ((l.id*th_per_locale + t)..(target.size-1) by npes){
            if(target[i].read() > -1){ offset[my_rank]+=1;}
          }
          
        }
      }
    
    offset = (+ scan offset);

    // redistribute the darts in a CYCLIC fashion to create the final permutation
    coforall l in Locales do on l{
        var th_per_locale = dataParTasksPerLocale;
        coforall t in (0..#th_per_locale){
          var my_rank = l.id*th_per_locale + t;
          var npes = numLocales*th_per_locale;
          
          var pos = 0;
          if(my_rank) {pos = offset[my_rank - 1];}
        
          for i in ((l.id*th_per_locale + t)..(target.size-1) by npes){
            var p = target[i].read();
            if(p > -1){
              perm[pos] = p;
              pos += 1;
            }
          }
        }
      }
  }
  if (!is_perm(perm)){
    writeln("ERROR!!!!! not a perm!");
  }
  return(perm);
}

/***************************************************************/
proc invert_perm(perm:[] int){
  var inv_perm = newCyclicArr(0..#perm.size, int);
  [i in 0..#perm.size] inv_perm[perm[i]] = i;
  return(inv_perm);
}
/***************************************************************/
proc is_perm(perm) : bool{

  var N = perm.size;
  var flag = newCyclicArr(0..#N, int(8));
  flag = 0;
  var error = 0;
  forall i in perm with (+ reduce error){
    if((i >= N) || (i < 0)){
      error += 1;
    }else{
      flag[i] = 1;
    }
  }
  forall i in flag with (+ reduce error){
    if(i != 1){
      error +=1;
    }
  }
  if(error){
    return(false);
  }else{
    return(true);
  }
}

/***************************************************************/
// why can't I put ": Sparsemat" as a return type?
proc permute_matrix(A: Sparsemat, rperm: []int, cperm: []int){
  var th_per_locale = dataParTasksPerLocale;
  
  // get the rowcounts array for the original matrix
  var rowcounts_orig = A.getrowcounts();
  
  // get the permuted rowcounts
  var rowcounts = newCyclicArr(0..#A.nr, int);
  forall i in 0..#rperm.size{
    rowcounts[i] = rowcounts_orig[rperm[i]];
  }
  var lnnz: [0..#numLocales] int;
  coforall l in Locales do on l{
      lnnz[l.id] = (+ reduce rowcounts.localSlice[rowcounts.localSubdomain()]);
    }

  var Ap = new owned Sparsemat(A.nr, A.nc, lnnz);
  Ap.get_offsets_from_rowcounts(rowcounts);
  
  // shuffle rows
  
  /* coforall l in Locales do on l{ */
  /*     const local_nonzero_inds = Ap.nonzero.localSubdomain(); */
  /*     ref lnonzero = Ap.nonzero.localSlice[local_nonzero_inds]; */
  /*     coforall t in 0..#th_per_locale{ */
  /*       var my_rank = l.id*th_per_locale + t; */
  /*       for Aprow in Ap.my_task_rows(my_rank){ */
  /*         var Arow = rperm[Aprow]; */
  /*         var nzpos = Ap.offset[Aprow]; */
  /*         for col in A.row_nz(Arow){ */
  /*           lnonzero[nzpos] = col; */
  /*           nzpos += local_nonzero_inds.stride; */
  /*         } */
  /*       } */
  /*     } */
  /*   } */

  // this is a cleaner version of the above code
  // I think this zipped iteration makes it so that Aprow
  // is local (even though it looks like the follower...)
  forall (Arow, Aprow) in zip(rperm, 0..#Ap.nr){
    var nzpos = Ap.offset[Aprow];
    for col in A.row_nz(Arow){
      Ap.nonzero[nzpos] = col;
      nzpos += Ap.nonzero.domain.stride;
    }
  }
  
  // invert cperm
  var cperm_inv = invert_perm(cperm);
  
  // relabel column indicies
  forall c in Ap.nzs(){
    c = cperm_inv[c];
  }
  //forall i in 0..#Ap.nonzero.size{
  //Ap.nonzero[i] = cperm_inv[Ap.nonzero[i]];
  //}
  
  return(Ap);
}

proc triangles_pretty(L: Sparsemat) {
  var th_per_locale = dataParTasksPerLocale;
  var rp = newCyclicArr(0..#L.nr, int);
  var cp = newCyclicArr(0..#L.nc, int);
  var pos: atomic int;
  pos.write(L.nr-1);
  
  // triangle main
  // L is the lower triangular matrix submatrix of the adjacency matrix
  // really computing +/+/( L.&(L*U)) 
  var totcnt = newCyclicArr(0..#numLocales, int);
  coforall l in Locales do on l{
    var cnt: [0..#dataParTasksPerLocale] int;
    coforall t in 0..#th_per_locale{
      var my_rank = l.id*th_per_locale + t;
      for row in L.my_task_rows(my_rank){    
        for col in L.row_nz(row){
          // count the intersection of the columns in row(row) and row(col) of L
          for wc in L.row_nz(col) {
            for wr in L.row_nz(row){
              if (wc == wr) {
                writeln("tri: ",row, col, wc);
                cnt[t] += 1;
              }
              if (wr > wc ) {break;}
            }
          }
        }
      }
    }
    totcnt[l.id] = (+ reduce cnt);
  }
  var count = (+ reduce totcnt);
  return(count);
}

proc triangles(L: Sparsemat) {
  var th_per_locale = dataParTasksPerLocale;
  var rp = newCyclicArr(0..#L.nr, int);
  var cp = newCyclicArr(0..#L.nc, int);
  var pos: atomic int;
  pos.write(L.nr-1);
  
  // triangle main
  // L is the lower triangular matrix submatrix of the adjacency matrix
  // really computing +/+/( L.&(L*U)) 
  var totcnt = newCyclicArr(0..#numLocales, int);
  coforall l in Locales do on l{
    var cnt: [0..#dataParTasksPerLocale] int;
    coforall t in 0..#th_per_locale{
      var my_rank = l.id*th_per_locale + t;
      for row in L.my_task_rows(my_rank){    
        for col in L.row_nz(row){
          // count the intersection of the columns in row(row) and row(col) of L
          var start = L.row_start(row);
          for wc in L.row_nz(col) {
            for i in start..(L.row_end(row)-1) {
              var wr = L.nonzero[i];
              if (wc == wr) {
                writeln("row= ", row, " col= ", col, " w= ", wc);
                cnt[t] += 1;
                start = i+1;
                break;
              }
              if (wr > wc ) {
                start = i;
                break;
              }
            }
          }
        }
      }
    }
    totcnt[l.id] = (+ reduce cnt);
  }
  var count = (+ reduce totcnt);
  return(count);
}

/*****************************************************************/
/*                           MAIN                                */
/*****************************************************************/

config const N = 12;
config const p = 0.1;
config const seed = 1038;
config const detail = 3;

var t: Timer;
t.start();

var A = generate_erdos_renyi_half(N, p, seed, false, false);
A.print(detail);

var numtris = triangles(A);

t.stop();
writeln("Found ", numtris, " triangles");

writeln("Matrix generation time ", t.elapsed()," secs");t.clear();


