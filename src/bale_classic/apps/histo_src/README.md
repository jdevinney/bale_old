# histo (or histogram)

This is the simplest application in bale. It is important because it
represents a pattern of communication that is frequently used in
parallel applications. This is the pattern where PEs are
asynchronously sending lots of small data to other PEs. The histogram
app in bale is a special case of this, where the data being tranferred
is an index into a distributed table and the action we wish to perform
is to increment that table's value at that index.

In SHMEM, this looks like

    for(i = 0; i < N; i++)
      shmem_atomic_add(&table[index[i], 1);

where table is a distributed array and index is a local array of indices into table.

Clearly, it does not matter what order the updates are done in the
histogram application, in fact there are no dependencies at all. All
that matters is that we complete all the updates. This makes it an
obvious target for aggregation. The histogram application when written
with exstack looks a little more complicated...

    while( exstack_proceed(ex, (i==T)) ){
      for( ; i < T; i++){
        col = pckindx[i] >> 16;
        pe  = pckindx[i] & 0xffff;
        if( !exstack_push(ex, &col, pe) )
          break;
      }
      
      exstack_exchange(ex);
    
      while((colp = exstack_pull(ex, NULL)))
        lcounts[*colp]++;
    }

See apps/histo_src/histo_exstack.upc for the full implementation.