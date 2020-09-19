# histo (or histogram)

## Description

In the histogram app, each PE generates a list of uniform random indices into a distributed table and then for each index, we increment the table's value at that index. 

In SHMEM, this looks like

```c
for(i = 0; i < N; i++)
  shmem_atomic_add(&table[index[i/THREADS], 1, index[i] % THREADS);
```

where table is a distributed array with M total elements and index is a local array of global indices into the table.

## Discussion
### Why it is in bale
Histogram is the simplest application in bale, yet it is worthy of our attention because it represents a pattern of communication that is frequently used in parallel applications.
The performance of this simple loop is key to the performance of most of the other apps in bale.

### Performance requirements and implications
 This is the pattern where PEs are asynchronously sending lots of small and easy-to-perform updates to other PEs. The histogram pattern is rather easy to write in plain UPC and SHMEM and even in our exstack/conveyor aggregation libraries. We should make the distinction between an app like histogram, where the updates are simple and can be done easily using atomics, or even just puts, versus an app like Single Source Shortst Path ([SSSP](../sssp_src/README.md)) where the update is complicated and would not be simple (or in some cases possible) to achieve with atomics and puts.

Clearly, it does not matter what order the updates are done in the
histogram application, in fact there are no dependencies at all. All
that matters is that we complete all the updates. This makes it an
obvious target for aggregation. The histogram application when written
with exstack looks a little more complicated than the SHMEM implementation...

```c
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
```
### AGIness
The shmem loop looks pretty good to me.

See apps/histo_src/histo_exstack.upc for the full implementation.

