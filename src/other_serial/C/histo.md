## histo
### Definition
We form the histogram of a large number of `int64_t`'s into a large table.
The loop is as simple as:
```
foreach idx in index[]
  counts[idx]++
```
On a serial thread this shows the difference between random stores
and streaming stores.  On a large parallel machine this is the simplest
case of managing the latency and bandwidth of the interconnection network.

### Algorithms(s)
We have the generic algorithm (as written above).

We also have a buffered version where we sort the indices
into buffers, based on their bits. When a buffer get full
we write all the updates from the buffers contents all at once.

A third version first sorts the indices before running the loop.
This ought to be closer to the streams performance as it 
is streaming with holes in it.

### Discussion
In serial, we don't have the problem of atomic updates to `histo`.
Running a version of histo on a node, with node level threads,
might reveal something about atomic updates to memory,
without the performance being dominated by the interconnection network.

### References
https://www.cs.virginia.edu/stream
