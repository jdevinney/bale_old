## ig
### Definition
We do an index_gather of a large number of entries from a large table.
The loop is simply:
```
for i in 0,...n-1
  target[i] = source[ index[i] ]
```
This is complement of [histo](histo.md).
### Algorithms(s)
We have the generic implementation and a buffered implementation.  
In the buffered version, we collect the index[i] values into
buffers based on their high bits. When a buffer is full we 
do all the gathers for the indices in that buffer before continuing.

### Discussion
This is surprising complicated in the parallel case. 


### References
