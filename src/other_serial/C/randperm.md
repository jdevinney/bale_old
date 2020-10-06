## randperm
### Definition
We fill an array of int64_t's with a flat random permutation.

### Algorithms(s)
#### the Fisher-Yates algorithm
```
   fill the array, rand[ ], with indices 0 thru n-1. 
   for l=n-1, n-2, ... ,1
     swap rand[l] with a random entry in 0,1,...,l
```
By definition this picks "n balls from a urn without replacement".
This is a standard serial algorithm that is in fact a serial algorithm.
You have to process the entries from right to left one at a time.

#### the "dart throwing" algorithm
This is here to shadow the version in bale_classic. It is in 
bale_classic because it is fun.

Pick a dart board (the array) that is bigger than the permutation needed.
Then randomly throw darts at the dart board, 
re-throwing any dart that hits a entry that is already occupied. 
Then we squeeze out the holes.  
Note we pick the dartboard to be twice the size of the array 
so that even the last dart has a 50/50 chance of hitting an open entry.

### the sorting algorithm
We form an array of (index, key) pairs. Then we randomly fill the keys
and sort on the keys. Then we read the permutation from the indices.

### Discussion
This app is also only interesting in bale_classic and only because
the dart throwing algorithm is fun and exercises 
our parallel programming models in interesting ways.

### References
ref for fisher-yates (Knuth probably)
