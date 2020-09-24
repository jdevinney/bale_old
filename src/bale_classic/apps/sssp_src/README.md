# Single Source Shortest Path  (sssp)

## Definition
Given a directed, weighted graph and an given vertex, *v_0*,
find the lightest weighted path from *v_0* to all other vertices,
where the weight of a path is the sum of the weights of the 
edges in the path.

## Discussion
This is well studied problem.
It is used in classes on algorithms, data structures, and graph theory.
It is also used as a benchmark in Graph500 benchmark suite.
Good references for these uses are easily found.

This is not an academic exercise nor are we interested in benchmarking.
We are using this problem to evaluate what we have learned about our programming models.
Are approach is that of someone given the assignment to get a SSSP application
running reasonable well, at scale, using our programming models.

We consider three algorithms: Dijsktra's, Delta-Stepping, and Bellman-Ford.
These algorithm work by assign a tentative weight, *tent(v)*, to each vertex. 
Initially the weight of the given vertex is zero and all other weights are set to infinity.
The algorithms proceed by changing the state of the vertices from *unreached* (with infinite weight)
to *unsettled* to *settled* (the weight to that vertex is then known). 
The vertices change state via the process of *relaxing edges*.
One relaxes and edge by setting the tentative weight of the head of the edge to the 
minimum of itself and the tentative weight of the tail plus the weight of the edge.
The algorithms vary in how much parallelism they enjoy and how and when vertices are *settled*.

We have implementations of all three algorithms in the other_serial/C "cousin" directory 
and have implementation of Delta-Stepping and Bellman-Ford in this directory.

#### Parallel Considerations
Dijsktra's algorithm is a serial algorithm. It moves vertices from be *unsettled* to *settled* 
one at a time and relaxes each edge exactly once.  It does the minimal amount of work, but
there is no opportunity for parallelism.
Bellman-Ford has a lot of parallelism, but it does a lot of redundant work.
In fact, if one had a processor for each edge, Bellman-Ford could be described 
simply as repeatedly relaxing all the edges until none of the tentative weights improve.
The Delta-Stepping algorithm is in between the two.
It is like Dijsktra's in that most of the edges are relaxed only once, but remaining edges
can be relaxed in parallel.


```c
while( convey_advance(conveyor, (i==T)) ){
  for( ; i < T; i++){
    col = pckindx[i] >> 16;
    pe  = pckindx[i] & 0xffff;
    if( !convey_push(ex, &col, pe) )
      break;
  }

  while(convey_pull(conveyor, &col, NULL) == convey_OK )
    lcounts[col]++;
}
```

#### Why it is in bale?
It is 


#### From the Book?


See apps/sssp_src/ for the all implementations.

