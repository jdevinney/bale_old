#ifndef SPMAT_ENUMS_H
#define SPMAT_ENUMS_H

typedef enum graph_model {FLAT, GEOMETRIC, KRONECKER} graph_model;
typedef enum edge_type {DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED} edge_type;
typedef enum self_loops {NOLOOPS, LOOPS} self_loops;
#endif
