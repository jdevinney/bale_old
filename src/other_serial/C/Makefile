
#CFLAGS= -O3
#CFLAGS= -std=c99
LFLAGS= -lm -largp
CC= gcc -Wall ${CFLAGS}

SRC= spmat_utils.h spmat_utils.c 
  SRC+=histo.c 
  SRC+=ig.c 
  SRC+=randperm.c 
  SRC+=toposort.c 
  SRC+=triangle.c 
  SRC+=unionfind.c 
  SRC+=transpose_matrix.c
  SRC+=permute_matrix.c
  SRC+=sssp.c
  SRC+=sssp_dijsktra.c
  SRC+=sssp_bellmanford.c
  SRC+=demo_spmat.c

EXE= histo 
  EXE+=ig 
  EXE+=randperm 
  EXE+=topo 
  EXE+=triangles 
  EXE+=unionfind
  EXE+=transpose_matrix
  EXE+=permute_matrix
  EXE+=sssp
  EXE+=demo_spmat

DOC= mainpage.h Doxyfile

OBJ= spmat_utils.o

GRAPHS= demo.mm 10_ring.mm 2comps.mm

all: std_options.o spmat_utils.o ${EXE}

std_options.c: std_options.c std_options.h
	${CC} -c std_options.c

spmat_utils.o: spmat_utils.c spmat_utils.h
	${CC} -c spmat_utils.c -lm

demo_spmat: demo_spmat.c spmat_utils.o 
	${CC} demo_spmat.c spmat_utils.o ${LFLAGS} -o demo_spmat

histo: histo.c spmat_utils.o 
	${CC} histo.c spmat_utils.o ${LFLAGS} -o histo

ig: ig.c spmat_utils.o 
	${CC} ig.c spmat_utils.o ${LFLAGS} -o ig

randperm: randperm.c spmat_utils.o
	${CC} randperm.c spmat_utils.o ${LFLAGS} -o randperm

topo: toposort.c spmat_utils.o
	${CC} toposort.c spmat_utils.o ${LFLAGS} -o topo

triangles: triangle.c spmat_utils.o
	${CC} triangle.c spmat_utils.o ${LFLAGS} -o triangles

unionfind: unionfind.c spmat_utils.o
	${CC} unionfind.c spmat_utils.o ${LFLAGS} -o unionfind

permute_matrix: permute_matrix.c spmat_utils.o
	${CC} permute_matrix.c spmat_utils.o ${LFLAGS} -o permute_matrix

transpose_matrix: transpose_matrix.c spmat_utils.o
	${CC} transpose_matrix.c spmat_utils.o std_options.o ${LFLAGS} -o transpose_matrix

sssp: sssp.c sssp_dijsktra.c sssp_delta.c sssp_bellmanford.c spmat_utils.o
	${CC} sssp.c sssp_dijsktra.c sssp_delta.c sssp_bellmanford.c spmat_utils.o ${LFLAGS} -o sssp

doc:
	doxygen

ci:
	ci Makefile  ${SRC} ${GRAPHS} ${DOC} runall.sh
co: 
	co ${SRC} ${GRAPHS} ${DOC} runall.sh

tar: 
	mkdir -p Bale_C
	cp Makefile ${SRC} ${GRAPHS} ${DOC} runall.sh Bale_C
	tar cfv Bale_C.tar Bale_C
	rm -r Bale_C

clean:
	rm -f *.out ${EXE} ${OBJ}
