config const N=10; // number of updates
config const M=10; // size of table
// allocate main table and array of random ints
const Mspace = {0..M-1};
const D = Mspace dmapped Cyclic(startIdx=Mspace.low);
var A: [D] atomic int;
const Nspace = {0..(N*numLocales - 1)};
const D2 = Nspace dmapped Block(Nspace);
var rindex: [D2] int;

/* set up loop */
fillRandom(rindex, 208); // the 208 is a seed
forall r in rindex{
  r = mod(r, M);
}

/* main loop */
forall r in rindex{
  A[r].add(1); //atomic add
}

//We can write the main loop in a more node-centric way though. 

//coforall loc in Locales do on loc do{
//	forall i in rindex.localSubdomain(){
//		A[rindex[i]].add(1);
//}
//}

//Most economical of all, we can also write the main loop in this “vector” way:
//A[rindex].add(1);
