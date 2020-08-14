config const M = 20; //table size
config const N = 10; //num reads

const Nspace = {0..N*numLocales-1};
const Mspace = {0..M-1};

const D = Mspace dmapped Cyclic(startIdx=Mspace.low);
var T: [D] int;

const D2 = Nspace dmapped Block(Nspace);
var rindex: [D2] int;
var tmp: [D2] int;

fillRandom(rindex, 208); //208 is just the seed
rindex = mod(rindex, M);

T = {0..#M}; // identity permutation

/* MAIN LOOP */
forall i in rindex.domain{
  tmp[i] = T[rindex[i]];
}

//We could switch the main loop to this…

//coforall loc in Locales do on loc do{
//	var inds = D2.localSubdomain();
//	forall i in inds{
//		tmp[i] = T[rindex[i]];
//}
//}

//We can check for success with this line:
if(!tmp.equals(rindex)){
  writeln("Error!");
}
