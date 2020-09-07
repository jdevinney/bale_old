use sparsemat::wall_seconds;
use sparsemat::SparseMat;
use std::fs::File;
use std::io::Write;
use std::io::Error;
use std::ops::Range;
use std::path::Path;

// A helper function for dumping only part of a big data structure.
// This should really go somewhere else than the delta_stepper lib.
pub fn display_ranges(max_disp: usize, num_items: usize) -> Vec<Range<usize>> {
        let mut ranges: Vec<Range<usize>> = Vec::new();
        if max_disp <= num_items && max_disp > 0 {
            ranges.push(0..max_disp/2);
            ranges.push(num_items-max_disp/2..num_items);
        } else {
            ranges.push(0..num_items);
        }
        //ranges.iter().enumerate()
        ranges
}

// Output structure for single-source shortest path
#[derive(Debug, Clone)]
pub struct SsspInfo {
    distance: Vec<f64>,
    source: usize,
    pub laptime: f64,
}

impl SsspInfo {
    // Dump output distances to a file
    pub fn dump(&self, max_disp: usize, filename: &str) -> Result<(),Error> { //2 outs but perm has 1??
        let path = Path::new(&filename);
        let mut file = File::create(path)?;
        write!(file, "vtx: dist\n")?;
        // for (i, r) in display_ranges(max_disp, self.distance.len()) {
        //     if i > 0 {
        //         writeln!("...");
        //     }
        for r in display_ranges(max_disp, self.distance.len()) {
            for v in r {  // should use zip 0-0
                write!(file, "{}: {}\n", v, self.distance[v])?;
            }
        }
        Ok(())
    }
}

// A potential edge relaxation to be examined
struct Request {
    w: usize,    // head of edge being relaxed
    dist: f64,   // new distance from source to w using that edge
//  v: usize,    // could include tail of edge (v,w) in request to build shortest path tree
}

// Struct for a single bucket, a circular doubly-linked list, and its iterator.
// The buckets live inside the BucketSearcher struct below.
struct Bucket {
    header: usize,   // which list elt is this bucket's header? (answer: nv+this bucket index)
    size: usize,     // number of vertices in this bucket
}
impl Bucket { 
    fn new(header: usize) -> Bucket { 
        let mut size = 0;
        Bucket{
            header,
            size,
        }
    }
}

// A struct and methods for all the data structures in delta stepping
// Each bucket is a circular doubly-linked list, linked by prev_elt and next_elt,
// indexed by an "elt" that is either a vertex in the bucket, or a "bucket header".
// A vertex's "elt" number is its vertex number.
// The header of bucket i (with reuse) is "elt" number num_buckets + i.

struct BucketSearcher<'a> {
    graph: &'a SparseMat,           // the graph being searched
    delta: f64,                     // width of a bucket
    num_buckets: usize,             // number of actual buckets, taking into account reuse
    tentative_dist: Vec<f64>,       // current tentative distance from source to this vtx
    prev_elt: Vec<usize>,           // back link in list of vertices in each bucket (including bucket header)
    next_elt: Vec<usize>,           // forward link in list of vertices in each bucket (including header)
    activated: Vec<bool>,           // has this vtx ever been activated? 
    vtx_bucket: Vec<Option<usize>>, // what bucket if any is this vtx in?
    buckets: Vec<Bucket>,           // the buckets themselves

    // iterator notes:
    // BucketSearcher should also have an iterator that updates the active bucket, skipping to the
    // next nonempty bucket (% num_buckets) and stopping when all buckets are empty.
    // (But I'm not going to put active_bucket in the struct.)
    // Also there should be an iterator that yields the vertices in a bucket,
    // and keeps working when its current vertex is removed from the bucket.
      
    // parallel notes: 
    // Vertices belong to PEs in round-robin order. 
    // Each PE has its own copy of every bucket, but only containing vertices it owns.
    // (Say nv >= 100 * num_buckets * THREADS. Then storage for buckets is relatively small.)
    // Arrays bucket_size and bucket_header are local to each PE.
    // Arrays indexed only by vtx are split among PEs (but don't have to be shared, I think): 
    //     tentative_dist, activated, bucket.
    // Arrays prev_elt and next_elt are also split, linking together only the vtxs on the local PE.
    // The bucket-header nodes at the end of prev_elt and next_elt have copies on each PE.
    // The parallelism all happens in relax_requests, where a request to relax an edge with head w
    // gets conveyed to the PE that owns w.
}

impl<'a> BucketSearcher<'a> { 

    // Create a bucket structure for a weighted graph
    fn new(graph: &SparseMat, delta: f64) -> BucketSearcher { // I guess this uses lifetime elision
        let nv = graph.numrows;
        let mut max_edge_len: f64 = 0.0;
        if let Some(edge_len) = &graph.value { // 0-0 yuck! and max_edge_len = edge_len.max(); fails
            for c in edge_len {
                if *c > max_edge_len { // 0-0 the compiler wants *c here, which I don't understand
                    max_edge_len = *c;
                }
            }
        }
        println!("max_edge_len = {}", max_edge_len);
        // upper bound on number of buckets we will ever need at the same time
        let num_buckets = (max_edge_len/delta).ceil() as usize + 1;
        // tentative distances all start out infinite, including source
        let mut tentative_dist = vec![f64::INFINITY; nv];    
        // circular linked lists have room for the bucket headers,
        // and initially every list is empty (every element points to itself).
        let mut prev_elt: Vec<usize> = (0..nv+num_buckets).collect();
        let mut next_elt: Vec<usize> = (0..nv+num_buckets).collect();
        // initially no vtx has ever been activated.
        let mut activated = vec![false; nv];
        // initially every vtx is in no bucket.
        let mut vtx_bucket = vec![Option::<usize>::None; nv];
        // initially every bucket contains no vertices.
        let mut buckets: Vec<Bucket> = Vec::new();
        for header in nv..num_buckets+nv {
            buckets.push(Bucket::new(header));
        }
        buckets[0].size += 1;  // just for test
        
        //barrier here

        BucketSearcher {
            graph,
            delta,
            num_buckets,
            tentative_dist,
            prev_elt,
            next_elt,
            activated,
            vtx_bucket,
            buckets,
        }
    }

    // Dump bucket structure state to a file
    fn dump(&self, max_disp: usize, filename: &str) -> Result<(),Error> { // why 2 outs when perm has 1?
        let path = Path::new(&filename);
        let mut file = File::create(path)?;
        let nv = self.graph.numrows;
        writeln!(file, "==========================================================")?;
        writeln!(file, "BucketSearcher: nv={}, num_buckets={}, delta={}", nv, self.num_buckets, self.delta)?;
        writeln!(file, "elt: prev_elt next_elt vtx_bucket activated tentative_dist")?;
        // for (i,r) in display_ranges(max_disp, nv) {
        //     if i > 0 {
        //         writeln!("...");
        //     }
        for r in display_ranges(max_disp, nv) {
            for v in r {  // should be firstlist(r).zip(secondlist(r)) 0-0
                writeln!(
                    file, 
                    "{}: {} {} {} {} {}", 
                    v,
                    self.prev_elt[v], 
                    self.next_elt[v], 
                    if let Some(b) = self.vtx_bucket[v] {
                        b.to_string()
                    } else { 
                        "N".to_string()
                    },
                    self.activated[v], 
                    self.tentative_dist[v]
                )?;
            }
        }
        for v in nv..nv+self.num_buckets {  // should use zip 0-0
            writeln!(file, "{}: {} {}", v, self.prev_elt[v], self.next_elt[v])?;
        }
        writeln!(file, "bucket (bucket_size): elt elt ...")?;
        // for elt in display_ranges(max_disp, self.num_buckets) {
        //     if i > 0 {
        //         writeln!("...")?;
        //     }
        for r in display_ranges(max_disp, self.num_buckets) {
            for b in r {  
                write!(file, "{} ({}):", b, self.buckets[b].size)?;
                // for e in self.buckets[b] {
                //     write!(file, " {}", e)?;
                // }
                write!(file, "\n")?;
            }
        }
        writeln!(file, "==========================================================")?;
        Ok(())
    }   

    // what bucket does a vtx with this tentative distance go in?
    fn home_bucket(&self, dist: f64) -> usize {
        ((dist/self.delta).floor() as usize) % self.num_buckets
    }

    // remove a vtx from the bucket it's in. (harmless if not in a bucket)
    fn remove_from_bucket(&mut self, w: usize) {
        if let Some(b) = self.vtx_bucket[w] {
            self.prev_elt[self.next_elt[w]] = self.prev_elt[w];
            self.next_elt[self.prev_elt[w]] = self.next_elt[w];
            self.prev_elt[w] = w;
            self.next_elt[w] = w;
            self.vtx_bucket[w] = None;
            self.buckets[b].size -= 1;
        }
    }

    // insert a vtx into a bucket it's not in. (vtx must not be in a bucket)
    fn place_in_bucket(&mut self, w:usize, new_bucket: usize) {
        assert!(self.vtx_bucket[w] == None);
        assert!(self.prev_elt[w] == w);
        assert!(self.next_elt[w] == w);
        self.prev_elt[w] = self.buckets[new_bucket].header;
        self.next_elt[w] = self.next_elt[self.buckets[new_bucket].header];
        self.prev_elt[self.next_elt[w]] = w;
        self.next_elt[self.prev_elt[w]] = w;
        self.vtx_bucket[w] = Some(new_bucket);
        self.buckets[new_bucket].size += 1;
    }

    // make a list of all relaxation requests from light edges with tails in bucket
    // we need an iterator over vertices in a bucket (that lets you remove a vtx in flight)
    fn find_light_requests(&self, _bucket: &Bucket) -> Vec<Request> { 
        let mut requests: Vec<Request> = Vec::new();
        requests.push(Request{w: 0, dist: 0.0}); //0-0
        requests
    }

    // make a list of all relaxation requests from heavy edges with tails on vtxlist
    fn find_heavy_requests(&self, _vtxlist: Vec<usize>) -> Vec<Request> { 
        let mut requests: Vec<Request> = Vec::new();
        requests.push(Request{w: 0, dist: 0.0}); //0-0
        requests
    }

    // return vertices in bucket that have not been activated(removed from an active bucket) before
    fn newly_active_vertices(&self, _bucket: &Bucket) -> Vec<usize> { 
        let mut new_vtxs: Vec<usize> = Vec::new(); 
        new_vtxs.push(0); // 0-0
        new_vtxs
    }

    // remove all vertices from the bucket and mark them activated
    fn empty_bucket(&mut self, _bucket: &Bucket) { 
        // 0-0
    }

    // relax all the requests from this phase (could be parallel)
    fn relax_requests(&mut self, requests: Vec<Request>) {
        // convey the request r=(w,d) to the PE that owns vtx w here, and have it call relax
        for r in requests {
            self.relax(r);
        }
        // barrier here (maybe also barrier at beginning of relax_requests?)
    }

    // relax an incoming edge to vtx r.w with new source distance r.dist, and rebucket r.w if necessary
    // this will be called by r.w's PE, so there is no race on tent[r.w] 
    fn relax(&mut self, r: Request) {
        if r.dist < self.tentative_dist[r.w] {
            let new_bucket = self.home_bucket(r.dist);
            if let Some(old_bucket) = self.vtx_bucket[r.w] {
                // r.w was in a bucket,
                if old_bucket != new_bucket {
                    // move r.w from old to new bucket
                    self.remove_from_bucket(r.w);
                    self.place_in_bucket(r.w, new_bucket);
                }
            } else {
                // r.w was not in a bucket, put it in new bucket
                self.place_in_bucket(r.w, new_bucket);
            }
            self.tentative_dist[r.w] = r.dist;
        }
    }
}

pub trait DeltaStepping {
    fn delta_stepping(&self, source: usize) -> SsspInfo;
    fn check_result(&self, info: &SsspInfo, dump_files: bool) -> bool;
}

impl DeltaStepping for SparseMat {

    /// This implements the sequential AGI variant of delta stepping.
    /// # Argument: source vertex
    fn delta_stepping(&self, source: usize) -> SsspInfo {
        assert!(self.numrows == self.numcols);
        assert!(source < self.numrows);

        let t1 = wall_seconds().expect("wall second error");

        // choose a value for delta, the bucket width
        let delta = 1.0; // 0-0 need to fix this, probably 1/(max degree)
        
        // initialize buckets, activated flags, etc.
        let mut searcher = BucketSearcher::new(&self, delta);

        // use relax to set tent(source) to 0, which also puts it in bucket 0
        searcher.relax(Request{w: source, dist: 0.0});

        searcher.dump(20, "buckets.out").expect("bucket dump failed");

        // outer loop: for each nonempty bucket in order ...
        // for active_bucket in searcher  
        let active_bucket = &searcher.buckets[0]; // just so this will compile without the iterator
        {
            let mut removed: Vec<usize> = Vec::new(); // vertices removed from active bucket, R in paper
        
            // middle loop: while active bucket (B[i] in paper) is not empty ...
            while active_bucket.size > 0 { 

                // find light edges with tails in active bucket;
                // empty active bucket, keeping a set "removed" of unique vtxs removed from this bucket;
                // relax light edges, which may put some removed and other vtxs into the active bucket.

                let requests = searcher.find_light_requests(active_bucket); 
                removed.append(&mut searcher.newly_active_vertices(active_bucket));
                searcher.empty_bucket(active_bucket);
                searcher.relax_requests(requests);
                
                break; // 0-0 so as not to have an infinite loop while debugging
            } // end of middle looop
            
            // relax heavy edges with tails in removed set, which cannot add vtxs to active bucket
            let requests = searcher.find_heavy_requests(removed);
            searcher.relax_requests(requests);

        } // end of outer loop

        // return the info struct, which will now own the distance array
        SsspInfo {
            distance: searcher.tentative_dist,
            source: source,
            laptime: wall_seconds().expect("wall second error") - t1,
        }
    }

    /// check the result of delta stepping
    ///
    /// # Arguments
    /// * info data from the run to check
    /// * dump_files debugging flag
    fn check_result(&self, info: &SsspInfo, dump_files: bool) -> bool {
        println!("check_result: source is {}, dump_files is {}", info.source, dump_files);
        if dump_files {
            info.dump(20, "dist.out").expect("cannot write dist.out");
        }
        for v in 0..self.numrows {
            if v == info.source && info.distance[v] != 0.0 {
                return false;
            } 
        }
        true
    }
}
