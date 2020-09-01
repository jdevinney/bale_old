use sparsemat::wall_seconds;
use sparsemat::SparseMat;
use std::fs::File;
use std::io::Write;
use std::io::Error;
use std::ops::Range;
use std::path::Path;

// Output structure for single-source shortest path
#[derive(Debug, Clone)]
pub struct SsspInfo {
    distance: Vec<f64>,
    source: usize,
    pub laptime: f64,
}

impl SsspInfo {
    // Create an output structure with all distances infinite
    fn new(source: usize, nv: usize) -> SsspInfo {
        let distance = vec![f64::INFINITY; nv];
        SsspInfo {
            distance,
            source: source,
            laptime: 0.0,
        }
    }
    
    // Dump output distances to a file
    pub fn dump(&self, maxdisp: usize, filename: &str) -> Result<(),Error> { //jg: why 2 outs when perm has 1?
        let path = Path::new(&filename);
        let mut file = File::create(path)?;

        write!(file, "vtx: dist\n")?;

        let mut ranges: Vec<Range<usize>> = Vec::new();
        if maxdisp <= self.distance.len() && maxdisp > 0 {
            ranges.push(0..maxdisp/2);
            ranges.push(self.distance.len()-maxdisp/2..self.distance.len());
        } else {
            ranges.push(0..self.distance.len());
        }
        for r in ranges {
            for v in r { 
                write!(file, "{}: {}\n", v, self.distance[v])?;
            }
        }
        Ok(())
    }

}

pub trait DeltaStepping {
    fn delta_stepping(&self, source: usize) -> SsspInfo;
    fn check_result(&self, info: &SsspInfo, dump_files: bool) -> bool;
}

impl DeltaStepping for SparseMat {
    /// This routine implements the agi variant of delta stepping.
    /// # Arguments: source vertex
    fn delta_stepping(&self, source: usize) -> SsspInfo {
        assert!(self.numrows == self.numcols);
        let nv = self.numrows;
        assert!(source < nv);
        let mut ret = SsspInfo::new(source, nv);
        ret.distance[source] = 0.0;
        let t1 = wall_seconds().expect("wall second error");
        ret.laptime = wall_seconds().expect("wall second error") - t1;
        ret
    }

    /// check the result of delta stepping
    ///
    /// For now, just check the initial state of the tentative distances.
    /// # Arguments
    /// * info data from the run to check
    /// * dump_files debugging flag
    fn check_result(&self, info: &SsspInfo, dump_files: bool) -> bool {
        if dump_files {
            info.dump(20, "dist.out").expect("cannot write dist.out");
        }
        for v in 0..self.numrows {
            if v == info.source {
                assert!(info.distance[v] == 0.0);
            } else {
                assert!(f64::is_infinite(info.distance[v]));
            }
        }
        println!("check_result assertions passed, source is {}, dump_files is {}", info.source, dump_files);
        true
    }
}
