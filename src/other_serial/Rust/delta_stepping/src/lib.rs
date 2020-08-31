use sparsemat::wall_seconds;
use sparsemat::SparseMat;

#[derive(Debug, Clone)]
pub struct SsspInfo {
    distance: Vec<f64>,
    source: usize,
    pub laptime: f64,
}

impl SsspInfo {
    fn new(source: usize, nv: usize) -> SsspInfo {
        let distance = vec![f64::INFINITY; nv];
        SsspInfo {
            distance,
            source: source,
            laptime: 0.0,
        }
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
        for v in 0..self.numrows {
            if v == info.source {
                assert!(info.distance[v] == 0.0);
            } else {
                assert!(f64::is_infinite(info.distance[v]));
            }
        }
        println!("check_result assertions passed, dump_files is {}", dump_files);
        true
    }
}
