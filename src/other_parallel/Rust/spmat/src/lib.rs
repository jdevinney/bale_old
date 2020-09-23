#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Sparse Matrix Library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For licence information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use crate::perm::Perm;
use convey_hpc::collect::{IVal, PType, ValueCollect};
use convey_hpc::session::ConveySession;
use convey_hpc::Convey;
use serde::de::DeserializeOwned;
use serde::ser::Serialize;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use std::time::{SystemTime, UNIX_EPOCH};
/// A routine to give access to the wall clock timer on most UNIX-like systems.
///    Uses rust's SystemTime.
pub fn wall_seconds() -> f64 {
    let n = SystemTime::now().duration_since(UNIX_EPOCH).unwrap();
    (n.as_secs() as f64) + (n.as_micros() as f64) * 1.0e-6
}

/// A sparse matrix.  Data is stored as a num_rows_this_rank size
/// vector of offsets in the nonzero vector, which has column numbers
/// for each row.
pub struct SparseMat {
    numrows: usize, // the total number of rows in the matrix
    numrows_this_rank: usize,
    numcols: usize, // the nonzeros have values between 0 and numcols
    nnz: usize,     // total number of nonzeros in the matrix
    nnz_this_rank: usize,
    /// the row offsets into the array of nonzeros, size is nrows+1,
    /// offsets[nrows] is nnz
    offset: Vec<usize>,
    nonzero: Vec<usize>,    // the global array of nonzeros
    convey: Option<Convey>, // a conveyor for our use
}

impl std::fmt::Debug for SparseMat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let debug_rows = self.numrows_this_rank.min(10);
        let mut rows: Vec<(usize, Vec<usize>)> = Vec::new();
        for i in 0..debug_rows {
            let start = self.offset[i];
            let end = self.offset[i + 1];
            let debug_cols = (end - start).min(10);
            let mut row: Vec<usize> = Vec::new();
            for j in start..start + debug_cols {
                row.push(self.nonzero[j]);
            }
            rows.push((end - start, row));
        }
        f.debug_struct("SparseMat")
            .field("rank", &self.my_rank())
            .field("numrows", &self.numrows)
            .field("numcols", &self.numcols)
            .field("nnz", &self.nnz)
            .field("nnz.len()", &self.nonzero.len())
            .field("offset.len()", &self.offset.len())
            .field("first_part", &rows)
            .finish()
    }
}

impl SparseMat {
    /// new creates a distributed sparse matrix across all ranks using
    /// number of rows, number of columns, and number of rows on
    /// current rank
    pub fn new(numrows: usize, numcols: usize, nnz_this_rank: usize) -> Self {
        let convey = Convey::new().expect("fixme");
        let total_nnz = convey.reduce_sum(nnz_this_rank);
        SparseMat {
            numrows: numrows,
            numrows_this_rank: convey.per_my_rank(numrows),
            numcols: numcols,
            nnz: total_nnz,
            nnz_this_rank: nnz_this_rank,
            offset: vec![0; convey.per_my_rank(numrows) + 1],
            nonzero: vec![0; nnz_this_rank],
            convey: Some(convey),
        }
    }

    /// new_local creates a local sparse matrix on each rank using
    /// number of rows, number of columns, and number of rows on
    /// current rank
    pub fn new_local(numrows: usize, numcols: usize, nnz_this_rank: usize) -> Self {
        SparseMat {
            numrows: numrows,
            numrows_this_rank: numrows,
            numcols: numcols,
            nnz: nnz_this_rank,
            nnz_this_rank: nnz_this_rank,
            offset: vec![0; numrows + 1],
            nonzero: vec![0; nnz_this_rank],
            convey: None,
        }
    }

    /// Convenience function returning my_rank
    pub fn my_rank(&self) -> usize {
        if let Some(convey) = &self.convey {
            convey.my_rank
        } else {
            0
        }
    }

    /// Convenience function returning num_ranks
    pub fn num_ranks(&self) -> usize {
        if let Some(convey) = &self.convey {
            convey.num_ranks
        } else {
            0
        }
    }

    /// Convenience function returning the number of elements on
    /// my_rank for a total number of elements
    pub fn per_my_rank(&self, n: usize) -> usize {
        if let Some(convey) = &self.convey {
            convey.per_my_rank(n)
        } else {
            n
        }
    }

    /// Convenience function returning the offset and rank in a
    /// distributed sparsemat for a given index
    pub fn offset_rank(&self, n: usize) -> (usize, usize) {
        if let Some(convey) = &self.convey {
            convey.offset_rank(n)
        } else {
            (n, 0)
        }
    }

    /// create a session without a pull_fn
    pub fn barrier(&self) {
        if let Some(convey) = &self.convey {
            convey.barrier()
        }
    }

    /// create a session without a pull_fn
    pub fn begin<'a, T: Copy + Serialize + DeserializeOwned>(
        &'a self,
        pull_fn: impl FnMut(T, usize) + 'a,
    ) -> ConveySession<'a, T> {
        if let Some(convey) = &self.convey {
            convey.begin(pull_fn)
        } else {
            unimplemented!("Conveyors not implemented for local SparseMat");
        }
    }

    /// convenience function to create a session
    pub fn session<'a, T: Copy + Serialize + DeserializeOwned>(&'a self) -> ConveySession<'a, T> {
        if let Some(convey) = &self.convey {
            convey.session()
        } else {
            unimplemented!("Conveyors not implemented for local SparseMat");
        }
    }

    /// execute a simple function
    pub fn simple<T, I>(&self, map: I, pull_fn: impl FnMut(T, usize)) -> ()
    where
        T: Copy + Serialize + DeserializeOwned,
        I: Iterator<Item = (T, usize)>,
    {
        if let Some(convey) = &self.convey {
            convey.simple(map, pull_fn)
        } else {
            unimplemented!("Conveyors not implemented for local SparseMat");
        }
    }

    /// execute a simple function with return
    pub fn simple_return<T1, T2, I>(
        &self,
        map: I,
        pull_fn: impl FnMut(T1) -> T2,
        ret_fn: impl FnMut(usize, T2),
    ) -> ()
    where
        T1: Copy + Serialize + DeserializeOwned,
        T2: Copy + Serialize + DeserializeOwned,
        I: Iterator<Item = ((usize, T1), usize)>,
    {
        if let Some(convey) = &self.convey {
            convey.simple_return(map, pull_fn, ret_fn)
        } else {
            unimplemented!("Conveyors not implemented for local SparseMat");
        }
    }
    /// dumps a sparse matrix to a file in a ASCII format
    /// # Arguments
    /// * maxrows the number of rows that are written, 0 means everything,
    ///           otherwise write the first and last maxrows/2 rows
    /// * filename the filename to written to
    pub fn dump(&self, maxrows: usize, filename: &str) -> Result<(), std::io::Error> {
        let path = Path::new(&filename);
        let mut file = File::create(path)?;

        let use_maxrow: bool = if maxrows < self.numrows && maxrows > 0 {
            true
        } else {
            false
        };
        let start_row = if use_maxrow {
            self.numrows - maxrows / 2
        } else {
            self.numrows
        };
        let stop_row = if use_maxrow {
            maxrows / 2
        } else {
            self.numrows
        };

        writeln!(file, "\n--------- offsets:")?;
        for off in &self.offset[0..stop_row] {
            write!(file, "{} ", off)?;
        }
        if use_maxrow {
            write!(file, " ... ")?;
            for off in &self.offset[start_row..self.numrows] {
                write!(file, "{}", off)?;
            }
        }

        writeln!(file, "\n--------- nonzeros:")?;
        for i in 0..stop_row {
            for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                writeln!(file, "{} {}", i, nz)?;
            }
        }
        if use_maxrow {
            write!(file, ".\n.\n.\n")?;
            for i in start_row..self.numrows {
                for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                    writeln!(file, "{} {}", i, nz)?;
                }
            }
        }
        Ok(())
    }

    /// writes a sparse matrix to a file in a MatrixMarket ASCII formats
    /// # Arguments
    /// * filename the filename to written to
    pub fn write_mm_file(&self, filename: &str) -> Result<(), std::io::Error> {
        let path = Path::new(&filename);
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        self.write_mm(&mut writer)
    }

    /// write a matrix in matrix_market format
    ///  TODO: make this work in parallel
    pub fn write_mm<W>(&self, writer: &mut W) -> Result<(), std::io::Error>
    where
        W: Write,
    {
        writeln!(writer, "%%MasterMarket matrix coordinate position")?;
        writeln!(writer, "{} {} {}", self.numrows, self.numcols, self.nnz)?;

        for rank in 0..self.num_ranks() {
            if rank == self.my_rank() {
                for i in 0..self.numrows {
                    for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                        writeln!(writer, "{} {}", i + 1, nz + 1)?;
                    }
                }
            }
            self.barrier();
        }
        Ok(())
    }

    /// perdicate: is the SparseMat lower triangular
    pub fn is_lower_triangular(&self, unit_diagonal: bool) -> bool {
        let mut lower_cnt = 0;
        let mut diag_missing_cnt = 0;
        for i in 0..self.numrows_this_rank {
            let global_row = i * self.num_ranks() + self.my_rank();
            let mut pivot = false;
            for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                if *nz < global_row {
                    lower_cnt += 1;
                } else if *nz == global_row {
                    pivot = true;
                }
            }
            if !pivot {
                diag_missing_cnt += 1;
            }
        }
        let total_lower = self.reduce_sum(lower_cnt);
        let total_diag_missing = if unit_diagonal {
            self.reduce_sum(diag_missing_cnt)
        } else {
            0
        };
        total_lower != 0 || total_diag_missing != 0
    }

    /// perdicate: is the SparseMat upper triangular
    pub fn is_upper_triangular(&self, unit_diagonal: bool) -> bool {
        let mut lower_cnt = 0;
        let mut diag_missing_cnt = 0;
        for i in 0..self.numrows_this_rank {
            let global_row = i * self.num_ranks() + self.my_rank();
            let mut pivot = false;
            for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                if *nz > global_row {
                    lower_cnt += 1;
                } else if *nz == global_row {
                    pivot = true;
                }
            }
            if !pivot {
                diag_missing_cnt += 1;
            }
        }
        let total_lower = self.reduce_sum(lower_cnt);
        let total_diag_missing = if unit_diagonal {
            self.reduce_sum(diag_missing_cnt)
        } else {
            0
        };
        total_lower != 0 || total_diag_missing != 0
    }

    /// createa a new SparseMat which is the permutation of rows and
    /// columns as defined by given permutations
    pub fn permute(&self, rperminv: &Perm, cperminv: &Perm) -> Self {
        let mut rowcounts = vec![0_usize; self.numrows_this_rank];
        assert_eq!(rperminv.len(), self.numrows_this_rank);
        assert_eq!(cperminv.len(), self.numrows_this_rank);
        assert!(rperminv.is_perm());
        assert!(cperminv.is_perm());
        let mut nnz = 0_usize;
        self.simple(
            (0..rperminv.len()).map(|idx| {
                let (offset, rank) = self.offset_rank(rperminv.entry(idx));
                let cnt = self.offset[idx + 1] - self.offset[idx];
                ((offset, cnt), rank)
            }),
            |item: (usize, usize), _from_rank| {
                rowcounts[item.0] = item.1;
                nnz += item.1;
            },
        );

        assert_eq!(self.nnz, self.reduce_sum(nnz));
        let mut permuted = SparseMat::new(self.numrows, self.numcols, nnz);
        permuted.offset[0] = 0;

        // step 2: set up new offsets
        for i in 1..permuted.numrows_this_rank + 1 {
            permuted.offset[i] = permuted.offset[i - 1] + rowcounts[i - 1];
        }
        //println!("{:?}{:?}", permuted, permuted.offset);
        assert_eq!(
            permuted.offset[permuted.numrows_this_rank],
            permuted.nnz_this_rank
        );
        // step 3: distrbute nonzeros
        let mut wrkoff = vec![0_usize; self.numrows_this_rank];
        {
            let mut session = self.begin(|item: (usize, usize), _from_rank| {
                let index = permuted.offset[item.0] + wrkoff[item.0];
                permuted.nonzero[index] = item.1;
                wrkoff[item.0] += 1;
            });
            let mut row = 0;
            for i in 0..self.nnz_this_rank {
                while i == self.offset[row + 1] {
                    row += 1;
                }
                let (offset, rank) = self.offset_rank(rperminv.entry(row));
                session.push((offset, self.nonzero[i]), rank);
            }
            session.finish();
        }

        // step 4: do column permutation (essentailly indexgather)
        {
            let my_nnz = permuted.nnz_this_rank;
            let cloned_nonzero = permuted.nonzero.clone();
            self.simple_return(
                (0..my_nnz).map(|i| {
                    let (offset, rank) = self.offset_rank(cloned_nonzero[i]);
                    ((i, offset), rank)
                }),
                |item: usize| cperminv.entry(item),
                |index: usize, val: usize| {
                    permuted.nonzero[index] = val;
                },
            );
        }
        permuted
    }

    /// Create a new SparseMat which is the transpose of the given SparseMat
    pub fn transpose(&self) -> Self {
        let mut colcnt = vec![0_usize; self.per_my_rank(self.numrows)];
        let mut nnz = 0_usize;

        // distributed calculation of column counts, with resulting local nonzeros
        self.simple(
            self.nonzero.iter().map(|nz| self.offset_rank(*nz)),
            |item: usize, _from_rank| {
                colcnt[item] += 1;
                nnz += 1;
            },
        );

        let mut trans = SparseMat::new(self.numcols, self.numrows, nnz);
        if self.nnz != trans.nnz {
            println!("self: {:?} colcnt: {:?}", self, colcnt);
        }
        assert_eq!(self.nnz, trans.nnz);

        trans.offset[0] = 0;
        for idx in 1..self.offset.len() {
            trans.offset[idx] = trans.offset[idx - 1] + colcnt[idx - 1];
        }

        let mut wrkoff = vec![0_usize; trans.numrows];
        {
            let mut session = self.begin(|item: (usize, usize), _from_rank| {
                let index = trans.offset[item.0] + wrkoff[item.0];
                trans.nonzero[index] = item.1;
                wrkoff[item.0] += 1;
            });
            let mut row = 0;
            for i in 0..self.nonzero.len() {
                while i == self.offset[row + 1] {
                    row += 1;
                }
                let (offset, rank) = self.offset_rank(self.nonzero[i]);
                let tcol = row * self.num_ranks() + self.my_rank();
                session.push((offset, tcol), rank);
            }
            session.finish();
        }
        trans
    }

    /// Predicate: are two SparseMats equal?
    pub fn compare(&self, other: &SparseMat) -> bool {
        if (self.numcols != other.numcols)
            || (self.numrows != other.numrows)
            || (self.numrows_this_rank != other.numrows_this_rank)
            || (self.nnz != other.nnz)
            || (self.nnz_this_rank != other.nnz_this_rank)
        {
            println!("counts differ {:?} {:?}", self, other);
            return false;
        }
        if self.offset != other.offset {
            println!("offsets differ {:?} {:?}", self, other);
            return false;
        }
        // need to sort nonzeros before compare
        let mut self_sorted = self.nonzero.clone();
        let mut other_sorted = other.nonzero.clone();
        self_sorted.sort();
        other_sorted.sort();
        if self_sorted != other_sorted {
            println!("nonzero differ {:?} {:?}", self, other);
            return false;
        }
        true
    }

    /// Create a new SparseMat which is the sum of 2 given
    pub fn add(&self, other: &SparseMat) -> Self {
        // need to check for errors
        let mut sum = SparseMat::new(
            self.numrows,
            self.numcols,
            self.nnz_this_rank + other.nnz_this_rank,
        );
        let mut nnz = 0;
        for row in 0..self.numrows_this_rank {
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                sum.nonzero[nnz] = *nz;
                nnz += 1;
            }
            for nz in &other.nonzero[other.offset[row]..other.offset[row + 1]] {
                sum.nonzero[nnz] = *nz;
                nnz += 1;
            }
            sum.offset[row + 1] = nnz;
        }
        sum
    }

    /// Generate an erdos_renyi random sparse matrix
    pub fn gen_erdos_renyi_graph(
        num_vert: usize,
        prob: f64,
        unit_diag: bool,
        mode: u8,
        seed: i64,
    ) -> Self {
        if mode == 0 {
            let upper = SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, false, seed);
            let lower = upper.transpose();
            upper.add(&lower)
        } else if mode == 1 {
            SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, true, seed)
        } else if mode == 2 {
            SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, false, seed)
        } else {
            let upper = SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, false, seed);
            let lower = SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, true, seed + 1);
            upper.add(&lower)
        }
    }

    /// Worker function: generate an erdos_renyi triangular random sparse matrix
    fn erdos_renyi_tri(
        num_vert: usize,
        prob: f64,
        unit_diag: bool,
        lower: bool,
        _seed: i64,
    ) -> Self {
        let mut tri = SparseMat::new(num_vert, num_vert, 0);
        let l_max = (std::u32::MAX as f64).ln();
        let d = (1.0 - prob).ln();
        let numrows = num_vert;
        let mut row = tri.my_rank();
        let mut col = if lower { 0 } else { tri.my_rank() + 1 };
        let num_ranks = tri.num_ranks();

        // if we are upper triangular and unit_diag, we need to set
        // the first row
        if unit_diag && !lower {
            tri.nonzero.push(tri.my_rank());
        }
        while row < num_vert {
            let mut r = rand::random::<u32>();
            if r == std::u32::MAX {
                r -= 1;
            }

            col += (1.0 + ((((std::u32::MAX - r) as f64).ln() - l_max) / d).floor()) as usize;
            while (col >= (if lower { row } else { numrows })) && row < numrows {
                if lower {
                    if unit_diag {
                        tri.nonzero.push(row);
                    }
                    col = col - row;
                    row += num_ranks;
                    tri.offset[row / num_ranks] = tri.nonzero.len();
                } else {
                    row += num_ranks;
                    col = row + 1 + col - numrows;
                    tri.offset[row / num_ranks] = tri.nonzero.len();
                    if (row < numrows) && unit_diag {
                        tri.nonzero.push(row);
                    }
                }
            }
            if row < numrows {
                tri.nonzero.push(col);
            }
        }
        tri.offset[row / num_ranks] = tri.nonzero.len();
        tri.nnz_this_rank = tri.nonzero.len();
        tri.nnz = tri.reduce_sum(tri.nnz_this_rank);
        tri
    }
}

/// implement value colectives accessed through the Sparse matrix
macro_rules! decl {
    ($ty: ident) => {
        impl ValueCollect<$ty> for SparseMat {
            type T = $ty;
            fn reduce(
                &self,
                value: $ty,
                combine_fn: impl Fn($ty, $ty) -> Self::T,
                style: PType,
                init: IVal,
            ) -> Self::T {
                if let Some(convey) = &self.convey {
                    convey.reduce(value, combine_fn, style, init)
                } else {
                    value
                }
            }
        }
    };
}

// For all standard types
decl!(u128);
decl!(u64);
decl!(u32);
decl!(u16);
decl!(u8);
decl!(usize);

decl!(i128);
decl!(i64);
decl!(i32);
decl!(i16);
decl!(i8);
decl!(isize);

decl!(f32);
decl!(f64);

/// Bring in the permutation library
pub mod perm;

#[cfg(test)]
mod tests {
    use super::SparseMat;
    use crate::Perm;
    use convey_hpc::testing_support::TestingMutex;
    #[test]
    fn rand_mat_tri_upper() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, false, false, 0);
        assert_eq!(mat.is_lower_triangular(false), false);
        assert_eq!(mat.is_upper_triangular(false), true);
        assert_eq!(mat.my_rank(), mutex.convey.my_rank);
    }
    #[test]
    fn rand_mat_tri_lower() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, false, true, 0);
        assert_eq!(mat.is_lower_triangular(false), true);
        assert_eq!(mat.is_upper_triangular(false), false);
        assert_eq!(mat.my_rank(), mutex.convey.my_rank);
    }
    #[test]
    fn rand_mat_tri_upper_unit() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, true, false, 0);
        assert_eq!(mat.is_lower_triangular(true), false);
        assert_eq!(mat.is_upper_triangular(true), true);
        assert_eq!(mat.my_rank(), mutex.convey.my_rank);
    }
    #[test]
    fn rand_mat_tri_lower_unit() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, true, true, 0);
        assert_eq!(mat.is_lower_triangular(true), true);
        assert_eq!(mat.is_upper_triangular(true), false);
        assert_eq!(mat.my_rank(), mutex.convey.my_rank);
    }
    #[test]
    fn rand_mat_perm() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(1000, 0.05, true, false, 0);
        let rperm = Perm::random(1000, 0);
        let cperm = Perm::random(1000, 0);
        let _permuted = mat.permute(&rperm, &cperm);
        assert_eq!(mat.is_lower_triangular(true), false);
        assert_eq!(mat.is_upper_triangular(true), true);
        //assert_eq!(permuted.is_lower_triangular(true), false);
        //assert_eq!(permuted.is_upper_triangular(true), false);
        assert_eq!(mat.my_rank(), mutex.convey.my_rank);
    }
    #[test]
    fn transpose1() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(1000, 0.05, true, false, 0);
        let tmat = mat.transpose();
        let ttmat = tmat.transpose();
        println!("mat:{:?}", mat);
        println!("tmat:{:?}", tmat);
        println!("ttmat:{:?}", ttmat);
        assert_eq!(mat.compare(&tmat), false);
        assert_eq!(mat.compare(&ttmat), true);
        assert_eq!(mat.my_rank(), mutex.convey.my_rank);
    }
}
