/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
//
//
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
//
*****************************************************************/

use clap::{App, Arg};
use convey_hpc::Convey;
use toposort::generate_toposort_input;
use toposort::TopoSort;

/*
 * Demo application that finds an upper triangular form for a matrix.
 * That is, we are given a matrix that is a random row and column permutation
 * of a an upper triangular matrix (with ones on the diagonal).
 * This algorithm finds a row and column permutation that would return it
 * to an upper triangular form.
 */

/* toposort_page Topologically sort a morally upper triangular matrix.

   First we generate the problem by generating an upper triangular matrix
   and applying row and column permutations.

   The output of toposort is a row and a column permutation that, if applied,
   would result in an upper triangular matrix.

   We set the row and column permutations,  rperm and cperm, one pivot at a time.

   N = number of rows
   for( pos=N-1; pos > 0; pos-- ) {
     pick a row, r, with a single nonzero, c.
     say (r,c) is the pivot and set rperm[pos] = r and cprem[pos] = c
     Note: a pivot always exists because the matrix is morally upper tri.

     cross out that row r and col c
   }

   Meaning of cross out:
   Rather than changing the matrix by deleting rows and column and then searching the
   new matrix for the next pivot.  We do the obvious thing of keeping row counts, where
   rowcnt[i] is the number of non-zeros in row i and we use a really cool trick
   of keeping the sum of the live column indices for the non-zeros in each row.
   That is, rowsum[i] is the sum of the column indices, not the sum of the non-zero elements,
   for the non-zeros in row i.  To "delete a column" one decrements the rowcnt by one and
   the rowsum by the corrsponding column index.
   The cool trick is that, when the rowcnt gets to one, the rowsum is the column that is left.
*/

fn main() {
    let matches = App::new("TopoSort")
        .version("0.1.0")
        .about("Implements a test of TopoSort")
        .arg(
            Arg::with_name("numrows_per_rank")
                .short("n")
                .long("numrows_per_rank")
                .takes_value(true)
                .help("The number of rows per rank"),
        )
        .arg(
            Arg::with_name("nz_per_row")
                .short("Z")
                .long("nz_per_row")
                .takes_value(true)
                .help("Number of non-zeros per row"),
        )
        .arg(
            Arg::with_name("er_prob")
                .short("e")
                .long("erdos_renyi_prob")
                .takes_value(true)
                .help("Probability of an edge in the erdos renyi graph"),
        )
        .arg(
            Arg::with_name("dump_files")
                .short("d")
                .long("dump_files")
                .takes_value(false)
                .help("Produce short dumps as the algorithm progresses"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .long("quiet")
                .takes_value(false)
                .help("produce less chatty output"),
        )
        .get_matches();

    // input args, just constants for now
    let numrows_per_rank: usize = matches
        .value_of("numrows_per_rank")
        .unwrap_or("20000")
        .parse()
        .expect("numrows_per_rank: not an integer");
    let mut erdos_renyi_prob: f64 = matches
        .value_of("er_prob")
        .unwrap_or("0.0")
        .parse()
        .expect("er_prob: not a float");
    let mut nz_per_row: f64 = matches
        .value_of("nz_per_row")
        .unwrap_or("35")
        .parse()
        .expect("nz_per_row: not a float");

    let convey = Convey::new().expect("Conveyor system initializiotn failed");
    let seed = 12346;
    let quiet = matches.is_present("quiet") || convey.my_rank > 0;
    let dump_files = matches.is_present("dump_files");

    let numrows = numrows_per_rank * convey.num_ranks;

    if erdos_renyi_prob == 0.0 {
        // use nz_per_row to get erdos_renyi_prob
        erdos_renyi_prob = (2.0 * nz_per_row) / (numrows - 1) as f64;
        if erdos_renyi_prob > 1.0 {
            erdos_renyi_prob = 1.0;
        }
    } else {
        // use erdos_renyi_prob to get nz_per_row
        nz_per_row = erdos_renyi_prob * numrows as f64
    }

    if !quiet {
        println!("Running toposort on {} ranks", convey.num_ranks);
        println!("Number of rows per rank     (-n) {}", numrows_per_rank);
        println!("Avg # non-eros per row      (-Z) {}", nz_per_row);
        println!("Erdos-Renyi edge probabilty (-e) {}", erdos_renyi_prob);

        println!("creating input matrix for toposort");
    }

    let mat = generate_toposort_input(numrows, erdos_renyi_prob, seed, dump_files);

    if !quiet {
        println!(
            "input matrix: {} rows, {} nonzeros",
            mat.numrows(),
            mat.nnz()
        );
    }

    if dump_files {
        mat.dump(20, "mat.out").expect("could not write mat.out");
        mat.write_mm_file("topo_mat.mm")
            .expect("could not write topo_mat.mm");
    }

    let tmat = mat.transpose();

    if dump_files {
        tmat.dump(20, "trans.out")
            .expect("could not write trans.out");
        tmat.write_mm_file("topo_tmat.mm")
            .expect("could not write topo_tmat.mm");
    }

    if !quiet {
        println!("Running toposort on mat (and tmat) ...");
    }

    let mut matret;
    for i in 0..2 {
        if i == 0 {
            if !quiet {
                print!("   using queue toposort: ");
            }
            matret = mat.toposort_queue(&tmat);
        } else {
            if !quiet {
                print!("   using queue2 toposort: ");
            }
            matret = mat.toposort_queue2(&tmat);
        }
        if !mat.check_result(&matret, dump_files) {
            println!("ERROR: check_result failed");
        }
        if !quiet {
            println!(" {} seconds", matret.laptime)
        }
    }
}
