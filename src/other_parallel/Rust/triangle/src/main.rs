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
use convey_hpc::collect::ValueCollect;
use convey_hpc::Convey;
use triangle::Triangle;

/*
  \page triangles_page Triangles

This uses matrix algebra approach to counting triangles in a graph.

The adjacency matrix, <b>A</b>, for the graph is a {0,1} matrix
where the rows and cols correspond to the vertices
and \f$a_{ij} \f$ = <b>A[i][j]</b> is 1 exactly when there is a edge between
vertices <i>v_i</i> and <i>v_j</i>.

The triangle with vertices <i>{v_i, v_j, v_k}</i> has associated
edges <i>{v_i, v_j}</i>, <i>{v_j, v_k}</i> and <i>{v_k, v_i}</i>
which correspond to non-zero entries
\f$a_{ij}\f$,
\f$a_{jk}\f$, and
\f$a_{ki}\f$
in the adjacency matrix.
Hence the sum
\f$ \sum_{i,j,k} a_{ij}a_{jk}a_{ki} \f$ counts the triangles in the graph.
However, it counts each triangle 6 times according to the 6 symmetries of a triangle
and the 6 symmetric ways to choose the three nonzeros in <b>A</b>.
To count each triangle once, we compute the sum
\f[ \sum_{i=1}^{n}\sum_{j=1}^{i-1}\sum_{k=1}^{j-1} a_{ij}a_{jk}a_{ik} =
    \sum_{i=1}^{n}\sum_{j=1}^{i-1} a_{ij} \sum_{k=1}^{j-1} a_{jk}a_{ik}. \f]

This picks out a unique labelling from the 6 possible and it means that
all the information we need about edges is contained in the lower triangular
part of symmetric adjacency matrix.  We call this matrix <b>L</b>.

The mathematical expression:
for each nonzero \f$ a_{ij} \f$ compute the dot product
of row \f$ i\f$ and row \f$ j \f$ becomes
\verbatim
  For each non-zero L[i][j]
     compute the size of the intersection of the nonzeros in row i and row j
\endverbatim

Usage:
- -a = 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
- -e = p: specify the Erdos-Renyi probability p
- -K = str: Generate a Kronecker product graph with specified parameters. See below
- -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate
- -N = n: Specify the number of rows_per_thread in the matrix (if using the Erdos-Renyi generator).
- -r "file" : Specify a filename containing a matrix in MatrixMarket format to read as input.
- -b = count: Specify the number of packages in an exstack(2) stack

Explanation of -K option. Using a special string format you must specify a mode,
and a sequence of numbers. For example:
"0 3 4 5 9"
The first integer is always the mode and the valid modes are 0, 1, or 2
Mode 0 graphs have no triangles, mode 1 graphs have lots of triangles and mode 2 graphs
have few triangles.
After the first number, the next numbers are the parameters for the two kronecker product graphs. We split
the sequence of numbers in half to get two sequences.
In our example above we would produce the product of K(3,4) and K(5,9).

See "Design, Generation, and Validation of Extreme Scale Power-Law Graphs" by Kepner et. al.
 */

fn main() {
    let matches = App::new("Triangle")
        .version("0.1.0")
        .about("Implements a test of Triangle Counting")
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
            Arg::with_name("gen_kron")
                .short("K")
                .long("gen_kron")
                .takes_value(true)
                .help("Kroniker product graph size"),
        )
        .arg(
            Arg::with_name("alg")
                .short("a")
                .long("alg")
                .takes_value(true)
                .help("algorithm variant, 0 is L & L * U, 1 is L & U * L"),
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
    let alg: u16 = matches
        .value_of("alg")
        .unwrap_or("0")
        .parse()
        .expect("alg: not an integer");
    let gen_kron = matches.is_present("gen_kron");
    let gen_kron_str = matches.value_of("gen_kron").unwrap_or("");

    let convey = Convey::new().expect("Conveyor system initializtion failed");
    let seed = 12346;
    let quiet = matches.is_present("quiet") || convey.my_rank > 0;
    let _dump_files = matches.is_present("dump_files");

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
        println!("Running triangle on {} ranks", convey.num_ranks);
        println!("Number of rows per rank     (-n) {}", numrows_per_rank);
        println!("Algorithm                   (-a) {}", alg);
        if gen_kron {
            println!("Avg # non-zeros per row     (-Z) {}", nz_per_row);
            println!("Erdos-Renyi edge probabilty (-e) {}", erdos_renyi_prob);
        } else {
            println!("Generating Kronicker Prod   (-K) {}", gen_kron_str);
        }
        println!("creating input matrix for triangle");
    }

    let (mat, correct_answer) = if gen_kron {
        let mut iter = gen_kron_str.split(" ");
        let mode: u16 = iter.next().unwrap_or("0").parse().expect("bad mode arg");
        let mut args: Vec<u16> = Vec::new();
        for item in iter {
            args.push(item.parse().expect("bad kron value"));
        }
        let ans = match mode {
            0 => 0.0,
            1 => {
                let mut ans = 1.0;
                for item in &args {
                    ans *= 3 as f64 * *item as f64 + 1.0;
                }
                ans *= 1.0 / 6.0;
                let mut x = 1.0;
                for item in &args {
                    x *= *item as f64 + 1.0;
                }
                ans - 0.5 * x + 1.0 / 3.0
            }
            2 => {
                (1.0 / 6.0) * 4.0f64.powf(args.len() as f64) - 2.0f64.powf((args.len() - 1) as f64)
                    + 1.0 / 3.0
            }
            _ => panic!("unsupported kronecker mode"),
        };

        let half = args.len() / 2;
        (
            triangle::generate_kronecker_graph(&args[..half], &args[half..], mode),
            ans.round(),
        )
    } else {
        (
            spmat::SparseMat::gen_erdos_renyi_graph(numrows, erdos_renyi_prob, false, 1, seed),
            0.0,
        )
    };
    if !quiet {
        println!(
            "input matrix: {} rows, {} nonzeros, expected count: {}",
            mat.numrows(),
            mat.nnz(),
            correct_answer
        );
    }

    let tmat = mat.transpose();
    let upper = if alg == 1 { Some(&tmat) } else { None };

    if !quiet {
        println!("Running triangle on mat (and tmat) ...");
    }

    let mut triret;
    for i in 0..1 {
        if i == 0 {
            if !quiet {
                print!("   using triangle_push: ");
            }
            triret = mat.triangle_push(upper);
        } else {
            if !quiet {
                print!("   using triangle_pull: ");
            }
            triret = mat.triangle_pull(upper);
        }
        let tot_shref = mat.reduce_sum(triret.sh_refs);
        let tot_tri_cnt = mat.reduce_sum(triret.tri_cnt);
        if !quiet {
            println!(
                " {:.4} seconds, {} Shrefs, {} Triangles",
                triret.laptime, tot_shref, tot_tri_cnt
            )
        }
    }
}
