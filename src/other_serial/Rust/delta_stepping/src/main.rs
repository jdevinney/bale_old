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
use sparsemat::SparseMat;
use delta_stepping::DeltaStepping;

/*
 * Application that finds shortest path lengths from a single source in 
 * a directed graph, using the Meyer/Sanders delta-stepping algorithm. 
 * This serial Rust implementation is based on the serial C in bale,
 * and is intended as a first step toward a Rust/Conveyors version.
 */

/* Find single-source shortest path lengths in a directed graph.

   First we generate the problem by making a random directed graph with edge costs c(v,w).

   Every vertex has a "tentative distance" during the algorithm. The source has tentative
   distance 0, and all other vertices intially have have tentative distance inf. We
   proceed by "relaxing" edges (in a clever order): relaxing edge (v,w) changes the
   tentative distance of w, tent(w), to min(tent(w), tent(v) + c(v,w)). Eventually each
   vertex's tent() distance becomes final, or "settled"; initially only the source is
   settled.

   Unsettled vertices are kept in "buckets" by tent() value; bucket i contains vertices 
   with tent() at least i*\Delta and less than (i+1)*\Delta, where \Delta is a parameter.
   The algorithm has three nested loops. 
   
   The outer (serial) loop is over buckets; an iteration processes vertices in the lowest 
   nonempty bucket until it is empty. 
   
   The middle (serial) loop is over "phases"; a phase consists of removing all the vertices 
   in the bucket and relaxing all the "light" edges out of them (an edge is "light" if it 
   has cost at most \Delta, "heavy" otherwise). The edge relaxations in a phase may cause 
   vertices to enter the active bucket; phases continue until the bucket is empty. At that 
   point all the vertices that were in that bucket are settled.  Following the light-edge 
   phases, one more phase relaxes all the heavy edges from vertices deleted from the active 
   bucket; this cannot cause any vertices to go into that bucket. 
   
   The inner (potentially parallel) loop implements the edge relaxations in a single phase.
   Those relaxations can be done in any order, provided they are done atomically.
*/

fn main() {
    let matches = App::new("DeltaStepping")
        .version("0.1.0")
        .about("Implements a test of DeltaStepping")
        .arg(
            Arg::with_name("numrows")
                .short("n")
                .long("numrows")
                .takes_value(true)
                .help("The number of rows in test matrix"),
        )
        .arg(
            Arg::with_name("source")
                .short("s")
                .long("source")
                .takes_value(true)
                .help("The number of the source vertex"),
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
    let numrows: usize = matches
        .value_of("numrows")
        .unwrap_or("10")
        .parse()
        .expect("numrows: not an integer");
    let source: usize = matches
        .value_of("source")
        .unwrap_or("2")
        .parse()
        .expect("source: not an integer");
    let erdos_renyi_prob: f64 = matches
        .value_of("er_prob")
        .unwrap_or("0.3")
        .parse()
        .expect("er_prob: not a float");

    let seed = 12346;
    let quiet = matches.is_present("quiet");
    let dump_files = matches.is_present("dump_files");

    if !quiet {
        println!("creating input matrix for delta_stepping");
    }

    let mut mat = SparseMat::erdos_renyi_graph(numrows, erdos_renyi_prob, false, seed);
    mat.randomize_values();

    if !quiet {
        println!("input matrix stats:");
        mat.stats();
    }

    if dump_files {
        mat.dump(20, "mat.out").expect("could not write mat.out");
    }

    mat.write_mm_file("sssp_mat.mm")
        .expect("could not write sssp_mat.mm");

    if !quiet {
        println!("Running delta_stepping on mat from source {} ...", source);
    }
    let matret = mat.delta_stepping(source);

    if !mat.check_result(&matret, dump_files) {
        println!("ERROR: check_result failed");
    }
    if !quiet {
        println!(" {} seconds", matret.laptime)
    }
}
