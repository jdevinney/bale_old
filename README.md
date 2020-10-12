# bale

bale is a collection of code meant as a vehicle for discussion for parallel programming models.

The bale effort attempts to:

- demonstrate some challenges of implementing interesting (i.e. irregular) scalable distributed parallel applications.
- demonstrate an approach to achieving high performance for the internode communication in such applications
- explore concepts that make it easier to write, maintain, and get top performance from such applications

The bale repo includes:

- custom communication layers for HPC systems built on top of PGAS programming models (UPC, SHMEM) that offer performance
  and/or programmability benefits relative to hand-coding against lower level communication models.
- example implementations of benchmarks on a variety of programming models and platforms:
  - standard sequential programming langauges (e.g. C, Rust) under src/other_serial
  - standard parallel programming languages/models (e.g. Chapel, Rust) under src/other_parallel
  - custom bale parallel libraries (e.g. exstack, conveyors) under src/bale_classic

Each variant of bale (sequential, parallel, and custom parallel) is a separate project and comes with its own build and run instructions. All source code is found in the src directory.

**bale_classic**: (src/bale_classic)

* [README](src/bale_classic/README.md)
* [Detailed Install instructions](src/bale_classic/INSTALL.md)
* [Build and Run Demos](src/bale_classic/DEMO.md)

**C**: (src/other_serial/C)

* [README](src/other_serial/C/README.md)

**Serial Rust**: (src/other_serial/Rust)

* [README](src/other_serial/Rust/README.md)

**Parallel Rust**: (src/other_parallel/Rust)

* [README](src/other_parallel/Rust/README.md)

**Chapel**: (src/other_parallel/Chapel)

* [README](src/other_parallel/Chapel/README.md)

