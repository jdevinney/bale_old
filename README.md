# bale

Bale is a collection of code meant as a vehicle for discussion for parallel programming models.

The bale effort attempts to:

* paint a picture of the challenges of implementing a scalable distributed parallel application.
* demonstrate an approach to achieving high performance for the inter-node-communication in such applications
* explore concepts that make it easier to write, maintain, and get top performance from these applications

In addition to [bale_classic](src/bale_classic/README.md), there are now several other instances of bale in other languages. These "cousins" of bale imitate the original in one way or another (usually they implement the "apps"
directory in bale_classic). Each of these cousin directories is a separate project and comes with its own build and run instructions. All source code is found in the src directory.

**bale_classic**: (src/bale_classic)

* [README](src/bale_classic/README.md)
* [Detailed Install instructions](src/bale_classic/INSTALL.md)
* [Install Demos](src/bale_classic/INST_DEMO.md)

**C**: (src/other_serial/C)

* [README](src/other_serial/C/README.md)

**Serial Rust**: (src/other_serial/Rust)

* [README](src/other_serial/Rust/README.md)

**Parallel Rust**: (src/other_parallel/Rust)

* [README](src/other_parallel/Rust/README.md)

**Chapel**: (src/other_parallel/Chapel)

* [README](src/other_parallel/Chapel/README.md)

