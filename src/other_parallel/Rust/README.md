# Bale Applications written in Rust

This is the implementation of the Bale applications, using the Rust
programming language.

It builds on the rust conveyor library [available
here](https://github.com/wwc559/convey_private). 

The build is setup to run correctly if you place convey_private
in the same directory where you placed bale.  Modify Cargo.toml,
spmat/Cargo.toml, triangle/Cargo.toml, and toposort/Cargo.toml
if you have placed it in a different location.

To build this and say `cargo build --release --workspace` and then you can 
run any of the bale apps with:

```
oshrun -n 4 target/release/histo_convey
oshrun -n 4 target/release/ig_convey
oshrun -n 4 target/release/permute_convey
oshrun -n 4 target/release/randperm_convey
oshrun -n 4 target/release/toposort
oshrun -n 4 target/release/triangle
```



