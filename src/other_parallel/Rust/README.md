# Bale Applications written in Rust

This is the implementation of the Bale applications, using the Rust
programming language.

It builds on the rust conveyor library [available
here](https://github.com/wwc559/convey_private). The
build is setup to run correctly if you place convey_private
in the same directory where you placed bale.  Modify Cargo.toml
if you have placed it in a different location.

To build this and say `cargo build --release` and then you can 
run it with:

```
oshrun -n 4 target/release/histo_convey
```



