# Bale Applications written in Rust

This project explores the implementation of the Bale applications,
[available here](https://github.com/jdevinney/bale), using the Rust
programming language.

It builds on the rust conveyor library [available
here](https://github.com/wwc559/convey_private).

**Note: This repository is currently private, please do not redistribute.**
Open an 'issue' to get access for others.

To build this just git clone both this repo and the conveyor rust
library in the same parent directory and say `cargo build --release`
and then you can run it with:

```
oshrun -n 4 target/release/histo_convey
```



