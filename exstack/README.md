# exstack

This library contains both exstack and exstack2.

exstack was originally written in 2011 and was our first attempt at an
aggregation library. The main functions in its API (push and pop)
remain in its descendants (exstack2 and
[conveyors](../convey/README.md)). In a typical exstack loop there are
three phases. First, each PE pushes items onto its local
out-buffers. Once a PE sees that one of its out-buffers has become
full, that PE goes to the "exchange" phase where it waits for all
other PEs to join it. Once all PEs are in the exchange phase, all
out-buffers are sent to their destination where they land in
in-buffers. PEs then enter into the pop phase where they pop items off
of their in-buffers and do whatever computation is required. See
[histogram](../apps/histo_src/README.md) for a simple example. While
exstack is naive compared to conveyors, it still acheives very good
performance on most of the bale apps.

exstack2 was our attempt to improve both the performance and
ease-of-use of exstack. In fact, bale grew out of our attempt to write
exstack2.
