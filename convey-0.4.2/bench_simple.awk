# Copyright (c) 2018, Institute for Defense Analyses,
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
#
# This material may be reproduced by or for the U.S. Government 
# pursuant to the copyright license under the clauses at DFARS 
# 252.227-7013 and 252.227-7014.

BEGIN {
    if (buf == 0) {
        print "set approximate buffer bytes with -vbuf=..."
        exit
    }

    bytes[0] = 8; bytes[1] = 16; bytes[2] = 32; bytes[3] = 128
    print "#!/bin/sh -e"
    print "export MPP_INIT_QUIET=1"
    nodes = 32
    for (i = 0; i < 13; i++) {
        print "$NLAUNCH", 32 * nodes, "./alltoall? -- <<EOF"
        for (j = 0; j < 4; j++) {
            size = bytes[j]
# each process wants to send 256 MiB, so load = 2^28 / (32 * nodes * size)
            load = 8388608 / (nodes * size)
            cap = buf / size
# make sure buffers occupy at most 1GB
            while (cap * size * nodes * 32 * 2 > 1000 * 1000 * 1000) {
                cap = cap / 2
            }
            for (k = 0; k < 2; k++) {
                option = (k == 0) ? "" : "-x "
                printf("-c%d -n3 %ssimple %d %d\n", cap, option, load, size)
            }
        }
        print "EOF"
        if (i % 2 == 0) nodes = (45 * nodes) / 32
        if (i % 2 == 1) nodes = (64 * nodes) / 45
    }
    print "exit 0"
}
