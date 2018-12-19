# Copyright (c) 2018, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# This material may be reproduced by or for the U.S. Government 
# pursuant to the copyright license under the clauses at DFARS 
# 252.227-7013 and 252.227-7014.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the copyright holder nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.

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
            cap = buf / (size + 4)
            if (nodes <= 64) {
                printf("-b2 -c%d -n3 -t16 matrix %d %d\n", cap, load, size);
            }
            printf("-b2 -c%d -n3 -t32 matrix %d %d\n", cap, load, size);
            printf("-b2 -c%d -n3 -t16 tensor %d %d\n", cap, load, size);
            if (nodes >= 1024) {
                printf("-b2 -c%d -n3 -t32 tensor %d %d\n", cap, load, size);
            }
        }
        print "EOF"
        if (i % 2 == 0) nodes = (45 * nodes) / 32
        if (i % 2 == 1) nodes = (64 * nodes) / 45
    }
    print "exit 0"
}
