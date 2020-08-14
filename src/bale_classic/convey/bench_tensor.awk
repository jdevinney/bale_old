# Copyright (c) 2019, Institute for Defense Analyses
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
    if (ppn <= 0 || ppn % 2 == 1) {
        print "set procs/node to an even number with -v ppn=..."
        exit
    }
    if (max <= 0) {
        print "set maximum nodes with -v max=..."
        exit
    }

    bytes[0] = 8; bytes[1] = 16; bytes[2] = 32; bytes[3] = 128
    print "#!/bin/sh -e"
    print "export MPP_INIT_QUIET=1"
    nodes = 2; power = 2
    for (i = 0; nodes <= max; i++) {
        print "$NLAUNCH", ppn * nodes, "./alltoall? -- <<EOF"
        for (j = 0; j < 4; j++) {
            size = bytes[j]
# each process wants to send 256 MiB
            load = 268435456 / (nodes * ppn * size)
            if (nodes <= 64) {
                printf("-b2 -n3 vector %d %d\n", load, size);
            }
            if (nodes <= 128) {
                printf("-b2 -n3 -t%d matrix %d %d\n", ppn / 2, load, size);
            }
            printf("-b2 -n3 -t%d matrix %d %d\n", ppn, load, size);
            if (nodes > ppn) {
                printf("-b2 -n3 -t%d tensor %d %d\n", ppn / 2, load, size);
            }
            if (nodes >= 1024) {
                printf("-b2 -n3 -t%d tensor %d %d\n", ppn, load, size);
            }
        }
        print "EOF"
        if (i % 2 == 1) {
            power *= 2
            nodes = power
        } else {
            if (power >= 32) nodes = 45 * power / 32
            else if (power >= 8) nodes = 11 * power / 8
            else nodes = 3 * power / 2
        }
    }
    print "exit 0"
}
