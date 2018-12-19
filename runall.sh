#!/bin/sh
#
#
#  Copyright(C) 2018, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
#  This material may be reproduced by or for the US Government
#  pursuant to the copyright license under the clauses at DFARS
#  252.227-7013 and 252.227-7014.
# 
#
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the copyright holder nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
# 
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
#  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
#  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
#  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
#  OF THE POSSIBILITY OF SUCH DAMAGE.
# 

# this script runs all of the bale apps with some reasonable parameters on 1, 4, and 16 nodes
# It exits if any application exits abnormally.
# You specify cores_per_node on the command line as the first argument: for example
# ./runall 32
#
cores_per_node=""
quick_run=0
mask=31
while getopts ":c:qM:" opt; do
    case $opt in
        c ) cores_per_node=$OPTARG;;
        q ) quick_run=1;;
        M ) mask=$OPTARG;;
        \? ) echo 'usage: runall [-c cores_per_node] [-M models_mask] [-q]'
             exit 1
    esac
done

if [ -z $cores_per_node ]; then
    echo "Please specify cores per node"
    exit 1
fi

LAUNCHER=""
if [ -x "$(command -v srun)" ]; then
   LAUNCHER='srun'
fi
if [ -x "$(command -v aprun)" ]; then
   LAUNCHER='aprun'
fi
if [ -z $LAUNCHER ]; then
    echo "Can't find srun or aprun!"
    exit 1
fi

options="-M $mask"
if [ $quick_run == 1 ]; then
    # this makes the tests run a little quicker, if you want to run a longer
    # set of tests, use a second argument (can be anything!)
    options+=" -n 1000"
fi

for app in histo ig topo randperm permute_matrix transpose_matrix triangles write_sparse_matrix
#for app in triangles
do
    # just run the command with -h
    echo;echo; echo XXXXXXXXXXXX $app XXXXXXXXXXXXXXX
    echo;
    $LAUNCHER -n1 $HOME/bale/build_$PLATFORM/apps/$app -h
    echo;
    
    for i in `seq 0 2 4`
    do
        echo;
        nodes=$((2**i))
        cores=$(($nodes*$cores_per_node))

        if [ $app != 'triangles' ]; then
            cmd="$LAUNCHER -n $cores $HOME/bale/build_$PLATFORM/apps/$app -c $cores_per_node $options"
            eval "$cmd"
            if [ $? -ne 0 ]; then
                echo "ERROR! $cmd"
                exit 1
            fi
        else
            for a in `seq 0 1`
            do
                cmd="$LAUNCHER -n $cores $HOME/bale/build_$PLATFORM/apps/$app -a $a -c $cores_per_node $options"
                eval "$cmd"
                if [ $? -ne 0 ]; then
                    echo "ERROR! $cmd"
                    exit 1
                fi
                $LAUNCHER -n $cores $HOME/bale/build_$PLATFORM/apps/$app -a $a -c $cores_per_node $options -K '"1 4 9 16 25"'
                if [ $? -ne 0 ]; then
                    echo "ERROR! in triangles run"
                    exit 1
                fi
            done
        fi
    done
    
done
