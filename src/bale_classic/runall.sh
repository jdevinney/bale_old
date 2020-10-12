#!/bin/sh
#
#
#  Copyright(C) 2019-2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
#
#  All rights reserved.
#  
#  This file is a part of Bale.  For license information see the
#  LICENSE file in the top level directory of the distribution.
#  

# this script runs all of the bale apps with some reasonable parameters on 1, 4, and 16 nodes
# It exits if any application exits abnormally.
# You specify cores_per_node on the command line using the -c option: for example
# ./runall.sh -c 32
#
if [ -z ${PLATFORM+x} ]; then
    PLATFORM="unknown"
fi
cores_per_node=""
quick_run=0
mask=31
path=$HOME/bale/build_$PLATFORM/bin
LAUNCHER=""
UPC_OPTS=""
options=""
while getopts ":c:qM:p:e:o:l:" opt; do
    case $opt in
        c ) cores_per_node=$OPTARG;;
        q ) quick_run=1;;
        M ) mask=$OPTARG;;
        p ) path=$OPTARG;;
	o ) options=$OPTARG;;
	l ) LAUNCHER=$OPTARG;;
        \? ) echo 'usage: runall -c cores_per_node [-M models_mask] [-p PATH_TO_BINARIES] [-q] [-l LAUNCHER] [-o options]'
             exit 1
    esac
done

if [ -z $cores_per_node ]; then
    echo "Please specify cores per node"
    exit 1
fi

if [ -z $LAUNCHER ]; then
     if [ -x "$(command -v srun)" ]; then
	 LAUNCHER='srun'
     fi
     if [ -x "$(command -v aprun)" ]; then
	 LAUNCHER='aprun'
     fi
     if [ -x "$(command -v oshrun)" ]; then
	 LAUNCHER='oshrun'
     fi
     if [ -x "$(command -v upcrun)" ]; then
         LAUNCHER='upcrun'
     fi
fi
if [ -z $LAUNCHER ]; then
    echo "Can't find oshrun, srun, upcrun, or aprun!"
    echo "Assuming this is GNU or CLANG UPC..."
    export UPC_OPTS=" -n 1"
    unset PREAMBLE
    echo "$UPC_OPTS"
fi

: ${BALEDIR:=$PWD}
echo $mask
options=$options" -M $mask"
if [ ${quick_run} -eq 1 ]; then
    # this makes the tests run a little quicker, if you want to run a longer
    # set of tests, use a second argument (can be anything!)
    options+=" -n 1000"
fi

for app in histo ig topo randperm permute_matrix transpose_matrix triangles write_sparse_matrix
#for app in triangles
do
    if [ ! -z "$UPC_OPTS" ]; then
        UPC_OPTS=' -n 1'        
    else
        PREAMBLE="$LAUNCHER -n 1"
    fi
    # just run the command with -h
    echo;echo; echo XXXXXXXXXXXX $app XXXXXXXXXXXXXXX
    echo;
    cmd=$PREAMBLE $BALEDIR/build_$PLATFORM/apps/$app $UPC_OPTS -h
    echo $cmd
    eval "$cmd"
    echo;
    
    for i in `seq 0 2`
    do
        echo;
        nodes=$((2**i))
        cores=$(($nodes*$cores_per_node))
        if [ ! -z "$UPC_OPTS" ]; then
            UPC_OPTS=' -n $cores'            
        else
            PREAMBLE="$LAUNCHER -n $cores"
        fi
        if [ $app != 'triangles' ]; then
            cmd="$PREAMBLE $BALEDIR/build_$PLATFORM/apps/$app ${UPC_OPTS} -c $cores_per_node $options"
            eval "$cmd"
            if [ $? -ne 0 ]; then
                echo "ERROR! $cmd"
                exit 1
            fi
        else
            for a in `seq 0 1`
            do
                cmd="$PREAMBLE $BALEDIR/build_$PLATFORM/apps/$app ${UPC_OPTS} -a $a -c $cores_per_node $options"
                echo "$cmd"
                eval "$cmd"
                if [ $? -ne 0 ]; then
                    echo "ERROR! $cmd"
                    exit 1
                fi
                triopt="-K '1 4 9 16'"
                cmd2="$PREAMBLE $BALEDIR/build_$PLATFORM/apps/$app ${UPC_OPTS} -a $a -c $cores_per_node $options ${triopt}"
                echo "$cmd2"
                eval "$cmd2"
                if [ $? -ne 0 ]; then
                    echo "ERROR! in triangles run"
                    exit 1
                fi
            done
        fi
    done
    
done
