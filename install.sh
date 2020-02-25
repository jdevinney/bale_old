#!/bin/bash 
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

#/******************************************************************
# *****************************************************************/ 

# NOTES on using this script:
# You should define an environment variable called PLATFORM that describes your machine.
# In the bale source directory type
#   sh install.sh [options]
# That will create a build directory called build_PLATFORM in the bale source directory.
#
# The installed libraries and binaries will go in build_PLATFORM/lib and build_PLATFORM/bin respectively.
#
# If you are on a machine with Open Shmem you need to add -c "LIBS=-loshmem" as an option.
# If you want to turn on debugging, you can add -c "CFLAGS=-g".


export HERE=$PWD
if [ -z ${PLATFORM+x} ]; then
    PLATFORM="unknown"
fi
export BUILDDIR=$HERE/build_$PLATFORM
export INSTALLDIR=$BUILDDIR
mkdir -p $BUILDDIR

packages=("libgetput" "exstack" "convey" "spmat" "apps")

fromscratch=0
justmake=0
option="--with-upc"
PES=1
config_opts=""
while getopts ":p:usfmj:c:" opt; do
    case $opt in
        p ) INSTALLDIR=$(readlink -f $OPTARG);;
        u ) option="--with-upc";;
        s ) option="--with-shmem";;
        f ) fromscratch=1;;
        m ) justmake=1;;
        c ) config_opts=$OPTARG;;
        j ) PES=$OPTARG;;
        \? ) echo 'usage: install.sh [-p install_dir] [-u] [-s] [-f] [-m] [-j make_threads] [-c configure_opts]'
            exit 1
    esac
done
shift $(($OPTIND - 1))

if [ $justmake -eq 1 ]; then
    fromscratch=0
fi

export PKG_CONFIG_PATH=$INSTALLDIR/lib/pkgconfig/:$PKG_CONFIG_PATH

for i in ${packages[@]};
do  
    echo ""
    echo "*****************************************************"
    echo $i
    echo "*****************************************************"    
    if [ $fromscratch -eq 1 ]; then
        cd $HERE/$i
        cmd="autoreconf -fi"
        echo $cmd
        eval "$cmd"
        if [ $? -ne 0 ]; then
            echo "autoreconf of $i failed!"
            exit 1
        fi
    fi

    mkdir -p $BUILDDIR/$i
    cd $BUILDDIR/$i
    
    if [ $justmake -eq 0 ]; then
        cmd="../../$i/configure --prefix=$INSTALLDIR $option $config_opts"
        echo $cmd
        eval "$cmd"
        if [ $? -ne 0 ]; then
            echo "configure of $i failed!"
            exit 1
        fi
    fi
    #if  [ $fromscratch -eq 1 ]; then
    make clean
    #fi
    make -j $PES
    if [ $? -ne 0 ]; then
        echo "build of $i failed!"
        exit 1
    fi

    make install
    if [ $? -ne 0 ]; then
        echo "install of $i failed!"
        exit 1
    fi
done
