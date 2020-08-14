#!/bin/sh

if [ "$#" -lt 3 ]; then
    echo "usage: $0 -n THREADS ..."
    exit 1
fi

if [ "$1" != "-n" ]; then
    echo "usage: $0 -n THREADS ..."
    exit 1
fi

shift
threads="$1"
shift
program="$1"
shift
$program -fupc-threads-$threads -- "$@"
