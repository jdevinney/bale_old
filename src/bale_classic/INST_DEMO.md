# Building and Running bale on some common platforms

Assuming you have cloned bale in a directory called $BALEDIR.

## ... on a Linux SMP with OpenMPI/OSHMEM

You should have a reasonably modern C complier (gcc 8.3.0 seems to be working well for us) and 'oshcc' in your path. You also need OpenMPI 4.x.x or higher.

```bash
cd $BALEDIR
export PLATFORM=linux_oshmem
export CC=oshcc
./bootstrap.sh
./install.sh -s
```

This builds everything in $BALEDIR/build_linux_oshmem. Binaries appear in $BALEDIR/build_linux_oshmem/bin.

**Note**: Due to a bug in OpenMPI/OpenSHMEM, exstack2 and conveyor apps sometimes hang. We recommend you avoid running them under oshmem by running application with the implementation mask set to 3 (-M 3).

## ... on a Linux SMP with Sandia OpenSHMEM

You should have a reasonably modern C complier (gcc 8.3.0 seems to be working well for us) and 'oshcc' in your path. You also need SOS 1.4.5 or higher, preferrably built on XPMEM. See note below.

Otherwise, exactly the same instructions as above substituting "linux_sos" for PLATFORM.

**Note**: Due to a bug in SOS built using CMA transport layer, exstack2 and conveyor apps sometimes hang. We recommend you avoid running them under SOS by running application with the implementation mask set to 3 (-M 3).

## ... on a Linux SMP with GNU UPC (GUPC) or Clang UPC (CUPC)

You should have a reasonably modern C complier (gcc 8.3.0 seems to be working well for us) and gupc/upc in your path for GUPC or CUPC. You should also have GUPC 5.2.0.1 or Clang-UPC 3.9.1.1.

```bash
cd $BALEDIR
export PLATFORM=linux_gupc
export UPC=gupc (for GUPC)
export UPC=upc (for CUPC)
unset CC
./bootstrap.sh
./install.sh
```

This builds everything in $BALEDIR/build_linux_oshmem. Binaries appear in $BALEDIR/build_linux_oshmem/bin.

## ... on a Cray XC30

- If you want to use Cray UPC:
    - Use the PrgEnvCray module to access the Cray UPC compiler.
    - `export UPC="cc -hupc"`
- If you want to use Cray SHMEM: Use PrgEnvCray and load the cray-shmem module.
- If you want to use Cray OpenSHMEMX: Use PrgEnvgnu and load cray-openshmemx. 

```bash
cd $BALEDIR
export PLATFORM=xc30
./bootstrap.sh
./install.sh
```



# Run a test
Try running a simple test (remember to use -M 3 with OpenMPI/oshmem or SOS)
### with oshrun
```bash
oshrun -n 4 $BALEDIR/src/bale_classic/build_$PLATFORM/bin/histo
```
### with slurm
```bash
srun -n 4 $BALEDIR/src/bale_classic/build_$PLATFORM/bin/histo
```
### with gupc
```bash
$BALEDIR/src/bale_classic/build_$PLATFORM/bin/histo -n 4
```

You should see something like this...

```bash

***************************************************************
Bale Version 3.00 (UPC 201311): 2110-08-17.12:50
Running command on 4 PEs: ../build_ucs3_gupc/bin/histo
***************************************************************

num_updates_per_pe: 100000
table_size_per_pe: 1000
Standard options:
----------------------------------------------------
buf_cnt (buffer size)    (-b): 1024
seed                     (-s): 122222
cores_per_node           (-c): 0
Models Mask              (-M): 15

       AGI:    0.009 seconds     0.000 GB/s injection bandwidth
   Exstack:    0.003 seconds     0.000 GB/s injection bandwidth
  Exstack2:    0.004 seconds     0.000 GB/s injection bandwidth
  Conveyor:    0.012 seconds     0.000 GB/s injection bandwidth
```