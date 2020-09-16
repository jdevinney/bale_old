# Building bale on some common platforms

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

## ... on a Cray XC30 on UPC

Use the PrgEnvCray module to access the Cray UPC compiler

```bash
cd $BALEDIR
export PLATFORM=xc30_upc
export UPC="cc -hupc"
./bootstrap.sh
./install.sh
```

## ... on a Cray XC30 with Cray-shmem



