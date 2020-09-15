Start with the [README](README.md) or Doxygen documentation for more information about bale.

# Environment Configuration

For this document, let BALEDIR = the directory containing this file.

You can build bale on top of UPC or SHMEM. 

- UPC:   You should set the $UPC environment variable to point to the UPC compiler and 
        the UPCFLAGS variable to whatever options you pass to the compiler (for example:
        UPCFLAGS="-gupc -network=ibv".
- SHMEM: If you are using openshmem, set CC=oshcc (the openshmem compiler).

# Other Dependencies
Bale also depends on:
- argp (part of glibc)
- autoconf (2.69)
- automake (1.13.4)


# Build and Install

To make building and installing easier, we have included a few scripts.

First, in the bale_classic directory, run the **run_autoreconf** script.

Next, run the **install.sh** script (see below for options). This script visits each subpackage in bale and runs configure, make, and make install. You could of course do this yourself by hand; the install script automates the process, makes it easy to keep speparte build and install directories for multiple platforms and keeps these directories separate from the source directory. The default build and install directory is $BALEDIR/build_$PLATFORM. If you don't set the $PLATFORM variable, your platform will be set to "unknown". 

There are several important options for the install.sh script...

- -u OR -s : specify you would like to build libgetput on top of UPC or SHMEM. default: UPC
- -p : specify an alternate install dir
- -f : run autoreconf steps also
- -c : specify options to configure step
- -m : just make, don't configure


The configure process creates architecture dependent files (e.g. Makefile) and creates symbolic links back to the source. This allows us to rename files to satisfy various compilers (e.g. 
some compilers require that UPC programs end in .upc)

If you are interested in optimizing the conveyor code, after running install.sh
once (or otherwise building the subpackages) go into build_$PLATFORM/convey
and run 'make tune' with the LAUNCH variable set to a command prefix for
launching a parallel job of a size you care about.  For example:

    make tune LAUNCH="srun -n512"

The tuning should take a few minutes.  After it finishes, you can return to
this directory and run ./install.sh with the -m option (in addition to any
previously used options) to rebuild everything without reconfiguring.  See
also the AUTOTUNING section of convey/INSTALL.

# Installing on Mac OS X

This is a bit tricky :)
1. install osx compiling environment, bring up the mac app store and install xcode
   After installation you need to accept the license:

    % sudo xcodebuild -license

2. to get autotools, best to use brew (https://brew.sh):

    % brew install autoconf
    
    % brew install automake
    
    % brew install libtool

3. to get clang-upc, which seems to run well, go to:
    https://clangupc.github.io/clang-upc/install.html
   and follow the directions. You will need to do

    % brew install cmake

   to get this to work.  On recent OS X versions (10.15, maybe 10.14) you need
   to pass another argument to cmake: -DDEFAULT_SYSROOT:STRING="$(xcrun --show-sdk-path)"

4. to use shmem, SOS openShmem seems best for now.  Go to:
    https://github.com/Sandia-OpenSHMEM/SOS/wiki/OFI-Build-Instructions
   and follow the directions.  For recent versions OS X (10.15, maybe 10.14) you
   need to set an environment variable before runnig anything will work.
    export SHMEM_OFI_DOMAIN=lo0
   performance is horrid because it only has a socket based implementation

5. when using the runit.sh script with clang, you need to say
    ./runit.sh -l $PWD/clang_upc_run.sh -c 2
   Note that -c should be small but not less than 2, otherwise some things hang

# Documentation 

Bale is documented using Doxygen. To generate the html documentation you must have 
doxygen on your system (I recommend version 1.8.6). To generate the docs, simply change
to the $BALEDIR directory and type:

   doxygen

Then navigate your browser to $BALEDIR/html/index.html.
