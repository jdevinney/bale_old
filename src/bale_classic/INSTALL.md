Start with the [README](README.md). Also see the [DEMO.md](DEMO.md) for instructions for some common platforms.

# Environment Configuration

For this document, let BALEDIR = the directory containing this file.

You can build bale on top of UPC or SHMEM. 

- UPC:   You should set the $UPC environment variable to point to the UPC compiler and 
        the UPCFLAGS variable to whatever options you pass to the compiler (for example:
        UPCFLAGS="-gupc -network=ibv".
- SHMEM: If you are using openshmem, set CC=oshcc (the openshmem compiler).


# Build and Install

To make building and installing easier, we have included a few scripts. The main components of bale are each independent autotools projects and most people would be familiar with how to build and install these. But to make life easier...

1. First, in the bale_classic directory, run the **bootstrap.sh** script.

2. Next, run the **make_bale** script (see below for options). This script visits each subpackage in bale and runs configure, make, and make install. The make_bale script automates the usual process, makes it easy to keep separate build and install directories for multiple platforms, and keeps these directories separate from the source directory. The default build and install directory is $BALEDIR/build_$PLATFORM. If you don't set the $PLATFORM variable, your platform will be set to "unknown". 

   There is one required and several important options for the make_bale script...

	- -u OR -s : You **must** specify UPC or SHMEM.
        - -b : specify an alternate build directory
	- -i : specify an alternate install directory    
	- -c : specify options to configure step
	- -f : force configure to run in all packages again (once you run configure once, the script assumes you don't need to run it again)
        - -j : specify a parallel build
        
The configure process creates architecture dependent files (e.g. Makefile) and creates symbolic links back to the source. This allows us to rename files to satisfy various compilers (e.g. 
some compilers require that UPC programs end in .upc)

##### conveyor autotuning

If you are interested in optimizing the conveyor code, after running make_bale
once (or otherwise building the subpackages) go into build_$PLATFORM/convey
and run 'make tune' with the LAUNCH variable set to a command prefix for
launching a parallel job of a size you care about.  For example:

```bash
make tune LAUNCH="srun -n512"
```

The tuning should take a few minutes.  After it finishes, you can return to
this directory and run make_bale with the -m option (in addition to any
previously used options) to rebuild everything without reconfiguring.  See
also the AUTOTUNING section of convey/INSTALL.