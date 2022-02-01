# Installation

## Prerequisites

* an implementation of the [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) standard such as [mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/) or [Intel® MPI Library](https://software.intel.com/en-us/mpi-library)
* [zlib](https://zlib.net/)
* [autoconf 2.69](https://www.gnu.org/software/autoconf/)
* [automake 1.13](https://www.gnu.org/software/automake/)
* [make](https://www.gnu.org/software/make/)
* [htslib](https://github.com/samtools/htslib)(optionnal, needed to produce BAM output)


The MPI compiler but be available in your PATH or set with the CC environment variable.

Custom options can be used with `configure` such as `--prefix` to set the destination installation path or `CC` for the MPI compiler, for example:

`./configure CC=mpi_bin_path --prefix=${HOME}/local`

This should be only what you need to know about how to use `./configure` but if you are interested, more details are available in the [README-configure](README-configure) and on the command line `./configure --help`.

If `automake` and `autoconf` have been installed in custom directories, be sure their are available in your PATH:

`export PATH=path_to_automake/automake-1.13/bin:path_to_autoconf/autoconf-2.69/bin:${PATH}`

If needed, you can set your PATH according to your configuration directly in your `${HOME}/.bashrc` file.

If you have automake 1.15 you can edit in the `configure.ac` file and change the line `AM_INIT_AUTOMAKE([1.13 foreign -Wall])` with `AM_INIT_AUTOMAKE([1.15 foreign -Wall])`.

## Build from the git repository

```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
# Checkout the branch of the version you want to install, for example:
# git checkout version-1.0
aclocal
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
autoconf
./configure --prefix=${HOME}/local/mpiSORT
# finally 
make
make install
```

## Build with htslib support

To link mpiSORT with htslib you need to build the htslib.a and import some headers.
So to do that build htslib somewhere like this

```
git clone https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
# only curl is needed
autoreconf -i 
./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$INSTALL_PATH
make
make install

```


Then configure mpiSORT with the htslib install dir.
By default compiler flag for lzma, bz2, curl are not present.

If you have installed htslib with curl or lzma, bz2 support
use the flags --enable-lzma, --enable-curl, --enable-bz2 
in during configuration. Ignore the options if not needed.

```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
# Checkout the branch of the version you want to install, for example:
# git checkout version-1.0
aclocal
automake --add-missing
autoreconf -i
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
# by default all libraries (crypto (s3), lzma, bz2, curl) are disable
# from the compilation of htslib we enable curl 
./configure --prefix=${HOME}/local/mpiSORT --enable-curl --with-libhts=${INSTALL_PATH}/htslib (other options: --enable-lzma --enable-bz2 --enable-crypto )
# finally
make
make install
``` 





## Build from a tar.gz archive

Download  the source code archive of the version you want to install from [https://github.com/bioinfo-pf-curie/mpiSORT/releases](https://github.com/bioinfo-pf-curie/mpiSORT/releases).

```
tar xzf mpisort-1.0.tar.gz
cd mpisort-1.0
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiSORT
make
make install
```


## Build from container recipes

[singularity](https://sylabs.io/docs/) recipes are provided in the [containers](../containers) folder. At least [singularity](https://sylabs.io/docs/) version 3.2 is required to build the image. We provided two recipes, the first one using [CentOS](https://www.centos.org/), the second one using [ubuntu](https://ubuntu.com/) as it gives the details on how to install the software such that you can reproduce the installation process if you want to install it locally on your computer.


`singularity` must be available in your PATH.

`sudo singularity build mpiSORT mpiSORT-ubuntu.def`.

This will output the `mpiSORT` executable (that is actually a singularity image) that you can launch as any executable.

## Package the source code into tar.gz archive


```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
# Checkout the branch of the version you want to package, for example:
# git checkout version-1.0
aclocal
autoconf
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiSORT
make dist
```
