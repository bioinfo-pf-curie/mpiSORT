# Installation

## Prerequisites

* an implementation of the Message Passing Interface (MPI) standard such as [mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/) or [Intel® MPI Library](https://software.intel.com/en-us/mpi-library)
* [zlib](https://zlib.net/)
* [autoconf](https://www.gnu.org/software/autoconf/)
* [automake 1.15](https://www.gnu.org/software/automake/)
* [make](https://www.gnu.org/software/make/)

The MPI compiler but be available in your PATH or set with the CC environment variable.

If you don't have automake 1.15 but a former version (such as 1.13), you can edit in the `configure.ac` file and change the line `AM_INIT_AUTOMAKE([1.15 foreign -Wall])` with `AM_INIT_AUTOMAKE([1.13 foreign -Wall])`.

If automake` and `autoconf` have been installed in custom directories, be sure your their are available in your PATH:

`export PATH=path_to_automake/automake-1.15/bin:path_to_autoconf/autoconf-2.69/bin:${PATH}`

If needed, you can set your PATH according to your configuration directly in your `${HOME}/.bashrc`.

Custom options can be used with `configure` such as `--prefix` to set the destination installation path or `CC` for the MPI compiler, for example:
./configure CC=mpi_bin_path --prefix

## Build from the git repository

```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
aclocal
autoconf
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiSORT
make
make install
```

## Build a tar.gz archive


```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
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


## Build from a tar.gz archive

Download  the source code archive from [https://github.com/bioinfo-pf-curie/mpiSORT/releases](https://github.com/bioinfo-pf-curie/mpiSORT/releases).

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

[singularity](https://sylabs.io/docs/) recipes are provided in the `containers` folder. At least [singularity](https://sylabs.io/docs/) version 3.2 is required to build the image. We provided two recipes, the first one using [CentOS](https://www.centos.org/), the second one using [ubuntu](https://ubuntu.com/) as it gives the details on how to install the software that you can reproduc if you want to install it locally on your computer.


`singularity` must be available in your PATH.

`sudo singularity build mpiSORT mpiSORT-ubuntu.def`.

This will output `mpiSORT` executable (that is acty a singularity image) that you can launch as any executable.

