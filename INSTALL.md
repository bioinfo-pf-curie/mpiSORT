# Installation

## Prerequisites

* an implementation of the Message Passing Interface (MPI) standard ([mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/), [Intel® MPI Library](https://software.intel.com/en-us/mpi-library)
* [zlib](https://zlib.net/)
* [autoconf](https://www.gnu.org/software/autoconf/)
* [automake 1.15](https://www.gnu.org/software/automake/)
* [make](https://www.gnu.org/software/make/)

The MPI compiler `mpicc` but be available in your PATH.

If you don't have automake 1.15 but a former version (such as 1.13), you can edit in the `configure.ac` file and change the line `AM_INIT_AUTOMAKE([1.15 foreign -Wall])` with `AM_INIT_AUTOMAKE([1.13 foreign -Wall])`.

If `mpicc`, `automake` and `autoconf` have been installed in custom directories, be sure your their are available in your PATH:

`export PATH=path_to_automake/automake-1.15/bin:path_to_autoconf/autoconf-2.69/bin:path_to_mpicc/bin:$PATH`

If needed, you can set your PATH in your `$HOME/.bashrc` file to set your PATH according to your configuration.


## Build from the git repository

```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
aclocal
autoconf
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc` at the configure stage
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure
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
# If not yet in your PATH, you can provide the PATH to `mpicc` at the configure stage
# ./configure CC=/usr/lib64/mpich/bin/mpicc
./configure
make dist
```


## Build from a tar.gz archive

or for distribution: 
make dist  
tar xzf .tar.gz  
cd mpisort-1.0  
./configure && make install && make

for passing mpi path: 
./configure CC=mpi_bin_path  
add --prefix in configure if you need  

## Build from container recipes



## Package the source code
