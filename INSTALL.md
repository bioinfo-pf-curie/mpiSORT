# Installation

## Prerequisites

* an implementation of the Message Passing Interface (MPI) standard [mpich](https://www.mpich.org/)
* [zlib](https://zlib.net/)
* [autoconf](https://www.gnu.org/software/autoconf/)
* [automake 1.15](https://www.gnu.org/software/automake/)
* [make](https://www.gnu.org/software/make/)

If you don't have automake 1.15 but a former version (such as 1.13), you can edit in the `configure.ac` file and change the line `AM_INIT_AUTOMAKE([1.15 foreign -Wall])` with `AM_INIT_AUTOMAKE([1.13 foreign -Wall])`.

If automake and autoconf have been installed in custom directories, be sure your their are available in your PATH:

`export PATH=path_to_automake/automake-1.15/bin:path_to_autoconf/autoconf-2.69/bin:$PATH`

If needed, you can set your PATH in your `$HOME/.bashrc` file.

## Build from the git repository



```
git clone https://github.com/bioinfo-pf-curie/mpiSORT.git
cd mpiSORT
aclocal
autoconf
automake --add-missing
./configure CC=/usr/lib64/mpich/bin/mpicc
make
make install
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


## Package the source code
