
You need to install zlib librairies <br />
You need automake 1.15 for the installation. <br />
You can install automake and autoconf in differents directories and export the path like this: <br />
export PATH=path_to_automake/automake-1.15/bin:path_to_autoconf/autoconf-2.69/bin:$PATH <br />
(The best is to add it in your .bashrc and to source it).

Download from git. In the folder mpiSORT type: <br />
./configure && make install && make <br />

or for distribution: <br />
make dist  <br />
tar xzf .tar.gz  <br />
cd mpisort-1.0  <br />
./configure && make install && make<br />

for passing mpi path: <br />
./configure CC=mpi_bin_path  <br />
add --prefix in configure if you need  <br />

example of command line in a terminal:  <br />

mpirun -n 5 psort INPUT_FILE OUTPUT_DIR -q 0  <br />

-q is for quality filtering of the reads.

