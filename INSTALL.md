
You need to install zlib librairies <br />
You need automake 1.15 for the installation. <br />

You can install automake and autoconf in differents directories and export the path like this: <br />

export PATH=path_to_automake/automake-1.15/bin:path_to_autoconf/autoconf-2.69/bin:$PATH <br />

(The best is to add it in your .bashrc and to source it).

Download from git. In the folder mpiSORT type: <br />

git clone https://github.com/bioinfo-pf-curie/mpiSORT.git <br />
cd mpiSORT <br />
git checkout master <br />
git pull <br />
aclocal  <br />
autoconf  <br />
automake --add-missing  <br />
./configure CC=mpi_bin_path   <br />
 make <br />

or for distribution: <br />
make dist  <br />
tar xzf .tar.gz  <br />
cd mpisort-1.0  <br />
./configure && make install && make<br />

for passing mpi path: <br />
./configure CC=mpi_bin_path  <br />
add --prefix in configure if you need  <br />

If you don't have automake 1.15 change in the configure.ac the line
AM_INIT_AUTOMAKE([1.15 foreign -Wall]) with AM_INIT_AUTOMAKE([1.13 foreign -Wall]) 

