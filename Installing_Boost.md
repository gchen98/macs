# Introduction #

Some mac users have stated that Boost is not included in their OS distribution.  These steps will outline what to do in case you don't have root access.

# Details #

1) Visit http://www.boost.org/users/download/ and click Download under the Current Release

2) Download any of the links (different compressed formats).  We'll select the one with the .tar.gz extension.

3) Save the file in a desired location. Here, we'll choose /Users/gary

4) Decompress the file.  For example:
tar -xzvf boost\_1\_51\_0.tar.gz

5) change to the directory.  For example:
cd boost\_1\_51\_0

6) Build the boost libraries.  First configure the builder and tell it to build the output into the current directory.  For example:
./bootstrap.sh --prefix=$PWD

6) Run the builder.  Example:
./b2 install

7) In the macs project directory, have the makefile point to the proper locations for boost.  In our example let the LIB variable equal -I /Users/gary/boost\_1\_51\_0/include and let the LINKFLAGS variable equal -L /Users/gary/boost\_1\_51\_0/lib

8) Run make in the macs project.  If unsuccessful, please post the console output in Issues.