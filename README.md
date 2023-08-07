# SACfe: Secure Access Control in Functional Encryption with Unbounded Data

This repository conatins the code for the schemes UZP-IPFE and UAB-IPFE.

Steps to reproduce the experiments.

1. Install dependencies of CiFEr
2. Build our modified CiFEr library
3. Compile and run an example

### Install dependencies of CiFEr

The requirements have to be installed manually (via package manager or building 
the source code).
- [CMake](https://cmake.org/download/) (version 3.11+)
- [GMP](https://gmplib.org/)
- [libsodium](https://download.libsodium.org/doc/)
- [AMCL](https://github.com/miracl/amcl)
- [Protobuf](https://github.com/protocolbuffers/protobuf)

To be able to build CiFEr as described below, AMCL must be compiled with BN254
curve. This can be done manually, but for convenience, we provide a Bash script
that runs a modified AMCL setup (a Python script) and installs a minimal version
 of AMCL in the standard directory `/usr/local/lib` and header files in
`/usr/local/include`. These default values can be changed in
`external/amcl/setup_amcl.sh`. To use the script, run:

```
cd external/amcl
sudo ./setup_amcl.sh
cd ../..
```

### Build our modified CiFEr library

To build and install, first clone the repo, then run the following commands in the 
source code directory:

```
mkdir build
cd build
cmake ..
make
sudo make install
```

This builds the shared library (`libcifer.so`) and installs it.
By default, it places the shared library in `/usr/local/lib` and the header 
files in `/usr/local/include`

### Compile an example

To compile and run an example go to the example folder and run make.
There is an example for each proposed scheme e.g. uzpipfe / uabipfe :

```
cd example
make uzpipfe
./uzpipfe.out
```