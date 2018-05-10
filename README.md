OnlineRLBWT:
===============
Author: Tomohiro I

### Download

The source codes in 'module' directory are maintained in different repositories.
So, to download all the necessary source codes, do the following:
```sh
git clone https://github.com/itomomoti/OnlineRlbwt.git
cd OnlineRlbwt
git submodule init
git submodule update
```

### Compile

Compilation may require cmake version no less than 3.1, and a compiler supporting some features of C++14.
It has been tested with "Apple LLVM version 7.3.0 (clang-703.0.31)" and "g++6.3.0".

The following commands creates the executable in the build directory (default build type is release).
```sh
mkdir build
cd build
cmake ..
make
```


### Usage

Executables (without option shows help).

```sh
./OnlineRlbwt
./OnlineRindex
./OnlineRindex_Demo
./OnlineLz77ViaRlbwt
./DecompressLz77
```

