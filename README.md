OnlineRLBWT:
===============
Author: Tomohiro I

### Download

The source codes in 'module' directory are maintained in different repositories.
So, to download all the necessary source codes, do the following:
```sh
git clone --recursive https://github.com/itomomoti/OnlineRLBWT.git
```

### Compile

Compilation may require cmake version no less than 3.1, and a compiler supporting some features of C++14.
It has been tested with "Apple LLVM version 7.3.0 (clang-703.0.31)" and "g++6.3.0".

The following commands creates the executable in the build directory (default build type is release).
```sh
mkdir build && cd build
cmake ..
make
```


### Usage

Executing OnlineRLBWT without option shows help.

```sh
./OnlineRLBWT
```

Currently OnlineRLBWT gets input file by -i option, 
computes its RLBWT online while printing statistics to std::out,
then writes the "decompressed" string to the file specified by -o option.

Example.

```sh
./OnlineRLBWT -i in_filename -o out_filename
```

Then, in_filename and out_filename should be the same. Check it by cmp command.

```sh
cmp in_filename out_filename
```
