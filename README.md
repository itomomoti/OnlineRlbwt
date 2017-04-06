OnlineRLBWT:
===============
Author: Tomohiro I

### Compile

Compilation may require cmake version no less than version 3.1, and a compiler supporting some features of C++14.

```
mkdir build
cd build
```

Then, launch cmake as (default build type is release):

```
cmake ..
```

Finally, build the executable:

```
make
```

The above command creates the executables in the build directory. 

### Usage

Executing OnlineRLBWT without option shows help.

```
OnlineRLBWT
```

Currently OnlineRLBWT gets input file by -i option, computes its RLBWT online,
then output "decompressed" string to the file specified by -o option.

Example.

```
./OnlineRLBWT -i ../testfile -o testoutput
```

Then, ../testfile and ./testoutput should be the same. Check it by diff command.

```
diff ../testfile ./testoutput
```
