# Microbenchmarks for Dinocchio

This is a repository directly cloned from the implementation of [Rinocchio](https://github.com/zkFHE/ringSNARK) with minor changes for our microbenchmarks of the underlying encoding scheme used.

To get the unit time of operations for our encoding schemes:
```bash
sudo apt-get update && sudo apt-get install build-essential # if needed
sudo apt-get install autoconf # if no autoconf
sudo apt-get install cmake # if no cmake
sudo apt install libboost-all-dev

cd dringSNARK
git submodule init && git submodule update --recursive
mkdir build && cd build && cmake ..
make
```
And then run, which should print out numbers similar to Figure 8 in our submission:
```bash
./bin/bench
```

To generate the numbers in Figure 5 and Figure 6, run:
```bash
cd dringSNARK
python3 runtime.py
```

To generate the numbers in Figure 7, run:
```bash
cd dringSNARK
python3 runtime_N.py
```



