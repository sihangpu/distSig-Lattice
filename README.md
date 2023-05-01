# Description
This is a proof-of-concept implementation of the lattice-based distributed signature (i.e., n-out-of-n threshold signature) in C++.

## Instruction
Make sure you have installed `NTL`, `cmake`, and `ninja` , otherwise run the following command to install.
```bash
brew install cmake ninja ntl
```


Build the project:
```bash
 cd <ROOT_PATH_OF_THE_REPO>
 mkdir build && cd build
 cmake -GNinja ..
 ninja
```

By default, the exectuable files in `<ROOT_PATH_OF_THE_REPO>/build/src/` 

