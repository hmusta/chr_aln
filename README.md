# chr_aln
Rearrangement-aware whole-genome global alignment

## Usage instructions
### Compilation instructions
#### Cloning the repository
```bash
git clone --recursive https://github.com/hmusta/chr_aln
```

#### Building instructions (Release mode by default)
```bash
mkdir build
cd build
cmake ..
make -j$NTHREADS
```
where `$NTHREADS` can be set as an environment variable, or replaced in this command by however many threads you would like to use for compilation.

Other available modes include `RelDebInfo` (release mode with debugging info), `Debug`, and `Profile`.
