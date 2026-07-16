# chr_aln
Rearrangement-aware whole-genome global alignment

## Usage instructions
### Compilation instructions
#### Cloning the repository
```bash
git clone --recursive https://github.com/hmusta/chr_aln
```

#### Building instructions (Release mode by default)
First, build MUMMer4 in the `mummer` directory. Then compile `chr_aln` as follows:
```bash
mkdir build
cd build
cmake ..
make -j$NTHREADS
```
where `$NTHREADS` can be set as an environment variable, or replaced in this command by however many threads you would like to use for compilation.

Other available modes include `RelDebInfo` (release mode with debugging info), `Debug`, and `Profile`.

### Running chr_aln
#### Compute maximal unique matches (MUMs) with MUMMer4
```bash
mummer -mum -l $k -b -n -c REF.fa QRY.fa > REF_QRY.out
```

#### Run chr_aln
```bash
chr_aln REF.fa QRY.fa REF_QRY.out $NTHREADS $CHECK_INVERSIONS REF_QRY.chain.out > REF_QRY.out
```
