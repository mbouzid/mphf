# mphf

## Dependencies

- sdsl-lite (git clone https://github.com/simongog/sdsl-lite.git)

## Acknowledgment
-  BBHash (github.com/rizkg/BBHash)

## Install

` make && make clean `

## Usage

- indexing abundance k-mer matrix (w/compress)
  ` ./bin/mphf index [counts.tsv] [hash.bin] [counts.bin] y`

- file format for [counts.tsv] 
` tag  sample1 ... sampleM
  <kmer1>  abundance11  ... abundance1M
    ...
  <kmerN> abundanceN1 ... abundanceNM`
  
  e.g. 
  `tag  sample1 sample2 sample3 sample4
  AAAAAAAAAAAAAAAAAAAAAAAAAAATTAT	12  78	11  91
  CAAAAAAAAAAAAAAAAAAAAAAAAAAAATA	59  22	14  78
  CCCTAAAAAAAAAAAAAAAAAAAAAAAAAAT	31  52	11  0`

- querying k-mer 
  ` ./bin/mphf query [counts.bin] [hash.bin] <number of samples> <k-mer>`

## Test

- index
  `cd test && sh index.sh`
- query
  `cd test && sh query.sh`
