# mphf-index

## Dependencies

- sdsl-lite (git clone https://github.com/simongog/sdsl-lite.git)


## Install

` make `

## Usage

- indexing abundance k-mer matrix
  ` ./bin/sdsl-index index <input counts matrix> <output hash filename> <output counts matrix> `

- querying k-mer
  ` ./bin/sdsl-index query <input counts matrix> <input hash> <k-mer> <number of samples> `



