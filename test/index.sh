#!/bin/sh

# compressed matrix
../bin/mphf index counts.tsv hash.bin counts_compressed.bin y

# non-compressed matrix
../bin/mphf index counts.tsv hash.bin counts.bin n


