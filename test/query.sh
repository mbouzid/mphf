#!/bin/sh

# non-compressed matrix
echo "non compressed matrix"
## single query
echo "single query"
../bin/mphf query_disk counts.bin hash.bin 4 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC

## multiple queries
echo "multiple queries"
../bin/mphf queries_disk counts.bin hash.bin queries 4


# compressed matrix
echo "compressed matrix"
## single query
echo "single query"
../bin/mphf query counts_compressed.bin hash.bin 4 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC

## multiple queries
echo "multiple queries"
../bin/mphf queries counts_compressed.bin hash.bin queries 4

