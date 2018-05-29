#ifndef tools_h
#define tools_h
#define kmer_length 31

#include "BooPHF-internal.h"
#include "hash.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>

// Rizk, Chikhi
typedef boomphf::mphf <mer_hasher_uint64 > bphf64_t ; 

// Johnathan Boccara's blog : https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
std::vector <std::string> split (const std::string &);

// Jérôme Audoux
uint64_t str_to_int (const char *, size_t);
std::string int_to_str (uint64_t);


bphf64_t * build_mphf (const std::vector <uint64_t> & , size_t );
void save_mphf (const char * , bphf64_t *);
bphf64_t * load_mphf(const char *);

double ticking();

#endif
