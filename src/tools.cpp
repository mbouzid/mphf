#include "tools.h"


std::vector <std::string> split (const std::string & text)
{
  std::istringstream iss(text);
  std::vector<std::string> results((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
  return (results);
}

uint64_t str_to_int (const char* str, size_t l)
{
  uint64_t strint = 0;
  for (size_t i = 0; i < l; i++) {
    uint8_t curr = 0;
    switch (str[i]) {
      case 'A': { curr = 0; break; }
      case 'T': { curr = 3; break; }
      case 'C': { curr = 1; break; }
      case 'G': { curr = 2; break; }
    }
    strint = strint << 2;
    strint = strint | curr;
  }
  return strint;
}

std::string int_to_str(uint64_t kmer)
{
  uint8_t base;
  std::string str;
  for (int i=kmer_length; i>0; i--) {
    base = (kmer >> (i*2-2)) & 3ULL;
    char chr;
    switch (base) {
      case 0: { chr = 'A'; break; }
      case 3: { chr = 'T'; break; }
      case 1: { chr = 'C'; break; }
      case 2: { chr = 'G'; break; }
    }		
    str.push_back(chr);
  }
  return str;
}

bphf64_t * build_mphf (const std::vector <uint64_t> & data, size_t nthreads)
{	
  double gammaFactor = 2.0; 
  bphf64_t * bphf = new bphf64_t (data.size(),data,nthreads,gammaFactor);

  return bphf;
}

void save_mphf (const char * output, bphf64_t * bphf)
{
  std::ofstream foutput (output);
  bphf->save(foutput);
  foutput.close(); 
}

bphf64_t * load_mphf(const char * input)
{
  std::ifstream finput (input);
  bphf64_t * bphf = new bphf64_t();
  bphf->load(finput);
  finput.close();

  return bphf;
}


double ticking()
{
  struct timeval timet;
  gettimeofday(&timet, NULL); 
  return timet.tv_sec +(timet.tv_usec/1000000.0);
}

