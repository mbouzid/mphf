#include <cstdlib>
#include <map>
#include <algorithm>
#include <omp.h>

#include "tools.h"

#include <sdsl/int_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/coder.hpp>

const int nthreads = 20;
const long unsigned int trace_sz = 1000000;


bphf64_t * building_saving_mphf(char * output_hash, const std::vector <uint64_t> & kmers)
{
  double btime, etime ;
 
  std::cout << "building and saving hash" << std::endl;
  btime = ticking();
  
  bphf64_t * bphf = build_mphf(kmers,nthreads);
  save_mphf (output_hash,bphf);
  
  etime = ticking();
  std::cout << "elapsed: " << etime-btime << " s" << std::endl;
  return bphf;
}

int sorting_by_hash (std::vector <uint64_t> & kmers, bphf64_t * bphf)
{
  double btime, etime ;
 
  auto by_hashkey = [&bphf] (uint64_t a, uint64_t b) { return (bphf->lookup(a) < bphf->lookup(b)); };
  std::cout << "sorting kmers by hash key" << std::endl; 

  btime = ticking();
  std::sort (kmers.begin(), kmers.end(), by_hashkey);
  etime = ticking();

  std::cout << "elapsed: " << etime-btime << " s" << std::endl;
  
  return true;
}




int index (int argc, char * argv [])
{
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " <input_counts> <output_hash> <output_counts> " << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv [1];	    // input matrix counts
  char * output_hash = argv [2];	    // output mphf
  char * output_counts = argv[3];	    // output matrix counts

  std::vector <uint64_t> kmers;		    // vector of kmers for building the mphf
  size_t nsamples;
  double btime, etime ;
 
  std::ifstream finput_counts (input_counts);
  std::string line;

  std::getline(finput_counts,line);	    // get header
  nsamples = split(line).size() -1;    // number of samples = number of elements in header -1 (tag)
  std::cout << "nsamples=" << nsamples << std::endl;
  std::cout << "loading kmers & offsets, saving counts" << std::endl;
  btime = ticking();
  
  std::string tmp("tmp_"+std::string(output_counts));
  sdsl::int_vector_buffer <> counts(tmp,std::ios::out);
  
  
  long unsigned int i = 0;
  std::map <uint64_t, uint64_t> idx;
  while ( std::getline(finput_counts,line) )
  { 
    std::vector <std::string> row (split(line)); 
    uint64_t i_kmer = str_to_int(row.at(0).c_str(),kmer_length);
    
    for (size_t j = 0; j < nsamples; ++j)
    {
      counts[(i*nsamples)+j] = std::stoll(row.at(j+1));
    }
        
    kmers.push_back(i_kmer);
    idx.emplace(i_kmer,i);

    ++i;
  }
    
  finput_counts.close();
  finput_counts.clear(); 
  counts.close(); 

  etime = ticking();
  std::cout << "elapsed: " << etime-btime << " s" << std::endl;

  bphf64_t * bphf = building_saving_mphf(output_hash,kmers);
  sorting_by_hash(kmers,bphf); 
  delete bphf;

  std::cout << "loading counts into vlc_vector" << std::endl;
  btime = ticking();
  
   
  sdsl::int_vector_buffer <> input_tmp(tmp,std::ios::in); 
  sdsl::vlc_vector <sdsl::coder::fibonacci> input_tmp_vlc (input_tmp);
  input_tmp.close();
  etime = ticking();

  std::cout << "elapsed: " << etime-btime << "s" << std::endl;

  std::cout << "permutates counts" << std::endl;

  std::string tmp1 ("tmp1_"+std::string(output_counts));
  sdsl::int_vector_buffer <> output_tmp (tmp1,std::ios::out); 
  
  btime = ticking();
  long unsigned int sz = kmers.size();
  for (size_t l = 0; l < sz; ++l)
  {
    size_t m = idx.at(kmers.at(l));
    idx.erase(kmers.at(l));
    for (size_t j = 0; j < nsamples; ++j)
    {
      output_tmp[(l*nsamples)+j] = input_tmp_vlc[(m*nsamples)+j]; 
    } 

    if ( (l+1) % trace_sz == 0 )
    {
      etime = ticking();
     
      std::cout << (l+1) << "/" << sz << std::endl;
      std::cout << "elapsed:" << etime-btime << "s" << std::endl;
      
      btime = ticking();
    }

  }
  output_tmp.close();

  etime = ticking();
  std::cout << "elapsed: " << etime-btime << "s" << std::endl;

  std::cout << "writing as vlc_vector" << std::endl;
  btime = ticking();

  sdsl::int_vector_buffer <> input_output_tmp (tmp1,std::ios::in);
  std::ofstream foutput_counts (output_counts);
  sdsl::vlc_vector < sdsl::coder::fibonacci > encoded_counts(input_output_tmp); 
  input_output_tmp.close();

  encoded_counts.serialize(foutput_counts);

  foutput_counts.close();
    
  etime = ticking();

 

  std::cout << "elapsed: " << etime-btime << " s" << std::endl;
  return (EXIT_SUCCESS);

}
int query_buffer(int argc, char * argv [])
{ 
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " <input_counts_buffer> <input_hash> <kmer> <nsamples>" << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv[1];
  char * input_hash = argv[2];
  char * kmer = argv[3];
  size_t nsamples = atoi(argv[4]);

  uint64_t i_kmer = str_to_int(kmer,kmer_length);

  bphf64_t * bphf = load_mphf(input_hash);

  uint64_t key = bphf->lookup(i_kmer);
  if ( key >= bphf->size() )
  {
    std::cerr << kmer << " not_found" << std::endl;
    exit(EXIT_FAILURE);
  }
  delete bphf;
  
  sdsl::int_vector_buffer <> counts(std::string(input_counts),std::ios::in);

  std::cout << "key=" << key << std::endl; 
  for (size_t j = 0; j < nsamples; ++j)
  {
    std::cout << counts[(key*nsamples)+j] << " ";
  }

  std::cout << std::endl;
  counts.close();

  return (EXIT_SUCCESS);
}




int query (int argc, char * argv [])
{ 
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " <input_counts> <input_hash> <kmer> <nsamples>" << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv[1];
  char * input_hash = argv[2];
  char * kmer = argv[3];
  size_t nsamples = atoi(argv[4]);

  uint64_t i_kmer = str_to_int(kmer,kmer_length);

  bphf64_t * bphf = load_mphf(input_hash);

  uint64_t key = bphf->lookup(i_kmer);
  if ( key >= bphf->size() )
  {
    std::cerr << kmer << " not_found" << std::endl;
    exit(EXIT_FAILURE);
  }
  delete bphf;
  
  sdsl::vlc_vector <sdsl::coder::fibonacci> counts;
  std::ifstream finput_counts (input_counts);
  counts.load(finput_counts);
  finput_counts.close();

  std::cout << "key=" << key << std::endl; 
  for (size_t j = 0; j < nsamples; ++j)
  {
    std::cout << counts[(key*nsamples)+j] << " ";
  }

  std::cout << std::endl;

  return (EXIT_SUCCESS);
}

int queries (int argc, char * argv [])
{ 
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " <input_counts> <input_hash> <kmer> <nsamples>" << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv[1];
  char * input_hash = argv[2];
  char * input_kmers = argv[3];
  size_t nsamples = atoi(argv[4]);

  double btime, etime;
  

  std::vector <uint64_t> kmers ;
  
  std::cout << "loading kmers" << std::endl;
  btime=  ticking();
  std::ifstream finput_kmers (input_kmers);
  std::string line;
  
  while ( finput_kmers >> line )
  { 
    uint64_t i_kmer = str_to_int(line.c_str(),kmer_length);
    kmers.push_back(i_kmer);
  }
  finput_kmers.close();
  etime = ticking();

  std::cout << "elapsed: " << etime-btime << "s" << std::endl;
  
  std::cout << "loading hash" << std::endl;
  btime = ticking();
  bphf64_t * bphf = load_mphf(input_hash);
  etime = ticking();

  std::cout << "elapsed: " << etime-btime << "s" << std::endl;
  
  std::cout << "loading counts" << std::endl;
  btime = ticking();
  sdsl::vlc_vector <sdsl::coder::fibonacci> counts;
  std::ifstream finput_counts (input_counts);
  counts.load(finput_counts);
  finput_counts.close();
  etime = ticking();
  
  std::cout << "elapsed: " << etime-btime << "s" << std::endl;

  std::cout << "find counts " << std::endl;
  btime = ticking();
  for (uint64_t i_kmer : kmers )
  {
    uint64_t key = bphf->lookup(i_kmer);
    std::cout << int_to_str(i_kmer) ;
    if ( key >= bphf->size() )
    {
     std::cout << " not_found" << std::endl;
    }
    else
    {
      for (size_t j = 0; j < nsamples; ++j)
      {
	std::cout << " " << counts[(key*nsamples)+j];
      }
      std::cout << std::endl;
    }
  }
  etime = ticking();
  std::cout << "elapsed:" << etime-btime << "s" << std::endl;

  return (EXIT_SUCCESS);
}

int as_vlc (int argc, char * argv [])
{
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " <input counts buffer> <output vlc > " << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_buffer = argv[1];
  char * output_counts = argv[2];


  double btime, etime ;

  std::cout << "vlc_vector from int_vector_buffer" << std::endl;
  btime = ticking();
  sdsl::int_vector_buffer <> finput_buffer (std::string(input_buffer),std::ios::in);
  sdsl::vlc_vector < sdsl::coder::fibonacci > encoded_counts(finput_buffer); 
  etime = ticking();
  std::cout << "elapsed: " << etime-btime << "s" << std::endl;
  
  std::cout << "serialization " << std::endl;
  std::ofstream foutput_counts (output_counts,std::ios::binary);
  encoded_counts.serialize(foutput_counts);
  foutput_counts.close();
    
  etime = ticking(); 

  std::cout << "elapsed: " << etime-btime << " s" << std::endl;
 
  

}


int main (int argc, char * argv [])
{

  if ( strcmp(argv[1],"index") == 0)
    return index(argc-1,argv+1);
  else if ( strcmp(argv[1],"query") == 0)
    return query(argc-1,argv+1);
  else if ( strcmp(argv[1],"queries") == 0)
    return queries(argc-1,argv+1);
  else if ( strcmp(argv[1],"as_vlc") == 0)
    return as_vlc (argc-1,argv+1);
  else if ( strcmp(argv[1],"query_buffer") == 0)
    return query_buffer(argc-1,argv+1);
  return (EXIT_SUCCESS);
}
