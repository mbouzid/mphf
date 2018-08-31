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
 
  std::cerr << "building and saving hash" << std::endl;
  btime = ticking();
  
  bphf64_t * bphf = build_mphf(kmers,nthreads);
  save_mphf (output_hash,bphf);
  
  etime = ticking();
  std::cerr << "elapsed: " << etime-btime << " s" << std::endl;
  return bphf;
}

int sorting_by_hash (std::vector <uint64_t> & kmers, bphf64_t * bphf)
{
  double btime, etime ;
 
  auto by_hashkey = [&bphf] (uint64_t a, uint64_t b) { return (bphf->lookup(a) < bphf->lookup(b)); };
  std::cerr << "sorting kmers by hash key" << std::endl; 

  btime = ticking();
  std::sort (kmers.begin(), kmers.end(), by_hashkey);
  etime = ticking();

  std::cerr << "elapsed: " << etime-btime << " s" << std::endl;
  
  return true;
}

int save_kmers (const std::vector <uint64_t> & kmers, const char * fname)
{ 
  double btime, etime;

  long unsigned int nelem = kmers.size();
  sdsl::int_vector<64> mers(nelem);

  for (long unsigned int i = 0; i < nelem; ++i)
  {
    mers[i] = kmers.at(i);
  }

  sdsl::vlc_vector <sdsl::coder::fibonacci> encoded (mers); 

  std::ofstream foutput_kmers (fname);
  encoded.serialize(foutput_kmers);

  foutput_kmers.close(); 
  
  return true;
}



int index (int argc, char * argv [])
{
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " <input_counts> <output_hash> <output_counts> <compress:y/n>" << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv [1];	    // input matrix counts
  char * output_hash = argv [2];	    // output mphf
  char * output_counts = argv[3];	    // output matrix counts
  bool compress = ( (strcmp(argv[4],"y")==0) ? true : false);
  

  std::vector <uint64_t> kmers;		    // vector of kmers for building the mphf
  size_t nsamples;
  double btime, etime ;
  long unsigned int nelem; 

  std::ifstream finput_counts (input_counts);
  std::string line;

  std::getline(finput_counts,line);	    // get header
  nsamples = split(line).size() -1;    // number of samples = number of elements in header -1 (tag)
  std::cerr << "nsamples=" << nsamples << std::endl;
  std::cerr << "loading kmers & offsets, saving counts" << std::endl;
  btime = ticking();
  
  std::string tmp("tmp_"+std::string(output_counts));
  std::string kmers_output("kmers_"+std::string(output_counts));
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
  std::cerr << "elapsed: " << etime-btime << " s" << std::endl;

  nelem = kmers.size();

  bphf64_t * bphf = building_saving_mphf(output_hash,kmers);
  sorting_by_hash(kmers,bphf); 
  delete bphf;
  
  
  std::cerr << "saving kmers " << std::endl;
  save_kmers(kmers,kmers_output.c_str());

  std::cerr << "loading counts into vlc_vector" << std::endl;
  btime = ticking();
  
   
  sdsl::int_vector_buffer <> input_tmp(tmp,std::ios::in); 
  sdsl::vlc_vector <sdsl::coder::fibonacci> input_tmp_vlc (input_tmp);
  input_tmp.close();
  etime = ticking();

  std::cerr << "elapsed: " << etime-btime << "s" << std::endl;

  std::cerr << "permutates counts" << std::endl;

  std::string tmp1 ("tmp1_"+std::string(output_counts));
  sdsl::int_vector_buffer <> output_tmp (tmp1,std::ios::out); 
  
  btime = ticking();

  
  
  for (size_t l = 0; l < nelem; ++l)
  {
    long unsigned int m = idx.at(kmers.at(l)); 
    for (size_t j = 0; j < nsamples; ++j)
    {
      output_tmp[(l*nsamples)+j] = input_tmp_vlc[(m*nsamples)+j]; 
    } 
  

  }
  output_tmp.close();
    


  etime = ticking();
  std::cerr << "elapsed: " << etime-btime << "s" << std::endl;

  input_tmp_vlc = sdsl::vlc_vector<sdsl::coder::fibonacci>();

  if (compress)
  {
    std::cerr << "writing as vlc_vector" << std::endl;
    btime = ticking();

    sdsl::int_vector_buffer <> input_output_tmp (tmp1,std::ios::in);
    std::ofstream foutput_counts (output_counts);
    sdsl::vlc_vector < sdsl::coder::fibonacci > encoded_counts (input_output_tmp); 
    input_output_tmp.close();

    encoded_counts.serialize(foutput_counts);

    foutput_counts.close();
    
    etime = ticking();

    std::cerr << "elapsed: " << etime-btime << " s" << std::endl;
  }
  else
  {
    std::cerr << "writing as uint32_t" << std::endl;

    btime = ticking();

    std::ofstream foutput_counts (output_counts,std::ios::binary);
    sdsl::int_vector_buffer <> input_output_tmp (tmp1,std::ios::in);
   
    for(long unsigned int i = 0; i < nelem; ++i)
    { 
      uint32_t * cnts = new uint32_t [nsamples](); 
      for(size_t j = 0; j < nsamples; ++j)
      {
	cnts[j] = (uint32_t)input_output_tmp[(i*nsamples)+j];
      }
      

      foutput_counts.write(reinterpret_cast<const char *>(cnts),sizeof(uint32_t)*nsamples);  
      delete [](cnts);
    }
   
    foutput_counts.close();
    input_output_tmp.close();

    etime = ticking();
  
    std::cerr << "elapsed: " << etime-btime << " s" << std::endl;
  }
 
  remove(tmp.c_str());
  remove(tmp1.c_str());
  return (EXIT_SUCCESS);

}



int find (int argc, char * argv [])
{
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " <input_kmers> <kmer>" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  char * input_kmers = argv[1];
  char * kmer = argv[2];
  double btime, etime;

  std::cout << "loading vlc_vector" << std::endl;
  btime = ticking();
  sdsl::vlc_vector <sdsl::coder::fibonacci> kmers;
  std::ifstream finput_kmers (input_kmers);
  kmers.load(finput_kmers);
  finput_kmers.close();
  etime = ticking();
  std::cout << "elapsed: " << etime-btime << "s"<< std::endl;
  bool found = ( std::find(kmers.begin(),kmers.end(),str_to_int(kmer,kmer_length)) != kmers.end() );
  std::cout << kmer << "\t" << ( found ? "found" : "not_found" ) << std::endl;
  
  return (EXIT_SUCCESS);
}


int find_some(int argc, char * argv [])
{
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " <input_kmers> <input_kmers_list>" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  char * input_kmers = argv[1];
  char * input_kmers_list = argv[2];
  double btime, etime;

  std::cout << "loading vlc_vector" << std::endl;
  
  btime = ticking();
  sdsl::vlc_vector <sdsl::coder::fibonacci> kmers;
  std::ifstream finput_kmers (input_kmers);
  kmers.load(finput_kmers);
  finput_kmers.close();
  etime = ticking();
  
  std::ifstream finput_kmers_list (input_kmers_list);
  std::string line;
  std::vector <std::string> kmers_list;
  while (std::getline(finput_kmers_list,line))
  {
    kmers_list.push_back(line);
  }
  finput_kmers_list.close();

  for (const auto & kmer : kmers_list)
  {
    bool found = ( std::find(kmers.begin(),kmers.end(),str_to_int(kmer.c_str(),kmer_length)) != kmers.end() );
    std::cout << kmer << "\t" << ( found ? "found" : "not_found" ) << std::endl;
  }

  return (EXIT_SUCCESS);
}

int query (int argc, char * argv [])
{ 
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " <input_counts> <input_hash> <nsamples> <kmer>" << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv[1];
  char * input_hash = argv[2];
  size_t nsamples = atoi(argv[3]);
  char * kmer = argv[4];

  double btime, etime;
  uint64_t i_kmer = str_to_int(kmer,kmer_length);

 // std::cout << "loading hash" << std::endl;
  btime = ticking();
  bphf64_t * bphf = load_mphf(input_hash);
  etime = ticking();

 // std::cout << "elapsed: " << etime-btime << "s" << std::endl;

  uint64_t key = bphf->lookup(i_kmer);
  if ( key >= bphf->size() )
  {
    std::cout << kmer << " not_found" << std::endl;
    exit(EXIT_FAILURE);
  }
  delete bphf;
  
  //std::cout << "loading vlc_vector" << std::endl;
  btime = ticking();
  sdsl::vlc_vector <sdsl::coder::fibonacci> counts;
  std::ifstream finput_counts (input_counts);
  counts.load(finput_counts);
  finput_counts.close();
  etime = ticking();

  std::cout << kmer ; 
  for (size_t j = 0; j < nsamples; ++j)
  {
    std::cout << "\t" << counts[(key*nsamples)+j] ;
  }

  std::cout << std::endl;
  //std::cout << "elapsed: " << etime-btime << "s" << std::endl;

  return (EXIT_SUCCESS);
}

int queries (int argc, char * argv [])
{ 
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " <input_counts> <input_hash> <kmers> <nsamples>" << std::endl;
    exit(EXIT_FAILURE);
  }

  char * input_counts = argv[1];
  char * input_hash = argv[2];
  char * input_kmers = argv[3];
  size_t nsamples = atoi(argv[4]);

  double btime, etime;
  

  std::vector <uint64_t> kmers ;
  
 // std::cout << "loading kmers" << std::endl;
  btime=  ticking();
  std::ifstream finput_kmers (input_kmers);
  std::string line;
  
  while ( std::getline(finput_kmers,line) )
  { 
    uint64_t i_kmer = str_to_int(line.c_str(),kmer_length);
    kmers.push_back(i_kmer);
  }
  finput_kmers.close();
  etime = ticking();
  
  //std::random_shuffle(kmers.begin(),kmers.end());

 // std::cout << "elapsed: " << etime-btime << "s" << std::endl;
  
  //std::cout << "loading hash" << std::endl;
  btime = ticking();
  bphf64_t * bphf = load_mphf(input_hash);
  etime = ticking();

//  std::cout << "elapsed: " << etime-btime << "s" << std::endl;
  
 // std::cout << "loading counts" << std::endl;
  btime = ticking();
  sdsl::vlc_vector <sdsl::coder::fibonacci> counts;
  std::ifstream finput_counts (input_counts);
  counts.load(finput_counts);
  finput_counts.close();
  etime = ticking();
  
  //std::cout << "elapsed: " << etime-btime << "s" << std::endl;

 // std::cout << "find counts " << std::endl;
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
	std::cout << "\t" << counts[(key*nsamples)+j];
      }
      std::cout << std::endl;
    }
  }
  etime = ticking();
  //std::cout << "elapsed:" << etime-btime << "s" << std::endl;

  return (EXIT_SUCCESS);
}

int query_disk (int argc, char * argv [])
{
  if (argc < 2)
  {
    std::cerr << argv[0] << " <input_count> <input_hash> <ncol> <k-mer>" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  char * input_counts = argv[1];
  char * input_hash = argv[2];
  size_t ncol = atoi(argv[3]);
  char * kmer = argv[4];
  

  double btime, etime;

  //std::cout << "loading hash" << std::endl;
  btime = ticking();

  bphf64_t * bphf = load_mphf(input_hash); 
  
  etime = ticking();
  //std::cout << "elapsed: " << etime-btime << "s" << std::endl; 

  uint64_t key = bphf->lookup(str_to_int(kmer,kmer_length));
  delete bphf;

  std::cout << kmer ;
  if (key < ULLONG_MAX)
  {

    std::ifstream finput_counts (input_counts,std::ios::binary);  
    off_t sz = sizeof(uint32_t)*ncol;
    finput_counts.seekg(key*sz);
    uint32_t * counts = new uint32_t[ncol]();
    finput_counts.read (reinterpret_cast <char *> (counts),sizeof(uint32_t)*ncol);
    finput_counts.close();
    
    for (size_t i = 0; i < ncol; ++i)
    {
      std::cout << "\t" << counts[i] ;
    }
    std::cout << std::endl;
    
    delete [](counts);
  }
  else
  {
    std::cout << "\tnot_found"<< std::endl;
  }
  
  return (EXIT_SUCCESS);

}


int queries_disk (int argc, char * argv [])
{
  if (argc < 4)
  {
    std::cerr << argv[0] << " <input_count> <input_hash> <input_tags> <ncol>" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  char * input_counts = argv[1];
  char * input_hash = argv[2];
  char * input_tags = argv[3];
  
  size_t ncol = atoi(argv[4]);
  

  double btime, etime;

  //std::cout << "loading hash" << std::endl;
  btime = ticking();

  bphf64_t * bphf = load_mphf(input_hash); 
  
  etime = ticking();
  //std::cout << "elapsed: " << etime-btime << "s" << std::endl; 

  //std::cout << "loading tags" << std::endl;
  btime = ticking();
  std::vector <uint64_t> kmers;

  std::ifstream finput_tags (input_tags);
  std::string line;
  while ( std::getline(finput_tags,line) )
  {
    uint64_t i_kmer = str_to_int(line.c_str(),kmer_length);
    kmers.push_back(i_kmer);
  }

  finput_tags.close();
  etime = ticking();
  
  //std::cout << "elapsed: " << etime-btime << "s" << std::endl; 


  std::ifstream finput_counts (input_counts,std::ios::binary);  
  
  off_t sz = sizeof(uint32_t)*ncol;
  long unsigned int bphfsz = bphf->size();
  
  btime = ticking();
  for (uint64_t i_kmer : kmers)
  {
    uint64_t key = bphf->lookup(i_kmer);
    std::cout << int_to_str(i_kmer);
    if ( key < ULLONG_MAX )
    {
      finput_counts.seekg(key*sz);
      uint32_t * counts = new uint32_t[ncol]();
      finput_counts.read (reinterpret_cast <char *> (counts),sz);

      for (size_t i = 0; i < ncol; ++i)
      {
	std::cout << "\t" << (uint32_t)counts[i] ;
      }
      std::cout << std::endl;
    
      delete[](counts);
    }
    else
    {
      std::cout << "\tnot_found" << std::endl;
    }
  }
  finput_counts.close();
  
  etime = ticking();
  //std::cout << "elapsed: " << etime-btime << "s" << std::endl; 



  delete bphf;
  return (EXIT_SUCCESS);

}



int main (int argc, char * argv [])
{
  srand(time(NULL));
  
  if (argc < 2)
  {
    std::cerr << "**-" << argv[0] << "-**"<< std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << "\tindex\t<input_counts> <output_hash> <output_counts> <compress:y/n>" << std::endl;
    std::cerr << "\tquery/query_disk\t <input_counts> <input_hash> <kmer> <nsamples>" << std::endl;
    std::cerr << "\tqueries/queries_disk\t <input_counts> <input_hash> <kmers> <nsamples>" << std::endl;
    std::cerr << "\tfind\t<input_kmers> <kmer> " << std::endl;
    std::cerr << "\tfind_some\t<input_kmers> <kmers_list> " << std::endl;
    exit(EXIT_FAILURE);
  }

  if ( strcmp(argv[1],"index") == 0)
    return index(argc-1,argv+1);
  else if ( strcmp(argv[1],"query") == 0)
    return query(argc-1,argv+1);
  else if ( strcmp(argv[1],"queries") == 0)
    return queries(argc-1,argv+1);
  else if ( strcmp(argv[1],"query_disk") == 0)	
    return query_disk(argc-1,argv+1);
  else if ( strcmp (argv[1],"queries_disk") == 0)
    return queries_disk(argc-1,argv+1);
  else if ( strcmp(argv[1],"find") == 0)
    return find(argc-1,argv+1);
  else if ( strcmp(argv[1],"find_some") == 0)
    return find_some(argc-1,argv+1);
  return (EXIT_SUCCESS);
}
