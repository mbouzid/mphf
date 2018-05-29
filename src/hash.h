#ifndef hash_h
#define hash_h

#include <string>

class mer_hasher
{
  
  public:
  typedef std::string Item;
  typedef std::pair<uint64_t,uint64_t> hash_pair_t;

  public:
  hash_pair_t operator ()  (const Item& key) const  {  
  hash_pair_t result;
  result.first =  singleHasher (key, 0xAAAAAAAA55555555ULL);;
  result.second =  singleHasher (key, 0x33333333CCCCCCCCULL);;

  return result;
  }

  uint64_t singleHasher (std::string key, uint64_t seed=0) const
  {


  uint64_t hash  =  hash_fn(key);
  hash ^= seed;

  return hash;
  }

   std::hash<std::string> hash_fn;
  
};



class mer_hasher_uint64
{
  
  public:
  typedef uint64_t Item;
  typedef std::pair<uint64_t,uint64_t> hash_pair_t;

  public:

  hash_pair_t operator ()  (const Item& key) const  {   
    hash_pair_t result;
    result.first =  singleHasher (key, 0xAAAAAAAA55555555ULL);;
    result.second =  singleHasher (key, 0x33333333CCCCCCCCULL);;

  return result;
  }

  uint64_t singleHasher (uint64_t key, uint64_t seed=0) const
  {
    
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccd;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53;
    key ^= key >> 33;
			
    key ^= seed;
				
    return key;
  }

};

#endif
