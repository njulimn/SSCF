#ifndef PUBLICSTRUCT
#define PUBLICSTRUCT
#include "../util/key.h"
#define HASH_NUM 3
namespace SeesawCF
{

struct HABFilterKey{
    Slice * data_;   //key的string和cost大小
    uint8_t hash_map_[HASH_NUM];  //哈希函数的数组
};

struct Tuple_V{
    bool singleflag_;
    HABFilterKey * keyid;
    Tuple_V(): singleflag_(true), keyid(NULL){};
};

struct HashExpressorCell{
    bool endbit;
    uint8_t hash_idx_;
};

typedef struct {
    size_t bytes;
    //int    fd; //by xp: what is it?
    char  *array;
} bitmap_t;


typedef struct {
    uint64_t id;
    uint32_t count; //by xp: number of elements in this bloom filter
    uint32_t _pad;
} counting_bloom_header_t;










} // namespace dynamichabf





#endif