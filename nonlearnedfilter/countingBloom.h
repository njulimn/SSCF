#ifndef COUNTINGBLOOM
#define COUNTINGBLOOM

#include <iostream>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <limits.h>
#include <exception> 
#include <string.h>
#include <chrono> 
#include <fstream>
#include <sys/mman.h>
#include "../util/key.h"
#include "../util/hashutil.h"
#include <algorithm>

#define HASH_NUM 3
#define ALLOCATION_RATIO 4 // Allocation ratio between the size of bloom and hashexpressor
using namespace std;

namespace countingBloom{

const int hashset_size_ = 7;

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

class CountingBloom{
    
    counting_bloom_header_t *header;
    uint32_t *hashes; //by xp: the point to hash function point lists
    uint32_t *offset_per_range;
    size_t nfuncs; //by xp: number of hash functions
    size_t num_bytes; //by xp: space occupied by this counting bloom filter
    bitmap_t *bitmap;
    size_t bitPerCounter;
    size_t offset;
    size_t numCounter;
    public:
    CountingBloom(size_t numCounter,size_t bitPerCounter, size_t nfunc);
    ~CountingBloom();
    //counting_bloom_t *new_counting_bloom_from_file(unsigned int capacity, double error_rate, const char *filename);
    int counting_bloom_add(Slice &key);
    int counting_bloom_remove(Slice &key);
    int counting_bloom_check(Slice &key);

    size_t GetHaveInsertNum(){return header->count;};//已经插入到bloom中的个数

    //bitmap operation
    bitmap_t *bitmap_resize(bitmap_t *bitmap, size_t old_size, size_t new_size);
    bitmap_t *new_bitmap(size_t bytes);
    int bitmap_increment_adapt(bitmap_t *bitmap, unsigned int index, long offset);
    int bitmap_decrement(bitmap_t *bitmap, unsigned int index, long offset);
    int bitmap_check_adapt(bitmap_t *bitmap, unsigned int index, long offset);
    int bitmap_flush(bitmap_t *bitmap);
    void free_bitmap();
    void printBloom();
    void printHaveInsertNum(){
        cout<<"After delete, the number of item in CBF: "<<this->header->count<<endl;
    }
    
};

//能够根据要求每个counter的比特数bitPerCounter进行初始化
CountingBloom::CountingBloom(size_t _numCounter,size_t _bitPerCounter, size_t nfunc)
{
    
    numCounter=_numCounter;
    // cout<<"numCounter: "<<numCounter<<endl;
    bitPerCounter = _bitPerCounter;
    //整个bloom可使用的空间sumSize, byte
    size_t bloomSize=(numCounter*bitPerCounter+7)/8;
    bitmap = NULL;
    offset = 0 + sizeof(counting_bloom_header_t); //by xp: because all filters are sorted
    this->nfuncs = nfunc; //by xp: return the smallest integral value that is not less than x
    /* rounding-up integer divide by 2 of bloom->size */
    this->num_bytes = (numCounter*bitPerCounter+7)/8 + sizeof(counting_bloom_header_t); //by xp: the size of (the header + the bitmap) of this filter
    this->hashes = (uint32_t *)calloc(this->nfuncs, sizeof(uint32_t));
    this->bitmap = new_bitmap(this->num_bytes);
    //cout<<"have new bitmap\n";
    this->header = (counting_bloom_header_t *)(this->bitmap->array);
    

}

CountingBloom::~CountingBloom()
{
        if(hashes != NULL){
            free(hashes);
            hashes = NULL;
        }else{
            std::cout<<"error in free\n";
        }
        // if(bitmap->array != NULL){
        //     free(bitmap->array);
        //     bitmap->array = NULL;
        // }else{
        //     std::cout<<"error in free\n";
        // }


        // if(bitmap != NULL){
        //     free(bitmap);
        //     bitmap=NULL;
        // }else{
        //     std::cout<<"error in free\n";
        // }
        // if(header!=NULL){
        //     free(header);
        //     header=NULL;
        // }else{
        //     std::cout<<"error in free\n";
        // }

        free_bitmap();
        // if(bitmap->array != NULL){
        //     // std::cout<<"aa\n"<<endl;
        //     free(bitmap->array);
        //     // bitmap->array=NULL;
        //     // std::cout<<"bb\n"<<endl;
        // }else{
        //     std::cout<<"error in free\n";
        // }

        // if(bitmap != NULL){
        //     // std::cout<<"aa\n"<<endl;
        //     free(bitmap);
        //     // std::cout<<"bb\n"<<endl;
        // }

}
int CountingBloom::counting_bloom_add(Slice &key){
    
    
    for(uint8_t j=0; j<nfuncs; j++){         
        unsigned int bitpos = XXH3_128bits_withSeed(&key.str[0], key.str.size(), j+1).low64 % numCounter;

        bitmap_increment_adapt(this->bitmap, bitpos, offset);
    }
    this->header->count++;
    return 0;
}



int CountingBloom::counting_bloom_check(Slice &key){


    for(int i=0 ; i<nfuncs; i++){
        uint32_t bitpos = XXH3_128bits_withSeed(&key.str[0], key.str.size(), i+1).low64 % numCounter;
        if (!(bitmap_check_adapt(this->bitmap, bitpos, offset))) {
            return 0;
        }
    }
    return 1;
}

int CountingBloom::counting_bloom_remove(Slice &key)
{
    // if(!counting_bloom_check(key)) return 0;
    
    for (size_t i = 0; i < nfuncs; i++) {
        uint32_t bitpos = XXH3_128bits_withSeed(&key.str[0], key.str.size(), i+1).low64 % numCounter;
        bitmap_decrement(this->bitmap, bitpos, this->offset);
    }
    this->header->count--;
    
    return 1;
}



void CountingBloom::free_bitmap()
{
    /*if ((munmap(bitmap->array, bitmap->bytes)) < 0) {
        perror("Error, unmapping memory");
    }
    close(bitmap->fd);
    */
    free(bitmap->array);
    free(bitmap);
}



bitmap_t *CountingBloom::bitmap_resize(bitmap_t *bitmap, size_t old_size, size_t new_size)
{
/*    int fd = bitmap->fd; by xp: why using fd to transfer data?
	struct stat fileStat;
    
    fstat(fd, &fileStat); by xp: fd to be stat-ed by struct stat
    size_t size = fileStat.st_size;
    */
    /* grow file if necessary */
/*    if (size < new_size) {
        if (ftruncate(fd, new_size) < 0) {  by xp: truncate fd to new_size length
            perror("Error increasing file size with ftruncate");
            free_bitmap(bitmap);
            close(fd);
            return NULL;
        }
    }
    lseek(fd, 0, SEEK_SET); //by xp: repositions the offset 0 of fd
    */
    /* resize if mmap exists and possible on this os, else new mmap */
    if (bitmap->array != NULL) {
#if __linux
        //bitmap->array = mremap(bitmap->array, old_size, new_size, MREMAP_MAYMOVE);
   	  bitmap->array = (char *)realloc(bitmap->array, new_size); //注意。。。。。。。。。。。。。。。。。。。。。。。。
      
        if (bitmap->array == NULL) {
            perror("Error resizing realloc");
            free_bitmap();
           // close(fd);
            return NULL;
        }
#else
/*        if (munmap(bitmap->array, bitmap->bytes) < 0) {
            perror("Error unmapping memory");
            free_bitmap(bitmap);
           // close(fd);
            return NULL;
        }
        bitmap->array = NULL; */
#endif
    }
    if (bitmap->array == NULL) {
//        bitmap->array = mmap(0, new_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	    bitmap->array = (char *) malloc(new_size*sizeof(char));
        // bitmap->array=new char[new_size*sizeof(char)];
        memset(bitmap->array,0,new_size*sizeof(char));
        if (bitmap->array == NULL) {
            perror("Error init bitmap array, can't alloc mem for it.\n");
            free_bitmap();
            //close(fd);
            return NULL;
        }
    }
    
    bitmap->bytes = new_size;
    return bitmap;
}

bitmap_t *CountingBloom::new_bitmap(size_t bytes)
{
    bitmap_t *bitmap;
    
    if ((bitmap = (bitmap_t *)malloc(sizeof(bitmap_t))) == NULL) {
        std::cout<<"not bitmap_resize\n"<<std::endl;
        return NULL;
    }
    bitmap->bytes = bytes;
//    bitmap->fd = fd; //by xp: file handle for bloom_file
    bitmap->array = NULL;
    
    if ((bitmap = bitmap_resize(bitmap, 0, bytes)) == NULL) {
        
        return NULL;
    }
    memset(bitmap->array, 0, bitmap->bytes);
    
    return bitmap;
}




int CountingBloom::bitmap_increment_adapt(bitmap_t *bitmap, unsigned int index, long offset)
{

    
    long access = (index*bitPerCounter)/8 + offset; //by xp: ???
    if(bitPerCounter==4){
     
        uint8_t temp;
        uint8_t n = bitmap->array[access];
        if (index % 2 != 0) { //by xp: counter increasement
            temp = (n & 0x0f);
            n = (n & 0xf0) + ((n & 0x0f) + 0x01);
        } else {
            temp = (n & 0xf0) >> 4;
            n = (n & 0x0f) + ((n & 0xf0) + 0x10);
        }
    
        if (temp == 0x0f) {
            fprintf(stderr, "Error, 4 bit int Overflow\n");
            return -1;
        }
        bitmap->array[access] = n;
        

     }else if(bitPerCounter==2){
        
        size_t j=index%4;
        uint8_t temp;
        uint8_t n = bitmap->array[access];
        temp=(n>>j*2)& 0x03;
        if(temp == 0x03){  //如果已经满了
            fprintf(stderr, "Error, 4 bit int Overflow\n");
            return -1;          
        }
        uint8_t addNum=1ULL << bitPerCounter*j;
        n=n+addNum;
        bitmap->array[access] = n;
    }
    return 0;
}


/* increments the four bit counter */
int CountingBloom::bitmap_decrement(bitmap_t *bitmap, unsigned int index, long offset)
{
    long access = index / 2 + offset;
    uint8_t temp;

    uint8_t n = bitmap->array[access];
    
    if (index % 2 != 0) {
        temp = (n & 0x0f);
        n = (n & 0xf0) + ((n & 0x0f) - 0x01);
    } else {
        temp = (n & 0xf0) >> 4;
        n = (n & 0x0f) + ((n & 0xf0) - 0x10);
    }
    
    if (temp == 0x00) {
        fprintf(stderr, "Error, Decrementing zero\n");
        return -1;
    }
    
    bitmap->array[access] = n;
    return 0;
}


/* decrements the four bit counter */
int CountingBloom::bitmap_check_adapt(bitmap_t *bitmap, unsigned int index, long offset)
{
    if(bitPerCounter==4){
        long access = index / 2 + offset;
        if (index % 2 != 0 ) {
           return bitmap->array[access] & 0x0f;
        } else {
           return bitmap->array[access] & 0xf0;
        }
    }else if(bitPerCounter==2){
        long access = (index*bitPerCounter)/8 + offset;
        size_t j=index%4;
        uint8_t temp;
        uint8_t n = bitmap->array[access];
        return temp=(n>>(j*2))& 0x03;        
     }
    return 0;

}

void CountingBloom::printBloom(){
    std::cout<<"bitmap->bytes"<<bitmap->bytes<<std::endl;
    for( size_t i=0;i<bitmap->bytes;i++){
        for(int j=0;j<4;j++){
          if(((bitmap->array[i]>>j*2)& 0x03)!=0x00){
               std::cout<<int((bitmap->array[i]>>j*2)& 0x03)<<std::endl;
          }
          
        }


    }
}

/* by xp: this function is deprecated by our demands.*/
int CountingBloom::bitmap_flush(bitmap_t *bitmap)
{
    if ((msync(bitmap->array, bitmap->bytes, MS_SYNC) < 0)) {
        perror("Error, flushing bitmap to disk");
        return -1;
    } else {
        return 0;
    }
}












}
#endif