#ifndef HASHMODULATOR
#define HASHEXPRESSOR_ONE
#define HASH_NUM 3
#include <iostream>
#include <assert.h>
#include "../util/xxhash/xxhash.h"
#include "publicStruct.h"
#include <limits.h>
#include "../util/murmurhash/MurMurHash.h"
using namespace std;
namespace SeesawCF
{

class HashModulator{
    public:
        HashModulator(uint8_t kfuncs, uint32_t bit_size_);
        ~HashModulator();
        bool Add(HABFilterKey &key);
        void Add(HABFilterKey &key, uint8_t *insert_sequence_);
        void Add(Slice &key, size_t cellpos, size_t bitpos_he, uint64_t newHashIdx);
        void AddCounterByOne(Slice &key, size_t cellpos, size_t bitpos_he, size_t cell_counter);
        bool GetInsertSequence(HABFilterKey &key);
        //bool GetInsertSequence(HABFilterKey &key, uint64_t hashValue);
        //size_t Query(HABFilterKey &key, size_t notMapIdx);
        bool Query(Slice &key, uint8_t &deltaIdx);
        void QueryIsEmpty(Slice &key, size_t &cell_counter, size_t &newHashIdx, size_t &cellpos);
        void Clear(size_t cell_counter, size_t cellpos);
        // int getSize(){return 8 * data_.size();};
        int getCellNum(){return cell_num_;};
        
        size_t addOneTime = 0;

        //bool GetInsertSequence(Slice & key, uint64_t hashValue);
        bool GetInsertSequence(Slice & key, size_t &newHashIdx, size_t &counter_cell, size_t &cellpos, size_t & bitpos_he);

        void AssignDefaultHashset_(HABFilterKey *keys_, uint32_t &size_);
        void AssignDefaultHashset_(HABFilterKey &key);

        size_t GetHaveInsertNum(){
            return haveInsertNum;
        }
        size_t GetCellNum(){
            return cell_num_;
        }

        //vector<uint> testHash;
        int hashCount=0;
    private:
        void setbit(const int &bitinx);
        void clearBit(const int &bitinx);
        bool getbit(const int &bitinx);
        void getCell(const int &cellinx, uint8_t &hash_idx);
        void getCell_oldIdx(const int &cellinx, uint8_t &oldIdx);
        bool getCell_deltaIdx(const int &cellinx, uint8_t &oldIdx, uint8_t &deltaIdx);
        void setCell(const int &cellinx, uint8_t hash_idx_);
        void setCell_oldIdx(const int &cellinx, uint8_t oldIdx, uint8_t deltaIdx);

        void FindInsertSequence(uint64_t &hash_key_, std::vector<uint8_t> &hash_map_, std::vector<uint8_t> &insert_sequence_, int left_, int right_, int init_cell_count_, int &min_cell_count_);

        

        uint32_t cell_num_;
        uint8_t cell_size_;
        uint8_t k_;
        uint32_t t;
        // std::string data_;
        char * data_;
        // char * bitArray;
        size_t sizeCounter_;
        size_t sizeOldIndex_;
        size_t sizeDeltaIndex_;
        size_t haveInsertNum = 0;

        

};

uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;   //取末六位
    return (n << c) | ( n >> ((-c) & mask));
}

uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

// for scene of few good hash functions or caculate quickly 
uint32_t GetHashFromHash(uint64_t hash, uint32_t index, uint32_t blockLength) {
    uint32_t r = rotl64(hash, index * 21);
    return (uint32_t) reduce(r, blockLength);
}

bool compare_func(Slice * key1, Slice * key2){
    return key1->cost > key2->cost;
}

HashModulator::HashModulator(uint8_t kfuncs, uint32_t bit_size_){
    t = 0;
    k_ = kfuncs;
    sizeCounter_ = 3;
    sizeDeltaIndex_ = 1;
    cell_size_ = sizeCounter_ + sizeDeltaIndex_;
    cell_num_ = bit_size_ / cell_size_;    
    // std::cout<< "hashexpressor's cell_num: " << cell_num_ << std::endl;

    // data_.resize((cell_num_ * cell_size_ +7 )/ 8 , 0);
    // data_ = (char*)malloc((cell_num_ * cell_size_ +7 )/ 8 );
    data_ = new char[(cell_num_ * cell_size_ +7 )/ 8 ];
    
    memset(data_, 0 , (cell_num_ * cell_size_ +7 )/ 8 );
    // std::cout<<"aa:"<<(cell_num_ * cell_size_ +7 )/ 8<<endl;
    // for(int i =0; i<(cell_num_ * cell_size_ +7 )/ 8;i++){
    //     if(data_[i]!=0)
    //     {
    //         cout<<"yes"<<endl;
    //     }
    // }
    // bitArray = & data_[0];
    //testHash.resize(cell_num_+1);
}
HashModulator::~HashModulator(){
    // std::cout<<"yes\n"<<endl;
    delete []data_;
}



bool HashModulator::Query(Slice &key, uint8_t &deltaIdx){
    
    //uint64_t bitpos =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_) *cell_size_;
    // uint64_t bitpos =  XXH3_128bits_withSeed(&key.str[0], key.str.size(), 0).low64 % cell_num_ * cell_size_;
    //uint64_t bitpos = MurmurHash64B(&key.str[0], key.str.size(), 1) % cell_num_ * cell_size_;

    uint64_t cellpos =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_);
    size_t counter_cell = ((data_[cellpos>>1] >> ( (cellpos & 0x1)<<2)) & 0x0E) >> sizeDeltaIndex_; 
    if(!counter_cell)
    {
        
        return false;
    }else
    {
        deltaIdx +=  ((data_[cellpos>>1] >> ( (cellpos & 0x1)<<2)) & 0x01) ; 
        return true;
    }


    // size_t w = 1;
    // size_t cell_counter=0;

    // for(size_t i = sizeDeltaIndex_; i< cell_size_; i++)
    // { 
    //     cell_counter += w * getbit(bitpos + i);
    //     w = w<<1; 
    // }
    // w=1;
    // deltaIdx += (1 + getbit(bitpos));

    // return cell_counter;
  
}

void HashModulator::QueryIsEmpty(Slice &key, size_t &cell_counter, size_t &newHashIdx, size_t &cellpos){
    
    cellpos =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_);
    size_t temp = (data_[cellpos>>1] >> ((cellpos & 0x1)<<2));
    cell_counter = (temp & 0x0E) >> 1; 


    // if(cell_counter == 0)
    // {

    //     return true;
    // }
    newHashIdx += (temp & 0x01) ; 
    return;
    






    // size_t w = 1;
    // size_t i = sizeDeltaIndex_;

    // bitpos_he =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_) * cell_size_;
    // //uint64_t bitpos =  XXH3_128bits_withSeed(&key.str[0], key.str.size(), 0).low64 % cell_num_ * cell_size_;
    // //uint64_t bitpos = MurmurHash64B(&key.str[0], key.str.size(), 1) % cell_num_ * cell_size_;

    // for(; i < cell_size_; i++)
    // {
    //     cell_counter += w * getbit(bitpos_he + i);
    //     w = w<<1; 
    // }
    // if(!cell_counter){
    //     return true;
    // }else{
    //     newHashIdx += (k_ + 1 + getbit(bitpos_he));
    //     return false;
    // }
    
}

void HashModulator::Clear(size_t cell_counter, size_t cellpos){

    //uint64_t temp =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_);
    // uint64_t bitpos =  XXH3_128bits_withSeed(&key.str[0], key.str.size(), 0).low64 % cell_num_ * cell_size_;
    //size_t temp = MurmurHash64B(&key.str[0], key.str.size(), 1) % cell_num_;
    // testHash[temp]++;

    uint8_t n = data_[cellpos>>1], temp = (cellpos & 0x01)<<2;
    n = n & ~(0x0e << temp) ;
    data_[cellpos>>1] =  n | ( (cell_counter << 1) << temp );

    // bitpos_he =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_) * cell_size_;

    // size_t i = sizeDeltaIndex_;


    // while(cell_counter > 0){
    //     assert(i<cell_size_);
    //     if(cell_counter % 2 > 0)
    //     {
    //         setbit(bitpos_he + i);
    //     }else {
    //         clearBit(bitpos_he + i);
    //     }
    //     i++;
    //     cell_counter = cell_counter >> 1;
    // }
    
}

void HashModulator::Add(HABFilterKey &key, uint8_t *insert_sequence_){
 
    size_t prevHashIdx=0;
    
    uint8_t hash_idx_;
    
    uint64_t hash_key_ = XXH3_64bits(&key.data_->str[0], key.data_->str.size());    
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);

    
    hash_idx_ =  insert_sequence_[k_-1] - k_;
    setCell(h, hash_idx_);

    haveInsertNum++;
}

// void HashExpressor_One::Add(Slice &key, uint64_t oldHashValue, uint8_t newHashIdx){
 

//     uint32_t h = oldHashValue % cell_num_;
//     uint8_t hash_idx_ =  newHashIdx - k_;
//     setCell(h, hash_idx_);

//     haveInsertNum++;
// }

void HashModulator::Add(Slice &key, size_t cellpos, size_t bitpos_he, uint64_t newHashIdx){

    //uint8_t hash_idx_;
    
    //uint32_t bitpos =  GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_) *cell_size_;
    //uint64_t bitpos =  XXH3_128bits_withSeed(&key.str[0], key.str.size(), 0).low64 % cell_num_ * cell_size_;
    //uint64_t bitpos = MurmurHash64B(&key.str[0], key.str.size(), 1) % cell_num_ * cell_size_;

    newHashIdx = newHashIdx - k_ - 1;
    if(newHashIdx)
    {
        data_[cellpos>>1] |= ( 0x01 << ((cellpos & 0x01)<<2) );
    }else{
        data_[cellpos>>1] &= ~( 0x01 << ((cellpos & 0x01)<<2) );         
    }
    data_[cellpos>>1] |= ( 0x02 << ((cellpos & 0x01)<<2) );
    haveInsertNum++;
    //setbit(bitpos_he + sizeDeltaIndex_);

    // newHashIdx = newHashIdx - k_ - 1;

    // if(newHashIdx)
    // {
    //     //assert(bitpos<=data_.size()*8);
    //     char *array = &(data_[0]);
    //     int charIdx = (bitpos_he) >> 3; 
    //     int b = (bitpos_he) % 8;
    //     array[charIdx] |= (1 << b);   
    // }else{
    //     //assert(bitpos<=data_.size()*8);
    //     char *array = &(data_[0]);
    //     int charIdx = (bitpos_he) >> 3; 
    //     int b = (bitpos_he) % 8;
    //     array[charIdx] &= ~(1 << b);           
    // }
    // setbit(bitpos_he + sizeDeltaIndex_);
    // haveInsertNum++;
}

void HashModulator::AddCounterByOne(Slice &key, size_t cellpos, size_t bitpos_he, size_t cell_counter)
{
    if(cell_counter >7)
    {
        std::cout<<"counter overflow!"<<cellpos<<"\n";
        return;
    }
    //assert(cell_counter <=7);
    uint8_t n = data_[cellpos>>1];
    n = n & ~(0x0e << ((cellpos & 0x01)<<2)) ;
    data_[cellpos>>1] =  n | ( (cell_counter << 1) << ((cellpos & 0x01)<<2) );

    // size_t i = sizeDeltaIndex_;
    // while(cell_counter > 0){
    //     assert(i<cell_size_);
    //     if((cell_counter & 0x1) > 0)
    //     {
    //         setbit(bitpos_he + i);
    //     }else {
    //         clearBit(bitpos_he + i);
    //     }
    //     i++;
    //     cell_counter = cell_counter >> 1;
    // }
    
}

bool HashModulator::GetInsertSequence(HABFilterKey &key){

    
    
    uint8_t hash_idx_ = 0;

    uint64_t hash_key_ = XXH3_64bits(&key.data_->str[0], key.data_->str.size()); 
 
    uint32_t h =  GetHashFromHash(hash_key_, 0, cell_num_);

    getCell(h, hash_idx_);

    if(hash_idx_ != 0){
        return false;
    }else{
        return true;
    }
    
    
}

bool HashModulator::GetInsertSequence(Slice & key, size_t &newHashIdx, size_t &counter_cell, size_t &cellpos, size_t & bitpos_he){

    size_t w = 1;
    size_t i = sizeDeltaIndex_;
    cellpos = GetHashFromHash(XXH3_64bits(&key.str[0], key.str.size()), 0, cell_num_);
    //bitpos_he =  cellpos * cell_size_;

    //uint64_t bitpos = MurmurHash64B(&key.str[0], key.str.size(), 1) % cell_num_ * cell_size_;
    counter_cell = ((data_[cellpos>>1] >> ( (cellpos & 0x1)<<2)) & 0x0E) >> sizeDeltaIndex_; 
    if(!counter_cell)
    {
        return true;
    }else
    {
        newHashIdx += (k_ + 1 + ((data_[cellpos>>1] >> ( (cellpos & 0x1)<<2)) & 0x01) ); 
        return false;
    }



    // for(; i < cell_size_; i++)
    // {
    //     counter_cell += w * getbit(bitpos_he + i);
    //     w = w<<1; 
    // }

    // if(!counter_cell)
    // {
    //     return true;
    // }else{
    //     newHashIdx += (k_ + 1 + getbit(bitpos_he ));
    //     // w = 1;
    //     // i = 0;
    //     // for(; i < sizeDeltaIndex_ ; i++)
    //     // {
    //     //     newHashIdx += w * getbit(bitpos + i);
    //     //     w = w<<1; 
    //     // }
    //     return false;
    // }

}

// bool HashExpressor_One::GetInsertSequence(Slice & key, uint64_t hashValue){

    
    
//     uint8_t hash_idx_ = 0;
//     //uint32_t h = XXH3_128bits_withSeed(&key.str[0], key.str.size(), adpHashIdx).low64 % cell_num_;
//     uint32_t h = hashValue% cell_num_;
//     getCell(h, hash_idx_);

//     if(hash_idx_ != 0){
//         return false;
//     }else{
//         return true;
//     }
    
    
// }

void HashModulator::setbit(const int &bitinx){
    //assert(bitinx<=data_.size()*8);
    char *array = &(data_[0]);
    int charIdx = (bitinx) / 8; 
    int b = (bitinx) % 8;
    array[charIdx] |= (1 << b);
}

void HashModulator::clearBit(const int &bitinx){
    //assert(bitinx<=data_.size()<<3);
    char *array = &(data_[0]);
    int charIdx = (bitinx) >> 3; 
    int b = (bitinx) % 8;
    array[charIdx] &= ~(1 << b);
}



bool HashModulator::getbit(const int &bitinx){

    // assert(bitinx<=data_.size()<<3);

    char *array = &(data_[0]);
    int charIdx = (bitinx) >>3; 
    int b = (bitinx) % 8;
    
    
    return (array[charIdx] & (1 << b)) > 0;
}

void HashModulator::getCell(const int &cellinx, uint8_t &hash_idx_)
{   
    uint8_t w = 1;
    
    for(uint8_t i=cell_size_; i>0; i--)
    {
        
        hash_idx_ += w * getbit(cellinx*cell_size_ + i-1);
        w = w<<1; 
    }
}

void HashModulator::getCell_oldIdx(const int &cellinx, uint8_t &oldIdx)
{   
    uint8_t w = 1;
    
    for(uint8_t i = (cell_size_ >> 1); i>0; i--)
    {
        
        oldIdx += w * getbit(cellinx*cell_size_ + i-1);
        w = w<<1; 
    }
}

bool HashModulator::getCell_deltaIdx(const int &cellinx, uint8_t &oldIdx, uint8_t &deltaIdx)
{   
    uint8_t w = 1;
    uint8_t cell_oldIdx=0;
    for(uint8_t i = (cell_size_ >> 1); i>0; i--)
    {
        
        cell_oldIdx += w * getbit(cellinx*cell_size_ + i-1);
        w = w<<1; 
    }
    if(oldIdx != cell_oldIdx){
        return false;
    }else{
        w=1;
        
        for(uint8_t i = cell_size_; i>2; i--)
        {
            
            deltaIdx += w * getbit(cellinx*cell_size_ + i-1);
            w = w<<1; 
        }
        return true;        
    }
}

void HashModulator::setCell(const int &cellinx, uint8_t hash_idx_){

    uint8_t i = cell_size_;
    
    while(hash_idx_ > 0){
        assert(i>0);
        if(hash_idx_ % 2 > 0) setbit(cellinx * cell_size_ + i-1);
        i--;
        hash_idx_ = hash_idx_ >> 1;
    }

}

void HashModulator::setCell_oldIdx(const int &cellinx, uint8_t oldIdx, uint8_t deltaIdx){

    uint8_t i = cell_size_;
    
    while(deltaIdx > 0){
        assert(i>0);
        if(deltaIdx % 2 > 0) setbit(cellinx * cell_size_ + i-1);
        i--;
        deltaIdx = deltaIdx >> 1;
    }
    i=2;
    while(oldIdx > 0){
        assert(i>0);
        if(oldIdx % 2 > 0) setbit(cellinx * cell_size_ + i-1);
        i--;
        oldIdx = oldIdx >> 1;
    }
}

void HashModulator::AssignDefaultHashset_(HABFilterKey *keys_, uint32_t &size_){
    for(int i=0; i<size_; i++)
        for(int j=0; j<k_; j++)
            keys_[i].hash_map_[j] = j + 1;
}

void HashModulator::AssignDefaultHashset_(HABFilterKey &key){
    for(int j=0; j<k_; j++)
            key.hash_map_[j] = j + 1;
}


}

#endif