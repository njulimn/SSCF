#ifndef DYNAMIXHABF_ONE
#define DYNAMIXHABF_ONE

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
#include "./publicStruct.h"
#include "HashMudulator.h"
#define HASH_NUM 3
#define ALLOCATION_RATIO 4 // Allocation ratio between the size of bloom and hashexpressor
using namespace std;

namespace SeesawCF{

const int hashset_size_ = 7;

class CountingBloom{
    
    counting_bloom_header_t *header;
    
    uint32_t *hashes; //by xp: the point to hash function point lists
    
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
    int counting_bloom_add(HABFilterKey &key);
    int counting_bloom_add_default(HABFilterKey &key);
    int counting_bloom_add(HABFilterKey &key, size_t adaptHashIdx, size_t altAdaptIdx);
    //int counting_bloom_remove(HABFilterKey &key);
    int counting_bloom_check(HABFilterKey &key);
    int counting_bloom_check(HABFilterKey &key, size_t &notMapNum, size_t notMapIdx);
    int counting_bloom_check_default(HABFilterKey &key, size_t &notMapNum, size_t &notMapIdx);
    int counting_bloom_check(HABFilterKey &key, size_t newHashIdx);
    int counting_bloom_check_opt(HABFilterKey &key, size_t hashValue);
    size_t GetHaveInsertNum(){return header->count;};
    int counting_bloom_remove(HABFilterKey &key);

    //bitmap operation
    bitmap_t *bitmap_resize(bitmap_t *bitmap, size_t old_size, size_t new_size);
    bitmap_t *new_bitmap(size_t bytes);
    int bitmap_increment_adapt(bitmap_t *bitmap, unsigned int index, long offset);
    int bitmap_increment_adapt_opt(unsigned int index);
    int bitmap_decrement(unsigned int index);
    int bitmap_check_adapt(bitmap_t *bitmap, unsigned int index, long offset);
    int bitmap_check_adapt_opt(unsigned int index);
    void free_bitmap();

    size_t belowBloom = 0;
    
};

class SeesawCF{
    public:
        SeesawCF(double isOneRatio, float bits_per_key_, int pos_count_, int neg_count_, std::vector<Slice *> &neg_keys_, double cbfSpaceRatio);
        SeesawCF(double isOneRatio, float bits_per_key_, int pos_count_, int neg_count_, std::vector<Slice *> &neg_keys_, bool flag);

        ~SeesawCF();
        void Add(Slice &key);
        bool Contain(Slice *key);
        void Deletion(Slice *key);
        bool ContainTest(Slice &key);
        int m=0;
        int getHashExpressorCellNum(){return hashModulator_->getCellNum();}

        uint64_t adaptNum=0;
        size_t onlyBloomNum=0;
        size_t belowBloomNum=0;

        void PrintHE_haveInsertNum(){
            cout<<"haveInsertNum in hashexpressor: "<<hashModulator_->GetHaveInsertNum()<<endl;
        }
        size_t GetHaveInsertNum(){
            return hashModulator_->GetHaveInsertNum();
        }

        void PrintBelowBloom1(){
            std::cout<<"after shuff delete error num : "<<belowBloomNum<<endl;
        }
        void PrintBelowBloom(){
            std::cout<<"after shuff delete error num below : "<<countBloom_->belowBloom<<endl;
        }
        size_t GetCellNum(){
            return hashModulator_->getCellNum();
        }        
        void PrintCounterAboveOne(){
            size_t counterAboveOne=0;
            size_t sumCounter = 0;
            for (size_t i = 0; i < numCounter_cb; i++)
            {
                size_t temp =countBloom_->bitmap_check_adapt_opt(i);
                if(temp>=1){
                    counterAboveOne++; 
                }
            }
            std::cout<<"After delete, the number of cell which not be cleard: "<<counterAboveOne<<endl; 
        }

    private:
        size_t adaptCost = 14;//22
        size_t sumSize;  //the sum space of entire data structure 
        size_t cbf_Size;  //the size ofCounting bloom filter and cost filter 
        size_t hashEx_Size;
        size_t bit_per_counter;
        size_t nHashFunc;
        size_t numCounter_cb;  //the number of the bloom counters
        double maxCost = 0.0;  //the max cost field
        CountingBloom * countBloom_;
        HashModulator * hashModulator_;
        double total_cost = 0.0;
        double *cost_eachCounter; 
        int *cost_eachCounter_tmp;
        size_t altHashNum;
       
     
        uint32_t averCost;
        uint8_t *init_eachCost;
    
};


SeesawCF::SeesawCF(double isOneRatio, float bits_per_key_, int pos_count_, int neg_count_, std::vector<Slice *> &neg_keys_, bool flag){
    
    sumSize = bits_per_key_ * pos_count_;
    cbf_Size = sumSize * 90.0/100;
    hashEx_Size = sumSize - cbf_Size;

    bit_per_counter = 4;
    // nHashFunc = 3;
    nHashFunc = cbf_Size / 4.0/ pos_count_ * 0.69;
    // std::cout<<"The number of hashFunc: "<<nHashFunc<<endl;
    numCounter_cb = cbf_Size / (double)(bit_per_counter + 1.0); 
    // std::cout<<"numCounter_cb: "<<numCounter_cb<<endl;
    countBloom_ = new CountingBloom(numCounter_cb, bit_per_counter, nHashFunc);
    hashModulator_ =  new HashModulator(nHashFunc, hashEx_Size);
    
    cost_eachCounter_tmp = new int[numCounter_cb];
    memset(cost_eachCounter_tmp, 0, sizeof(int)* numCounter_cb);


    //cout<<"sumSize: "<<sumSize<<endl;
    size_t sumAboveAdaptCost = 0;
    size_t aboveAdaptCostCounter = 0;
    size_t isOneNum = pos_count_ * isOneRatio;
    for (size_t i = 0; i < neg_count_; i++)
    {
        total_cost += neg_keys_[i]->cost;
    }
    
    for (int i=0; i<isOneNum; i++)
    {
        for (int j=1; j<=nHashFunc; j++){
            uint32_t h = XXH3_128bits_withSeed(&neg_keys_[i]->str[0], neg_keys_[i]->str.size(), j).low64 % numCounter_cb;
            cost_eachCounter_tmp[h]=1;
        }
        sumAboveAdaptCost += neg_keys_[i]->cost;

    }
    // std::cout<<"sumAboveCost: "<<sumAboveAdaptCost<<" total_cost: "<<total_cost<<" Rate: "<<sumAboveAdaptCost/(double)total_cost<<endl;
    
}

SeesawCF::SeesawCF(double isOneRatio, float bits_per_key_, int pos_count_, int neg_count_, std::vector<Slice *> &neg_keys_, double cbfSpaceRatio){
    
    sumSize = bits_per_key_ * pos_count_;
    cbf_Size = sumSize * cbfSpaceRatio;
    hashEx_Size = sumSize - cbf_Size;

    bit_per_counter = 4;
    nHashFunc = cbf_Size / 4.0/ pos_count_ * 0.69;
    numCounter_cb = cbf_Size / (double)(bit_per_counter + 1.0); 
    countBloom_ = new CountingBloom(numCounter_cb, bit_per_counter, nHashFunc);
    hashModulator_ =  new HashModulator(nHashFunc, hashEx_Size);
    
    cost_eachCounter_tmp = new int[numCounter_cb];
    memset(cost_eachCounter_tmp, 0, sizeof(int)* numCounter_cb);


    //cout<<"sumSize: "<<sumSize<<endl;
    size_t sumAboveAdaptCost = 0;
    size_t aboveAdaptCostCounter = 0;
    size_t isOneNum = pos_count_ * isOneRatio;
    for (size_t i = 0; i < neg_count_; i++)
    {
        total_cost += neg_keys_[i]->cost;
    }
    
    for (int i=0; i<isOneNum; i++)
    {
        for (int j=1; j<=nHashFunc; j++){
            uint32_t h = XXH3_128bits_withSeed(&neg_keys_[i]->str[0], neg_keys_[i]->str.size(), j).low64 % numCounter_cb;
            cost_eachCounter_tmp[h]=1;
        }
        sumAboveAdaptCost += neg_keys_[i]->cost;

    }
    std::cout<<"sumAboveCost: "<<sumAboveAdaptCost<<" total_cost: "<<total_cost<<" Rate: "<<sumAboveAdaptCost/(double)total_cost<<endl;
    
}


SeesawCF::~SeesawCF(){
    delete countBloom_;
    delete hashModulator_;
    delete []cost_eachCounter_tmp;
}



// SeesawCF::SeesawCF(double isOneRatio, float bits_per_key_, int pos_count_, int neg_count_, std::vector<Slice *> &neg_keys_){
    
//     sumSize = bits_per_key_ * pos_count_;
//     cbf_Size = sumSize * 90.0/100;
//     hashEx_Size = (sumSize - cbf_Size);

//     bit_per_counter = 4;
//     nHashFunc = 3;
//     numCounter_cb = cbf_Size / (double)(bit_per_counter + 1.0); 
//     //std::cout<<"numCounter_cb: "<<numCounter_cb<<endl;
//     countBloom_ = new CountingBloom(numCounter_cb, bit_per_counter, nHashFunc);
//     hashModulator_ =  new HashModulator(nHashFunc, hashEx_Size);

//     cost_eachCounter = new double[numCounter_cb];
//     memset(cost_eachCounter, 0, sizeof(cost_eachCounter));
//     //cout<<"numcounter"<<numCounter_cb<<endl;
//     cost_eachCounter_tmp = new int[numCounter_cb];
//     memset(cost_eachCounter_tmp, 0, sizeof(cost_eachCounter_tmp));
//     init_eachCost = new uint8_t[nHashFunc];
//     memset(init_eachCost, 0, sizeof(init_eachCost));

//     //cout<<"sumSize: "<<sumSize<<endl;
//     size_t sumAboveAdaptCost = 0;
//     size_t aboveAdaptCostCounter = 0;
//     size_t isOneNum = pos_count_ * isOneRatio;
//     for (size_t i = 0; i < neg_count_; i++)
//     {
//         total_cost += neg_keys_[i]->cost;
//     }
    
//     for (int i=0; i<isOneNum; i++)
//     {
//         for (int j=1; j<=3; j++){
//             uint32_t h = XXH3_128bits_withSeed(&neg_keys_[i]->str[0], neg_keys_[i]->str.size(), j).low64 % numCounter_cb;
//             cost_eachCounter_tmp[h]=1;
//         }
//         // if( neg_keys_[i]->cost > adaptCost)
//         // {
//         //     aboveAdaptCostCounter++;
//         //     sumAboveAdaptCost += neg_keys_[i]->cost;
//         // }
//         sumAboveAdaptCost += neg_keys_[i]->cost;

//     }
//     std::cout<<"sumAboveCost: "<<sumAboveAdaptCost<<" total_cost: "<<total_cost<<" Rate: "<<sumAboveAdaptCost/(double)total_cost<<endl;
    
//     // size_t idOne = 0;
//     // for (int i=0; i< numCounter_cb; i++)
//     // {
//     //     // if(cost_eachCounter[i] >= 5){
//     //     //     above1000 ++;
//     //     // }

//     //     if (cost_eachCounter[i] < adaptCost)
//     //     {
//     //         cost_eachCounter_tmp[i] = 0;            
//     //     }else 
//     //     {
//     //         cost_eachCounter_tmp[i] = 1;
//     //         idOne++;
//     //     }
//     // } 
//     // cout<<"isOne Rate: "<<idOne/(double)numCounter_cb*100<<"%"<<endl;

// }



void SeesawCF::Add(Slice &key){

    uint32_t h;
    size_t bitpos_he, cellpos;
    for (size_t oldHashIdx=1; oldHashIdx <= nHashFunc; oldHashIdx++)
    {
        h = XXH3_128bits_withSeed(&key.str[0], key.str.size(), oldHashIdx).low64 % numCounter_cb;
     
        size_t h1,  cell_newHashIdx = 0, cell_counter = 0;

        if(cost_eachCounter_tmp[h])
        {
            if (hashModulator_->GetInsertSequence(key,  cell_newHashIdx, cell_counter, cellpos, bitpos_he)) //if it is empty
            {
                for (size_t newHashIdx = nHashFunc + 1; newHashIdx <= nHashFunc + 2; newHashIdx++)
                {
                        h1 = XXH3_128bits_withSeed(&key.str[0], key.str.size(), newHashIdx).low64 % numCounter_cb;
                        if (cost_eachCounter_tmp[h1]==0)  
                        {   
                            //hashModulator_->AddCounterByOne(key,bitpos_he, 1); 
                            hashModulator_->Add(key, cellpos, bitpos_he, newHashIdx);   
                            countBloom_->bitmap_increment_adapt_opt(h1);

                            for(oldHashIdx+=1 ; oldHashIdx <= nHashFunc; oldHashIdx++)
                            {         
                                h1 = XXH3_128bits_withSeed(&key.str[0], key.str.size(), oldHashIdx).low64 % numCounter_cb;
                                countBloom_->bitmap_increment_adapt_opt(h1);
                            }   
                            return;
                        } 
                }
                hashModulator_->AddCounterByOne(key,cellpos, bitpos_he, 1);
                countBloom_->bitmap_increment_adapt_opt(h);  
                for(oldHashIdx+=1 ; oldHashIdx <= nHashFunc; oldHashIdx++)
                {         
                    h1 = XXH3_128bits_withSeed(&key.str[0], key.str.size(), oldHashIdx).low64 % numCounter_cb;
                    countBloom_->bitmap_increment_adapt_opt(h1);
                }    
                return;             
            }else  //if the cell is not empty
            {
                //TODO:
                hashModulator_->AddCounterByOne(key, cellpos, bitpos_he, cell_counter+1);  


                h1 = XXH3_128bits_withSeed(&key.str[0], key.str.size(), cell_newHashIdx).low64 % numCounter_cb;
                if (cost_eachCounter_tmp[h1]==0)
                {
                    countBloom_->bitmap_increment_adapt_opt(h1);
                    
                    for(oldHashIdx+=1 ; oldHashIdx <= nHashFunc; oldHashIdx++)
                    {         
                        h1 = XXH3_128bits_withSeed(&key.str[0], key.str.size(), oldHashIdx).low64 % numCounter_cb;
                        countBloom_->bitmap_increment_adapt_opt(h1);
                    } 
                    return;
                }else
                {
                    countBloom_->bitmap_increment_adapt_opt(h);
                    //cell_counter+1
                    for(oldHashIdx+=1 ; oldHashIdx <= nHashFunc; oldHashIdx++)
                    {         
                        h1 = XXH3_128bits_withSeed(&key.str[0], key.str.size(), oldHashIdx).low64 % numCounter_cb;
                        countBloom_->bitmap_increment_adapt_opt(h1);
                    } 
                    return;
                }

            }

        }else{

            countBloom_->bitmap_increment_adapt_opt(h);
        }     
    } 

    return ;
}


bool SeesawCF::Contain(Slice *key)
{
    size_t notMapNum = 0;
    uint64_t cellPos;
    
    for (uint8_t i=1 ; i<=nHashFunc; i++){
        cellPos = XXH3_128bits_withSeed(&key->str[0], key->str.size(), i).low64 % numCounter_cb;
        if (!(countBloom_->bitmap_check_adapt_opt(cellPos))) 
        {
            notMapNum ++;
            if (cost_eachCounter_tmp[cellPos] !=1 || notMapNum >=2)
            {
                return false;
            }
        }
    }
    
    if (notMapNum){
        uint8_t deltaIdx = 1;
        if (!hashModulator_->Query(*key, deltaIdx)){
            return false;
        }else
        {
            cellPos = XXH3_128bits_withSeed(&key->str[0], key->str.size(), deltaIdx + nHashFunc).low64 % numCounter_cb;
            if (countBloom_->bitmap_check_adapt_opt(cellPos)  )//&& cost_eachCounter_tmp[bitpos] == 0
            {
                return true;
            }
        }
        return false;
        
    }else{
        return true;
    }
    
}

void SeesawCF::Deletion(Slice *key)
{
    size_t cell_counter , newIdx = nHashFunc+1, cellpos, hOld, initHashIdx ;
    hashModulator_->QueryIsEmpty(*key, cell_counter, newIdx, cellpos);
    if (cell_counter == 0)
    {
        for (initHashIdx = 1; initHashIdx <= nHashFunc; initHashIdx++)
        {
            //hOld = XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb;
            countBloom_->bitmap_decrement( XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb ) ;
        }
        return;
    }else
    {
        for (initHashIdx = 1; initHashIdx <= nHashFunc; initHashIdx++)
        {
            hOld = XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb;

            if (!cost_eachCounter_tmp[hOld])
            {
        
                countBloom_->bitmap_decrement(hOld) ;

            }else
            {
                hashModulator_->Clear(cell_counter-1, cellpos); 
                size_t hNew = XXH3_128bits_withSeed(&key->str[0], key->str.size(), newIdx).low64 % numCounter_cb;
                if( cost_eachCounter_tmp[hNew])
                {  
                    countBloom_->bitmap_decrement(hOld);
                    for (initHashIdx += 1; initHashIdx <= nHashFunc; initHashIdx++)
                    {
                        //hOld = XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb;
                        countBloom_->bitmap_decrement( XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb ) ;
                    }    
                    return;   
                    
                }else 
                {
                    countBloom_->bitmap_decrement(hNew);
                    for (initHashIdx += 1; initHashIdx <= nHashFunc; initHashIdx++)
                    {
                        //hOld = XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb;
                        countBloom_->bitmap_decrement( XXH3_128bits_withSeed(&key->str[0], key->str.size(), initHashIdx).low64 % numCounter_cb ) ;
                    }
                    return;
                }
            }
        }

        
    }
      

    return;

    




}

CountingBloom::CountingBloom(size_t _numCounter,size_t _bitPerCounter, size_t _nfunc)
{
    
    numCounter=_numCounter;
    //std::cout<<"numCounter: "<<numCounter<<endl;
    bitPerCounter = _bitPerCounter;
    // size_t bloomSize = (numCounter * bitPerCounter + 7) / 8;
    bitmap = NULL;
    offset = 0 + sizeof(counting_bloom_header_t); //by xp: because all filters are sorted
    this->nfuncs = _nfunc; //by xp: return the smallest integral value that is not less than x
    /* rounding-up integer divide by 2 of bloom->size */
    this->num_bytes = (_numCounter * _bitPerCounter + 7) / 8 + sizeof(counting_bloom_header_t); //by xp: the size of (the header + the bitmap) of this filter

    this->hashes = (uint32_t *)calloc(this->nfuncs, sizeof(uint32_t));

    this->bitmap = new_bitmap(this->num_bytes);
    // for(int i =0; i<num_bytes;i++){
    //     if(bitmap->array[i]!=0)
    //     {
    //         cout<<"yes"<<endl;
    //     }
    // }
    this->header = (counting_bloom_header_t *)(this->bitmap->array);

}


CountingBloom::~CountingBloom()
{
        // free(this->hashes);
        // this->hashes = NULL;
        // free(this->bitmap);
        // this->bitmap=NULL;
        // free(this->header);
        // this->header=NULL;
        if(hashes != NULL){
            free(hashes);
            // hashes = NULL;
        }else{
            std::cout<<"error in free\n";
        }
        // if(bitmap->array != NULL){
        //     free(bitmap->array);
        //     bitmap->array = NULL;
        // }else{
        //     std::cout<<"error in free\n";
        // }
        if(bitmap->array != NULL){
            // std::cout<<"aa\n"<<endl;
            free(bitmap->array);
            // bitmap->array=NULL;
            // std::cout<<"bb\n"<<endl;
        }else{
            std::cout<<"error in free\n";
        }

        if(bitmap != NULL){
            // std::cout<<"aa\n"<<endl;
            free(bitmap);
            // std::cout<<"bb\n"<<endl;
        }
        bitmap=NULL;
        // if(header!=NULL){
        //     std::cout<<"aa\n"<<endl;
        //     free(header);
        //     std::cout<<"bb\n"<<endl;
        // }else{
        //     std::cout<<"error in free\n";
        // }
}

int CountingBloom::counting_bloom_add(HABFilterKey &key){
    
    
    for(uint8_t j=0; j<nfuncs; j++){         
        unsigned int bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), key.hash_map_[j]).low64 % numCounter;

        bitmap_increment_adapt(this->bitmap, bitpos, offset);

    }
    this->header->count++;
    return 0;
}

int CountingBloom::counting_bloom_add_default(HABFilterKey &key){
    
    
    for(uint8_t j=0; j<nfuncs; j++){         
        unsigned int bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), j+1).low64 % numCounter;

        bitmap_increment_adapt(this->bitmap, bitpos, offset);

    }
    this->header->count++;
    return 0;
}

int CountingBloom::counting_bloom_add(HABFilterKey &key, size_t adaptHashIdx, size_t altAdaptIdx){
    
    uint8_t j = adaptHashIdx;
    unsigned int bitpos;
    // for(; j<adaptHashIdx; j++){         
    //     bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), j).low64 % numCounter;

    //     bitmap_increment_adapt(this->bitmap, bitpos, offset);

    // }
    
    bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), altAdaptIdx).low64 % numCounter;
    bitmap_increment_adapt(this->bitmap, bitpos, offset);  

    j++;

    for(; j<=nfuncs; j++){         
        bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), j).low64 % numCounter;
        bitmap_increment_adapt(this->bitmap, bitpos, offset);
    } 

    //this->header->count++;
    return 0;
}

int CountingBloom::counting_bloom_check(HABFilterKey &key, size_t &notMapNum, size_t notMapIdx){


    for(int i=0 ; i<nfuncs; i++){
        uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), key.hash_map_[i]).low64 % numCounter;
        if (!(bitmap_check_adapt(this->bitmap, bitpos, offset))) {
            notMapNum++;
            notMapIdx = i;
            if(notMapNum >=2)
                return 0;
        }
    }
    if(notMapNum) return 0;
    return 1;
}



int CountingBloom::counting_bloom_check_default(HABFilterKey &key, size_t &notMapNum, size_t &notMapIdx){


    for(int i=0 ; i<nfuncs; i++){
        uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), i + 1).low64 % numCounter;
        if (!(bitmap_check_adapt(this->bitmap, bitpos, offset))) {
            notMapNum++;
            notMapIdx = i + 1;
            if(notMapNum >=2)
                return 0;
        }
    }
    if(notMapNum) return 0;
    return 1;
}

int CountingBloom::counting_bloom_check(HABFilterKey &key, size_t newHashIdx){

    uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), newHashIdx).low64 % numCounter;
    if (!(bitmap_check_adapt(this->bitmap, bitpos, offset))) {
        return 0;
    }
    return 1;
  
}


int CountingBloom::counting_bloom_check_opt(HABFilterKey &key, size_t hashValue){

    //uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), newHashIdx).low64 % numCounter;
    if (!(bitmap_check_adapt(this->bitmap, hashValue, offset))) {
        return 0;
    }
    return 1;
  
}

int CountingBloom::counting_bloom_check(HABFilterKey &key){


    for(int i=0 ; i<nfuncs; i++){
        uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), key.hash_map_[i]).low64 % numCounter;
        if (!(bitmap_check_adapt(this->bitmap, bitpos, offset))) {
            return 0;
        }
    }
    return 1;
}


// int CountingBloom::counting_bloom_remove(HABFilterKey &key)
// {
//     if(!counting_bloom_check(key)) return 0;
    
//     for (size_t i = 0; i < nfuncs; i++) {
//         uint32_t bitpos = XXH3_128bits_withSeed(&key.data_->str[0], key.data_->str.size(), key.hash_map_[i]).low64 % numCounter;
//         bitmap_decrement(this->bitmap, bitpos, this->offset);
//     }
//     this->header->count--;
//     return 1;
// }




void CountingBloom::free_bitmap()
{

    free(bitmap->array);
    free(bitmap);
}



bitmap_t *CountingBloom::bitmap_resize(bitmap_t *bitmap, size_t old_size, size_t new_size)
{

    if (bitmap->array != NULL) {
#if __linux
        //bitmap->array = mremap(bitmap->array, old_size, new_size, MREMAP_MAYMOVE);
   	  bitmap->array = (char *)realloc(bitmap->array, new_size); 
      
        if (bitmap->array == NULL) {
            perror("Error resizing realloc");
            free_bitmap();
           // close(fd);
            return NULL;
        }
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
        std::cout<<"1"<<std::endl;
        return NULL;
    }
   
    bitmap->bytes = bytes;
//    bitmap->fd = fd; //by xp: file handle for bloom_file
    bitmap->array = NULL;
    
    if ((bitmap = bitmap_resize(bitmap, 0, bytes)) == NULL) {
        
        return NULL;
    }

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
            fprintf(stderr, "Error, 2 bit int Overflow\n");
            return -1;          
        }
        uint8_t addNum=1ULL << bitPerCounter*j;
        n=n+addNum;
        bitmap->array[access] = n;
    }
    return 0;
}

int CountingBloom::bitmap_increment_adapt_opt(unsigned int index)
{

    long access = (index>>1) + offset;
    uint8_t temp;
    uint8_t n = bitmap->array[access];
    if ((index & 0x01) != 0) { //by xp: counter increasement
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

    return 0;
}

/* increments the four bit counter */
int CountingBloom::bitmap_decrement(unsigned int index)
{


    // long access = (index >>1) + offset;
    // uint8_t temp, bitOffset = ((index & 0x01)<<2);

    // uint8_t n = bitmap->array[access];
    // if((n & (0xf0>>bitOffset))==0x00)
    // {
    //     fprintf(stderr, "Error, Decrementing zero\n");
    //     return -1;        
    // }
    
    // if (bitOffset != 0) {

    //     bitmap->array[access] = n & 0xf0 | (n & (0x0f)) - 0x01;
    // } else {

    //     bitmap->array[access] = n & 0x0f | (n & 0xf0) - 0x10;
    // }

    // return 0;




    long access = (index >>1) + offset;
    uint8_t temp, n = bitmap->array[access];
    
    if ((index & 0x01)) {
        temp = (n & 0x0f);
        n = n & 0xf0 | (n & 0x0f) - 0x01;
    } else {
        temp = (n & 0xf0) ;
        n = n & 0x0f | (n & 0xf0) - 0x10;
    }
    
    if (temp == 0x00) {
        // belowBloom ++;
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

/* decrements the four bit counter */
int CountingBloom::bitmap_check_adapt_opt(unsigned int index)
{
    //if(bitPerCounter==4){
        long access = (index >> 1) + offset;
        if ((index &0x1) != 0 ) {
           return bitmap->array[access] & 0x0f;
        } else {
           return bitmap->array[access] & 0xf0;
        }
    // }else if(bitPerCounter==2){
    //     long access = (index*bitPerCounter)/8 + offset;
    //     size_t j=index%4;
    //     uint8_t temp;
    //     uint8_t n = bitmap->array[access];
    //     return temp=(n>>(j*2))& 0x03;        
    //  }
    return 0;

}








}
#endif
