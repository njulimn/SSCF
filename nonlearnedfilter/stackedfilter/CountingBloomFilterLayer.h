#pragma once

#include <cstddef>
#include <vector>
#include "InterfaceAMQ.h"
#include "InterfaceElement.h"

template<typename element_type>
class CountingBloomFilterLayer : public InterfaceAMQ<CountingBloomFilterLayer, element_type> {
private:
public:
    std::vector<int> filter_;
    unsigned int num_hashes_;
    size_t numCounter;
    int seed_;

    size_t bit_Per_Counter;
    CountingBloomFilterLayer(size_t _numCounter, int num_hashes, int seed, size_t _bit_Per_Counter);

    CountingBloomFilterLayer() : CountingBloomFilterLayer(0, 0, 0, 0) {};

    ~CountingBloomFilterLayer();

    size_t getHash1(const element_type element);

    size_t getHash2(const element_type element);

    size_t getNthHash(const size_t hash1, const size_t hash2, const size_t hash_num);

    bool LookupElement(const element_type element);

    void InsertElement(const element_type element);

    void DeleteElement(const element_type element); 

    double GetLoadFactor();

    static size_t SizeFunctionImplementation(double fpr,
                                             size_t num_expected_elements) {
        size_t num_hashes = std::max<size_t>(ceil(-std::log2(fpr)), 1);
        size_t theoretical_size = -1. / (std::pow(1. - std::pow(fpr, 1. / num_hashes),
                                                  1. / (num_hashes * num_expected_elements)) - 1.);
        return theoretical_size; // Adjusted to keep load factors <=.5
    }
};


template<typename element_type>
CountingBloomFilterLayer<element_type>::CountingBloomFilterLayer(size_t _numCounter, int num_hashes,
                                                 int seed, size_t _bit_Per_Counter) {
    seed_ = seed;
    bit_Per_Counter = _bit_Per_Counter;
    num_hashes_ = num_hashes;
    //cout<<"numHash: "<<num_hashes<<endl;
    numCounter = _numCounter;
    this->total_size_ = numCounter * bit_Per_Counter ;
   
    filter_.resize(numCounter, 0);
    // if (filter_size <= 2) {
    //     filter_size_ = 2;
    //     filter_.resize(filter_size_, true);
    // }
}

template<typename element_type>
CountingBloomFilterLayer<element_type>::~CountingBloomFilterLayer() {}

_GLIBCXX17_INLINE
template<typename element_type>
bool CountingBloomFilterLayer<element_type>::LookupElement(const element_type element) {
    // size_t hash1 = getHash1(element);
    // size_t hash2 = getHash2(element);
    // for (unsigned int i = 0; i < num_hashes_; i++) {
    //     if (filter_[getNthHash(hash1, hash2, i)] == 0 ) {
    //         return false;
    //     }
    // }
    for (unsigned int i = 0; i < num_hashes_; i++) {
        unsigned int bitpos = XXH3_128bits_withSeed(element.get_value(), element.size(), i+1).low64 % numCounter;
        if(filter_[bitpos] ==0){
            return false;
        }
            
    }
    return true;
}

_GLIBCXX17_INLINE
template<typename element_type>
void CountingBloomFilterLayer<element_type>::InsertElement(const element_type element) {
    this->num_elements_++;
    for (unsigned int i = 0; i < num_hashes_; i++) {
        unsigned int bitpos = XXH3_128bits_withSeed(element.get_value(), element.size(), i+1).low64 % numCounter;
        if(filter_[bitpos] < ((1<<bit_Per_Counter) - 1 )){
            filter_[bitpos] += 1;
        }else{
            printf("error! overflow! \n");
        }
            
    }
}

_GLIBCXX17_INLINE
template<typename element_type>
void CountingBloomFilterLayer<element_type>::DeleteElement(const element_type element) {
    this->num_elements_--;
    for (unsigned int i = 0; i < num_hashes_; i++) {
        unsigned int bitpos = XXH3_128bits_withSeed(element.get_value(), element.size(), i+1).low64 % numCounter;
        if(filter_[bitpos] > 0){
            filter_[bitpos] -= 1;
        }else{
            printf("error! overflow! \n");
        }
            
    }
}

_GLIBCXX17_INLINE
template<typename element_type>
size_t CountingBloomFilterLayer<element_type>::getHash1(const element_type x) {
    return CityHash64WithSeed(x.get_value(), x.size(), 123456789 + seed_);
}

_GLIBCXX17_INLINE
template<typename element_type>
size_t CountingBloomFilterLayer<element_type>::getHash2(const element_type x) {
    return CityHash64WithSeed(x.get_value(), x.size(), 987654321 + seed_);
}

_GLIBCXX17_INLINE
template<typename element_type>
// Double hashing strategy recommended in Mitzenmacher Paper.
size_t CountingBloomFilterLayer<element_type>::getNthHash(const size_t hash1, const size_t hash2,
                                                  const size_t hash_num) {
    return (hash1 + hash_num * hash2 + hash_num * hash_num * hash_num) % numCounter;
}

template<typename element_type>
double CountingBloomFilterLayer<element_type>::GetLoadFactor() {
    size_t load = 0;
    for (size_t i = 0; i < numCounter; i++)
        if (filter_[i] == true) load++;
    return (double) load / (double) numCounter;
}

template
class CountingBloomFilterLayer<IntElement>;

template
class CountingBloomFilterLayer<StringElement>;

template
class CountingBloomFilterLayer<BigIntElement>;
