#ifndef WEIGHTCOUNTINGBLOOM
#define WEIGHTCOUNTINGBLOOM

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
#include "../util/xxhash/xxhash.h"
#include <algorithm>

#define HASH_NUM 3
#define ALLOCATION_RATIO 4 // Allocation ratio between the size of bloom and hashexpressor
using namespace std;



namespace weightedCountingBloom{

const int hashset_size_ = 7;

double getEr(double cost_e, double xe, double s){
    return (1-xe)*cost_e / (xe*cost_e + s);
}

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

struct WCBFKey{
    Slice * data_;
    uint8_t k_;
};


bool compare_func(WCBFKey &key1, WCBFKey &key2){
    return key1.data_->cost > key2.data_->cost;
}

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
    int counting_bloom_add(Slice &key, size_t hashNum);
    int counting_bloom_remove(Slice &key,size_t ke);
    int counting_bloom_check(Slice &key, size_t k);

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
        cout<<"After delete, the number of item in WCBF: "<<this->header->count<<endl;
    }
    
};

//能够根据要求每个counter的比特数bitPerCounter进行初始化
CountingBloom::CountingBloom(size_t _numCounter,size_t _bitPerCounter, size_t nfunc)
{
    
    numCounter=_numCounter;
    // cout<<"numCounter: "<<numCounter<<endl;
    bitPerCounter = _bitPerCounter;
    //整个bloom可使用的空间sumSize, byte
    size_t bloomSize = (numCounter*bitPerCounter+7)/8;
    bitmap = NULL;
    offset = 0 + sizeof(counting_bloom_header_t); //by xp: because all filters are sorted
    this->nfuncs = nfunc; //by xp: return the smallest integral value that is not less than x
    /* rounding-up integer divide by 2 of bloom->size */
    this->num_bytes = (numCounter*bitPerCounter+7)/8 + sizeof(counting_bloom_header_t); //by xp: the size of (the header + the bitmap) of this filter
    this->hashes = (uint32_t *)calloc(this->nfuncs, sizeof(uint32_t));
    this->bitmap = new_bitmap(this->num_bytes);
    this->header = (counting_bloom_header_t *)(this->bitmap->array);
    

}
int CountingBloom::counting_bloom_add(Slice &key, size_t hashNum){
    
    
    for(uint8_t j=1; j<=hashNum; j++){         
        unsigned int bitpos = XXH3_128bits_withSeed(&key.str[0], key.str.size(), j).low64 % numCounter;
        bitmap_increment_adapt(this->bitmap, bitpos, offset);
    }
    this->header->count++;
    return 0;
}



int CountingBloom::counting_bloom_check(Slice &key, size_t k){

    for(int i=1 ; i<=k; i++){
        uint32_t bitpos = XXH3_128bits_withSeed(&key.str[0], key.str.size(), i).low64 % numCounter;
        if (!(bitmap_check_adapt(this->bitmap, bitpos, offset))) {
            return 0;
        }
    }
    return 1;
}

int CountingBloom::counting_bloom_remove(Slice &key, size_t ke)
{

    for (size_t i = 1; i <= ke; i++) {
        uint32_t bitpos = XXH3_128bits_withSeed(&key.str[0], key.str.size(), i).low64 % numCounter;
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
        fprintf(stderr, "Error, Decrementing zero!\n");
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
    if ((msync(bitmap->array, bitmap->bytes, 0) < 0)) {
        perror("Error, flushing bitmap to disk");
        return -1;
    } else {
        return 0;
    }
}








CountingBloom::~CountingBloom()
{

        // free(this->hashes);
        // this->hashes = NULL;
        // free(this->bitmap);
        // this->bitmap = NULL;
        // free(this->header);
        // this->header = NULL;
        if(hashes != NULL){
            free(hashes);
            hashes = NULL;
        }else{
            std::cout<<"error in free\n";
        }
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


class WCBF{
    public:
        WCBF(float bits_per_key, size_t bits_Per_Counter, std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_ );
        ~WCBF();
        void AddAll(std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_); 
        bool Contain(Slice &key);
        double testFPR();
        void AddAll_latency(WCBFKey* wbf_pos_keys_, WCBFKey* wbf_neg_keys_ , size_t pos_num, size_t neg_num);
        void AddAll(std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_, double wfprNum[], double minWFPR[], double maxWFPR[], int n);
        void DeleteAll(std::vector<Slice *> &pos_keys_);
        void PrintHaveInsertNum(){
            countingBloom_->printHaveInsertNum();
        }
    private:
        float bits_per_key_;
        CountingBloom *countingBloom_;
        std::vector<Slice *> cost_list_;
        set<Slice> setCost_list;
        std::vector<Slice *> no_cost_list_;
        double s;
        double s2;
        uint32_t ke;
        uint32_t pos_n;
        uint32_t N;
        size_t bits_Per_Counter_;
        size_t sumSize;
        size_t cost_storing_num;
        size_t numCounter;
};

WCBF::WCBF(float bits_per_key, size_t bits_Per_Counter, std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_) : bits_per_key_(bits_per_key){
    bits_Per_Counter_ = bits_Per_Counter;
    sumSize = bits_per_key * pos_keys_.size();
    // cout<<"sumSize:"<<sumSize<<endl;
    //double wcbf_cost_radio = 0.0001;
    double wcbf_cost_radio = 0.00005;
    cost_storing_num = neg_keys_.size() * wcbf_cost_radio;

    size_t cost_bits = cost_storing_num * 8 * sizeof(Slice); 
    size_t countingbloom_bits = sumSize - cost_bits;

    if(countingbloom_bits < 0){
        cout<<"error , have no space for countingbloom1\n";
    }
    numCounter = countingbloom_bits / bits_Per_Counter;
    countingBloom_ = new CountingBloom(numCounter, bits_Per_Counter_, 0);
    
}

WCBF::~WCBF(){
    // std::cout<<"DELETE!"<<endl;
    delete countingBloom_;
    
}



void WCBF::AddAll(std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_, double wfprNum[], double minWFPR[], double maxWFPR[], int n){


    pos_n = pos_keys_.size();
    uint32_t neg_n = neg_keys_.size();
    WCBFKey* wbf_pos_keys_ = new WCBFKey[pos_n];
    WCBFKey* wbf_neg_keys_ = new WCBFKey[neg_n];
    

    for(int i=0; i<pos_n; ++i)
        wbf_pos_keys_[i].data_ = pos_keys_[i];
    for(int i=0; i<neg_n; ++i)
        wbf_neg_keys_[i].data_ = neg_keys_[i];

    std::sort(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_func);

    cost_list_.resize(cost_storing_num);
    for(int i=0; i<cost_storing_num; ++i)
        cost_list_[i] = wbf_neg_keys_[i].data_;
    for(int i=0; i<cost_storing_num; ++i)
    {
        setCost_list.insert(*wbf_neg_keys_[i].data_);
    }
    // for(int i=cost_storing_num; i<neg_n; ++i)
    //     no_cost_list_.push_back(wbf_neg_keys_[i].data_);

    N = pos_n + neg_n;
    s = 0;

    for(int i=0; i< neg_n; i++)
        s += (1.0-(double)pos_n / (double)N) * (neg_keys_[i]->cost + 1.0);
    for(int i=0; i< pos_n; i++)
        s += (1.0-(double)pos_n / (double)N) * 1.0;

    // std::cout << s << std::endl;

    s2 = 0;
    for(int i=0; i< neg_n; i++){
        double eri = getEr(neg_keys_[i]->cost + 1.0, (double)pos_n / (double)N, s);
        s2 += log(eri) / (double)N;   /////////////////////有问题
    }
    for(int i=0; i< pos_n; i++){
        double eri = getEr(1.0, (double)pos_n / (double)N, s);
        s2 +=  log(eri) / (double)N;   /////////////////有问题
    }

    double Ere_Pos = getEr(1, (double)pos_n / (double)N, s);

    ke = ceil((double)numCounter / (double)pos_n *log(2) + (log(Ere_Pos) - s2)/log(2));
    
    //uint32_t op_ke = uint32_t((double)numCounter / (double)pos_n * log(2));
    //if(ke < op_ke) ke = op_ke;
    //std::cout << "op_ke: "<<op_ke << std::endl;



    size_t count = 0;
    size_t bits_per_item_ite = 50;
    for(auto v:pos_keys_){
        countingBloom_->counting_bloom_add(*v, ke);
        count++;
        double true_bits_per_item = sumSize / (double)count; 
        if( bits_per_item_ite >=10 && true_bits_per_item < bits_per_item_ite)
        {
        
            double fpr_c_temp = 0;
            double total_cost_temp = 0.0;
            for(int i=0; i<neg_keys_.size(); i++){
                if(Contain(*neg_keys_[i]))
                {
                    fpr_c_temp += neg_keys_[i]->cost;
                }
                total_cost_temp += neg_keys_[i]->cost;
            } 
            double tempV = (double)(fpr_c_temp*100) / total_cost_temp;
            wfprNum[bits_per_item_ite] += tempV;
            
            if (tempV > maxWFPR[bits_per_item_ite]){
                maxWFPR[bits_per_item_ite] = tempV;
            }
            if(tempV < minWFPR[bits_per_item_ite]){
                //cout<<"tempv: "<<tempV<<endl;
                minWFPR[bits_per_item_ite] = tempV;
            }
            //wfprNum[bits_per_item_ite] +=(double)(fpr_c_temp*100) / total_cost_temp; ;
            //std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  Weighted FPR:  "<< (double)(fpr_c_temp*100) / total_cost_temp << "%" <<endl;    
            bits_per_item_ite -=2;          
        }
    }

    delete [] wbf_pos_keys_;
    delete [] wbf_neg_keys_;
}



void WCBF::AddAll_latency(WCBFKey* wbf_pos_keys_, WCBFKey* wbf_neg_keys_ , size_t pos_num, size_t neg_num){

    pos_n = pos_num;
    uint32_t neg_n = neg_num;

    cost_list_.resize(cost_storing_num);
    for(int i=0; i<cost_storing_num; ++i)
        cost_list_[i] = wbf_neg_keys_[i].data_;
        
    for(int i=0; i<cost_storing_num; ++i)
    {
        setCost_list.insert(*wbf_neg_keys_[i].data_);
    }

    // for(int i=cost_storing_num; i<neg_n; ++i)
    //     no_cost_list_.push_back(wbf_neg_keys_[i].data_);

    N = pos_n + neg_n;
    s = 0;

    for(int i=0; i< neg_n; i++)
        s += (1.0-(double)pos_n / (double)N) * (wbf_neg_keys_[i].data_->cost + 1.0);
    for(int i=0; i< pos_n; i++)
        s += (1.0-(double)pos_n / (double)N) * 1.0;

    // std::cout << s << std::endl;

    s2 = 0;
    for(int i=0; i< neg_n; i++){
        double eri = getEr(wbf_neg_keys_[i].data_->cost + 1.0, (double)pos_n / (double)N, s);
        s2 += log(eri) / (double)N;  
    }
    for(int i=0; i< pos_n; i++){
        double eri = getEr(1.0, (double)pos_n / (double)N, s);
        s2 +=  log(eri) / (double)N;   
    }

    double Ere_Pos = getEr(1, (double)pos_n / (double)N, s);

    ke = ceil((double)numCounter / (double)pos_n *log(2) + (log(Ere_Pos) - s2)/log(2));
    

    for (size_t i = 0; i < pos_num; i++)
    {
        countingBloom_->counting_bloom_add(*wbf_pos_keys_[i].data_, ke);
    }

}

void WCBF::DeleteAll(std::vector<Slice *> &pos_keys_){
    size_t posNum = pos_keys_.size();
    for (size_t i = 0; i < posNum; i++)
    {
        countingBloom_->counting_bloom_remove(*pos_keys_[i], ke);
    }
}



void WCBF::AddAll(std::vector<Slice *> &pos_keys_, std::vector<Slice *> &neg_keys_){


    pos_n = pos_keys_.size();
    uint32_t neg_n = neg_keys_.size();
    WCBFKey* wbf_pos_keys_ = new WCBFKey[pos_n];
    WCBFKey* wbf_neg_keys_ = new WCBFKey[neg_n];
    

    for(int i=0; i<pos_n; ++i)
        wbf_pos_keys_[i].data_ = pos_keys_[i];
    for(int i=0; i<neg_n; ++i)
        wbf_neg_keys_[i].data_ = neg_keys_[i];

    std::sort(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_func);

    cost_list_.resize(cost_storing_num);
    for(int i=0; i<cost_storing_num; ++i)
        cost_list_[i] = wbf_neg_keys_[i].data_;

    for(int i=cost_storing_num; i<neg_n; ++i)
        no_cost_list_.push_back(wbf_neg_keys_[i].data_);

    N = pos_n + neg_n;
    s = 0;

    for(int i=0; i< neg_n; i++)
        s += (1.0-(double)pos_n / (double)N) * (neg_keys_[i]->cost + 1.0);
    for(int i=0; i< pos_n; i++)
        s += (1.0-(double)pos_n / (double)N) * 1.0;

    // std::cout << s << std::endl;

    s2 = 0;
    for(int i=0; i< neg_n; i++){
        double eri = getEr(neg_keys_[i]->cost + 1.0, (double)pos_n / (double)N, s);
        s2 += log(eri) / (double)N;   /////////////////////有问题
    }
    for(int i=0; i< pos_n; i++){
        double eri = getEr(1.0, (double)pos_n / (double)N, s);
        s2 +=  log(eri) / (double)N;   /////////////////有问题
    }

    double Ere_Pos = getEr(1, (double)pos_n / (double)N, s);

    ke = ceil((double)numCounter / (double)pos_n *log(2) + (log(Ere_Pos) - s2)/log(2));
    
    //uint32_t op_ke = uint32_t((double)numCounter / (double)pos_n * log(2));
    //if(ke < op_ke) ke = op_ke;
    //std::cout << "op_ke: "<<op_ke << std::endl;
    size_t count = 0;
    size_t bits_per_item_ite = 50;
    for(auto v:pos_keys_){
        countingBloom_->counting_bloom_add(*v, ke);
        count++;
        double true_bits_per_item = sumSize / (double)count; 
        if( bits_per_item_ite >=10 && true_bits_per_item < bits_per_item_ite)
        {
        
            double fpr_c_temp = 0;
            double total_cost_temp = 0.0;
            for(int i=0; i<neg_keys_.size(); i++){
                if(Contain(*neg_keys_[i]))
                {
                    fpr_c_temp += neg_keys_[i]->cost;
                }
                total_cost_temp += neg_keys_[i]->cost;
            } 
            std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  Weighted FPR:  "<< (double)(fpr_c_temp*100) / total_cost_temp << "%" <<endl;    
            bits_per_item_ite -=2;          
        }
    }

    delete [] wbf_pos_keys_;
    delete [] wbf_neg_keys_;
}


bool WCBF::Contain(Slice &key){
    uint32_t k = ke;
    // for(Slice * key2 : cost_list_){
    //     if(key.str == key2->str){
    //         double Ere = getEr(key2->cost+1.0, (double)pos_n / (double)N, s);
    //         k = ceil((double)numCounter / (double)pos_n *log(2) + (log(Ere) - s2)/log(2));
    //         // cout<< "k:"<<k<<endl;
    //         //if (k<ke)
    //           //  k = int((double)numCounter / (double)pos_n * log(2));
    //         // std::cout << k << std::endl;
    //         break;
    //     }
    // }
    //set<Slice>::iterator iter = setCost_list.find(key);
    set<Slice>::iterator iter = setCost_list.find(key);
    if(iter != setCost_list.end()){
        double Ere = getEr(iter->cost + 1.0, (double)pos_n / (double)N, s);
        k = ceil((double)numCounter / (double)pos_n *log(2) + (log(Ere) - s2)/log(2));          
    }
    return countingBloom_->counting_bloom_check(key, k);
}

double WCBF::testFPR(){
    uint32_t k = ke;
    double fpr_c = 0;
    double total_cost_ = 0.0;
    for(Slice * key2 : cost_list_){
        double Ere = getEr(key2->cost+1.0, (double)pos_n / (double)N, s);
        k = ceil((double)numCounter / (double)pos_n *log(2) + (log(Ere) - s2)/log(2));
        if (k<ke)
            k = int((double)numCounter/ (double)pos_n * log(2));
        if(countingBloom_->counting_bloom_check(*key2, k))
            fpr_c += key2->cost;
        total_cost_ += key2->cost;
    }
    for(Slice * key2 : no_cost_list_){
        if(countingBloom_->counting_bloom_check(*key2, ke))
            fpr_c += key2->cost;
        total_cost_ += key2->cost;
    }
    return (double)(fpr_c*100) / total_cost_;
}



}





#endif