#ifndef DATALOADER
#define DATALOADER

#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <algorithm>
#include "key.h"

// Dataset directory
//#define SHALLA_PATH  "../data/dataset_hosefire_long/hosefile0_25"

//#define SHALLA_PATH  "../data/dataset_hosefire/dataset_hosefire_0_5/hosefile"
#define SHALLA_PATH  "../data/dataset_hosefire/dataset_hosefire_1_5/hosefire1_5_data"
#define HOSEFIRE_PATH  "../data/dataset_hosefire/dataset_hosefire_1_0/hosefire1_0_data"

//#define SHALLA_PATH  "../data/dataset_hosefire_0_25/hosefile0_25_data1"
//#define SHALLA_PATH  "../data/shalla_cost"
//#define SHALLA_PATH  "../../hashbloom/HashAdaptiveBF/data/dataset_hosefire/hosefile0_5"
#define YCSB_PATH  "../data/dataset_ycsb/ycsbt_1_5/ycsbt_cost_1.5-"
// #define YCSB_PATH  "../data/ycsb_generator/ycsbt_cost_1.5-"
// #define YCSB_PATH  "../data/ycsb_generator/temp/ycsbt_cost_1.0-"

#define SHALLA_PATH_SK  "../data/dataset_shalla/shalla_2_0/shalla_cost_2.0-"
// #define SHALLA_PATH_SK  "../data/dataset_skewness/skewness1_0/shalla_1_0-"
//#define SHALLA_PATH_SK  "../data/dataset_skewness/tenDataset/shalla_cost_1.5-"
// #define SHALLA_PATH_SK "../data/dataset_skewness/temp/shalla_cost_2.0-"

#define YCSB_PATH_SK  "../data/dataset_shalla/shalla_2_5/shalla_cost2.0-"

#define HOSEFIRE_PATH_1  "../data/dataset_hosefire_0_25/hosefile0_25_data1"
#define HOSEFIRE_PATH_2  "../data/dataset_hosefire_0_25/hosefile0_25_data2"
#define HOSEFIRE_PATH_3  "../data/dataset_hosefire_0_25/hosefile0_25_data3"
#define HOSEFIRE_PATH_4  "../data/dataset_hosefire_0_25/hosefile0_25_data4"
#define HOSEFIRE_PATH_5  "../data/dataset_hosefire_0_25/hosefile0_25_data5"




#define RANDOM_KEYSTR_PATH "../util/randomKeyStr.txt"
#define RANDOM_COST_TYPE zipf

enum {uniform, hotcost, normal, zipf};

// data loader for datasets and random data
class dataloader{
    public:
        std::vector<Slice *> pos_keys_; // positive keys
        std::vector<Slice *> neg_keys_; // negative keys

        ~dataloader();
        bool load(std::string data_name_, bool using_cost_);
        bool load(std::string data_name_, bool using_cost_, std::string epcho);
        bool loadRandomKey(int positives_number, int negatives_number, bool using_cost_);
        bool loadYCSB(bool using_cost_, std::string epcho);
        bool loadShalla(bool using_cost_, std::string epcho);
        bool loadHosefire(bool using_cost_, std::string epcho);
        bool loadShalla_sk(bool using_cost_, std::string epcho);
        bool loadNegtive();
        bool loadShalla_vulNeg(std::string name, bool using_cost_, std::string epcho);
    private:
        
        
        bool loadHosefire_1(bool using_cost_, std::string epcho);
        bool loadHosefire_2(bool using_cost_, std::string epcho);
        bool loadHosefire_3(bool using_cost_, std::string epcho);
        bool loadHosefire_4(bool using_cost_, std::string epcho);
        bool loadHosefire_5(bool using_cost_, std::string epcho);
        

};

// generate random keys and costs
class KeyBuilder{
    public:
        KeyBuilder();
        std::string GetKeyStr();
        bool ReadKeys(std::vector<Slice *> &v, int start_position_);
        void GenKeyStrAndToFile();
        void GenKeysUniformCosts(std::vector<Slice *> &keys, int interval);
        void GenKeysHotCosts(std::vector<Slice *> &keys, double hotNumberpro, int hotcost, int coldcost);
        void GenKeysNormalCosts(std::vector<Slice *> &keys, int u, int d);
        void GenKeysZipfCosts(std::vector<Slice *> &keys, double a, double c);
    private:
        std::vector<std::string> key_strs;
};

bool dataloader::load(std::string data_name_, bool using_cost_){
    if(data_name_ == "shalla_0") return loadShalla_sk(using_cost_, "0");
    else if(data_name_ == "shalla_1") return loadShalla_sk(using_cost_, "1");
    else if(data_name_ == "shalla_2") return loadShalla_sk(using_cost_, "2");
    else if(data_name_ == "shalla_3") return loadShalla_sk(using_cost_, "3");
    else if(data_name_ == "shalla_4") return loadShalla_sk(using_cost_, "4");
    else if(data_name_ == "ycsb") return loadYCSB(using_cost_, "");
    else if(data_name_ == "hosefire_1") return loadHosefire_1(using_cost_, "");
    else if(data_name_ == "hosefire_2") return loadHosefire_2(using_cost_, "");
    else if(data_name_ == "hosefire_3") return loadHosefire_3(using_cost_, "");
    else if(data_name_ == "hosefire_4") return loadHosefire_4(using_cost_, "");
    else if(data_name_ == "hosefire_5") return loadHosefire_5(using_cost_, "");
    

    else return false;
}

bool dataloader::load(std::string data_name_, bool using_cost_, std::string epcho){
    if(data_name_ == "shalla") return loadShalla(using_cost_, epcho);
    else if(data_name_ == "ycsb") return loadYCSB(using_cost_, epcho);
    else return false;
}
bool dataloader::loadShalla_sk(bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(SHALLA_PATH_SK+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadShalla_vulNeg(std::string name, bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(name+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadNegtive(){
    std::cout << "negative item reading..."  << std::endl;
    std::ifstream is("../data/addNegtive.txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = cost;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadShalla(bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(SHALLA_PATH+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}
bool dataloader::loadHosefire(bool using_cost_, std::string epcho){
    std::cout << "hosefire reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_1(bool using_cost_, std::string epcho){
    std::cout << "hosefire1 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_1+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_2(bool using_cost_, std::string epcho){
    std::cout << "hosefire2 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_2+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_3(bool using_cost_, std::string epcho){
    std::cout << "hosefire3 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_3+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_4(bool using_cost_, std::string epcho){
    std::cout << "hosefire4 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_4+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_5(bool using_cost_, std::string epcho){
    std::cout << "hosefire5 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_5+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}



bool dataloader::loadYCSB(bool using_cost_, std::string epcho){
    std::cout << "ycsb reading..."  << std::endl;
    std::ifstream is(YCSB_PATH+epcho+".txt");
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "FILTERKEY" || optype == "1") pos_keys_.push_back(key);
            else if(optype == "OTHERKEY" || optype == "0"){          
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadRandomKey(int positives_number, int negatives_number, bool using_cost_){
    KeyBuilder kb;
    for(int i=0; i<positives_number; i++){
        Slice *key = new Slice();
        pos_keys_.push_back(key);
    }
    for(int i=0; i<negatives_number; i++){
        Slice *key = new Slice();
        neg_keys_.push_back(key);
    }
    if(!kb.ReadKeys(pos_keys_, 0)) return false; 
    if(!kb.ReadKeys(pos_keys_, 0)) return false; 
    if(using_cost_)
        switch (RANDOM_COST_TYPE)
        {
        case uniform:
            kb.GenKeysUniformCosts(neg_keys_, 5);
            break;
        case hotcost:
            kb.GenKeysHotCosts(neg_keys_, 0.01, 100, 1);
            break;
        case normal:
            kb.GenKeysNormalCosts(neg_keys_, 20, 10);
            break;
        case zipf:
            kb.GenKeysZipfCosts(neg_keys_, 1.25, 1.0);
            break;
        default:
            break;
        }
    return true;
}

dataloader::~dataloader(){
    for(Slice *key : pos_keys_)
        delete key;
    for(Slice *key : neg_keys_)
        delete key;
}

KeyBuilder::KeyBuilder(){
    std::fstream ifs(RANDOM_KEYSTR_PATH);
    if(ifs.is_open()){
        std::string s;
        while(std::getline(ifs,s))
            key_strs.push_back(s);
    }else{
        std::cout << "Keystr file not exists, generate again..." << std::endl;
        GenKeyStrAndToFile();
    }
    ifs.close();
}
std::string KeyBuilder::GetKeyStr(){
    int k=rand()%10+1;
    char arr[10];
    for(int i=1;i<=k;i++){
        int x,s;                         
        s=rand()%2;                     
        if(s==1) x=rand()%('Z'-'A'+1)+'A';        
        else x=rand()%('z'-'a'+1)+'a';      
        arr[i-1] = x;                  
    }
    return std::string(arr,k);
}

bool KeyBuilder::ReadKeys(std::vector<Slice *> &v, int start_position_){
    int size_ = v.size();
    if(start_position_ + size_ >= key_strs.size()) return false;
    for(int j=0; j<size_; j++)
        v[j]->str = key_strs[start_position_+j];
    return true;
}

void KeyBuilder::GenKeyStrAndToFile(){
    int gen_key_size_ = 200000;
    int i = key_strs.size();
    std::ofstream ofs(RANDOM_KEYSTR_PATH);
    while(i < gen_key_size_){
        if(0 == i%10000) std::cout << i << "keys have been created..." << std::endl;
        std::string str = GetKeyStr();
        if(std::find(key_strs.begin(), key_strs.end(), str) == key_strs.end()){
            key_strs.push_back(str);
            ofs << str << std::endl;
            i++;
        }
    }
    ofs.close();
}

void KeyBuilder::GenKeysUniformCosts(std::vector<Slice *> &keys, int interval){
    for(int i=0; i<keys.size(); i++)
        keys[i]->cost = 1+(i+1)*interval;
}
void KeyBuilder::GenKeysHotCosts(std::vector<Slice *> &keys, double hotNumberpro, int hotcost, int coldcost){
    int hotNumSize = hotNumberpro * keys.size();
    for(int i=0; i<keys.size(); i++){
        if(i <= hotNumSize) 
            keys[i]->cost = hotcost;
        else 
            keys[i]->cost = coldcost;
    }
}
void KeyBuilder::GenKeysNormalCosts(std::vector<Slice *> &keys, int u, int d){
    for(int i=0; i<keys.size(); i++){
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::normal_distribution<double> distribution(u, d);
        keys[i]->cost = distribution(generator);
    }
}

void KeyBuilder::GenKeysZipfCosts(std::vector<Slice *> &keys, double a, double c){
    int r = 10000;
    double pf[10000];
    double sum = 0.0;
    for (int i = 0; i < r; i++)        
        sum += c/pow((double)(i+2), a);  
    for (int i = 0; i < r; i++){ 
        if (i == 0)
            pf[i] = c/pow((double)(i+2), a)/sum;
        else
            pf[i] = pf[i-1] + c/pow((double)(i+2), a)/sum;
    }
     for (int i = 0; i < keys.size(); i++){
        int index = 0;
        double data = (double)rand()/RAND_MAX;  
        while (data > pf[index])  
            index++;
        keys[i]->cost = index;
    }
}

#endif