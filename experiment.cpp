#include <time.h>

#include "util/dataloader.h"

#include "habf/fasthabf.h"
#include "nonlearnedfilter/wbf.h"
#include "habf/dynamic_cbf_one.h"
#include "habf/countingBloom.h"
#include "habf/gitCountingBloom.h"
#include "nonlearnedfilter/wcbf.h"
#include "nonlearnedfilter/stackedfilter/StackedFilter.h"
#include <vector>

using namespace std;
using namespace std::chrono;

void TestDHABF(dataloader &dl, double bits_per_key);
void TestCountingBloom_XXH128(dataloader &dl, double bits_per_key);
void TestCountingBloom_murmur(dataloader &dl, double bits_per_key);
void TestWeightedCountingBloom(dataloader &dl, double bits_per_key);
void TestStackedFilter(dataloader &dl, double bits_per_key);

void ConvertDataset(dataloader &dl, std::vector<StringElement> &positives, std::vector<StringElement> &negatives, std::vector<double> &cdf);
bool compare_func(Slice * key1, Slice * key2);

int main(){

    const string dataset_ = "shalla"; // shalla or ycsb
    const double space_size_ = 1.2; // total space size of data structure (MB)
    dataloader dl;
    if(dl.load(dataset_, true)){
        cout << "keys size (positives, negatives): (" << dl.pos_keys_.size() << "," << dl.neg_keys_.size() << ")" << endl;
        cout << "The number of positive items: " << dl.pos_keys_.size() << "  |   " <<"The number of negtive items: " <<  dl.neg_keys_.size()<<endl;
        std::cout << "Space size allocated for each data structure: " << space_size_ << " MB" << endl;
        double bits_per_key = (double)(space_size_ * 1024 * 1024 * 8) / (double) dl.pos_keys_.size();
        std::cout << "bits_per_key: " << bits_per_key << endl;

        if(bits_per_key < 5){
            std::cout << "bits_per_key is too small" << endl;
            return 0;
        }

        std::cout << "---------------------------------------------------" << endl;

        cout <<"dynamicHABF testing......" << endl;
        TestDHABF(dl, bits_per_key);

        std::cout << "---------------------------------------------------" << endl;

        cout <<"countingBloom testing......" << endl;
        TestCountingBloom_XXH128(dl, bits_per_key);

        std::cout << "---------------------------------------------------" << endl;

        cout <<"WCBF testing......" << endl;
        TestWeightedCountingBloom(dl, bits_per_key);

        std::cout << "---------------------------------------------------" << endl;

        cout <<"stackedFilter testing......" << endl;
        TestStackedFilter(dl, bits_per_key);        


    }else{
        cout << "dataset reading failed"  << endl;
    }
}

void TestDHABF(dataloader &dl, double bits_per_key){

    // construction testing
    size_t sumSize = bits_per_key * dl.pos_keys_.size();
    
    dynamichabf::DynamicHABFilter_One dychabf(bits_per_key ,dl.pos_keys_.size(), dl.neg_keys_.size(), dl.neg_keys_ );

    auto t1 = steady_clock::now();

    size_t count=0;
    int bits_per_item_ite = 50;

    for (Slice *v : dl.pos_keys_ ){
        
        dychabf.Add( *v);
        count++;
        double true_bits_per_item = sumSize / (double)count; 
        if( bits_per_item_ite >=10 && true_bits_per_item < bits_per_item_ite)
        {
        
            double fpr_c_temp = 0;
            double total_cost_temp = 0.0;
            size_t fpNum_temp = 0;
            for(int i=0; i<dl.neg_keys_.size(); i++){
                if(dychabf.Contain(dl.neg_keys_[i]))
                {
                    fpr_c_temp += dl.neg_keys_[i]->cost;
                    fpNum_temp ++;
                }
                total_cost_temp += dl.neg_keys_[i]->cost;
            } 
            std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  Weighted FPR:  "<< (double)(fpr_c_temp*100) / total_cost_temp << "%" <<endl;    
            //std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  FPR:  "<< (double)(fpNum_temp*100) / dl.neg_keys_.size() << "%" <<endl;    
            bits_per_item_ite -=2;          
        }
    }

    auto t2 = steady_clock::now();
    double construct_time_ = duration<double>(t2 - t1).count();

    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;

    // query testing
    size_t fpr_Num = 0;
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;

    t1 = steady_clock::now();

    for(int i=0; i<dl.pos_keys_.size(); i++)
    {
        if(!dychabf.Contain(dl.pos_keys_[i]))
        {
            fnr_c++;
        } 
    }
        
    for(int i=0; i<dl.neg_keys_.size(); i++)
    {
        if(dychabf.Contain(dl.neg_keys_[i]))
        {
            fpr_c += dl.neg_keys_[i]->cost;
            fpr_Num ++;
        }

        total_cost_ += dl.neg_keys_[i]->cost;
    }

    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FPR: "<< (double)(fpr_Num*100) / (double)dl.neg_keys_.size() << "%" <<endl;
    std::cout<<"FNR: "<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
    std::cout<<"onlyBloom "<<dychabf.onlyBloomNum<<endl;
    std::cout<<"belowBloom "<<dychabf.belowBloomNum<<endl;
    //std::cout<<"average adapt hash: "<<dychabf.adaptNum/(double)dl.pos_keys_.size()<<endl;
    dychabf.PrintHE_haveInsertNum();

}

void TestCountingBloom_XXH128(dataloader &dl, double bits_per_key){

    auto t1 = steady_clock::now();
    size_t nfuncs = 4;
    size_t bitPerCounter =4;
    size_t sumSize = bits_per_key * dl.pos_keys_.size();

    cout<<"sumSize"<<sumSize<<endl;

    size_t counter_Num = sumSize / bitPerCounter;
    cout<<"counter_Num: "<<counter_Num<<endl;

    countingBloom::CountingBloom cbf(counter_Num, bitPerCounter, nfuncs);
    
    // construction testing
    size_t fpNum = 0;
    size_t count=0;
    
    int bits_per_item_ite = 50;
    
    for (Slice *v : dl.pos_keys_ ){
        cbf.counting_bloom_add( *v);
        count++;
        double true_bits_per_item = sumSize / (double)count; 
        if( bits_per_item_ite >=10 && true_bits_per_item < bits_per_item_ite)
        {
            
            double fpr_c_temp = 0;
            double total_cost_temp = 0.0;
            size_t fpNum_temp = 0;
            for(int i=0; i<dl.neg_keys_.size(); i++){
                if(cbf.counting_bloom_check(*dl.neg_keys_[i])){
                    fpr_c_temp += dl.neg_keys_[i]->cost;
                    fpNum_temp ++;
                }
                    
                total_cost_temp += dl.neg_keys_[i]->cost;
            }         
            std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  Weighted FPR:  "<< (double)(fpr_c_temp*100) / total_cost_temp << "%" <<endl;    
            //std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  FPR:  "<< (double)(fpNum_temp*100) / dl.neg_keys_.size() << "%" <<endl; 
            bits_per_item_ite -= 2;   
        }
    }

    auto t2 = steady_clock::now();
    double construct_time_ = duration<double>(t2 - t1).count();

    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;

    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;

    t1 = steady_clock::now();
    for (Slice *v : dl.pos_keys_ )
    {
       if(!cbf.counting_bloom_check(*v)){
           fnr_c++;
       } 
    }
    
    for (Slice *v : dl.neg_keys_ )
    {
        if(cbf.counting_bloom_check(*v))
        {
            fpr_c += v->cost;
            fpNum++;
        }  
        total_cost_ += v->cost;
    }
    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    
    std::cout<<"Query time/key (ns): "<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR: "<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<< "FPR: "<<(fpNum*100)/(double)dl.neg_keys_.size()<<"%"<<endl;;
    std::cout<<"FNR: "<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;
}

void TestWeightedCountingBloom(dataloader &dl, double bits_per_key){

    size_t bits_Per_Counter = 4;
    weightedCountingBloom::WCBF wcbf(bits_per_key, bits_Per_Counter, dl.pos_keys_, dl.neg_keys_);
    
    auto t1 = steady_clock::now();

    wcbf.AddAll(dl.pos_keys_, dl.neg_keys_);

    auto t2 = steady_clock::now();
    double construct_time_ = duration<double>(t2 - t1).count();
    std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;

    // query testing
    double fpr_c = 0;
    int fnr_c = 0;
    double total_cost_ = 0.0;
    size_t fpNum = 0;
    t1 = steady_clock::now();

    for(int i=0; i<dl.pos_keys_.size(); i++)
    {
        if(!wcbf.Contain(*dl.pos_keys_[i]))
        {
            fnr_c++;
        } 
    }
        
    for(int i=0; i<dl.neg_keys_.size(); i++)
    {
        if(wcbf.Contain(*dl.neg_keys_[i]))
        {
            fpr_c += dl.neg_keys_[i]->cost;
            fpNum++;
        }
        total_cost_ += dl.neg_keys_[i]->cost;
    }
        
    t2 = steady_clock::now();
    double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
    std::cout<<"FPR: "<<(fpNum*100)/(double)dl.neg_keys_.size()<<"%"<<endl;;
    std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
    std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
    std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;   
}

void TestStackedFilter(dataloader &dl, double bits_per_key){
        

        std::vector<StringElement> positives, negatives;
        std::vector<double> cdf;

        ConvertDataset(dl, positives, negatives, cdf);
        
        size_t bit_Per_Counter =4;
        size_t sumSize = bits_per_key * dl.pos_keys_.size();

        StackedFilter<CountingBloomFilterLayer, StringElement> filter1(bits_per_key * positives.size(), positives, negatives, cdf, bit_Per_Counter);
        
        auto t1 = steady_clock::now();
        size_t count= 0;
        size_t bits_per_item_ite = 50;

        for(StringElement &s : positives)
        {
            filter1.InsertPositiveElement(s);

            count++;
            double true_bits_per_item = sumSize / (double)count; 
            if( bits_per_item_ite >=10 && true_bits_per_item < bits_per_item_ite)
            {
                
                double fpr_c_temp = 0;
                double total_cost_temp = 0.0;
                for(int i=0; i<negatives.size(); ++i){
                    //false_positives += filter.LookupElement(negatives[i]) * cdf[i];
                    if(filter1.LookupElement(negatives[i])){
                        fpr_c_temp +=  cdf[i];
                        
                    }
                    total_cost_temp += cdf[i];
                }    
                std::cout<<"bits_per_item = "<< bits_per_item_ite<< "  Weighted FPR:  "<< (double)(fpr_c_temp*100) / total_cost_temp << "%" <<endl;    
                bits_per_item_ite -= 2;   
            }
        }


        auto t2 = steady_clock::now();
        double construct_time_ = duration<double>(t2 - t1).count();
        std::cout<<"construct time/key (ns):"<< 1000000000.0 * construct_time_ / (double)(dl.pos_keys_.size()) <<endl;


        size_t fnr_c = 0 ;
        size_t fpNum = 0;
       
        double total_cost_ = 0;
        double fpr_c = 0;

        t1 = steady_clock::now();

        for(int i=0; i<positives.size(); ++i){
            //false_positives += filter.LookupElement(negatives[i]) * cdf[i];
            if(!filter1.LookupElement(positives[i])){
                fnr_c++;
            }
        }

    
 
        
        for(int i=0; i<negatives.size(); ++i){
            //false_positives += filter.LookupElement(negatives[i]) * cdf[i];
            if(filter1.LookupElement(negatives[i])){
                fpr_c +=  cdf[i];
                fpNum++;
            }
            total_cost_ += cdf[i];
        }

        t2 = steady_clock::now();

        double query_perkey_time =  1000000000.0 * duration<double>(t2 - t1).count() / (double)(dl.pos_keys_.size()+dl.neg_keys_.size());
        std::cout<<"FPR: "<<(fpNum*100)/(double)dl.neg_keys_.size()<<"%"<<endl;;
        std::cout<<"query time/key (ns):"<< query_perkey_time <<endl;
        std::cout<<"Weighted FPR:"<< (double)(fpr_c*100) / total_cost_ << "%" <<endl;
        std::cout<<"FNR:"<< (double)(fnr_c*100) / (double)dl.pos_keys_.size() << "%" <<endl;  
}

void ConvertDataset(dataloader &dl, std::vector<StringElement> &positives, std::vector<StringElement> &negatives, std::vector<double> &cdf){
    positives.resize(dl.pos_keys_.size());
    negatives.resize(dl.neg_keys_.size());
    cdf.resize(dl.neg_keys_.size());
    sort(dl.neg_keys_.begin(), dl.neg_keys_.end(), compare_func);
    for(int i=0; i<positives.size(); ++i){
        positives[i] = dl.pos_keys_[i]->str;
    }
    for(int i=0; i<negatives.size(); ++i){
        negatives[i] = dl.neg_keys_[i]->str;
    }
    for(int i=0; i<cdf.size(); ++i){
        cdf[i] = dl.neg_keys_[i]->cost;
    }
}

bool compare_func(Slice * key1, Slice * key2){
    return key1->cost > key2->cost;
}
