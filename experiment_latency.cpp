//Description: Used to evaluate the latecy of various filter

#include <time.h>
#include "util/dataloader.h"
#include "SSCF/SeesawCF.h"
#include "nonlearnedfilter/countingBloom.h"
#include "nonlearnedfilter/wcbf.h"
#include "nonlearnedfilter/stackedfilter/StackedFilter.h"
#include <vector>

#define SHALLA_NAME "shalla_"
#define YCSBT_NAME "ycsbt_"
#define HOSEFIRE_NAME "hosefire_"
const size_t datasetNum = 5; 

using namespace std;
using namespace std::chrono;

void TestSSCF(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets);
void TestWeightedCountingBloom(dataloader * dlSet,  double * bits_per_key_sets);
void TestStackedFilter(dataloader * dlSet,  double * bits_per_key_sets, double vulNegKeyRatio);

void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf);
bool compare_func(Slice * key1, Slice * key2);
bool compare_WcbfKey_func(weightedCountingBloom::WCBFKey &key1, weightedCountingBloom::WCBFKey &key2);

int main(){
    const double space_size_ = 4.25; // total space size of data structure (MB) 4.25
    dataloader * dlSet = new dataloader[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        //load hosefire
        // if(!dlSet[i].loadHosefire(true, to_string(i+1)))
        // {
        //     std::cout<< "error!\n";
        // }
        
        //load shalla
        if(!dlSet[i].loadShalla_sk( true, to_string(i)))
        {
            std::cout<< "error!\n";
        }

        //load ycsbt
        // if(!dlSet[i].loadYCSB(true, to_string(i)))
        // {
        //     std::cout<< "error!\n";
        // }
    }
    double *bits_per_key_sets = new double[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        bits_per_key_sets[i] = (double)(space_size_ * 1024 * 1024 * 8) / dlSet[i].pos_keys_.size();
        std::cout<<"Dataset "<< i+1 << " :bits_Per_Item: "<<bits_per_key_sets[i]<<endl;
    }

    for(size_t i = 0; i<datasetNum; i++){
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);
        }
    }
    double vulNegKeyRatio =0.05;
    std::cout << "Space size allocated for each data structure: " << space_size_ << " MB" << endl;
    
    std::cout << "--------------SSCF testing----------------" << endl;
    TestSSCF(dlSet, bits_per_key_sets, vulNegKeyRatio);

    std::cout << "-------------countingBloom testing---------------" << endl;
    TestCountingBloom(dlSet, bits_per_key_sets);

    std::cout << "------------------WCBF testing-------------------" << endl;
    TestWeightedCountingBloom(dlSet, bits_per_key_sets);

    std::cout << "--------------stackedFilter testing-------------------" << endl;
    TestStackedFilter(dlSet, bits_per_key_sets, vulNegKeyRatio);   
}

void TestSSCF(dataloader * dlSet,  double * bits_per_key_sets, double vulNegKeyRatio)
{
    SeesawCF::SeesawCF *sscfArr[datasetNum];

    for (size_t i = 0; i < datasetNum; i++)
    {
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            std::cout<<"error!\n";
        }
        sscfArr[i] = new SeesawCF::SeesawCF(vulNegKeyRatio, bits_per_key_sets[i], dlSet[i].pos_keys_.size(), dlSet[i].neg_keys_.size(), dlSet[i].neg_keys_ , true);
    }
    
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now();

        for (Slice *v : dlSet[i].pos_keys_ )
        {          
            sscfArr[i]->Add( *v);
        }

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }
    
    std::cout<<"Average Insertion Latency (Time/Key ns): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;
        size_t j=0;
        for (j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            total_cost_ += dlSet[i].neg_keys_[j]->cost;
        }        

        auto t1 = steady_clock::now();

        for (j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if (!sscfArr[i]->Contain(dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
            
        for (j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            if (sscfArr[i]->Contain(dlSet[i].neg_keys_[j]))
            {
                wfpr_c += dlSet[i].neg_keys_[j]->cost;
                fpr_Num ++;
            }
        }
        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size()+dlSet[i].neg_keys_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;

        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / total_cost_;
        if((double)(wfpr_c*100) / total_cost_ > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        if((double)(wfpr_c*100) / total_cost_ < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Average Query Latency (Time/Key ns): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;


    std::cout<<"~~~~~~~SSCF Start Deleting~~~~~~~\n";
    double sumDeleteTime = 0;
    double minDeleteTime = MAXFLOAT;
    double maxDeleteTime = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        //random shuffle
        vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            deleteItem[j] = dlSet[i].pos_keys_[j]; 
        }
        random_shuffle(deleteItem.begin(),deleteItem.end());

        auto t1 = steady_clock::now();
        for( Slice * v : deleteItem)
        {
            sscfArr[i]->Deletion(v);
        }
        auto t2 = steady_clock::now();   

        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
        sumDeleteTime += tempV;
        minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
        maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
        sscfArr[i]->PrintCounterAboveOne();
    }    
    std::cout<<"Average Deletion Latency (Time/Key ns): "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    double sumFNR_afterDelete = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        
        size_t fnr_c = 0;
        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if(!sscfArr[i]->Contain(dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
        sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }    
    std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl;

    for (size_t i = 0; i < datasetNum; i++) delete sscfArr[i] ;
}

void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets){

    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[5];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        nfuncs = bits_per_key_sets[i] / 4 * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] * dlSet[0].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now(); 

        for (Slice *v : dlSet[i].pos_keys_ )
        {    
            cbfArr[i]->counting_bloom_add( *v);
        }   

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;

        for (Slice *v : dlSet[i].neg_keys_ )
        {
            total_cost_ += v->cost;
        }

        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if (!cbfArr[i]->counting_bloom_check(*v))
            {
                fnr_c++;
            } 
        }
        
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            if(cbfArr[i]->counting_bloom_check(*v))
            {
                wfpr_c += v->cost;
                fpr_Num ++;
            }  
        }
        auto t2 = steady_clock::now();   
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size()+dlSet[i].neg_keys_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / total_cost_;
        if((double)(wfpr_c*100) / total_cost_ > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        if((double)(wfpr_c*100) / total_cost_ < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    
    std::cout<<"~~~~~~~CBF Start Deleting~~~~~~~\n";
    double sumDeleteTime = 0;
    double minDeleteTime = MAXFLOAT;
    double maxDeleteTime = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        //random shuffle
        vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            deleteItem[j] = dlSet[i].pos_keys_[j]; 
        }
        random_shuffle(deleteItem.begin(),deleteItem.end());
        
        auto t1 = steady_clock::now();
        for (Slice *v : deleteItem )
        {
            cbfArr[i]->counting_bloom_remove(*v);
        }
        auto t2 = steady_clock::now();       
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
        sumDeleteTime += tempV;
        minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
        maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
        cbfArr[i]->printHaveInsertNum();
    }  
    std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    
    double sumFNR_afterDelete=0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if(!cbfArr[i]->counting_bloom_check(*v)){
                fnr_c++;
            } 
        }
    
        sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl; 
}

void TestWeightedCountingBloom(dataloader * dlSet, double * bits_per_key_sets)
{
    size_t bits_Per_Counter = 4;
    weightedCountingBloom::WCBF *wcbfArr[datasetNum];
    for (size_t i = 0; i < datasetNum; i++)
    {
        wcbfArr[i] =  new weightedCountingBloom::WCBF(bits_per_key_sets[i], bits_Per_Counter, dlSet[i].pos_keys_, dlSet[i].neg_keys_);
    }

    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t pos_n = dlSet[i].pos_keys_.size();
        size_t neg_n = dlSet[i].neg_keys_.size();
        weightedCountingBloom::WCBFKey* wbf_pos_keys_ = new weightedCountingBloom::WCBFKey[dlSet[i].pos_keys_.size()];
        weightedCountingBloom::WCBFKey* wbf_neg_keys_ = new weightedCountingBloom::WCBFKey[dlSet[i].neg_keys_.size()];
        
        for(int j=0; j<pos_n; ++j){
             wbf_pos_keys_[j].data_ = dlSet[i].pos_keys_[j];
        }
        for(int j=0; j<neg_n; ++j){
            wbf_neg_keys_[j].data_ = dlSet[i].neg_keys_[j];
        }
      
        if(!is_sorted(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_WcbfKey_func)){
            std::sort(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_WcbfKey_func);
        }            

        auto t1 = steady_clock::now();
        wcbfArr[i]->AddAll_latency(wbf_pos_keys_, wbf_neg_keys_, pos_n, neg_n);
        auto t2 = steady_clock::now();
        
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
        delete [] wbf_pos_keys_;
        delete [] wbf_neg_keys_;
    }
    std::cout<<"Average Insertion (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        double wfpr_c = 0;
        size_t fnr_c = 0;
        double total_cost_ = 0.0;
        size_t fpr_Num = 0;
        for (int j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            total_cost_ += dlSet[i].neg_keys_[j]->cost;
        }

        auto t1 = steady_clock::now();

        for (int j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if(!wcbfArr[i]->Contain(*dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
            
        for (int j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            if (wcbfArr[i]->Contain(*dlSet[i].neg_keys_[j]))
            {
                wfpr_c += dlSet[i].neg_keys_[j]->cost;
                fpr_Num++;
            }
        }
        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size()+dlSet[i].neg_keys_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;

        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / total_cost_;
        if((double)(wfpr_c*100) / total_cost_ > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        if((double)(wfpr_c*100) / total_cost_ < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }

    std::cout<<"Average Query Latency (Time(ns)/key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    
    std::cout<<"~~~~~~~WCBF Start Deleting~~~~~~~\n";
    double sumDeleteTime = 0;
    double minDeleteTime = MAXFLOAT;
    double maxDeleteTime = 0;
    double sum_DeleteTime=0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            deleteItem[j] = dlSet[i].pos_keys_[j]; 
        }
        random_shuffle(deleteItem.begin(),deleteItem.end());        
 
        auto t1 = steady_clock::now();

        wcbfArr[i]->DeleteAll(deleteItem);

        auto t2 = steady_clock::now();

        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
        sumDeleteTime += tempV;
        minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
        maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
        wcbfArr[i]->PrintHaveInsertNum();
    }
    std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;

    double sumFNR_afterDelete = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        int fnr_c = 0;
        for(int j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if(!wcbfArr[i]->Contain(*dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
        sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();   
    }
    std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl;
}

void TestStackedFilter(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio)
{     
    std::vector<std::vector<StringElement>> vec_positives;
    std::vector<std::vector<StringElement>> vec_negatives;
    std::vector<std::vector<double>> vec_cdf;
    StackedFilter<CountingBloomFilterLayer, StringElement> * filterArr[datasetNum];
    size_t bit_Per_Counter =4;
    vec_positives.resize(datasetNum);
    vec_negatives.resize(datasetNum);
    vec_cdf.resize(datasetNum);

    for (size_t i = 0; i < datasetNum; i++)
    {
        vec_negatives[i].resize(dlSet[i].neg_keys_.size());
        vec_positives[i].resize(dlSet[i].pos_keys_.size());
        vec_cdf[i].resize(dlSet[i].neg_keys_.size());
    }
    
    ConvertDataset(dlSet, vec_positives, vec_negatives, vec_cdf);

    for (size_t i = 0; i < datasetNum; i++)
    {
        filterArr[i] =  new StackedFilter<CountingBloomFilterLayer, StringElement>(bits_per_key_sets[i] * vec_positives[i].size(), vec_positives[i], vec_negatives[i], vec_cdf[i], bit_Per_Counter, vulNegKeyRatio);
    }

    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t haveAddCount= 0;
        size_t bits_per_item_ite = 50;
        size_t sumSize_temp = bits_per_key_sets[i] * vec_positives[i].size();

        auto t1 = steady_clock::now();

        for(StringElement &s : vec_positives[i])
        {
            filterArr[i]->InsertPositiveElement(s);
        }

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }
    std::cout<<"Average Insertion (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0 ;
        size_t fpr_Num = 0;
        double total_cost_ = 0;
        double wfpr_c = 0;
        for (int j=0; j<vec_negatives[i].size(); ++j)
        {
            total_cost_ += vec_cdf[i][j];
        }

        auto t1 = steady_clock::now();
        for (int j=0; j<vec_positives[i].size(); ++j)
        {
            if (!filterArr[i]->LookupElement(vec_positives[i][j]))
            {
                fnr_c++;
            }
        }

        for (int j=0; j<vec_negatives[i].size(); ++j)
        {
            if (filterArr[i]->LookupElement(vec_negatives[i][j]))
            {
                wfpr_c +=  vec_cdf[i][j];
                fpr_Num++;
            }
        }
        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size()+dlSet[i].neg_keys_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;

        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / total_cost_;
        if((double)(wfpr_c*100) / total_cost_ > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        if((double)(wfpr_c*100) / total_cost_ < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }

    std::cout<<"Average Query Latency (Time(ns)/key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR:"<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;

    std::cout<<"~~~~~~~SF Start Deleting~~~~~~~\n";
    double sumDeleteTime = 0;
    double minDeleteTime = MAXFLOAT;
    double maxDeleteTime = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {

        random_shuffle(vec_positives[i].begin(),vec_positives[i].end());
        auto t1 = steady_clock::now();
        for(int j=0; j<vec_positives[i].size(); ++j){
            filterArr[i]->DeleteElement(vec_positives[i][j]);
        }
        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
        sumDeleteTime += tempV;
        minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
        maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    }   
    std::cout<<"~~~~~~~SF Finish Deleting~~~~~~~\n";  
    std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;  
    double sumFNR_afterDelete = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0 ;
        for (int j=0; j<vec_positives[i].size(); ++j){
            if (!filterArr[i]->LookupElement(vec_positives[i][j]))
            {
                fnr_c++;
            }
        }  
        sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();        
    } 
    std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl;    
}

void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf){
    
    for (size_t i = 0; i < datasetNum; i++)
    {

        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);
        }
 
        for(int j=0; j<dlSet[i].pos_keys_.size(); ++j){
            vec_positives[i][j] = dlSet[i].pos_keys_[j]->str;
        }
        for(int j=0; j<dlSet[i].neg_keys_.size(); ++j){
            vec_negatives[i][j] = dlSet[i].neg_keys_[j]->str;
        }
        for(int j=0; j<dlSet[i].neg_keys_.size(); ++j){
            vec_cdf[i][j] = dlSet[i].neg_keys_[j]->cost;
        }        
    }
}

bool compare_func(Slice * key1, Slice * key2){
    return key1->cost > key2->cost;
}
bool compare_WcbfKey_func(weightedCountingBloom::WCBFKey &key1, weightedCountingBloom::WCBFKey &key2){
    return key1.data_->cost > key2.data_->cost;
}
