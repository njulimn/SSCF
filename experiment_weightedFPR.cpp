//Description: Used to evaluate a) Weighted FPR b) FPR c) FNR

#include <time.h>
#include "util/dataloader.h"
#include "SSCF/SeesawCF.h"
#include "nonlearnedfilter/countingBloom.h"
#include "nonlearnedfilter/wcbf.h"
#include "nonlearnedfilter/naiveFilter.h"
#include "nonlearnedfilter/stackedfilter/StackedFilter.h"
#include <vector>

#define SHALLA_NAME "shalla_"
#define YCSBT_NAME "ycsbt_"
const size_t datasetNum = 5; 

using namespace std;
using namespace std::chrono;

void TestSSCF(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets);
void TestNaiveFilter(dataloader * dlSet, double * bits_per_key_sets);
void TestWeightedCountingBloom(dataloader * dlSet,  double * bits_per_key_sets);
void TestStackedFilter(dataloader * dlSet,  double * bits_per_key_sets, double vulNegKeyRatio);

void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf);
bool compare_func(Slice * key1, Slice * key2);
bool compare_WcbfKey_func(weightedCountingBloom::WCBFKey &key1, weightedCountingBloom::WCBFKey &key2);
int main(){
    double bitsPerItem = 20;
    // const double space_size_ = 4.25; // total space size of data structure (MB)
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
        // };
    }
    

    // std::cout<<"pos item num: "<<dlSet[0].pos_keys_.size()<<endl;
    
    for(size_t i = 0; i<datasetNum; i++){
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);
        }
    }

    double *bits_per_key_sets = new double[5];
    double vulNegKeyRatio = 0.06;
    for (size_t round = 0; round<9; round++)
    {
        
        for (size_t i = 0; i < datasetNum; i++)
        {
            bits_per_key_sets[i] = bitsPerItem;

        }
        double spaceSize = bitsPerItem * dlSet[0].pos_keys_.size() / 8 / 1024 / 1024;
        std::cout<<"---------------bitsPerItem = "<<bitsPerItem<<"  SpaceSize = "<<spaceSize<<" MB----------------------------------\n";
        

        std::cout << "-------------------SSCF testing------------------" << endl;
        TestSSCF(dlSet, bits_per_key_sets, vulNegKeyRatio);
        
        std::cout << "-------------CountingBloom testing---------------" << endl;
        TestCountingBloom(dlSet, bits_per_key_sets);

        std::cout << "------------------WCBF testing-------------------" << endl;
        TestWeightedCountingBloom(dlSet, bits_per_key_sets);

        std::cout << "--------------StackedFilter testing--------------" << endl;
        TestStackedFilter(dlSet, bits_per_key_sets, vulNegKeyRatio);

        bitsPerItem += 2;
        
    }

    delete []bits_per_key_sets;

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
    
    for (size_t i = 0; i < datasetNum; i++)
    {
        for (Slice *v : dlSet[i].pos_keys_ )
        {
                sscfArr[i]->Add( *v);  
        }
    }

    double sumFPR = 0;
    double sumWeightedFPR = 0;
    double sumFNR = 0;
    double maxFPR = 0, minFPR = MAXFLOAT, maxWeightedFPR = 0, minWeightedFPR = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;

        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if(!sscfArr[i]->Contain(dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
            
        for(int j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            if(sscfArr[i]->Contain(dlSet[i].neg_keys_[j]))
            {
                wfpr_c += dlSet[i].neg_keys_[j]->cost;
                fpr_Num ++;
            }
            total_cost_ += dlSet[i].neg_keys_[j]->cost;
        }
        
        double tempFPR = (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();
        sumFPR +=tempFPR;
        maxFPR = tempFPR > maxFPR ? tempFPR : maxFPR;
        minFPR = tempFPR < minFPR ? tempFPR : minFPR; 

        double tempWeightedFPR = (double)(wfpr_c*100) / total_cost_;
        sumWeightedFPR += tempWeightedFPR;
        maxWeightedFPR = tempWeightedFPR > maxWeightedFPR ? tempWeightedFPR : maxWeightedFPR;
        minWeightedFPR = tempWeightedFPR < minWeightedFPR ? tempWeightedFPR : minWeightedFPR;  

        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - minWeightedFPR - maxWeightedFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)(sumFPR - minFPR - maxFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"Average FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;
    
    for (size_t i = 0; i < datasetNum; i++) delete sscfArr[i] ;
    
}

void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets)
{
    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[5];

    for (size_t i = 0; i < datasetNum; i++)
    {
        nfuncs = bits_per_key_sets[i] / 4 * log(2);
        // std::cout<<"nHashFunc: "<<nfuncs<<endl;
        cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[i] * dlSet[i].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
    }
    
    for (size_t i = 0; i < datasetNum; i++)
    {
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            cbfArr[i]->counting_bloom_add( *v);
        }     
    }

    double sumFPR = 0;
    double sumWeightedFPR = 0;
    double sumFNR = 0;
    double maxFPR = 0, minFPR = MAXFLOAT, maxWeightedFPR = 0, minWeightedFPR = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;

        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if(!cbfArr[i]->counting_bloom_check(*v)){
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

            total_cost_ += v->cost;
        }

        double tempFPR = (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();
        sumFPR +=tempFPR;
        maxFPR = tempFPR > maxFPR ? tempFPR : maxFPR;
        minFPR = tempFPR < minFPR ? tempFPR : minFPR; 

        double tempWeightedFPR = (double)(wfpr_c*100) / total_cost_;
        sumWeightedFPR += tempWeightedFPR;
        maxWeightedFPR = tempWeightedFPR > maxWeightedFPR ? tempWeightedFPR : maxWeightedFPR;
        minWeightedFPR = tempWeightedFPR < minWeightedFPR ? tempWeightedFPR : minWeightedFPR;  

        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - minWeightedFPR - maxWeightedFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)(sumFPR - minFPR - maxFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"Average FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;

    for (size_t i = 0; i < datasetNum; i++) delete cbfArr[i];

}


void TestWeightedCountingBloom(dataloader * dlSet, double * bits_per_key_sets)
{
    size_t bits_Per_Counter = 4;
    weightedCountingBloom::WCBF *wcbfArr[datasetNum];
    for (size_t i = 0; i < datasetNum; i++)
    {
        wcbfArr[i] =  new weightedCountingBloom::WCBF(bits_per_key_sets[i], bits_Per_Counter, dlSet[i].pos_keys_, dlSet[i].neg_keys_);
    }

    // double sumInsertTime = 0;
    // double maxInsertTime = 0;
    // double minInsertTime = MAXFLOAT;

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

        // auto t1 = steady_clock::now();
        wcbfArr[i]->AddAll_latency(wbf_pos_keys_, wbf_neg_keys_, pos_n, neg_n);
        // auto t2 = steady_clock::now();
        
        // double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        // sumInsertTime += tempV;
        // minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        // maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
        delete [] wbf_pos_keys_;
        delete [] wbf_neg_keys_;
    }
    // std::cout<<"Insert Time/Key(average): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
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

    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;

    for (size_t i = 0; i < datasetNum; i++) delete wcbfArr[i];
}

void TestStackedFilter(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio){
        
    std::vector<std::vector<StringElement>> vec_positives;
    std::vector<std::vector<StringElement>> vec_negatives;
    std::vector<std::vector<double>> vec_cdf;
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
    StackedFilter<CountingBloomFilterLayer, StringElement> * filterArr[datasetNum];

    size_t bit_Per_Counter =4;
    for (size_t i = 0; i < datasetNum; i++)
    {
        filterArr[i] =  new StackedFilter<CountingBloomFilterLayer, StringElement>(bits_per_key_sets[i] * vec_positives[i].size(), vec_positives[i], vec_negatives[i], vec_cdf[i], bit_Per_Counter,vulNegKeyRatio);
    }
    
    for (size_t i = 0; i < datasetNum; i++)
    {
        for (StringElement &s : vec_positives[i])
        {
            filterArr[i]->InsertPositiveElement(s);
        }       
    }

    double sumFPR = 0;
    double sumWeightedFPR = 0;
    double sumFNR = 0;
    double maxFPR = 0, minFPR = MAXFLOAT, maxWeightedFPR = 0, minWeightedFPR = MAXFLOAT;
    
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0 ;
        size_t fpr_Num = 0;
        double total_cost_ = 0;
        double wfpr_c = 0;
        
        for(int j=0; j<vec_positives[i].size(); ++j)
        {
            if(!filterArr[i]->LookupElement(vec_positives[i][j]))
            {
                fnr_c++;
            }
        }
           
        for(int j=0; j<vec_negatives[i].size(); ++j)
        {
            if(filterArr[i]->LookupElement(vec_negatives[i][j]))
            {
                wfpr_c += vec_cdf[i][j];
                fpr_Num++;
            }
            total_cost_ += vec_cdf[i][j];
        }

        double tempFPR = (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();
        sumFPR +=tempFPR;
        maxFPR = tempFPR > maxFPR ? tempFPR : maxFPR;
        minFPR = tempFPR < minFPR ? tempFPR : minFPR; 

        double tempWeightedFPR = (double)(wfpr_c*100) / total_cost_;
        sumWeightedFPR += tempWeightedFPR;
        maxWeightedFPR = tempWeightedFPR > maxWeightedFPR ? tempWeightedFPR : maxWeightedFPR;
        minWeightedFPR = tempWeightedFPR < minWeightedFPR ? tempWeightedFPR : minWeightedFPR;  

        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();    
    }

    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - minWeightedFPR - maxWeightedFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)(sumFPR - minFPR - maxFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"Average FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;

    for (size_t i = 0; i < datasetNum; i++) delete filterArr[i];
    
}


void TestNaiveFilter(dataloader * dlSet, double * bits_per_key_sets)
{
    naiveFilter::NaiveFilter * nfArr[datasetNum];
    double usedNegRatio = 0.12;
    for (size_t i = 0; i < datasetNum; i++)
    {
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            std::cout<<"error!\n";
        }
        nfArr[i] = new naiveFilter::NaiveFilter(bits_per_key_sets[i], dlSet[i].pos_keys_.size(), usedNegRatio, dlSet->neg_keys_);
    }


    for (size_t i = 0; i < datasetNum; i++)
    {
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            nfArr[i]->Add( *v);

        }     
    }

    double sumFPR = 0;
    double sumWeightedFPR = 0;
    double sumFNR = 0;
    double maxFPR = 0, minFPR = MAXFLOAT, maxWeightedFPR = 0, minWeightedFPR = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;

        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if(!nfArr[i]->Query(*v)){
                fnr_c++;
            } 
        }
        
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            if(nfArr[i]->Query(*v))
            {
                wfpr_c += v->cost;
                fpr_Num ++;
            }  
            total_cost_ += v->cost;
        }

        double tempFPR = (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();
        sumFPR +=tempFPR;
        maxFPR = tempFPR > maxFPR ? tempFPR : maxFPR;
        minFPR = tempFPR < minFPR ? tempFPR : minFPR; 

        double tempWeightedFPR = (double)(wfpr_c*100) / total_cost_;
        sumWeightedFPR += tempWeightedFPR;
        maxWeightedFPR = tempWeightedFPR > maxWeightedFPR ? tempWeightedFPR : maxWeightedFPR;
        minWeightedFPR = tempWeightedFPR < minWeightedFPR ? tempWeightedFPR : minWeightedFPR;  

        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Weighted FPR:"<< (double)(sumWeightedFPR - minWeightedFPR - maxWeightedFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"FPR: "<< (double)(sumFPR - minFPR - maxFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;

    for (size_t i = 0; i < datasetNum; i++)
    {

        delete nfArr[i];
    }
}

void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf)
{
    for (size_t i = 0; i < datasetNum; i++)
    {
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            std::cout<<"Don't have been sorted\n";
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
