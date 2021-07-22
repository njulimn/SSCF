//used to Evaluation the experiment (Cost-weighted FPR vs N(n)/N(p) with 1.0 skewness)
//Please select the Shalla dataSet with a slope of 1 before evaluate
#include <time.h>
#include "util/dataloader.h"
#include "nonlearnedfilter/wbf.h"
#include "SSCF/SeesawCF.h"
#include "nonlearnedfilter/countingBloom.h"
#include "nonlearnedfilter/wcbf.h"
#include "nonlearnedfilter/stackedfilter/StackedFilter.h"
#include <vector>

#define SHALLA_NAME "shalla_"
#define YCSBT_NAME "ycsbt_"
#define SHALLA_PATH_T  "../data/dataset_shalla/shalla_1_0/shalla_cost_1.0-"
const size_t datasetNum = 5; 

using namespace std;
using namespace std::chrono;

void TestSSCF(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets);
void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf);
bool compare_func(Slice * key1, Slice * key2);

int main(){
    double bitsPerItem = 28; 
    dataloader * dlSet = new dataloader[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        //load hosefire
        // if(!dlSet[i].loadHosefire(true, to_string(i+1)))
        // {
        //     std::cout<< "error!\n";
        // }
        
        //load shalla
        if(!dlSet[i].loadShalla_vulNeg(SHALLA_PATH_T, true, to_string(i)))
        {
            std::cout<< "error!\n";
        }

        //load ycsbt
        // if(!dlSet[i].loadYCSB(true, to_string(i)))
        // {
        //     std::cout<< "error!\n";
        // };
    }
    double *bits_per_key_sets = new double[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        bits_per_key_sets[i] = bitsPerItem;

    }

    for(size_t i = 0; i<datasetNum; i++){
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);
        }
    }
    
    double vulNegKeyRatio = 0.02;
    std::cout << "----------------- SeesawCF testing -------------------" << endl;
    for (size_t i = 0; i < 15; i++)
    {
        std::cout << "----- The Ratio of vulnerable negative keys: "<<vulNegKeyRatio<<" ------" << endl;
        TestSSCF(dlSet, bits_per_key_sets, vulNegKeyRatio);
        vulNegKeyRatio+=0.02;
    }
    
    std::cout << "----------------- countingBloom testing -------------------" << endl;
    TestCountingBloom(dlSet, bits_per_key_sets);



    
}

void TestSSCF(dataloader * dlSet,  double * bits_per_key_sets, double vulNegKeyRatio)
{
    SeesawCF::SeesawCF *sscfArr[datasetNum];
    double cbfSpaceRatio = 0.82;
    for (size_t i = 0; i < datasetNum; i++)
    {
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            std::cout<<"error!\n";
        }
        sscfArr[i] = new SeesawCF::SeesawCF(vulNegKeyRatio, bits_per_key_sets[i], dlSet[i].pos_keys_.size(), dlSet[i].neg_keys_.size(), dlSet[i].neg_keys_ , cbfSpaceRatio) ;
    }
    
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t sumSize = bits_per_key_sets[i] * dlSet[i].pos_keys_.size();
        size_t haveInsertNum = 0;
        size_t bits_per_item_ite = 50;

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
    
    std::cout<<"Weighted FPR: "<< (double)(sumWeightedFPR - minWeightedFPR - maxWeightedFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"FPR: "<< (double)(sumFPR - minFPR - maxFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;
    
    // double sumHaveInsertNumHash =0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     sscfArr[i]->PrintHE_haveInsertNum();
    //     sumHaveInsertNumHash += sscfArr[i]->GetHaveInsertNum() / (double)sscfArr[i]->GetCellNum();
    // }
    // std::cout<<sumHaveInsertNumHash / 5.0<<endl;

    for (size_t i = 0; i < datasetNum; i++) delete sscfArr[i] ;
}

void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets){

    size_t nfuncs ;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[5];
    size_t sumSize = bits_per_key_sets[0] * dlSet[0].pos_keys_.size();
    size_t numCounter = sumSize / bitPerCounter;
    for (size_t i = 0; i < datasetNum; i++)
    {
        nfuncs = bits_per_key_sets[i] / 4 * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] * dlSet[0].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
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
    std::cout<<"FPR: "<< (double)(sumFPR - minFPR - maxFPR) / (datasetNum - 2) << "%" <<endl;
    std::cout<<"FNR: "<< (double)sumFNR / datasetNum << "%" <<endl;
    for (size_t i = 0; i < datasetNum; i++) delete cbfArr[i];
}


void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf)
{
    for (size_t i = 0; i < datasetNum; i++)
    {
        sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);/////////////是有问题的
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
