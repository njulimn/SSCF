# Seesaw Counting Filter

# About this repo
This repo contains all the source codes of SSCF together with other comparison algorithms, which are detailed in the following table.

|Algorithm| Description|
|----|----|
|SSCF|[[Implemention]](./SSCF/SeesawCF.h)|
|Counting Bloom filter|Li Fan, “summary cache: a scalable wide-area web cache sharing protocol,” Transactions on Networking. [[Implementation]](/nonlearnedfilter/countingBloom.h)|
|WCBF|J. Bruck, J. Gao, and A. Jiang, “Weighted Bloom filter,” in Proceedings of International Symposium on Information Theory. IEEE, 2006. [[Implementation]] (/nonlearnedfilter/wcbf.h)|
|SF|K. Deeds, B. Hentschel, S. Idreos, “Stacked filters: learning to filter by structure,” in Proceedings of International Conference on Very Large Data Bases. 2020, [[Implementation]](/nonlearnedfilter/stackedfilter/)|

# Requirement 
   1. cmake@3+
   2. make
   3. recommended compiler: g++ or clang
# How to build
```bash
mkdir -p build && cd build && cmake .. && make
```
# How to run the benchmark
```bash
./experiment_weightedFPR
./experiment_latency
```


