# pygcta
Python implementation of GCTA

Original paper
[Jian et. al](http://www.cell.com/ajhg/abstract/S0002-9297(10)00598-7)

---

Input Data

* Genotypes (plink format)
* Phenotypes

Functions
- Preprocessing
  -  Filter Genotypes (missing values) [Genotype class]
  -  Normalize phenotypes [Phenotype Class (n x g)]
  -  Matching function
  -  Imputation for missing genotypes
- Kernel Cacluation
- Linear Model Function - constructor requires genotype and phenotype
  - Derivative of the likelihood function
  - Optimization function
  - GCTA class

Output Data

TODOs

* Which python modules to import? 

* D and K set up class framework with existing functions

