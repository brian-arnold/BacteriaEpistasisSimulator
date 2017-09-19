#ifndef ADDEPI_H
#define ADDEPI_H

#include <stdio.h>
#include <vector>
#include <iostream>
#include "MersenneTwister.h"
#include <fstream>
#include <limits>
#include <string>

//#include <boost/shared_ptr.hpp>
using namespace std;

/* ============
STRUCTURES
==============*/
struct Locus{
    // 2 bit signed integer, can only take values of 1, -1, and 0
    signed int Eff : 2 ;
};

struct Hap{
    vector<Locus> Loci ;
    unsigned int NumHap ;
};
struct XYZUfreqs{
    int x ;
    int y ;
    int z ;
    int u ;
};
struct TwoLocusAFS{
    int one ;
    int two ;
};

/* ============
FUNCTION PROTOTYPES
==============*/
void InitializeHaplotypes(int N, int L, Hap * H, int sim, string &FlPth) ;
void InitializePositions(int * Pos, int L, int sim, string &FlPth) ;
void InitializeEpistasisTable(double Vint, int L, double * ET) ;
double MarginalFitness(Hap * p, int &L, double &VaFit, double * ET) ;
void HomoRecomb(Hap * Recipient, Hap * Donor, vector<int> &RecPos) ;
void SexLikeRecomb(Hap * Recipient, Hap * Donor, int &St, int &St2, int &L) ;
void CalcTwoLocusHaplotypeFreqs(Hap * p[], int &I, int &L, vector<XYZUfreqs> * TLHFBS) ;
void CalculateLDstart(Hap * p[], int &I, int &L, vector<double> * DprimeBD) ;
void CalcChangeHaploFreqsDeltar(Hap * p[], int &I, int &L, vector<XYZUfreqs> * TLHFBS, vector<XYZUfreqs> * TLHFAS,  vector<double> * Deltar);
void CalculateLDmetrics(Hap * p[], int &I, int &L, vector<XYZUfreqs> * TLHF, vector<double> * DprimeBD, vector<double> * DQLEBD) ;
double CalculateThetaPi(Hap * p[], int &I, int &L) ;
double BinomCoeff(double n, double k) ;
fstream& GotoLine(std::fstream& file, unsigned int num) ;
void ReadParameterFile(int &Loc, int &Ind, int &Gen, double &SmpGn, double &RecRate, double &GeoTrLen, double &Vadd, double &Vint, int &GenSize, int &NmSm, string &FlPth) ;
void WrightFisherEvo(int &Loci, int &Individuals, int &Generations, double &SampGen, double &RecombinationRate, double &GeomTractLength, double &Va, double &Vi, int &GenomeSize, int &Sim, string &FilePath) ;
double MeanFitness(int &I, double * Warray) ;
double VarFitness(int &I, double * Warray, double &Wm) ;
double SkewFitness(int &I, double * Warray, double &Wm, double &Wv) ;
void CalculatePairwiseAFS(Hap * p[], int &I, int &L, vector<TwoLocusAFS> * PWAFS) ;


double poisdev(const double xm) ;
double gammln(const double xx) ;
double binldev(const double pp, const int n) ;
double gasdev() ;
int geometric(const double &p) ;

#endif 
