#ifndef ISM_H
#define ISM_H

#include <stdio.h>
#include <vector>
#include <iostream>
#include "MersenneTwister.h"
#include <fstream>
#include <limits>
#include <boost/shared_ptr.hpp>

//#include <boost/shared_ptr.hpp>
using namespace std;

/* ============
STRUCTURES
==============*/
struct Mutation
{
    double pos; // position of mutation
    int UnderSeln;
};
typedef boost::shared_ptr<Mutation> SharedMutPtr ;
struct Hap
{
    vector<SharedMutPtr> Muts;  	// neutral locus (infinite sites model, each value in the vector corresponds to one mutation)
    int NumHap;   // number of copies of this genome segment in the population
};
struct SNPAF{
    double pos ;
    int AF ;
};
/* ============
FUNCTION PROTOTYPES
==============*/

fstream& GotoLine(std::fstream& file, unsigned int num) ;
void ReadParameterFile(int &Loc, int &Ind, int &Gen, double &RecRate, double &GeoTrLen, double &MutRate, double &NScoeff, int &NmSm) ;
void WrightFisherEvo(int &Loci, int &Individuals, int &Generations, double &RecombinationRate, double &GeomTractLength, double &MutRate, double &NScoeff, int &Sim) ;
bool compareByPointerToPos(const SharedMutPtr a, const SharedMutPtr b);
bool comparePos (const SNPAF a, const SNPAF b) ;
void HomoRecomb(Hap &NewHaplo, Hap * Recipient, Hap * Donor, int &NumR, double &GeoP, int &L) ;
bool compareVectorOfVectors(const std::vector< double >& a, const std::vector< double >& b) ;
double Fitness(Hap *Ind, double &RNS) ;

double poisdev(const double xm) ;
double gammln(const double xx) ;
double binldev(const double pp, const int n) ;
double gasdev() ;
int geometric(const double &p) ;

#endif 
