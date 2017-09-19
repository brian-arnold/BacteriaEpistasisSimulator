#include "ISM.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/shared_ptr.hpp>

using namespace std;
extern MTRand rnd;

double Fitness(Hap *Ind, double &RNS)
{
    double w = 1.0 ;
    int size = (*(Ind)).Muts.size()-1 ;
    while( size > -1 ){
        if( (*(Ind)).Muts[size]->UnderSeln ){
            w*=RNS ;
        }
        size-- ;
    }
    return w ;
}

void HomoRecomb(Hap &NewHaplo, Hap *Recipient, Hap *Donor, int &NumR, double &GeoP, int &L)
{

    int i,j,k,size ;
    double Start, TractLen, End ;
    int NumBr ; // Number of breakpoints
    int RecipSize = (*(Recipient)).Muts.size();
    int DonSize = (*(Donor)).Muts.size();
    vector< vector<double> > HRs(NumR, vector<double>(2) ) ; // initiate a vector of size NumR, each element being a vecotr of 2 dbls
    vector<int> HRsToDelete ;
    vector<double> PosHR ;
    ///////////
    /*
    cout << "REC_EVENT " << NumR << " Recs\n" ;
    cout << "RECIPIENT\n" ;
    for(i=0; i<RecipSize; i++){
        cout << (*(Recipient)).Muts[i] -> pos << " " ;
    }
    cout << "\n" ;
    cout << "DONOR\n" ;
    for(i=0; i<DonSize; i++){
        cout << (*(Donor)).Muts[i] -> pos << " " ;
    }
    cout << "\n" ;
     */
    ////////////
    for(i=0; i<NumR; i++){
        //RANDOMLY SELECT POSITION + TRACT LENGTH
        Start = rnd.rand() ; //
        TractLen = double(geometric(GeoP))/(L) ; //
        if(TractLen > 0.0){
            End = Start+TractLen ;
            HRs[i][0] = Start ;
            HRs[i][1] = End ;
            //HRs[i].push_back( Start ) ;
            //HRs[i].push_back( End ) ;
        }
    }
    // Fuse overlapping HRs by either
    // (1) extending ones with partially overlapping breakpoints or
    // (2) removing HR completely nested within another
    
    //cout << "REC EVENTS ORIGINAL\n" ;
    //for(i=0; i<HRs.size(); i++){
    //    cout << HRs[i][0] << " " << HRs[i][1] << "\t" ;
    //}
    //cout << "\n" ;
    
    if(NumR>1){
        // HRs NEEDS TO BE SORTED by STARTS!!!! makes things simpler below
        std::sort(HRs.begin(), HRs.end(), compareVectorOfVectors );
        
        //cout << "REC EVENTS SORTED\n" ;
        //for(i=0; i<HRs.size(); i++){
        //    cout << HRs[i][0] << " " << HRs[i][1] << "\t" ;
        //}
        //cout << "\n" ;
        size = HRs.size() ; // num rec events
        for(i=0; i<size-1; i++){
            HRsToDelete.clear() ;
            for(j=i+1; j<size; j++){
                // (1) is 2nd brkpt of i nested within j
                if( (HRs[i][1] > HRs[j][0]) && (HRs[i][1] < HRs[j][1]) ){
                    // fuse, with start of i end of j, remove j LATER o.w. seg faults
                    HRs[i][1] = HRs[j][1] ;
                    HRsToDelete.push_back(j) ;
                // (2) is HR j completely nested within HR i
                }else if( (HRs[i][0] < HRs[j][0]) && (HRs[i][1] > HRs[j][1]) ){
                    HRsToDelete.push_back(j) ;
                }
            }
            // remove all interfering HRs here
            if( !HRsToDelete.empty() ){
                for(k=HRsToDelete.size()-1; k>=0; k--){// delete in reverse to avoid altering indexing as you go
                    HRs.erase(HRs.begin() + HRsToDelete[k] ) ;
                }
            }
            // update size here, after going through all j's
            size = HRs.size() ;
        }
    }
    //cout << "REC EVENTS AFTER DELETION\n" ;
    //for(i=0; i<HRs.size(); i++){
    //    cout << HRs[i][0] << " " << HRs[i][1] << "\t" ;
    //}
    //cout << "\n" ;

    for(i=0; i<HRs.size(); i++){
        PosHR.push_back( HRs[i][0] ) ;
        PosHR.push_back( HRs[i][1] ) ;
    }
    NumBr = PosHR.size() ;
    ///////////////
    //cout << "BREAKPOINTS\t" ;
    //for(i=0; i<NumBr; i++){
    //    cout << PosHR[i] << "\t" ;
    //}
    //cout << "\n" ;
    ///////////////
    
    
    int s1 = 0 ;
    int s2 = 0 ;
    for (i = 0; i < NumBr; i++){
        // % operator extracts remainder from a division, alternates between 0,1 starting with 0 when i=0
        // this switch regultes whether we're pushing on Recipient Muts or Donor muts
        if (i % 2 == 0){
            while ((s1 < RecipSize) && ( (*(Recipient)).Muts[s1]->pos < PosHR[i])) // starting from begin of Recip
            {
                NewHaplo.Muts.push_back( (*(Recipient)).Muts[s1] );
                s1++;
            }
            while ((s2 < DonSize) && ( (*(Donor)).Muts[s2]->pos < PosHR[i])){
                s2++; // since Muts from this region are from Recip, increment s2 to ignore them from Donor
            }
        }else{
            while ((s2 < DonSize) && ( (*(Donor)).Muts[s2]->pos < PosHR[i]))
            {
                NewHaplo.Muts.push_back( (*(Donor)).Muts[s2] );
                s2++;
            }
            while ((s1 < RecipSize) && ( (*(Recipient)).Muts[s1]->pos < PosHR[i])){
                s1++;
            }
        }
    }
    // filled up through last HR event, fill rest with recipient DNA
    while (s1 < RecipSize){
        NewHaplo.Muts.push_back( (*(Recipient)).Muts[s1] ) ;
        s1++;
    }
    
    //cout << "NEW RECOMBINANT HAPLOTYPE\n" ;
    //for(i=0; i<NewHaplo.Muts.size(); i++){
    //    cout << NewHaplo.Muts[i] -> pos << " " ;
    //}
    //cout << "\n" ;
    
}

ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


void ReadParameterFile(int &Loc, int &Ind, int &Gen, double &RecRate, double &GeoTrLen, double &MutRate, double &NScoeff, int &NmSm)
{
    ifstream fin("Parameters") ;
    GotoLine(fin, 2);

    fin >> Loc  ;
    fin >> Ind  ;
    fin >> Gen  ;
    fin >> RecRate  ;
    fin >> GeoTrLen  ;
    fin >> MutRate ;
    fin >> NScoeff ;
    fin >> NmSm  ;

    fin.close() ;
}

bool compareByPointerToPos(const SharedMutPtr a, const SharedMutPtr b)
{
    return a->pos < b->pos;
}
bool comparePos (const SNPAF a, const SNPAF b)
{
    return a.pos < b.pos;
}
bool compareVectorOfVectors(const std::vector< double >& a, const std::vector< double >& b)
{
    return a[0] < b[0];
}


