#include "AddEpi.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace std;
extern MTRand rnd;


void WrightFisherEvo(int &Loci, int &Individuals, int &Generations, double &SampGen, double &RecombinationRate, double &GeomTractLength, double &Va, double &Vi, int &GenomeSize, int &Sim, string &FilePath)
{
    
    int site,ind,i,j,gen,k ; // incrementors
    int RandInd, RandInd2, STOP ;
    double Wmax, ThetaPiStart, ThetaPi, Wmean, Wvar, Wskew, tmp ;
    // Recombination related
    int NumRecs, Start, Start2, TractLen, End, RecTotal, Free ;
    // Frequently used constants
    //double L = double(Loci) ;
    double AddFitConst = Va ;
    int twoN = 2*Individuals ;
    // selection coefficient sqrt(Va/L)
    //with additive fitness these need to be small to avoid negative absolute fitnesses
    int Positions[Loci] ;
    vector< vector<int> > RecombinantPositions ;
    vector<int> TempDonorSNPs ;
    
    vector <int> PrintingPoints ;
    int PPinc = int(Individuals*SampGen) ; // sampling intervals, in fractions of N generations
    int PP = PPinc ;
    //cout << PP << "\n" ;
    while (PP <= Generations){
        PrintingPoints.push_back(PP) ;
        PP+=PPinc ;
        //cout << PP << "\n" ;
    }
    int CurrentPrintingPoint = 0 ;
    int MaxPrintingPoint = (PrintingPoints.size() - 1) ;
    
    /*============
     DECLARE DYNAMIC VARIABLES
     MAKE SURE TO DELETE! :)
     ==============*/
    Hap * Haplotypes = new Hap [twoN] ; // array of Hap structs: vector<Locus> Loci ; unsigned int NumHap ;
    Hap ** Pop = new Hap *[Individuals] ; // array of pointers to Hap structs
    Hap ** NextGen = new Hap *[Individuals] ;  // array of pointers to Hap structs
    double * WAbs = new double [Individuals]; // "WAbs" will hold the Absolute fitness of each individual:
    double * EpistasisTable = new double [ ((Loci*(Loci-1))/2) ] ;  // array of pairwise epistatic selection coefficients
    vector<double> * DprimeByDist = new vector<double> [ (Loci-1) ] ; // array of vector<double>
    vector<double> * DqleByDist = new vector<double> [ (Loci-1) ] ;  // array of vector<double>
    vector<double> * DeltarByDist = new vector<double> [ (Loci-1) ] ;  // array of vector<double>
    vector<TwoLocusAFS> * PairwiseAFS = new vector<TwoLocusAFS> [ (Loci-1) ] ; // array of vectors of TwoLocusAFS: int one ; int two ;
    vector<XYZUfreqs> * TwoLocusHaploFreqsBS = new vector<XYZUfreqs> [ (Loci-1) ] ; // BS = before selection
    vector<XYZUfreqs> * TwoLocusHaploFreqsAS = new vector<XYZUfreqs> [ (Loci-1) ] ; // AS = after selection
    
    string TempSimFileName = "Results" + boost::lexical_cast<string>(Sim) ;
    ofstream out ;
    out.open(TempSimFileName.c_str()) ;
    
    out.precision(20) ;
    //out2.precision(20) ;
    
    // START CLOCK
    time_t begin, end;
    begin = time(0);
    
    
    /*============
     INITIALIZATION
     ==============*/
    STOP = 0 ;
    InitializeHaplotypes(Individuals, Loci, Haplotypes, Sim, FilePath) ; //NumHap initialized later
    InitializePositions(Positions, Loci, Sim, FilePath) ;
    InitializeEpistasisTable(Vi, Loci, EpistasisTable) ;
    
    
    out << "POSITIONS\n" ;
    for(i=0; i<(Loci); i++){
        out << Positions[i] << " " ;
    }
    out << "\n" ;
    tmp = (AddFitConst*Individuals*2) ;
    out << "AddFitConst_2Ns: " << tmp << "\n" ;
    
    out << "EpistasisTable_2Ns\n" ;
    // I'm pretty sure this print epistasis values according to distance, with the top two being adjacent loci
    for(i=0; i<(Loci-1); i++){
        out << (2*Individuals*EpistasisTable[i]) << " " ;
        k=0 ;
        for(j=(Loci-1); j>(i+1); j--){
            k+=j ;
            out << (2*Individuals*EpistasisTable[i+k]) << " " ;
        }
        out << "\n" ;
    }
    out << "\n" ;

    // Assign Individuals to Haplotypes, initiate number in population
    for(ind=0; ind<Individuals; ind++){
        Pop[ind] = &Haplotypes[ind] ;
        (*(Pop[ind])).NumHap = 1 ;
    }
    for(ind=Individuals; ind<twoN; ind++){
        Haplotypes[ind].NumHap = 0  ;
    }

    
    //CALCULATE STARTING METRICS
    CalculateLDstart(Pop, Individuals, Loci, DprimeByDist) ;
    out << "Dprime_BY_DIST_gen0" << "\n" ;
    for(i=0; i<(Loci-1); i++){
        for(j=0; j<DprimeByDist[i].size() ;j++){
            out << DprimeByDist[i][j] << " " ;
        }
        out << "\n" ;
    }
    out << "\n" ;
    
    ThetaPiStart = CalculateThetaPi(Pop, Individuals, Loci) ;
    out << "ThetaPi_gen0" << "\t" << ThetaPiStart << "\n" ;

    for(ind=0; ind<Individuals; ind++){
        // 1+MarginalFitness as in Piganeau et al 2001
        WAbs[ind] = 1.0 + MarginalFitness(Pop[ind], Loci, AddFitConst, EpistasisTable) ;
    }
    Wmean = MeanFitness(Individuals, WAbs) ;
    Wvar = VarFitness(Individuals, WAbs, Wmean) ;
    Wskew = SkewFitness(Individuals, WAbs, Wmean, Wvar) ;
    out << "Wmean_gen0" << "\t" << Wmean << "\n" ;
    out << "Wvar_gen0" << "\t" << Wvar << "\n" ;
    out << "Wskew_gen0" << "\t" << Wskew << "\n" ;
    
    out << "\n" ;
    //===========================
    // BEGIN WRIGHT-FISHER EVOLUTION
    // ===========================

    for (gen = 1; gen <= Generations; gen++){
        
        Wmax = 0.0 ;
        // CALCULATE FITNESS HERE
        for(ind=0; ind<Individuals; ind++){
            // 1+MarginalFitness as in Piganeau et al 2001
            WAbs[ind] = 1.0 + MarginalFitness(Pop[ind], Loci, AddFitConst, EpistasisTable) ;
            if (Wmax < WAbs[ind]){
                Wmax = WAbs[ind] ;
            }
        }
        
        // CREATE NEXT GENERATION
        for(ind=0; ind<Individuals; ind++){
            do{
                RandInd = int(rnd.randExc(Individuals)) ; // real number in [0,Nv), int -> includes 0 but not Nv
            } while (rnd.rand() > WAbs[RandInd]/Wmax);// Randomly sample individual based on it's relative
            // Need to sample from Pop to simulate drift, not the original population Haplotypes

        
            NumRecs = int(poisdev(RecombinationRate)) ; // this ind is a recipient of DNA
            //cout << "NumRecs: " << NumRecs << "\n" ;
            if(NumRecs>0){
                RecombinantPositions.clear() ;
                // check how many are effective recombinations
                for(i=0; i<NumRecs; i++){
                    TempDonorSNPs.clear() ;
                    Start = int(rnd.randExc(GenomeSize)) ;
                    do{
                        TractLen = geometric(GeomTractLength) ;
                    }while( TractLen >= GenomeSize-1 ) ;
                    End = (Start+TractLen) ;
                    if( End <= (GenomeSize-1)){
                        j=0 ;
                        for(j=0; j<Loci; j++){
                            // recombination happens to right of Start, so only consider Positions to right
                            if(Start <= Positions[j]){
                                while( (Positions[j] <= End) && (j<Loci) ){
                                    // push j into vector, j++
                                    TempDonorSNPs.push_back( j ) ;
                                    j++ ;
                                }
                                break ;
                            }
                        }
                    }else{
                        j=0 ;
                        for(j=0; j<Loci; j++){
                            // recombination happens to right of Start, so only consider Positions to right
                            if(Start <= Positions[j]){
                                while( (Positions[j] <= (GenomeSize-1)) && (j<Loci) ){ // Positions only include up to GenomeSize-1
                                    // push j into vector, j++
                                    TempDonorSNPs.push_back( j ) ; // record j, not Positions, b/c 0-(L-1) loci
                                    j++ ;
                                }
                                break ;
                            }
                        }
                        k=0 ;
                        while( (Positions[k] <= (End-(GenomeSize-1))) && (k<Loci)){
                            TempDonorSNPs.push_back( k ) ;
                            k++ ;
                        }
                        
                    }
                    if( !TempDonorSNPs.empty() ){
                        RecombinantPositions.push_back( TempDonorSNPs ) ;
                    }
                }
                
                if( !RecombinantPositions.empty() ){
                    // Find unoccupied haplotype
                    // replace it with this recomb recipient, set NumHap = 1
                    Free=0 ;
                    while( (Haplotypes[Free].NumHap > 0) && (Free < twoN)){
                        Free++ ;
                    }
                    
                    Haplotypes[Free].Loci = (*(Pop[RandInd])).Loci ; // Copy recipient haplotype to recombine
                    NextGen[ind] = &Haplotypes[Free] ; // reassign Pop pointer to recipient haplotype
                    Haplotypes[Free].NumHap = 1 ; // was previously 0
                    
                    
                    // RANDOMLY SELECT A DONOR, 2nd random ind
                    do{
                        RandInd2 = int(rnd.randExc(Individuals)) ; // real number in [0,Nv), int -> includes 0 but not Nv
                    } while (RandInd2 == RandInd) ; // Cant select same individual!!
                    
                    //As it stands, multiple recomb events can have same start position and interfere with eachother
                    for(i=0; i<RecombinantPositions.size(); i++){
                        HomoRecomb( NextGen[ind], Pop[RandInd2], RecombinantPositions[i]) ;
                    }
                }else{
                    NextGen[ind] = Pop[RandInd] ;
                }
            }else{
                NextGen[ind] = Pop[RandInd] ;
            }
        
        }
        
        // Reset haplotype counts
        for(ind=0; ind<twoN; ind++){
            Haplotypes[ind].NumHap = 0 ;
        }
        
        //Move population to next generation, pointer reassignment
        for(ind=0; ind<Individuals; ind++){
            Pop[ind] = NextGen[ind] ;
            (*(Pop[ind])).NumHap++ ;
        }
        ThetaPi = CalculateThetaPi(Pop, Individuals, Loci) ;
        // Calculate stats right after selection
        if(STOP == 1){ // starts as 0
            //Calculate statistics every SampGen
            //LD, Wbar, Wmax, diversity actual Va and Vi?
            
            CalcTwoLocusHaplotypeFreqs(Pop, Individuals, Loci, TwoLocusHaploFreqsAS) ;
            CalcChangeHaploFreqsDeltar(Pop, Individuals, Loci, TwoLocusHaploFreqsBS, TwoLocusHaploFreqsAS, DeltarByDist) ;
            
            Wmean = MeanFitness(Individuals, WAbs) ;
            Wvar = VarFitness(Individuals, WAbs, Wmean) ;
            Wskew = SkewFitness(Individuals, WAbs, Wmean, Wvar) ;
            ThetaPi = CalculateThetaPi(Pop, Individuals, Loci) ;
            out << "ThetaPi_gen" << gen << "\t" << ThetaPi << "\n" ;
            out << "Wmean_gen" << gen << "\t" << Wmean << "\n" ;
            out << "Wmax_gen" << gen << "\t" << Wmax << "\n" ;
            out << "Wvar_gen"<< gen << "\t" << Wvar << "\n" ;
            out << "Wskew_gen"<< gen << "\t" << Wskew << "\n" ;
           
            out << "\n" ;
            CalculatePairwiseAFS(Pop, Individuals, Loci, PairwiseAFS) ;
            
            CalculateLDmetrics(Pop, Individuals, Loci, TwoLocusHaploFreqsBS, DprimeByDist, DqleByDist) ;
            out << "Dprime_BY_DIST_gen"<< gen << "\n" ;
            for(i=0; i<(Loci-1); i++){
                for(j=0; j<DprimeByDist[i].size() ;j++){
                    out << DprimeByDist[i][j] << " " ;
                }
                out << "\n" ;
            }
            out << "\n" ;
             
            out << "PairwiseAFS_BY_DIST_gen"<< gen << "\n" ;
            for(i=0; i<(Loci-1); i++){
                for(j=0; j<PairwiseAFS[i].size() ;j++){
                    out << PairwiseAFS[i][j].one << "," << PairwiseAFS[i][j].two << " " ;
                }
                out << "\n" ;
            }
            out << "\n" ;
            if(gen > PrintingPoints[MaxPrintingPoint]){ // exit loop b/c you've printed everything you need
                break ;
            }else{
                STOP = 0 ; // set it back to 0 to print more generations of evolution
            }
        }
        if(gen == ((PrintingPoints[CurrentPrintingPoint])-1) ){// set STOP to 1 generation before you wanna print
            CalcTwoLocusHaplotypeFreqs(Pop, Individuals, Loci, TwoLocusHaploFreqsBS) ;
            STOP = 1 ;
            CurrentPrintingPoint++ ;
            
        }
        
            }
    //=========================
    //END WRIGHT-FISHER EVOLUTION
    //=========================
    
   
    
    end = time(0);
    int TimeTaken = int(difftime(end,begin)) ;
    out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    
    
    out.close() ;
    
    delete [] Haplotypes ;
    delete [] Pop ;
    delete [] NextGen ;
    delete [] WAbs ;
    delete [] EpistasisTable ;
    delete [] DprimeByDist ;
    delete [] DqleByDist ;
    delete [] DeltarByDist ;
    delete [] TwoLocusHaploFreqsBS ;
    delete [] PairwiseAFS ;
    
}
