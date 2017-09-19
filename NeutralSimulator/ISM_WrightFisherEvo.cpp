#include "ISM.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
extern MTRand rnd;


void WrightFisherEvo(int &Loci, int &Individuals, int &Generations, double &RecombinationRate, double &GeomTractLength, double &MutRate, double &NScoeff, int &Sim)
{
    /*============
     CONSTANTS
     ==============*/
    int site,ind,i,j,gen,k,last,size ; // incrementors
    int RandDonor, Par, Free, NumMuts, NumRecs ;
    double Wbar, Wmax, Wind ;
    int twoN = (2*Individuals) ;
    double RelNScoeff = (1.0-NScoeff) ;
    
    vector<SNPAF> SnpAlleleFreqs ;
    SNPAF tempSNP ;
    vector<SNPAF> SnpPolymorphic ;
    /*============
     DECLARE DYNAMIC VARIABLES
     MAKE SURE TO DELETE! :)
     ==============*/
    Hap * Haplotypes = new Hap [twoN] ;
    Hap ** Pop = new Hap *[Individuals] ;
    Hap ** NextGen = new Hap *[Individuals] ;
    Hap * PopSamp = new Hap [Individuals] ;
    int * Parents = new int [Individuals] ;
    double * Wtot = new double [Individuals]; // "Wtot" will hold the fitness of each individual:

    string TempSimFileName = "Results" + boost::lexical_cast<string>(Sim);
    ofstream out(TempSimFileName.c_str());

    // START CLOCK
    time_t begin, end;
    begin = time(0);
    
    /*============
     INITIALIZATION
     ==============*/
    
    // Initialize Haplotypes
    Haplotypes[0].NumHap = Individuals ;
    for(i=1; i < twoN; i++){
        Haplotypes[i].NumHap = 0 ;
    }
    //Initialize Population
    for(i=0; i < Individuals ; i++){
        Pop[i] = &Haplotypes[0] ;
    }
    
    
    //CALCULATE STARTING METRICS

 
   
    /*===========================
     BEGIN WRIGHT-FISHER EVOLUTION
     ===========================*/
    for (gen = 1; gen <= Generations; gen++){
        Wbar = 0.0 ;
        Wmax = 0.0 ;
        
        // Calculate fitness of each individual
        for(ind=0; ind<Individuals; ind++){
            Wind = Fitness(Pop[ind], RelNScoeff) ;
            Wtot[ind] = Wind ;
            
            Wbar += Wind ;
            if(Wmax < Wind){
                Wmax = Wind ;
            }
        }
        Wbar /= Individuals ;
        
        
        // CREATE NEXT GENERATION
        for(ind=0; ind<Individuals; ind++){
            // Neutral, no fitness costs
            do{
                Par = int(rnd.randExc(Individuals)) ; // real number in [0,Nv), int -> includes 0 but not Nv
            }while( rnd.rand() > Wtot[Par]/Wmax ) ;
            Parents[ind] = Par ;
        }
        // MUTATE + RECOMBINE HAPLOTYPES
        Free = 0 ;
        for(ind=0; ind<Individuals; ind++){ // foreach ind of next generation
            NumMuts = int(poisdev(MutRate)) ; // this ind is a recipient of DNA
            NumRecs = int(poisdev(RecombinationRate)) ; // this ind is a recipient of DNA
        
            if(NumMuts==0 && NumRecs==0){
                NextGen[ind] = Pop[Parents[ind]] ;
            }
            else{
            //else if( NumMuts>0 || NumRecs>0 ){
                while( (Haplotypes[Free].NumHap > 0) ){
                    Free++ ;
                }
                if(NumRecs==0){
                    // just assign haplo, mutate later
                    Haplotypes[Free].Muts = (*(Pop[Parents[ind]])).Muts ; // Pop[Parents[ind]] still points to same haplo tho
                }else{
                    // RECOMBINE
                    do{
                        RandDonor = int(rnd.randExc(Individuals)) ; // real number in [0,Nv), int -> includes 0 but not Nv
                    } while (RandDonor == Parents[ind]) ; // Cant select same individual!!
                    Haplotypes[Free].Muts.clear() ;
                    HomoRecomb( Haplotypes[Free], Pop[Parents[ind]], Pop[RandDonor], NumRecs, GeomTractLength, Loci) ;
                }
                // MUTATE
                for(i=0; i<NumMuts; i++){
                    // create new mutation
                    Haplotypes[Free].Muts.push_back( SharedMutPtr(new Mutation) ) ;
                    last =(Haplotypes[Free].Muts.size()-1) ;
                    // assign position
                    Haplotypes[Free].Muts[last]-> pos = rnd.rand() ;
                    if( int( (Haplotypes[Free].Muts[last]-> pos)*Loci )%3 == 0 ){
                        Haplotypes[Free].Muts[last]-> UnderSeln  = 0 ;
                    }else{
                        Haplotypes[Free].Muts[last]-> UnderSeln  = 1 ;
                    }
                }
                std::sort(Haplotypes[Free].Muts.begin(), Haplotypes[Free].Muts.end(), compareByPointerToPos) ;
                NextGen[ind] = &Haplotypes[Free] ;
                Free++ ;
            }
        }
        for(ind=0; ind<twoN; ind++){
            Haplotypes[ind].NumHap = 0 ;
        }
        //Move population to next generation, pointer reassignment
        for(ind=0; ind<Individuals; ind++){
            Pop[ind] = NextGen[ind] ;
            (*(Pop[ind])).NumHap++ ;
        }
        
        // PRINT HAPLOTYPE COMPOSITION
        /*
        cout << "HAPLOTYPE COMPOSITION\n" ;
        cout << "GENERATION:" << gen << "\n" ;
        for(i=1; i < twoN; i++){
            cout << "\t" << i << "\n" ;
            if(Haplotypes[i].NumHap > 0){
                for(j=0; j<Haplotypes[i].Muts.size(); j++){
                    cout << Haplotypes[i].Muts[j]->pos << " " ;
                }
                cout << "\n" ;
                cout << Haplotypes[i].NumHap << "\n" ;
            }
        }
        */
        //
        if(gen == Generations){
            //Calculate statistics every SampGen
            for (i = 0; i < Individuals; i++){
                PopSamp[i] = (*(Pop[ i ])) ; // PopSamp not a pointer!
            }
                //cout << gen << "\t" << SampCounter << "\n" ;            
        }
    }
    // Construct SnpAlleleFreqs, to see which mutations polymorphic/fixed
    for(ind=0; ind<Individuals; ind++){//individual
        size = SnpAlleleFreqs.size() ;
        // foreach mutation, go thru SnpAlleleFreqs to see if already present
        for(i=0; i<PopSamp[ind].Muts.size(); i++){//mutations
            //on 1st iteration, SnpAlleleFreqs empty and j==0, following for loop doesn't execute
            for(j=0; j<size ; j++){
                if( SnpAlleleFreqs[j].pos == PopSamp[ind].Muts[i]->pos ){
                    SnpAlleleFreqs[j].AF += 1 ;
                    break ; // if SNP found, breaks loop and next condition (j==size) is False
                }
            }
            //cycled through above loop, i.e. didn't find match for mut -> new mut
            if( j == size ){
                tempSNP.pos = PopSamp[ind].Muts[i]->pos ;
                tempSNP.AF = 1 ;
                SnpAlleleFreqs.push_back( tempSNP ) ;
            }
        }
    }
    // Find polymorphic SNPs
    int poly = 0 ;
    int fixed = 0 ;
    for(i=0; i<SnpAlleleFreqs.size(); i++){
        if( SnpAlleleFreqs[i].AF != Individuals ){
            SnpPolymorphic.push_back(SnpAlleleFreqs[i]) ;
            poly++ ;
        }else{
            fixed++ ;
        }
    }
    std::sort(SnpPolymorphic.begin(), SnpPolymorphic.end(), comparePos) ;
    //cout << "Polymorphic: " << poly << "\n" ;
    //cout << "Fixed: " << fixed << "\n" ;

    /*
    cout << "POLYMORPHIC SNPs\n" ;
    for(i=0; i<SnpPolymorphic.size(); i++){
        cout << SnpPolymorphic[i].pos << "\t" ;
    }
    cout << "\n" ;
     */
    out << "Loci " << Loci << "\n" ;
    out << "Individuals " << Individuals << "\n" ;
    out << "Generations " << Generations << "\n" ;
    out << "SampSizeTotal " << Individuals << "\n" ;
    out << "RecombinationRate " << RecombinationRate << "\n" ;
    out << "GeomTractLength " << GeomTractLength << "\n" ;
    out << "MutRate " << MutRate << "\n" ;
    out << "NScoeff " << NScoeff << "\n" ;
    out << "Wbar " << Wbar << "\n" ;
    out << "POSITIONS: " ;
    for(i=0; i<SnpPolymorphic.size(); i++){
        out << SnpPolymorphic[i].pos << " " ;
    }
    out << "\n" ;
    
    out << "HAPLOTYPES:\n" ;
    //Go though sample
    for(i=0; i<Individuals; i++){
        // Go through each polymorphic position in sample
        for(j=0; j<SnpPolymorphic.size(); j++){
            // does this polymorphic SNP exist in sample?
            size = PopSamp[i].Muts.size() ;
            for(k=0; k<size; k++){
                // if yes, print "1"
                if(PopSamp[i].Muts[k]->pos == SnpPolymorphic[j].pos){
                    out << "1" ;
                    break ;
                }
            }
            // cycled through all mutations, couldn't find in SnpPolymorphic
            if(k==size){
                // if no, print "0"
                out << "0" ;
            }
        }
        out << "\n" ;
    }
    
    /*
    //PRINT SAMPLE
    for (i=0; i<SampSizeTotal ;i++){
        cout << "Ind " << i << " " ;
        for(j=0; j<PopSamp[i].Muts.size(); j++){
            cout << PopSamp[i].Muts[j]->pos << " " ;
        }
        cout << "\n" ;
    }
    */
    
    /*=========================
     END WRIGHT-FISHER EVOLUTION
     =========================*/

    end = time(0);
    int TimeTaken = int(difftime(end,begin)) ;
    out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    
    
    out.close() ;
    
    delete [] Haplotypes ;
    delete [] Pop ;
    delete [] NextGen ;
    delete [] Parents ;
    delete [] PopSamp ;
    delete [] Wtot ;
}
