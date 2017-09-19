#include "AddEpi.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
MTRand rnd;


int main()
{
    // S appendix indicates starting variable to be passed to WF function
    int LociS, IndividualsS, GenerationsS, GenomeSizeS, SegSizeS, NumSim ;
    double SampGenS, UnscaledRecombRateS, RecombinationRateS, VaS, ViS ;
    string FilePath ;
    
    ReadParameterFile(LociS, IndividualsS, GenerationsS, SampGenS, UnscaledRecombRateS, VaS, ViS, GenomeSizeS, NumSim, FilePath) ;
    //RecombinationRateS = UnscaledRecombRateS*LociS ;
    RecombinationRateS = UnscaledRecombRateS ;
    int sim, rep ;
    for(sim=1; sim<=NumSim; sim++){
        WrightFisherEvo(LociS, IndividualsS, GenerationsS, SampGenS, RecombinationRateS, VaS, ViS, GenomeSizeS, sim, FilePath) ;
    }
    
    
    return 0 ;
}

