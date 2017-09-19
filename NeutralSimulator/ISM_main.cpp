#include "ISM.h"
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
    int LociS, IndividualsS, GenerationsS, SampIncS, NumSim  ;
    double UnscaledRecombRateS, RecombinationRateS, GeomTractLengthS, UnscalesMutationRateS, MutationRateS, NScoeffS ;
    
    ReadParameterFile(LociS, IndividualsS, GenerationsS, UnscaledRecombRateS, GeomTractLengthS, UnscalesMutationRateS, NScoeffS, NumSim) ;
    RecombinationRateS = UnscaledRecombRateS*LociS ;
    MutationRateS = UnscalesMutationRateS*LociS ;
   
    int sim ;
    for(sim=1; sim<=NumSim; sim++){
        WrightFisherEvo(LociS, IndividualsS, GenerationsS, RecombinationRateS, GeomTractLengthS, MutationRateS, NScoeffS, sim) ;
    }
    
    
    return 0 ;
}

