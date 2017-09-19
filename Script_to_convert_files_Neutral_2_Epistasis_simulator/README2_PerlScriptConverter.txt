This file, README2_PerlScriptConverter.txt, explains how to convert output from this
neutral simulator into input files that may be used for the forward-in-time epistasis
simulator.

This script is used to take output from the neutral simulator and create input files for 
the epistasis simulator. This way, for each neutral simulation, one can simulate many 
different scenarios of selection without having to run a fowrard-in-time simulator to
mutation-recombination-drift balance, which can take a long time depending on the 
populaiton size, mutation rate, and recombination rate.

Essentially, this script takes the output from the neutral simulator and divides the 
number of Loci (specified in the Parameters file for the neutral simultor) into windows.
Then, within each window, it randomly select a SNP, but you have some control over where
it selects SNPs within these windows. For instance, if you specify a window size of 100bp,
you can tell the program if you'd like it to randomly sample a SNP within those full 100bp,
or if you'd like it to sample SNPs only from the leftmost 50bp of the 100bp window. Having
this control can be useful if you'd like to guarantee that SNPs have a minimum amount of 
distance separating them.


This script takes 3 input arguments, and the script is used the following way:

perl Result_2_UshapedHaplotypeData_for_FwdSims.pl Arg1 Arg2 Arg3

where Arg1 - Arg3 are the input arguments.

Input argument 1: the number of results for which you would like to extract SNPs, in 
order to put them into the epistasis simulator. These results files are printed from the 
neutral simulator and are labelled "Results1", "Results2",...


Input argument 2: The window size. For example, if you've simulated a neutral 1kb segment
and specify 100 here, it will divide the segment into 10 windows, and search for SNPs
and randomly select ONE within each window.

Input argument 3: The interval within each window in which you'd like to sample SNPs. This
essentially is a sub-window. If this number is the same as that specified in argument 2, 
the program will randomly sample SNPs within the full window. For instance, using the 
example above where a 1kb segment was split into 10 100bp windows, one SNP would be 
randomly selected every 100bp. However, if you specify 50bp for this argument, then it will
only look for SNPs within the first 50bp of the 100bp window. This guarantees that each 
SNP is at least 50bp apart!


NOTE: If you've done neutral simulations with very little diversity, this program may 
behave unexpectedly and not output the correct number of SNPs, because there may be no 
SNPs within some genomic windows. Ways to circumvent this include 1.) conducting neutral 
simulations with a higher population size or higher mutation rate or 2.) increase the 
interval within which the program looks for SNPs (Argument 3).

