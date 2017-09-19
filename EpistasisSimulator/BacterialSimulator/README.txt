INSTALLATION/COMPILING/RUNNING
To compile this program, keep all the files contents that end in .cpp and .h within the 
same directory, and type:

g++ -O3 -I /Path/To/Boost/Library *.cpp

where the filepath to the boost C++ library (www.boost.org) is included after the –I flag.
This should create a program named “a.out”, which may be renamed.


3 Input Files:

“Parameters” Input File 1
This file must be in the same directory as the compiled program (above) as it contains 10 
input parameters for the simulator. These parameter values must be on the second line of 
the file, and I use the first line as a header to label the parameters. From left to 
right, these parameters include:
1.)	Number of loci (L)
2.)	Number of individuals or population size (N)
3.)	Number of generations to run the simulator
4.) The sampling interval to measure summary statistics, in terms of N generations, where 
	N is specified in (2) (e.g. a value of 0.1 corresponds to measuring summary statistics
	every 0.1N generations)
5.)	Total recombination rate, i.e. the physical per base pair rate multiplied by the 
	genomic segment simulated (see 8.))
6.)	Geometric distribution parameter to determine tract lengths transferred between 
	individuals, must be between 0 and 1, and the mean tract length is the reciprocal of 
	the value specified here
7.)	Additive effect of each locus. For simulating effects of size 2Ns, where N is 
	population size and s is the selection coefficient, the parameter value specified here
	would be s.
8.)	Epistatic effect of each locus pair. For simulating effects of size 2Ns, where N is 
	population size and s is the selection coefficient for each locus pair, the parameter 
	value specified here would be s.
9.)	Circular genome size (in bp)
10.)	Number of simulation replicates

On the third line of the Parameters file should be the full path to a directory (either 
the current directory or an external directory) that contains the HaplotypeMatrix and 
Positions files (detailed below).


"HaplotypeMatrix"  Input File 2
This file is the starting population, where each row is an individual with L alleles 
represented as +/-1. This way, starting conditions may be arbitrarily specified with 
different allele frequencies and linkage relationships between loci.


"Positions"  Input File 2
This file contains a single line, listing the positions (integers) of each locus,
separated by tabs. These positions should increase from left to right, and the right-most
position must not be larger than the genome size specified in the "Parameters" file.


Note that there is no mutation rate; each locus is represented as a +/- 1.


OUTPUT FILES

“Results” file
This is the output of the simulator, and these files are enumerated by simulation 
replicate (e.g. Results1, Results2...)
This file begins with some general information from the simulation. First printed are the 
positions of the loci, which should match those listed in the "Positons" input file.
Second, the additive fitness effect of each locus is printed (AddFitConst_2Ns), which 
should correspond to the input additive effect of each locus multiplied by two times the 
population size (N).
Third are the epistatic fitness effects for each pairwise interaction (EpistasisTable_2Ns)
also listed as values of 2Ns, where s is the selection coefficient specified in the 
"Parameters" file. There should be L-1 lines, where the first line represents the L-1 
interactions between nearest-neighbor (adjacent) loci, the second line represents the L-2 
next-nearest-neighbor interactions (with one intervening locus), and so forth.

Next in this file are some summary statistics from the beginning of the simulation at time
t=0. The first of these summary statistics is Dprime_BY_DIST_gen0, which represents all 
the D’ values for each locus pair, where D’ here is a popular measure of linkage 
disequilibrium (Lewontin 1964). D’ values are printed for each locus pair in the same way 
as the EpistasisTable_2Ns section.

More summary statistics are then printed, including ThetaPi_gen0 (average number of 
pairwise differences at generation 0), Wmean_gen0 (mean fitness), Wvar_gen0 (the variance 
in fitness), Wskew_gen0 (skew, or third moment of the fitness distribution). Also printed
is the pairwise allele frequency (PairwiseAFS) for each pair of loci, following the same 
scheme as the EpistasisTable_2Ns section. Here, allele frequencies are for the "1" allele
(as opposed to the "-1" allele).

These summary statistics are then printed at specific generations, and labels are followed
by “_gen#” where # is the number of generations that have elapsed since t=0. 

