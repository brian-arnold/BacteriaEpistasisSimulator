Neutral simulator:

Installation:
Have all files ending in '.cpp' and '.h' in the same directory. You must also have GNU
Compiler Collection (GCC) installed and also boost C++ libraries (version 1.58 or later).
Then, in this directory type:

g++ -O3 ./*.cpp -I /Directory/of/boost_1_58_0

where the filepath to the boost C++ library (www.boost.org) is included after the –I flag.
This should create a program named “a.out”, which may be renamed.

Input File:
The “Parameters” input file must contain two rows, the first row contains a header to label 
the input parameter values the second row contains the actual values. 
The following provides an example:

Loci	Ind	Gen		RecR		TrLen	MutR		Scoeff	NumSim
500		100	1000	0.00000575	0.002	0.00001345	0.0		1

where the columns contain the following information:
	1.) the number of loci, which may be interpreted as the number of contiguous sites/nucleotides
	2.) the number of individuals in the population (N)
	3.) the number of generation to run the simulator
	4.) the raw recombination rate per site, to simulate with 2Nr=0.01 per site (where N 
		is population size and r is recombination rate), simply divide the desired value 
		of 2Nr by 2N, where N is specified in column 2, and put this value of r as a 
		parameter value here. NOTE: this is the per site recombination rate
	5.) the mean recombination tract length, specified here as the parameter for a 
		geometric distribution. For instance, for a mean tract length of 500bp, put the 
		value of 1/500, or 0.002. Tract lengths for each recombination event will vary
		stochastically according to a geometric distribution.
	6.) the mutation rate per site, to simulate with 2Nu=0.01 per site (where N is 
		population size and u is mutation rate), simply divide the desired value of 2Nu by
		2N, where N is specified in column 2, and put this value of u as a parameter value
		here. NOTE: as with the recombination parameter specified above, this is the per 
		site mutation rate
	7.) an optional negative selection coefficient, where each individual's fitness is 
		penalized by 1-Scoeff for each mutation it has, such that all mutations are 
		deleterious. Here, fitness is modeled as multiplicative, such that an individual 
		with 10 mutations would have a fitness of 10*(1-Scoeff).
	8.) the number of simulations to perform per parameter set. The output files "Results"
		are appended with the simulation number.
		
Output File:
Output file(s) are labelled "Results" followed by the simulation replicate. The first 8
lines of the output file lists the parameter values used for this simulation (the same
values specified in the "Parameters" input file). This is then followed by "Wbar", the
mean fitness of the population, which should only be less than 1 if a value was specified 
for "Scoeff".
The 10th line contains the positions of the polymorphic sites, separated by spaces. These
positions are on the interval [0,1]. To obtain the nucleotide position of each 
polymorphic site, simply multiply these numbers by the number of loci (or length of the
DNA segment) specified as "Loci" in the "Parameters" input file and remove anything after
the decimal.
The 11th line, which says "HAPLOTYPES:", indicates that the population of haplotypes are 
printed on the following lines. Here, all N individuals are printed. Each row specifies 
the haplotype of an individual, represented as a series of 0's and 1's, where 0 is the 
ancestral allele and 1 is the derived allele. Only variable sites are printed (and kept 
track of during the simulation). If you wish to specify the the entire haplotype, you may
insert 0's at all positions in between the variable sites printed here (this program does
not do this for you). 











		
		
		