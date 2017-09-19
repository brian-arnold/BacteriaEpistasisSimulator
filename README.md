# BacteriaEpistasisSimulator

The files in this repository may be used as a pipeline to simulate a multilocus trait
with additive and epistatic effects. There are two parts to this pipeline:

1.) simulate neutral loci to mutation-recombination-drift balance
2.) use the output from these neutral simulations as input into simulations with selection
on a multilocus trait with additive and epistatic effects

This is done to simulate selection on standing genetic variation, where the variable loci
have J-shaped allele frequencies and exhibit equilibrium levels of linkage disequilibrium.

There are three subdirectories here for (1) the neutral forward-in-time simulator, (2) a 
script to convert the output of these neutral simulations so that they may be used in the 
downstream simulations with epistasis, and (3) the multilocus selection simulators with 
additive and epistatic effects, with two recombination models for Bacteria and Eukaryotes.
