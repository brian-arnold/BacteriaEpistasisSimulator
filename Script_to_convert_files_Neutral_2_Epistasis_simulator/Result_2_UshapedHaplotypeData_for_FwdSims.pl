#!usr/bin/perl
use warnings ;
use strict ;
use lib "/n/holyscratch/bomblies_lab/bjarnold/List-MoreUtils-0.33/lib" ;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;

#open TESTOUT, ">./testout.txt" ;


my $batches = $ARGV[0] ; # the number of results 
my $InterSnpDist = $ARGV[1] ; # the spacing of each window, or the distance between windows within which you randomly select SNPs
my $SnpRange = $ARGV[2] ; # the range of sites within which you randomly select SNPs

unless($batches && $InterSnpDist && $SnpRange){
	print "You didn't give the script 3 input arguments; please do this and try again, and refer to the README2_PerlScriptConverter.txt for help\n" ;
	exit ;
}

foreach my $batch (1..$batches){
		my $CurrentMin = 0 ; # increase by by ($InterSnpDist*$count) - $SnpRange, where $count is current slice
		my %SitesToSample ;
		my %HaploSitesToSample ;

		my %Random_Sites ;
		my %Random_HaploSites ;
		my @AFS ;

		my $Sample_size ;
		my $num_sites ;

		my $lower ;
		my $upper ;
		my $middle ;


		#############################################
		# DATA STRUCTURES
		#############################################
		my $sum ;
		my $sq_sum ;

		my %multi_allele_sites ; #  $multi_allele_sites{num hits} = count ;
		my %Seq_data ; #  $Seq_data{index} = scalar of 0's and 1's
		my %Segsites ; # $Segsites{index} = site ;
		my %Segsites_Bi_Unfiltered_Pos; ; # $Segsites{index} = site ;
		my %Segsites_Bi_Unfiltered_Pos_HaploPosition ; # for mapping back to positin in string of 0's and 1's
		my %Segsites_Bi_FilterAF ; # $Segsites_Bi_FilterAF{index} = site ;
		my @Segsites_Bi_FilterAF_indices ; # for mapping back to positin in string of 0's and 1's
		my %AFS ; # $AFS{1}{site} = AF ; 
		my %Pairwise_Dprime  ; # $Pairwise_Dprime{site1}{site2} = D' ;

		my %SumStat_Results ;
		my %PairWise_LD_vs_dist ;
		my %PairWise_PC_vs_dist ;
		my %PairWise_LD_dist_slice ;
		my %PairWise_PC_dist_slice ;
		my %PairWise_Rsquared_dist_slice ;

		#############################################
		# ASSEMBLE SEQUENCES PER GENE PER IND
		#############################################
		open IN, "<./Results${batch}" ;
		my $num_segsites_raw = 0 ;
		while(<IN>){
			chomp $_ ;
	
			if($_ =~ m/^Loci/){
				my @line = split(" ", $_) ;
				$num_sites = $line[1] ;
			}
			if($_ =~ m/^Individuals/){
				my @line = split(" ", $_) ;
				$Sample_size = $line[1] ;
			}
			if($_ =~ m/^POSITIONS/){
				my @line = split(" ", $_) ;
				foreach my $x ( 1..$#line ){
					$Segsites{($x-1)} = (int($line[$x]*$num_sites)) ;
				}	
			}
			if($_ =~ m/^HAPLOTYPES/){
				foreach my $ind (0..($Sample_size-1)){
					my $x = <IN> ; chomp $x ;
					$Seq_data{ $ind } = $x ;
				}
				$num_segsites_raw = length($Seq_data{ 0 }) ;

			}
		}
		close IN ;
		$SumStat_Results{"S_star"} = $num_segsites_raw/$num_sites ;



		#############################################
		### delete any seg site index which has same site, multiallelic
		#############################################
		foreach my $i ( 0..$num_segsites_raw-2 ){
			if(exists $Segsites{$i}){
				my $j = $i+1 ;
				if($Segsites{$j} == $Segsites{$i}){
					while($Segsites{$j} == $Segsites{$i}){
						$j++ ; #to capture potential runs of multiple hits at same site
						if(!exists $Segsites{$j}){ #reached end
							last ;
						}
					}
					## delete $i through $j-1 (-1 b.c. while loop added extra)
					foreach my $x ($i..($j-1)){
						delete $Segsites{$x} ;
					}
					$multi_allele_sites{($j-1-$i)+1} ++ ;
				}else{
					next ;
				}
			}
		}
		my $count = 0 ;
		foreach my $key (sort{$a<=>$b} keys %Segsites){
			$Segsites_Bi_Unfiltered_Pos_HaploPosition{$count} = $key ; # key is position of segsite in 0/1 string
			$Segsites_Bi_Unfiltered_Pos{$count} = $Segsites{$key} ; # this is position on chromo
			$count++ ;
		}


		$count = 1 ;
		foreach my $key (sort{$a <=> $b} keys %Segsites_Bi_Unfiltered_Pos){
			if( $Segsites_Bi_Unfiltered_Pos{$key} > ($CurrentMin + $SnpRange) ){
				$CurrentMin = ($InterSnpDist*$count) - $SnpRange ;
				$count++ ;
			}
			if($Segsites_Bi_Unfiltered_Pos{$key} > $CurrentMin && $Segsites_Bi_Unfiltered_Pos{$key} < ($CurrentMin + $SnpRange)){
				push @{$SitesToSample{$count}}, $Segsites_Bi_Unfiltered_Pos{$key} ;
				push @{$HaploSitesToSample{$count}}, $Segsites_Bi_Unfiltered_Pos_HaploPosition{$key} ;
			}
		}


		$count = 0 ;
		foreach my $key (sort{$a <=> $b} keys %SitesToSample){
			my $rand_index = int(rand(scalar @{$SitesToSample{$key}})) ;
			$Random_Sites{$count} = ${$SitesToSample{$key}}[$rand_index] ;
			$Random_HaploSites{$count} = ${$HaploSitesToSample{$key}}[$rand_index] ;
			$count++ ;
		}


		#############################################
		## CONSTRUCT NEW DATASET
		#############################################
		my %NewData ; # make new data struct to print individual for ea row of table
		foreach my $key (sort{$a<=>$b} keys %Random_HaploSites ){
			my $site = $Random_HaploSites{$key} ;
			my $af = 0 ;
			my $allele ;
			## MAKES DISTRIBUTION U-SHAPED
			if(rand(1) > 0.5){
				$allele=1 ;
			}else{
				$allele=0 ;
			}
			foreach my $ind (keys %Seq_data){
				if( substr($Seq_data{$ind}, $site, 1) == 1 ){
					$af++ ;
					if($allele==1){ # switches polarity -> u-shape AFS
						$NewData{$ind}{$site} = 1;
					}else{
						$NewData{$ind}{$site} = -1;
					}
				}else{
					if($allele==1){
						$NewData{$ind}{$site} = -1;
					}else{
						$NewData{$ind}{$site} = 1;
					}		
				}
			}
			push @AFS, $af ;
		}
		open OUT, ">./AlleleFreqs${batch}" ;
		foreach (@AFS){
			print OUT $_, "\t" ;
		}
		print OUT "\n" ;
		close OUT ;
		
		open OUT, ">./Positions${batch}" ;
		foreach my $key (sort{$a <=> $b} keys %Random_Sites){
			print OUT $Random_Sites{$key}, "\t" ;
		}
		print OUT "\n" ;
		close OUT ;

		open OUT, ">./HaplotypeMatrix${batch}" ;
		foreach my $ind (sort{$a <=> $b} keys %NewData){
			foreach my $site ( sort{$a <=> $b} keys %{$NewData{$ind}} ){
				print OUT $NewData{$ind}{$site}, "\t" ;
			}
			print OUT "\n" ;
		}
		close OUT ;
}
exit ;

