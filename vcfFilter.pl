#!/usr/bin/perl

use warnings;
use strict;

# filters poolSeq vcf files

# usage vcfFilter.pl infile.vcf
#
# change the marked variables below to adjust settings
#

#### stringency variables, edits as desired
## 10 samples * 50 = 500
my $minCoverage = 500; # minimum number of sequences; DP
my $bqrs = 4; # Mann-Whitney U-z test of Base Quality Bias
my $mqrs = 4; # Mann-Whitney U-z test of Mapping Quality Bias
my $rprs = 4; # Mann-Whitney U-z test of Read Position Bias
my $mq = 30; # minimum mapping quality; MQ

my @line;

open (OUT, "> filtered_cmacPool.vcf") or die "Could not write the outfile\n";

my $cnt = 0;
my $flag = 0;
my $in;
my $first = 1;

foreach $in (@ARGV){
	open (IN, $in) or die "Could not read the infile = $in\n";
	print "Working on $in\n";

	while (<IN>){
		chomp;
		$flag = 1;
		if (m/^\#/ and $first == 1){ ## header row, always write
			$flag = 1;
		}
		elsif (m/^ENA/){ ## this is a sequence line, you migh need to edit this reg. expr.
			$flag = 1;
			if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
				$flag = 0;
			}
			@line = split(/\s+/,$_);
			if(length($line[3]) > 1 or length($line[4]) > 1){
				$flag = 0;
			}
			m/DP=(\d+)/ or die "Syntax error, DP not found\n";
			if ($1 < $minCoverage){
				$flag = 0;
			}
			if(m/BQBZ=([0-9e\-\.]*)/){
				if (abs($1) > $bqrs){
					$flag = 0;
				}
			}
			if(m/MQBZ=([0-9e\-\.]*)/){
				if (abs($1) > $mqrs){
					$flag = 0;
				}
			}
			if(m/RPBZ=([0-9e\-\.]*)/){
				if (abs($1) > $rprs){
					$flag = 0;
				}
			}
			if(m/MQ=([0-9\.]+)/){
				if ($1 < $mq){
					$flag = 0;
#					print "fail MQ : ";
				}
			}
			if ($flag == 1){
				$cnt++; ## this is a good SNV
			}
		}
		else{ ## pops up on files after first
			##	print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
			$flag = 0;
		}
		if ($flag == 1){
			print OUT "$_\n";
		}
	}
	close (IN);
	$first = 0;
}
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
