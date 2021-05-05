#!/usr/bin/perl

##--- zwy 2014.12.1 ---##
##--- use this script in the directory which contain tophat results folder to detect differential expression genes by cufflinks ---##
use strict;
use Getopt::Std;
use Cwd;
use Cwd 'abs_path';

our ($opt_1,$opt_2,$opt_f,$opt_p,$opt_3,$opt_4);
getopts('1:2:f:p:3:4:');

my $lable1 = $opt_1;
my $lable2 = $opt_2;
my $bam1 = $opt_3;
my $bam2 = $opt_4;
my $genome_path = $opt_f;
my $cpu = (defined $opt_p)?$opt_p:6;

my $Usage = "\n###---use this script in the directory which contain tophat results folder---###
###---to detect differential expression genes by cufflinks.---###
#####------------------  zwy 2014.12.1  -----------------#####
#####
#################################################################
\n$0 -1 <LABLE1> -2 <LABLE2> -3 <bam1> -4 <bam2> -f <GENOME_FILE> -p <Threads>
		-1 <Chracter> pool1's name
		-2 <Chracter> pool2's name
		-3 <FILE>     bam1
		-4 <FILE>     bam2
		-f <FILE>     genome sequence
		-p <INT>      number of threads [6]
	\n";
die $Usage unless ($opt_1 && $opt_2 && $opt_f  && $opt_3 && $opt_4);
###############################################################
mkdir "cufflinks_results" or warn "Cant make 'cufflinks_results' directory!:$!";
############-- cufflinks --###########
my $counter = 0;
#my @file = map {abs_path($_)} @file;
my @cuffmerge;
		mkdir "cufflinks_results/$lable1" or warn "Cannot make $lable1 directory!:$!";
        	!system "/results/others/local/bin/cufflinks -p $cpu -o cufflinks_results/$lable1 --multi-read-correct -b $genome_path -L $lable1 $bam1" or die "Something wrong with cufflinks!";
       		push @cuffmerge, "cufflinks_results/$lable1/transcripts.gtf";
		
		mkdir "cufflinks_results/$lable2" or warn "Cannot make $lable2 directory!:$!";
       		!system "/results/others/local/bin/cufflinks -p $cpu -o cufflinks_results/$lable2 --multi-read-correct -b $genome_path -L $lable2 $bam2" or die "Something wrong with cufflinks!";
                push @cuffmerge, "cufflinks_results/$lable2/transcripts.gtf";
###########-- cuffmerge --############
open CUFMERGE,'>',"cuffmerge_list.txt" or die "Cannot!";
foreach (@cuffmerge) {
        print CUFMERGE "$_\n";
}
close CUFMERGE;
!system "/results/others/local/bin/cuffmerge -p $cpu -s $genome_path cuffmerge_list.txt" or die "Something wrong with cuffmerge!";
###########--  cuffdif  --#############
!system "/results/others/local/bin/cuffdiff -o cuffdiff_results --labels $lable1,$lable2 -p $cpu -b $genome_path -u merged_asm/merged.gtf $bam1 $bam2" or die "Something wrong with cuffdiff!";
