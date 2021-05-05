#!/usr/bin/perl

###--- zwy 2014.12.24 ---###
###---      v1.01     ---###
###--- This script allow multi-threading to predict CAPs Marker ---###
###--- This script require restriction_Site_for_targetscf.pl(v2.03) and ParaFly! ---###

use strict;
use Getopt::Std;
use Cwd 'abs_path';

###############################################################################################
#####----------- input file （require samtools vcf file and genome sequence)-------------######
### output file are target scffolds which have putative SNPs,Enzyme site list #########
our ($opt_v,$opt_s,$opt_e,$opt_r,$opt_D,$opt_p,$opt_f,$opt_a);
getopts("v:s:e:r:D:p:f:a");
my $vcf_file = $opt_v;                                  #### from samtools or GATK
my $genome_file = $opt_s;                               #### origin gennome sequence
my $seq_range_option = (defined $opt_e)?$opt_e:1000;    #### default 1Kb 
my $range = (defined $opt_r)?$opt_r:300;                #### default 300bp
my $Depth = (defined $opt_D)?$opt_D:13;                 #### Depth [13]
my $threads = (defined $opt_p)?$opt_p:6;                #### number of threads
my $split_line = (defined $opt_f)?$opt_f:50;            ####
my $alt_Enzyme = $opt_a;                                #### find the Enzyme which contain the alternative base such as 'B,V,D,H,N....'
###############################################################################################
###############################################################################################
my $Usage = "\n\tUsage:\n
	Notice: This script allow multi-threading to predict CAPs Marker within restriction_Site_for_targetscf.pl(v2.03) and ParaFly!
	You are using Enzyme_ParaFly v1.01!
	
        $0 -v <VCF FILE> -s <GENOME FILE> <OPTIONS>
                -v  <FILE>    from samtools or GATK (vcf format file)
                -s  <FILE>    genome sequence (fasta format)
                -e  <INT>     extract sequence length around the Enzyme site (bp) [1000]
                -r  <INT>     detect the Enzyme site which doesn't has the same Enzyme site in a range (bp) [300]
                -D  <INT>     Depth [13] (Not useful for GATK vcf file)
		-p  <INT>     number of threads [6]
		-f  <INT>     number of line in each split vcf file [50]
                -a  <LOGICAL> find the Enzyme which contain the alternative base such as 'B,V,D,H,N....'
                           Warning:this option will slow the progrem down!!!

                               ---Alternative base list---

                           K = G or T            B = C or G or T
                           M = A or C            D = A or G or T
                           R = A or G            H = A or C or T
                           S = C or G            V = A or C or G
                           W = A or T            N = A or C or G or T
                           Y = C or T
                \n";

die $Usage unless ($opt_v && $opt_s);
##############################################################################################
my $vcf_file_name = (split /\//,$vcf_file)[-1];
if (!-e $vcf_file_name) {
	system "cp $vcf_file ./";
	$vcf_file = $vcf_file_name;
}
##############################################################################################
my $time1 = time;
open GENOME,'<',"$genome_file" or die "Cannot open genome file:$!";
my ($seq_id,%seq);
while (<GENOME>) {
	chomp;
	if (/>/) {
		$seq_id = (split />/,(split)[0])[1];
		$seq{$seq_id} = '';
	}else{
		$seq{$seq_id} .= $_;
	}
}
close GENOME;                                                                        
##############################################################################################
open VCF,'<',"$vcf_file" or die "Cant open vcf file:$!";
my @vcf = <VCF>;
my $vcf_line = @vcf;
my ($i,@line,$chr,$pos,$tmp,$new_pos,$new_line,$fi_pos,$pre_pos,$seq_e,$pre_chr,$greater,@split_vcf_file);
foreach (@vcf) {
	chomp;
	if (/\A#/) {next};
	$chr = (split)[0];
	last;
}
my @vcf = ();
close VCF;
my $ii = 1;
mkdir 'mediate_enzyme_file' or warn "Cannot make mediate_enzyme_file directory:$!";
open my $VCFSPLIT,'>',"mediate_enzyme_file/$vcf_file-$chr-$ii" or die "Cannot creata $vcf_file-$chr-$ii";
push @split_vcf_file,"mediate_enzyme_file/$vcf_file-$chr-$ii";
open VCF,'<',"$vcf_file" or die "Cant open vcf file:$!";
while (<VCF>) {
	chomp;
	if (/\A#/) { next };
	@line = split;
	$chr = shift @line;
	$pos = shift @line;
	$i ++;
	if ($i == 1) {
		$fi_pos = (split)[1];
		$pre_chr = $chr;
	}
	if ($chr eq $pre_chr) {
		if ( $i <= $split_line*$ii ) {
			if ($pos > $seq_range_option) {
				$tmp = $pos - $fi_pos + $seq_range_option;
			}else{
				$tmp = $pos;
			}
			$new_pos = "$pos\_$tmp";
			$new_line = join("\t",@line);
			print $VCFSPLIT "$chr\t$new_pos\t$new_line\n";
			$pre_pos = $pos;
		} else {
			if ($fi_pos > $seq_range_option) {
				$seq_e = substr($seq{$chr},($fi_pos - $seq_range_option),($pre_pos - $fi_pos + 2*$seq_range_option));
			}else{
				$seq_e = substr($seq{$chr},0,($pre_pos + $seq_range_option));
			}
			open SEQ,'>',"mediate_enzyme_file/$vcf_file-$chr-$ii.fa";
			print SEQ ">$chr\n$seq_e\n";
			close SEQ;
			$ii ++;
			$fi_pos = $pos;
			close $VCFSPLIT;
			open $VCFSPLIT,'>',"mediate_enzyme_file/$vcf_file-$chr-$ii";
			push @split_vcf_file,"mediate_enzyme_file/$vcf_file-$chr-$ii";
               		@line = split;
                	$chr = shift @line;
                	$pos = shift @line;
                	if ($pos > $seq_range_option) {
                        	$tmp = $pos - $fi_pos + $seq_range_option;
                	}else{
                        	$tmp = $pos;
                	}
                	$new_pos = "$pos\_$tmp";
                	$new_line = join("\t",@line);
                	print $VCFSPLIT "$chr\t$new_pos\t$new_line\n";
			$pre_pos = $pos;
		}
	}else{
        	if ($fi_pos > $seq_range_option) {
               		$seq_e = substr($seq{$pre_chr},($fi_pos - $seq_range_option),($pre_pos - $fi_pos + 2*$seq_range_option));
                }else{
                        $seq_e = substr($seq{$pre_chr},0,($pre_pos + $seq_range_option));
                }
                open SEQ,'>',"mediate_enzyme_file/$vcf_file-$pre_chr-$ii.fa";
                print SEQ ">$pre_chr\n$seq_e\n";
                close SEQ;
                $ii = 1;
		$i = 1;
                $fi_pos = $pos;
                close $VCFSPLIT;
                open $VCFSPLIT,'>',"mediate_enzyme_file/$vcf_file-$chr-$ii";
		push @split_vcf_file,"mediate_enzyme_file/$vcf_file-$chr-$ii";
                @line = split;
                $chr = shift @line;
               	$pos = shift @line;
                if ($pos > $seq_range_option) {
                	$tmp = $pos - $fi_pos + $seq_range_option;
                }else{
                        $tmp = $pos;
                }
                $new_pos = "$pos\_$tmp";
                $new_line = join("\t",@line);
                print $VCFSPLIT "$chr\t$new_pos\t$new_line\n";
		$pre_chr = $chr;
		$pre_pos = $pos;
	}
}
if ($fi_pos > $seq_range_option) {
	$seq_e = substr($seq{$pre_chr},($fi_pos - $seq_range_option),($pre_pos - $fi_pos + 2*$seq_range_option));
}else{
        $seq_e = substr($seq{$pre_chr},0,($pre_pos + $seq_range_option));
}
open SEQ,'>',"mediate_enzyme_file/$vcf_file-$pre_chr-$ii.fa";
print SEQ ">$pre_chr\n$seq_e\n";
close SEQ;
close VCF;
close $VCFSPLIT;
###################################--- creating cmd file for ParaFly ---###
open CMD,'>',"parafly_for_enzyme_prediction.cmd";
foreach (@split_vcf_file) {
	if (defined $opt_a) {
		print CMD "restriction_Site_for_targetscf.pl -v $_ -s $_.fa -e $seq_range_option -r $range -D $Depth -a\n";
	}else{
		print CMD "restriction_Site_for_targetscf.pl -v $_ -s $_.fa -e $seq_range_option -r $range -D $Depth\n";
	}
}
##################################--- ParaFly
!system "ParaFly -c parafly_for_enzyme_prediction.cmd -CPU $threads -failed_cmds Correlation.cmd.failed -v" or die "Something wrong with ParaFly";
my @mediate_fasta = map { "$_.fa" } @split_vcf_file;
my @better_enzyme = map { "$_-Better_Enzyme_site.txt" } @split_vcf_file;
my @enzyme = map { "$_-Enzyme_site.txt" } @split_vcf_file;
unlink 'parafly_for_enzyme_prediction.cmd', 'parafly_for_enzyme_prediction.cmd.completed';
open BETTER,'>',"Better_enzyme_site.txt";
foreach (@better_enzyme) {
	chomp;
	open BE,'<',"$_" or die "Cant open $_:$!";
	while (<BE>) {
		chomp;
		print BETTER "$_\n";
	}
	close BE;
}
close BETTER;
open ENZYME,'>',"Enzyme_site.txt";
foreach (@enzyme) {
	chomp;
	open EN,'<',"$_" or die "Cant open $_:$!";
	while (<EN>) {
		chomp;
		print ENZYME "$_\n";
	}
}
close ENZYME;
unlink glob "mediate_enzyme_file/* mediate_enzyme_file/.*";
#unlink "$vcf_file_name";
rmdir 'mediate_enzyme_file';
my $time2 = time;
my $time = ($time2 - $time1) / 60;
printf "Total time is :%.0f min\n",$time;

