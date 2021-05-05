#!/usr/bin/perl -w
#
use strict;#
use warnings;#
use Getopt::Std;#
use File::Path qw(make_path remove_tree);#
#
our ($opt_z,$opt_o,$opt_D,$opt_d,$opt_c,$opt_f,$opt_3,$opt_4,$opt_5,$opt_6,$opt_a,$opt_b,$opt_q,$opt_Q,$opt_s,$opt_n);#
getopt ("z:o:D:d:c:f:3:4:5:6:a:b:q:Q:s:n:");#
########### snpMapper.pl Qun Hu and Weiyi Zhang
# Function: Get base counts from samtools' mpileup output
# the mpileup file should be produced with mutant.bam as first inputfile, wild_type.bam as second input file
# samtools command should be run as:
#     samtools mpile -ABf ref.fa mut.bam wt.bam >mut-wt_raw.mpileup
#
# 21:47 2014/3/24 v0.0.1 
# 19:47 2014/3/26 v1.0.1 modify the program structure
# 19:42 2014/3/30 v1.0.3 Incorporate positions that have 'skipped reference' meaning
# 16:05 2014/4/9  v1.04
# 22:41 2014/5/21 v1.05
# 22:18 2016/9/27 v1.07 (hisat2version) modify the tophat2 to hisat2 and this version is for web
# 17:57 2016/10/4 v1.08 Multi-thread the snpMapper
# 
# Future version should add algorithm to consider contamination in each pool -- but future will not come
# 
############################################################################# 
# this program takes mpileup file as input, and outoput several files regarding the SNP
# It can be modified to subject tetraploid species RAN-seq analysis

my $time1 = time(); 
my $version = 'v1.08';
my $Usage = "Usage:\n  $0 -OPTIONS VALUES

options:
--input options
	-f REFGENOME    reference genome for mapping
	-b hisat2index  hisat2 index
	 
--output options	 
	-o OUTDIR   output direcotry containg all output files [snpMappr_out]
	-3 FILE     bool1 left reads path
	-4 FILE     pool1 right reads path
	-5 FILE     pool2 left reads path
	-6 FILE	    pool2 right reads path
	-p FILE     mpileup format file produced by samtools mpileup command, if this option is not provided, then -c,-m,-f have to be provided
	-s STR      prefix string [sm]

--criteria options
	-D INT      minimum read coverage to be considered for valid SNP [10]
	-n FLOAT    pos whose wildtype allele frequency gt FLOAT will be ignored [0.95]
	-d FLOAT    minimum delta snp index, required [0.2]
	-q INT      skip alignments with mapQ smaller than INT [55]
	-Q INT      skip bases with baseQ/BAQ smaller than INT [30]
	
--perfromanat options
	-a INT      cpus cores used for the analysis [6]	
	\n";
	
die $Usage unless ($opt_f && $opt_b && $opt_3 && $opt_4 && $opt_5 && $opt_6);

# -------------- Start snpMapper --------------------
print "\n\nStarting $0 ...\n";
my $refgenome           = $opt_f; # reference genome
my $bowtie2_index       = $opt_b || 'bowtie2-index'; # intermediate outfile for bowtie2 index
my $fq1_left		= $opt_3;
my $fq1_right		= $opt_4;
my $fq2_left		= $opt_5;
my $fq2_right		= $opt_6;
my $user_defined_depth  = (defined $opt_D)?$opt_D:10;  # read coverage for a given position, default=5;
my $delta            = (defined $opt_d)?$opt_d:0; # delta of SNP-index of locus between mutant and wild type
my $nopoly           = (defined $opt_n)?$opt_n:0.95;
my $cpu              = (defined $opt_a)?$opt_a:6;
my $mapq             = (defined $opt_q)?$opt_q:55;
my $baseQ            = (defined $opt_Q)?$opt_Q:30;
my $prefix           = (defined $opt_s)?$opt_s:'sm';
my $outdir           = (defined $opt_o)?$opt_o:"${prefix}snpMapper_out"; # outfile
my $Putative_SNP     = "$outdir/${prefix}candidatesnps_D${user_defined_depth}d$delta.txt";
my $merge_reads      = (defined $opt_z)?$opt_z:'F';

###################### snpMapper.pl Main program starts from here ##################################

print "  Checking fastq files ...</br>";

make_path($outdir,{verbose=>1,mode=>0777}) unless (-e $outdir);
my $tophat_outdir1 = "$outdir/${prefix}hisat2_out1";
my $tophat_outdir2 = "$outdir/${prefix}hisat2_out2";
make_path($tophat_outdir1, {-verbose=>1,mode=>0777}) unless -e $tophat_outdir1 && -d $tophat_outdir1;
make_path($tophat_outdir2, {-verbose=>1,mode=>0777}) unless -e $tophat_outdir2 && -d $tophat_outdir2;


my @bowtie_index_files = ("$bowtie2_index.1.ht2","$bowtie2_index.2.ht2",
						  "$bowtie2_index.3.ht2","$bowtie2_index.4.ht2",
);

# checking whether genome has be indexed
my $flag_bowtie2 = 1;
foreach my $file (@bowtie_index_files){
	unless (-e $file){
	$flag_bowtie2 = 0;
	last;
	}
}
if ($flag_bowtie2 == 0){
	print "$bowtie2_index  -- Not Found! Creating...<.br>";
	&command_mkbowtie_index($refgenome,$bowtie2_index);
	print "done\n\n";
}else{
	print "$bowtie2_index -- Found! hisat2 index building skipped!</br>";
}

my ($mutbamfile, $wtbamfile) = ("NA","NA");
	print "  mutbamfile -- not found\n  Runing hisat2 for sample mutant ...";
	&command_tophat($bowtie2_index,$fq1_left,$fq1_right,$tophat_outdir1,$cpu);
	print "finished!\n\n";
	$mutbamfile = "$tophat_outdir1/acc.sorted.bam";
#=cut 1

	print "  wtbamfile -- not found\n  Runing hisat2 for sample wildtype ...";
	&command_tophat($bowtie2_index,$fq2_left,$fq2_right,$tophat_outdir2,$cpu);
	print "finished!</br></br>";
	$wtbamfile = "$tophat_outdir2/acc.sorted.bam";
#=cut 2
# ----  check whether all required files are present ----- #
check_file($mutbamfile,$wtbamfile,$refgenome);
#exit(0);

print "  -- Starting converting pileup format to an elite SNP file for simplicity ... </br>";
#open  my $fh_pileup, "/results/others/huqun/bin/samtools mpileup -q $mapq -Q $baseQ -ABf $refgenome $mutbamfile $wtbamfile|"  or die "Can not read input from samtools!\n";

my $header = '#chr	pos	refbase	mut.cov	mut.A	mut.C	mut.G	mut.T	wt.cov	wt.A	wt.C	wt.G	wt.T';

#####----- start muti-threads -----#####
#####----- checking chromosome length ----#####
open GENOME,'<',"$refgenome" or die "Cannot open refgenome:$!";
my ($seq_id,$seq_ido,%seq_len);
while (<GENOME>) {
		chomp;
		if (/>/) {
			$seq_ido = (split />/,(split)[0])[1];
			if ( $seq_ido =~ /(\d+)\z/ ) {
				$seq_id = $1;
			}
			if ($seq_id =~ /\A0(\d+)/){
				$seq_id = $1;
			}
			$seq_id = "$seq_ido\t$seq_id";
			$seq_len{$seq_id} = '';
		}else{
			$seq_len{$seq_id} += length($_);
		}
}
close GENOME;
### find max length of all chromosomes
### sort chromosome names
my @seq_name = keys %seq_len;
@seq_name = sort {(split /\t/,$a)[1] <=> (split /\t/,$b)[1]} @seq_name;
#################################################
if (-e "/data/output/med_bam"){
	system "rm /data/output/med_bam -rf";
}
mkdir "/data/output/med_bam";

open CMD1,'>',"samtools_mpileup_1.cmd";
open CMD2,'>',"samtools_mpileup_2.cmd";
open CMD3,'>',"samtools_mpileup_3.cmd";
my @deltafile;
my $med_count;
foreach (@seq_name) {
	my $sequence_name = $_;
	my ($seq_id,$seq_num) = (split /\t/,$_)[0,1];
	my $i = 0;
	$med_count ++;
	my $medout1 = "/data/output/med_bam/$seq_num\_$med_count\_1.bam";
	my $medout2 = "/data/output/med_bam/$seq_num\_$med_count\_2.bam";
	my $medout3 = "/data/output/med_bam/$seq_num\_$med_count\_x.part";
	push @deltafile, "/data/output/med_bam/$seq_num\_$med_count\_x.txt";
	print CMD1 "samtools view -b $mutbamfile $seq_id:0-$seq_len{$sequence_name} > $medout1\n";
	print CMD1 "samtools view -b $wtbamfile $seq_id:0-$seq_len{$sequence_name} > $medout2\n";
	print CMD2 "samtools mpileup -q $mapq -Q $baseQ -ABf $refgenome $medout1 $medout2 > $medout3\n";
	print CMD3 "perl snpMapper_sub.pl $medout3 $user_defined_depth\n";
}
!system "ParaFly -c samtools_mpileup_1.cmd -CPU $cpu -failed_cmds samtools_mpileup.cmd1.failed -v" or die "Something wrong with ParaFly";
!system "ParaFly -c samtools_mpileup_2.cmd -CPU $cpu -failed_cmds samtools_mpileup.cmd2.failed -v" or die "Something wrong with ParaFly";
#system "rm -rf med_bam/*.bam";
!system "ParaFly -c samtools_mpileup_3.cmd -CPU $cpu -failed_cmds samtools_mpileup.cmd3.failed -v" or die "Something wrong with ParaFly";
#system "rm -rf med_bam/*.part";

open my $snpout,'>',"$Putative_SNP";
print $snpout "$header\n";
foreach (@deltafile) {
	open IN2,'<',"$_" or die "$_";
	while (<IN2>) {
		chomp;
		print $snpout "$_\n";
	}
	close IN2;
}
close $snpout;
@deltafile = ();

#system "rm med_bam samtools_mpileup* -rf ";

my $time2 = time();
my $time  = ( $time2 - $time1 ) / 60;
print "$0 finished! Total time elapsed: $time</br>";
print "At this points, please find your candidate SNPs in $outdir! Wish you a good Luck!</br>";
exit(0);

############################ subroutine definition #######################
sub check_file{
  my @files = @_;
	
	my $flag = 1;
	foreach my $file (@files){
		next if $file =~ /^#/;
		next if $file =~ /^\s*$/;
		if (-e $file){
			print "  $file -- found</br>";
		}else{
			print "  $file -- not found</br>";
			$flag = 0;
		}
	}	
	if ($flag == 0){
		print ("some files does not exist,see above!</br>");
		exit(1);
	}
}

sub sum{
	my $sum = 0;
	for (@_){
		$sum += $_; 
	}
	return $sum;
}

sub command_mkbowtie_index {
	my ($refgenome,$bowtie2_index) = @_;
	!system "hisat2-build -p $cpu $refgenome $bowtie2_index 2> hisat2_build.log" or die "Error in hisat2 build index";
}

sub command_tophat{
	my ($bowtie_index,$s_left,$s_right,$outdir,$cpu)  = @_;
	print "hisat2 -x $bowtie_index -1 $s_left -2 $s_right -p $cpu -S $outdir/acc.sam";
	!system "hisat2 -x $bowtie_index -1 $s_left -2 $s_right -p $cpu -S $outdir/acc.sam 2> $outdir/hisat2_aln.out" or die "Error in hisat2";
	!system "samtools view -bhS $outdir/acc.sam > $outdir/acc.bam" or die "Error in samtools view (sam to bam)";
	system "samtools sort -O BAM -@ $cpu -o $outdir/acc.sorted.bam $outdir/acc.bam";
	system "samtools index $outdir/acc.sorted.bam";
	#unlink "$outdir/acc.sam"; unlink "$outdir/acc.bam";
}

