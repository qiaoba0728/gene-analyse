#!/usr/bin/perl

##---- zwy 2014.12.24 ----##!!!! Not gb version !!!!##
##---- version V2.03 ----##
use strict;
use Getopt::Std;
###############################################################################################
#####----------- input file （require samtools vcf file and genome sequence)-------------######
### output file are target scffolds which have putative SNPs,Enzyme site list #########
our ($opt_h,$opt_v,$opt_s,$opt_e,$opt_r,$opt_D,$opt_a);
getopts("h:v:s:e:r:D:a");
my $putativeSnp = $opt_h;                               #### from huqun progrem
my $vcf_file = $opt_v;                                  #### from samtools or GATK
my $genome_file = $opt_s;                               #### origin gennome sequence
my $seq_range_option = (defined $opt_e)?$opt_e:1000;    #### default 1Kb 
my $range = (defined $opt_r)?$opt_r:300;                #### default 300bp
my $Depth = (defined $opt_D)?$opt_D:13;                 #### Depth [13]
my $alt_Enzyme = $opt_a;                                #### find the Enzyme which contain the alternative base such as 'B,V,D,H,N....'
###############################################################################################
###############################################################################################
my $Usage = "\n\tUsage:\n
	$0 -h <HUQUN SNP FILE> -v <VCF FILE> -s <GENOME FILE> <OPTIONS>
		-h  <FILE>    from snpMapper progrem (Not nesseary!)
		-v  <FILE>    from samtools or GATK (vcf format file)
		-s  <FILE>    genome sequence (fasta format)
		-e  <INT>     extract sequence length around the Enzyme site (bp) [1000]
		-r  <INT>     detect the Enzyme site which doesn't has the same Enzyme site in a range (bp) [300]
		-D  <INT>     Depth [13] (Not useful for GATK vcf format)
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
###############################################################################################
my $vcf_line = `wc -l $vcf_file`;
my @vcf_putativesnps;
if (defined $opt_h) {
	print "\tSTEP1: extract the name of putative scf START!\n";
############################ STEP1 extract the name of putative scf
	open SNP,'<',"$putativeSnp" or die "Cant open putativeSnp file:$!";
	my(%huqun_chrom_name,@huqun_chrom,%origin_name,@origin_names);
	while(<SNP>){
	        chomp;
	        my $chrom_name = (split)[0]; ######## may need to change MARK1
	        $huqun_chrom_name{$chrom_name} = 1;
	}
	@huqun_chrom = (sort keys %huqun_chrom_name);
	my $number1 = @huqun_chrom;
	print "\t$number1 scffolds found!\n";
	close SNP;
	print "\tSTEP1 DONE!\n";
	print "\tSTEP2: merge samtools snps and huqun snps START!\n";
############################ STEP2 找到putative—snp对应的scf的上的samtools call出的SNP和INDEL
	my($chrom,$chrom2,@tmp,%number,$number2,$hh,$jishu);
	open OUTX ,'>','Un_finded_snps_In_Samtools_SNPs.log' or die "Cant create Un_finded_snps_In_Samtools_SNPs.log:$!";
	foreach $hh(@huqun_chrom){
	        open VCF,'<',"$vcf_file" or die "Cant open samtools vcf file:$!";
		while (<VCF>) {
			chomp;
			if (/#/){ next };
			if ((split)[0] =~ /$hh\z/) {  ### may need to change MARK1!!!
				push (@vcf_putativesnps,$_);
				$jishu = 1;
			}
		}
	        if ($jishu == 0) {
			print OUTX "$hh\n";
		}
	        close VCF;
	}
	close OUTX;
	#foreach (@vcf_putativesnps){
	#       $chrom = (split)[0];
	#       $chrom2 = (split /\|/,$chrom)[1];
	#       $number{$chrom2} = 1;
	#}
	#$number2 = (keys %number);
	#print "$number2\n";
	#$number2 = 0; %number =();
	print "\tSTEP2 DONE!\n";
}
print "\tGENOME FORMAT START!\n";
############################# STEP3 find the target scf and format it.
open GENOME,'<',"$genome_file" or die "Cant open genome sequence file(fasta format):$!";
#############  put genome sequence in a hash  ################
my ($counter,$scfname,$string,%genome);
$string = ();
while (<GENOME>) {
	chomp;
	if (/>/ && $counter == 0) { 
		$scfname = (split />/,(split)[0])[1];
		$counter = 1;
		next;
	}
	if (/>/ && $counter == 1) {
		$genome{$scfname} = $string;
		$scfname = (split />/,(split)[0])[1];
		$string = ();
		next;
	}
	if ($counter == 1) {
		$string .= $_;
	}
}
$genome{$scfname} = $string;
$scfname = ();
$string = ();
$counter = 0;
close GENOME;
#print "\tSTART OUTPUT TARGET SEQUENCE!\n";
#open OUT1,'>','Snp_target_scf.fasta';
#foreach (@huqun_chrom) {
#	print OUT1 ">$_\n$genome{$_}\n";
#}
#close OUT1;
print "\tGENOME FORMAT DONE!\n";
$counter = 0;
###########################################################################################################################

################################## Enzyme list (80 enzyme)####################################
my %Enzyme = (
	'AluI'	=>	'AGCT',
	'Bsh1236I'	=>	'CGCG',
	'BsuRI'	=>	'GGCC',
	'FspBI'	=>	'CTAG',
	'Hin1II'	=>	'CATG',
	'Hin6I'	=>	'GCGC',
	'HpaII'	=>	'CCGG',
	'MboI'	=>	'GATC',
	'RsaI'	=>	'GTAC',
	'TaiI'	=>	'ACGT',
	'TaqI'	=>	'TCGA',
	'TasI'	=>	'AATT',
	'Tru1I'	=>	'TTAA',
	'AanI'	=>	'TTATAA',
	'AatII'	=>	'GACGTC',
	'AjiI'	=>	'CACGTC',
	'Alw44I'	=>	'GTGCAC',
	'ApaI'	=>	'GGGCCC',
	'Acc65I'    =>    'GGTACC',
	'BamHI'	=>	'GGATCC',
	'BclI'	=>	'TGATCA',
	'BcuI'	=>	'ACTAGT',
	'BglII'	=>	'AGATCT',
	'BshTI'	=>	'ACCGGT',
	'Bsp119I'	=>	'TTCGAA',
	'Bsp1407I'	=>	'TGTACA',
	'Bsp68I'	=>	'TCGCGA',
	'BspTI'	=>	'CTTAAG',
	'Bst1107I'	=>	'GTATAC',
	'Bsu15I'	=>	'ATCGAT',
	'Cfr42I'	=>	'CCGCGG',
	'Cfr9I'	=>	'CCCGGG',
	'DraI'	=>	'TTTAAA',
	'Eco105I'	=>	'TACGTA',
	'Eco147I'	=>	'AGGCCT',
	'Eco32I'	=>	'GATATC',
	'Eco47III'	=>	'AGCGCT',
	'Eco52I'	=>	'CGGCCG',
	'Eco72I'	=>	'CACGTG',
	'EcoRI'	=>	'GAATTC',
	'EheI'	=>	'GGCGCC',
	'HindIII'	=>	'AAGCTT',
	'Kpn2I'	=>	'TCCGGA',
	'KpnI'	=>	'GGTACC',
	'KspAI'	=>	'GTTAAC',
	'MlsI'	=>	'TGGCCA',
	'MluI'	=>	'ACGCGT',
	'Mph1103I'	=>	'ATGCAT',
	'MunI'	=>	'CAATTG',
	'NcoI'	=>	'CCATGG',
	'NdeI'	=>	'CATATG',
	'NheI'	=>	'GCTAGC',
	'NsbI'	=>	'TGCGCA',
	'PaeI'	=>	'GCATGC',
	'PagI'	=>	'TCATGA',
	'PauI'	=>	'GCGCGC',
	'PdiI'	=>	'GCCGGC',
	'Pfl23II'	=>	'CGTACG',
	'PscI'	=>	'ACATGT',
	'Psp1406I'	=>	'AACGTT',
	'PstI'	=>	'CTGCAG',
	'PvuI'	=>	'CGATCG',
	'PvuII'	=>	'CAGCTG',
	'SacI'	=>	'GAGCTC',
	'SalI'	=>	'GTCGAC',
	'ScaI'	=>	'AGTACT',
	'SspI'	=>	'AATATT',
	'VspI'	=>	'ATTAAT',
	'XbaI'	=>	'TCTAGA',
	'XhoI'	=>	'CTCGAG',
	'XmaJI'	=>	'CCTAGG',
	'MreI'	=>	'CGCCGGCG',
	'MssI'	=>	'GTTTAAAC',
	'NotI'	=>	'GCGGCCGC',
	'MauBI'	=>	'CGCGCGCG',
	'PacI'	=>	'TTAATTAA',
	'SdaI'	=>	'CCTGCAGG',
	'SfaAI'	=>	'GCGATCGC',
	'SgrDI'	=>	'CGTCGACG',
	'SgsI'	=>	'GGCGCGCC',
	'SmiI'	=>	'ATTTAAAT',
	'AarI'      =>    'CACCTGC',
	'Alw26I'    =>    'GTCTC',
	'BauI'      =>    'CACGAG',
	'BfiI'      =>    'ACTGGG',
	'BfuI'      =>    'GTATCC',
	'Bpil'      =>    'GAAGAC',
	'BseGI'     =>    'GGATG',
	'BseMI'     =>    'GCAATG',
	'BseMII'    =>    'CTCAG',
	'BseNI'     =>    'ACTGG',
	'BspPI'     =>    'GGATC',
	'BveI'      =>    'AACTGC',
	'CseI'      =>    'GACGC',
	'Eam1104I'  =>    'CTCTTC',
	'Eco31I'    =>    'GGTCTC',
	'Eco57I'    =>    'CTGAAG',
	'Esp3I'     =>    'CGTCTC',
	'FaqI'      =>    'GGGAC',
	'GsuI'      =>    'CTGGAG',
	'HphI'      =>    'GGTGA',
	'LguI'      =>    'GCTCTTC',
	'lSP1109I'  =>    'GCAGC',
	'LweI'      =>    'GCATC',
	'MboII'     =>    'GAAGA',
	'MnII'      =>    'CCTC',
	'Mva1269I'  =>    'GAATGC',
	'SchI'      =>    'GAGTC',
	'SmuI'      =>    'CCCGC',
);
my %Alt_Enzyme = (
	'Hin4I'     =>    'GAYNNNNNVTC',
	'AasI'      =>    'GACNNNNNNGTC',
	'AdeI'      =>    'CACNNNGTG',
	'AjuI'      =>    'GAANNNNNNNTTGG',
	'AlfI'      =>    'GCANNNNNNTGC',
	'AloI'      =>    'GAACNNNNNNTCC',
	'Alw21I'    =>    'GWGCWC',
	'BcnI'      =>    'CCSGG',
	'BdaI'      =>    'TGANNNNNNTCA',
	'BfmI'      =>    'CTRYAG',
	'BglI'      =>    'GCCNNNNNGGC',
	'Bme1390I'  =>    'CCNGG',
	'BoxI'      =>    'GACNNNNGTC',
	'BpII'      =>    'GAGNNNNNCTC',
	'Bpu10I'    =>    'CCTNAGC',
	'Bpu1102I'  =>    'GCTNAGC',
	'BseDI'     =>    'CCNNGG',
	'BseJI'     =>    'GATNNNNATC',
	'BseLI'     =>    'CCNNNNNNNGG',
	'BseSI'     =>    'GKGCMC',
	'Bsh1285I'  =>    'CGRYCG',
	'BshNI'     =>    'GGYRCC',
	'BspLI'     =>    'GGNNCC',
	'BstXI'     =>    'CCANNNNNNTGG',
	'CaiI'      =>    'CAGNNNCTG',
	'CfrI'      =>    'YGGCCR',
	'Cfr10I'    =>    'RCCGGY',
	'Cfr13I'    =>    'GGNCC',
	'CpoI'      =>    'CGGWCCG',
	'Eam1105I'  =>    'GACNNNNNGTC',
	'Eco24I'    =>    'GRGCYC',
	'Eco47I'    =>    'GGWCC',
	'Eco57MI'   =>    'CTGRAG',
	'Eco81I'    =>    'CCTNAGG',
	'Eco88I'    =>    'CYCGRG',
	'Eco91I'    =>    'GGTNACC',
	'Eco130I'   =>    'CCWWGG',
	'EcoO109I'  =>    'RGGNCCY',
	'EcoRII'    =>    'CCWGG',
	'FspAI'     =>    'RTGCGCAY',
	'Hin1I'     =>    'GRCGYC',
	'HincII'    =>    'GTYRAC',
	'HinfI'     =>    'GANTC',
	'Hpy8I'     =>    'GTNNAC',
	'HpyF3I'    =>    'CTNAG',
	'HpyF10VI'  =>    'GCNNNNNNNGC',
	'NmuCI'     =>    'GTSAC',
	'OliI'      =>    'CACNNNNGTG',
	'PasI'      =>    'CCCWGGG',
	'PdmI'      =>    'GAANNNNTCC',
	'PfeI'      =>    'GAWTC',
	'PfoI'      =>    'TCCNGGA',
	'PpiI'      =>    'GAACNNNNNCTC',
	'Ppu21I'    =>    'YACGTR',
	'Psp51I'    =>    'RGGWCCY',
	'PsuI'      =>    'RGATCY',
	'PsyI'      =>    'GACNNNGTC',
	'RseI'      =>    'CAYNNNNRTC',
	'SatI'      =>    'GCNGC',
	'SduI'      =>    'GDGCHC',
	'SfiI'      =>    'GGCCNNNNNGGCC',
	'SmoI'      =>    'CTYRAG',
	'TaaI'      =>    'ACNGT',
	'TatI'      =>    'WGTACW',
	'TauI'      =>    'GCSGC',
	'TscAI'     =>    'NNCASTGNN',
	'TsoI'      =>    'TARCCA',
	'TstI'      =>    'CACNNNNNNTCC',
	'Van91I'    =>    'CCANNNNNTGG',
	'XagI'      =>    'CCTNNNNNAGG',
	'XapI'      =>    'RAATTY',
	'XceI'      =>    'RCATGY',
	'XmiI'      =>    'GTMKAC',
);
########################################################################################################
if (defined $opt_a) {
	my $Enzyme_number = (keys %Enzyme) + (keys %Alt_Enzyme);
	print "\tNow,create Enzyme list ($Enzyme_number enzyme)\n";
}else{
	my $Enzyme_number = (keys %Enzyme);
	print "\tNow,create Enzyme list ($Enzyme_number enzyme)\n";
}
############################ First step done ##########################
print "\tFinding the Enzyme recognition site\n";
open ENZYME,'>',"$vcf_file-Enzyme_site.txt";
open ENZYME2,'>',"$vcf_file-Better_Enzyme_site.txt";
my ($persentage,$iii,$dp1,$filter_parameter,$type,$qu,$dp4,$scf,$pos,$ref,$alt,$type1,$type2,$length,$string2,$sequence_for_filter,$sequence_for_filter2,$scflength,$EnzymeSite,$Enzyme_name,$altseq,$scflength,$seq1000,$pos2,$Enzyme_length,$match_number,$match_number2,$counter2,$true_pos,$mer_pos);
$persentage = 0.1;
my @Enzyme_list = ();
if (defined $opt_h) {
foreach (@vcf_putativesnps)  {
	chomp;
	if (/#/) { next };
        if ((split)[4] =~ /\,/) { next };
	if (/DP4=(\d+),(\d+),(\d+),(\d+)/) {                  ###--- if match DP4,the vcf come from samtools
		$dp4 = "DP4=$1,$2,$3,$4";
       	 	$type1 = (split /:/,(split)[9])[0];
        	$type2 = (split /:/,(split)[10])[0];
        	$qu = (split)[5];
        	if ( (($type1 =~ /1\/1/) && ($type2 =~ /1\/1/)) || $qu < 30 ) { next };
        	if (/DP4=(\d+),(\d+),(\d+),(\d+)/ && ((($1+$2) < $Depth) || (($3+$4) < $Depth))) {next};
	}else{                                                ###--- vcf come from GATK
		$type = (split /:/,(split)[9])[0];
		$filter_parameter = (split)[6];
		if ( ($type =~ /0\/1/) && ($filter_parameter =~ /PASS/) ) {
			$dp1 = (split /:/,(split)[9])[1];
			$dp4 = "DP=$dp1";
		}else{
			next;
		}
	}
	($scf,$mer_pos,$ref,$alt) = (split)[0,1,3,4];
	if ($mer_pos =~ /_/) {
		$pos = (split /_/,$mer_pos)[1];
		$true_pos = (split /_/,$mer_pos)[0];
	}else{
		$pos = $mer_pos;
		$true_pos = $pos;
	}
	if (!exists $genome{$scf}) { next };
        $scflength = length ($genome{$scf});
        if ($pos < 100 || (($scflength - $pos) < 100)) { next };
############# extract the alt position (+_3)
###############################################------ reference /# effective recognition sequence length = 4
	if ($pos < $seq_range_option && (($scflength - $pos) < $seq_range_option)) {
		$seq1000 = $genome{$scf};
		$pos2 = $pos;
	}elsif ($pos < $seq_range_option && (($scflength - $pos) > $seq_range_option)) {
		$seq1000 = substr($genome{$scf},0,($pos + $seq_range_option));
		$pos2 = $pos;
	}elsif ($pos > $seq_range_option && (($scflength - $pos) < $seq_range_option)) {
		$seq1000 = substr($genome{$scf},$pos - $seq_range_option);
		$pos2 = $seq_range_option;
	}else{
		$seq1000 = substr($genome{$scf},($pos - $seq_range_option),($seq_range_option * 2));
		$pos2 = $seq_range_option;
	}
##
	foreach (sort keys %Enzyme) {
		$EnzymeSite = $Enzyme{$_};
		$Enzyme_name = $_;
		$Enzyme_length = length($EnzymeSite);
		$string = substr($genome{$scf},($pos - $Enzyme_length),($Enzyme_length * 2 - 1));
		if ($string =~ /$EnzymeSite/) {
			$counter ++;
		}
		if ($counter == 0) {
			$string2 = &reverseCompliment($string);
			if ($string2 =~ /$EnzymeSite/) {
				$counter ++;
			}
		}
		$length = length($ref);
		$altseq = $genome{$scf};
		substr($altseq,($pos-1),$length) = $alt;
		substr($string,($Enzyme_length - 1),$length) = $alt;
		if ($string =~ /$EnzymeSite/) {
			$counter ++;
			$counter2 = 1
		}
		if ($counter2 == 0) {
			$string2 = &reverseCompliment($string);
			if ($string2 =~ /$_/) {
				$counter ++;
			}
		}
		if ($counter == 1) {
			print ENZYME ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
################################ try to find the better Enzyme site for developing marker! ###################
			$scflength = length ($genome{$scf});
			if (($pos > $range) && (($scflength - $pos) > $range)) {
				$sequence_for_filter = substr($genome{$scf},($pos - $range),($range * 2));
				$sequence_for_filter2 = substr($altseq,($pos - $range),($range * 2));
				$match_number = &match($sequence_for_filter,$EnzymeSite);
				$match_number2 = &match($sequence_for_filter2,$EnzymeSite);
				if ( ($match_number + $match_number2) < 2 ) {
					print ENZYME2 ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
				}
			}
		}
		$match_number = 0;
		$match_number2 = 0;
		$counter = 0;
		$counter2 = 0;
	}
	$counter = 0; ##!!!
####################################################################
	if (defined $opt_a) {
	foreach (sort keys %Alt_Enzyme) {
		$EnzymeSite = $Alt_Enzyme{$_};
		$Enzyme_name = $_;
		$Enzyme_length = length($EnzymeSite);
		$string = substr($genome{$scf},($pos - $Enzyme_length),($Enzyme_length * 2 - 1));
		@Enzyme_list = &change($EnzymeSite);
		foreach (@Enzyme_list) {
			if ($string =~ /$_/) {
				$counter ++;
			}
		}
		if ($counter == 0) {
			$string2 = &reverseCompliment($string);
			foreach (@Enzyme_list) {
				if ($string2 =~ /$_/) {
					$counter ++;
				}
			}
		}
		$length = length($ref);
		$altseq = $genome{$scf};
		substr($altseq,($pos-1),$length) = $alt;
		substr($string,($Enzyme_length - 1),$length) = $alt;
		foreach (@Enzyme_list) {
			if ($string =~ /$_/) {
				$counter ++;
				$counter2 = 1
			}
		}
		if ($counter2 == 0) {
			$string2 = &reverseCompliment($string);
			foreach (@Enzyme_list) {
				if ($string2 =~ /$_/) {
					$counter ++;
				}
			}
		}
		if ($counter == 1) {
			print ENZYME ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
################################ try to find the better Enzyme site for developing marker! ###################
			$scflength = length ($genome{$scf});
			if (($pos > $range) && (($scflength - $pos) > $range)) {
				$sequence_for_filter = substr($genome{$scf},($pos - $range),($range * 2));
				$sequence_for_filter2 = substr($altseq,($pos - $range),($range * 2));
				foreach (@Enzyme_list) {
					$match_number += &match($sequence_for_filter,$_);
					$match_number2 += &match($sequence_for_filter2,$_);
				}
				if ( ($match_number + $match_number2) < 2 ) {
					print ENZYME2 ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
				}
			}
		}
		$match_number = 0;
		$match_number2 = 0;
		$counter = 0;
		$counter2 = 0;
	}
	}
}
print "\tDONE!\n";
}else{
open VCF,'<',"$vcf_file" or die "Cant open samtools vcf file:$!";
while (<VCF>)  {
	$iii ++;
	if ( ($iii / $vcf_line) > $persentage ) {
		my $per = $persentage * 100;
		print "\t$per% of prediction comleted!\n";
		$persentage += 0.1;
	}
	chomp;
	if (/#/) { next };
        if ((split)[4] =~ /\,/) { next };
        if (/DP4=(\d+),(\d+),(\d+),(\d+)/) {                  ###--- if match DP4,the vcf is comefrom samtools
                $dp4 = "DP4=$1,$2,$3,$4";
                $type1 = (split /:/,(split)[9])[0];
                $type2 = (split /:/,(split)[10])[0];
                $qu = (split)[5];
                if ( (($type1 =~ /1\/1/) && ($type2 =~ /1\/1/)) || $qu < 30 ) { next };
                if (/DP4=(\d+),(\d+),(\d+),(\d+)/ && ((($1+$2) < $Depth) || (($3+$4) < $Depth))) {next};
        }else{                                                ###--- vcf comefrom GATK
                $type = (split /:/,(split)[9])[0];
                $filter_parameter = (split)[6];
                if ( ($type =~ /0\/1/) && ($filter_parameter =~ /PASS/) ) {
                        $dp1 = (split /:/,(split)[9])[1];
                        $dp4 = "DP=$dp1";
                }else{
                        next;
                }
        }
	($scf,$mer_pos,$ref,$alt) = (split)[0,1,3,4];
	if ($mer_pos =~ /_/) {
		$pos = (split /_/,$mer_pos)[1];
		$true_pos = (split /_/,$mer_pos)[0];
	}else{
		$pos = $mer_pos;
		$true_pos = $pos;
	}
	if (!exists $genome{$scf}) { next };
        $scflength = length ($genome{$scf});
        if ($pos < 100 || (($scflength - $pos) < 100)) { next };
############# extract the alt position (+_3)
###############################################------ reference /# effective recognition sequence length = 4
	if ($pos < $seq_range_option && (($scflength - $pos) < $seq_range_option)) {
		$seq1000 = $genome{$scf};
		$pos2 = $pos;
	}elsif ($pos < $seq_range_option && (($scflength - $pos) > $seq_range_option)) {
		$seq1000 = substr($genome{$scf},0,($pos + $seq_range_option));
		$pos2 = $pos;
	}elsif ($pos > $seq_range_option && (($scflength - $pos) < $seq_range_option)) {
		$seq1000 = substr($genome{$scf},$pos - $seq_range_option);
		$pos2 = $seq_range_option;
	}else{
		$seq1000 = substr($genome{$scf},($pos - $seq_range_option),($seq_range_option * 2));
		$pos2 = $seq_range_option;
	}
##
	foreach (sort keys %Enzyme) {
		$EnzymeSite = $Enzyme{$_};
		$Enzyme_name = $_;
		$Enzyme_length = length($EnzymeSite);
		$string = substr($genome{$scf},($pos - $Enzyme_length),($Enzyme_length * 2 - 1));
		if ($string =~ /$EnzymeSite/) {
			$counter ++;
		}
		if ($counter == 0) {
			$string2 = &reverseCompliment($string);
			if ($string2 =~ /$EnzymeSite/) {
				$counter ++;
			}
		}
		$length = length($ref);
		$altseq = $genome{$scf};
		substr($altseq,($pos-1),$length) = $alt;
		substr($string,($Enzyme_length - 1),$length) = $alt;
		if ($string =~ /$EnzymeSite/) {
			$counter ++;
			$counter2 = 1
		}
		if ($counter2 == 0) {
			$string2 = &reverseCompliment($string);
			if ($string2 =~ /$_/) {
				$counter ++;
			}
		}
		if ($counter == 1) {
			print ENZYME ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
################################ try to find the better Enzyme site for developing marker! ###################
			$scflength = length ($genome{$scf});
			if (($pos > $range) && (($scflength - $pos) > $range)) {
				$sequence_for_filter = substr($genome{$scf},($pos - $range),($range * 2));
				$sequence_for_filter2 = substr($altseq,($pos - $range),($range * 2));
				$match_number = &match($sequence_for_filter,$EnzymeSite);
				$match_number2 = &match($sequence_for_filter2,$EnzymeSite);
				if ( ($match_number + $match_number2) < 2 ) {
					print ENZYME2 ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
				}			
			}
		}
		$match_number = 0;
		$match_number2 = 0;
		$counter = 0;
		$counter2 = 0;
	}
	$counter = 0; ##!!!
##################################----- start alt enzyme site finding -----##############################
	if (defined $opt_a) {
	foreach (sort keys %Alt_Enzyme) {
		$EnzymeSite = $Alt_Enzyme{$_};
		$Enzyme_name = $_;
		$Enzyme_length = length($EnzymeSite);
		$string = substr($genome{$scf},($pos - $Enzyme_length),($Enzyme_length * 2 - 1));
		@Enzyme_list = &change($EnzymeSite);
		foreach (@Enzyme_list) {
			if ($string =~ /$_/) {
				$counter ++;
			}
		}
		if ($counter == 0) {
			$string2 = &reverseCompliment($string);
			foreach (@Enzyme_list) {
				if ($string2 =~ /$_/) {
					$counter ++;
				}
			}
		}
		$length = length($ref);
		$altseq = $genome{$scf};
		substr($altseq,($pos-1),$length) = $alt;
		substr($string,($Enzyme_length - 1),$length) = $alt;
		foreach (@Enzyme_list) {
			if ($string =~ /$_/) {
				$counter ++;
				$counter2 = 1
			}
		}
		if ($counter2 == 0) {
			$string2 = &reverseCompliment($string);
			foreach (@Enzyme_list) {
				if ($string2 =~ /$_/) {
					$counter ++;
				}
			}
		}
		if ($counter == 1) {
			print ENZYME ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
################################ try to find the better Enzyme site for developing marker! ###################
			$scflength = length ($genome{$scf});
			if (($pos > $range) && (($scflength - $pos) > $range)) {
				$sequence_for_filter = substr($genome{$scf},($pos - $range),($range * 2));
				$sequence_for_filter2 = substr($altseq,($pos - $range),($range * 2));
				foreach (@Enzyme_list) {
					$match_number += &match($sequence_for_filter,$_);
					$match_number2 += &match($sequence_for_filter2,$_);
				}
				if ( ($match_number + $match_number2) < 2 ) {
					print ENZYME2 ">$scf:$pos2:$true_pos:$Enzyme_name:$EnzymeSite:$dp4\n$seq1000\n";
				}			
			}
		}
		$match_number = 0;
		$match_number2 = 0;
		$counter = 0;
		$counter2 = 0;
	}
	}
}
print "\tDONE!\n";
close VCF;
}
####################### SUB PROGREM #####################
sub reverseCompliment{
	my ($origin,$reverse);
	$origin = shift @_;
	$reverse = reverse $origin;
	$reverse =~ s/A/Q/g;
	$reverse =~ s/T/W/g;
	$reverse =~ s/G/E/g;
	$reverse =~ s/C/R/g;
	$reverse =~ s/Q/T/g;
	$reverse =~ s/W/A/g;
	$reverse =~ s/E/C/g;
	$reverse =~ s/R/G/g;
	return $reverse;
}

sub match {
	my ($match,$sequence_match,$Enzymesite,$tt,$seq1,$seq2,$len,$end);
	$sequence_match = shift @_;
	$Enzymesite = shift @_;
	$seq1 = $sequence_match;
	while (1) {
		if ($seq1 =~ /$Enzymesite/g) {
			$len = length($Enzymesite);
			$end = pos($seq1);
			$tt = $end - $len;
			$match ++;
			$seq1 = substr($seq1,($tt + 1));
		}else{
			last;
		}
	}
	if ($Enzymesite eq (&reverseCompliment($Enzymesite))) {
		return $match;
	}
	$seq1 = $sequence_match;
	$seq2 = &reverseCompliment($Enzymesite);
	while (1) {
		if ($seq1 =~ /$seq2/g) {
			$len = length($seq2);
			$end = pos($seq1);
			$tt = $end - $len;
			$match ++;
			$seq1 = substr($seq1,($tt + 1));
		}else{
			last;
		}
	}
	return $match;
}

#################################---sub change---#########################################
sub change {
	my ($chN,@changing,@changedN,@changedY,@changedW,@changedS,@changedR,@changedM,@changedK,@changedV,@changedH,@changedD,@changedB) = ();
	foreach (@_) {
		if (/N/) {
			($chN = $_) =~ s/N/\./g;
			push @changedN,$chN;
		}else{
			push @changedN,$_;
		}
	}
###----- change Y ----###
	foreach (@changedN) {
		if (/Y/) {
			@changing = &changeY($_);
			push @changedY,@changing;
		}else{
			push @changedY,$_;
		}
	}
###---- change W ----###
	foreach (@changedY) {
		if (/W/) {
			@changing = &changeW($_);
			push @changedW,@changing;
		}else{
			push @changedW,$_;
		}
	}
###---- change S ----###
	foreach (@changedW) {
		if (/W/) {
			@changing = &changeS($_);
			push @changedS,@changing;
		}else{
			push @changedS,$_;
		}
	}
###---- change R ----###
	foreach (@changedS) {
		if (/R/) {
			@changing = &changeR($_);
			push @changedR,@changing;
		}else{
			push @changedR,$_;
		}
	}
###---- change M ----###
	foreach (@changedR) {
		if (/M/) {
			@changing = &changeM($_);
			push @changedM,@changing;
		}else{
			push @changedM,$_;
		}
	}
###---- change K ----###
	foreach (@changedM) {
		if (/K/) {
			@changing = &changeK($_);
			push @changedK,@changing;
		}else{
			push @changedK,$_;
		}
	}
###---- change V ----###
	foreach (@changedK) {
		if (/V/) {
			@changing = &changeV($_);
			push @changedV,@changing;
		}else{
			push @changedV,$_;
		}
	}
###---- change H ----###
	foreach (@changedV) {
		if (/H/) {
			@changing = &changeH($_);
			push @changedH,@changing;
		}else{
			push @changedH,$_;
		}
	}
###---- change D ----###
	foreach (@changedH) {
		if (/D/) {
			@changing = &changeD($_);
			push @changedD,@changing;
		}else{
			push @changedD,$_;
		}
	}
###---- change B ----###
	foreach (@changedD) {
		if (/B/) {
			@changing = &changeB($_);
			push @changedB,@changing;
		}else{
			push @changedB,$_;
		}
	}
	return @changedB;
}

sub changeY {
	my ($origin_site,@site,@y_site,$y1,$y2);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/Y/) {
				$origin_site = $_;
				($y1 = $origin_site) =~ s/Y/C/;
				($y2 = $origin_site) =~ s/Y/T/;
				push @y_site,$y1;
				push @y_site,$y2;
			}
		}
		@site = ();
		foreach (@y_site) {
			if (/Y/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@y_site = ();
		}
	}
	return @y_site;
}

sub changeW {
	my ($origin_site,@site,@w_site,$w1,$w2);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/W/) {
				$origin_site = $_;
				($w1 = $origin_site) =~ s/W/A/;
				($w2 = $origin_site) =~ s/W/T/;
				push @w_site,$w1;
				push @w_site,$w2;
			}
		}
		@site = ();
		foreach (@w_site) {
			if (/W/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@w_site = ();
		}
	}
	return @w_site;
}

sub changeS {
	my ($origin_site,@site,@s_site,$s1,$s2);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/S/) {
				$origin_site = $_;
				($s1 = $origin_site) =~ s/S/C/;
				($s2 = $origin_site) =~ s/S/G/;
				push @s_site,$s1;
				push @s_site,$s2;
			}
		}
		@site = ();
		foreach (@s_site) {
			if (/S/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@s_site = ();
		}
	}
	return @s_site;
}

sub changeR {
	my ($origin_site,@site,@r_site,$r1,$r2);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/R/) {
				$origin_site = $_;
				($r1 = $origin_site) =~ s/R/A/;
				($r2 = $origin_site) =~ s/R/G/;
				push @r_site,$r1;
				push @r_site,$r2;
			}
		}
		@site = ();
		foreach (@r_site) {
			if (/R/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@r_site = ();
		}
	}
	return @r_site;
}

sub changeM {
	my ($origin_site,@site,@m_site,$m1,$m2);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/M/) {
				$origin_site = $_;
				($m1 = $origin_site) =~ s/M/A/;
				($m2 = $origin_site) =~ s/M/C/;
				push @m_site,$m1;
				push @m_site,$m2;
			}
		}
		@site = ();
		foreach (@m_site) {
			if (/M/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@m_site = ();
		}
	}
	return @m_site;
}

sub changeK {
	my ($origin_site,@site,@k_site,$k1,$k2);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/K/) {
				$origin_site = $_;
				($k1 = $origin_site) =~ s/K/G/;
				($k2 = $origin_site) =~ s/K/T/;
				push @k_site,$k1;
				push @k_site,$k2;
			}
		}
		@site = ();
		foreach (@k_site) {
			if (/K/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@k_site = ();
		}
	}
	return @k_site;
}

sub changeV {
	my ($origin_site,@site,@v_site,$v1,$v2,$v3);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/V/) {
				$origin_site = $_;
				($v1 = $origin_site) =~ s/V/C/;
				($v2 = $origin_site) =~ s/V/G/;
				($v3 = $origin_site) =~ s/V/A/;
				push @v_site,$v1;
				push @v_site,$v2;
				push @v_site,$v3;
			}
		}
		@site = ();
		foreach (@v_site) {
			if (/V/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@v_site = ();
		}
	}
	return @v_site;
}

sub changeH {
	my ($origin_site,@site,@h_site,$h1,$h2,$h3);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/H/) {
				$origin_site = $_;
				($h1 = $origin_site) =~ s/H/C/;
				($h2 = $origin_site) =~ s/H/T/;
				($h3 = $origin_site) =~ s/H/A/;
				push @h_site,$h1;
				push @h_site,$h2;
				push @h_site,$h3;
			}
		}
		@site = ();
		foreach (@h_site) {
			if (/H/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@h_site = ();
		}
	}
	return @h_site;
}

sub changeD {
	my ($origin_site,@site,@d_site,$d1,$d2,$d3);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/D/) {
				$origin_site = $_;
				($d1 = $origin_site) =~ s/D/C/;
				($d2 = $origin_site) =~ s/D/G/;
				($d3 = $origin_site) =~ s/D/A/;
				push @d_site,$d1;
				push @d_site,$d2;
				push @d_site,$d3;
			}
		}
		@site = ();
		foreach (@d_site) {
			if (/D/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@d_site = ();
		}
	}
	return @d_site;
}

sub changeB {
	my ($origin_site,@site,@b_site,$b1,$b2,$b3);
	my $recycle = 0;
	@site = @_;
	while (1) {
		foreach (@site) {
			if (/B/) {
				$origin_site = $_;
				($b1 = $origin_site) =~ s/B/C/;
				($b2 = $origin_site) =~ s/B/G/;
				($b3 = $origin_site) =~ s/B/A/;
				push @b_site,$b1;
				push @b_site,$b2;
				push @b_site,$b3;
			}
		}
		@site = ();
		foreach (@b_site) {
			if (/B/) {
				push @site,$_;
				$recycle = 1;
			}else{
				$recycle = 0;
			}
		}
		if ($recycle == 0) {
			last;
		}else{
			@b_site = ();
		}
	}
	return @b_site;
}
