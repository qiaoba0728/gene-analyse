#!/usr/bin/perl

use strict;

my ($bam1,$bam2,$gtf) = @ARGV[0,1,2];
my $Usage = "\n\t$0 <bam1> <bam2> <gtf>
\n";
die $Usage unless (@ARGV == 3);

# transfer bam to sam #
!system "samtools view $bam1 > /data/output/sample1.sam" or die "ERROR with samtools:$!";
!system "gfold count -ann $gtf -tag /data/output/sample1.sam -o /data/output/sample1.read_cnt" or die "ERROR with gfold:$!";
#`rm /data/output/sample1.sam -rf`;


!system "samtools view $bam2 > /data/output/sample2.sam" or die "ERROR with samtools:$!";
!system "gfold count -ann $gtf -tag /data/output/sample2.sam -o /data/output/sample2.read_cnt" or die "ERROR with gfold:$!";
#`rm /data/output/sample2.sam -rf`;

!system "gfold diff -s1 /data/output/sample1 -s2 /data/output/sample2 -suf .read_cnt -o /data/output/sample1VSsample2.diff" or die "ERROR with gfold diff:$!" ;


