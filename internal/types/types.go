package types

import "context"

const (
	INPUT         = "/data/input/raw_data"
	OUTPUT        = "/data/output"
	LOG           = "/data/log"
	FASTP_OUT     = "/data/output/clean"
	REFERENCES    = "/data/input/references"
	GENOME_PREFIX = "/data/input/index/Species"
	GENOME_INDEX  = "/data/input/index"
	REPORT_OUT    = "/data/output/report_result"
	BSA_DNA_OUT   = "/data/output/bsa_result"
	HISAT2_OUT    = "/data/output/hisat2_result"
	SORTED_OUT    = "/data/output/sorted_result"
	SINGLE_OUT    = "/data/output/single_result"
	GATK_OUT      = "/data/output/gatk_result"
	GATK_G_OUT    = "/data/output/gatk_vcf_result"
	//DIFF_OUT = "/data/output/diff_result"
	KEGG_OUT = "/data/output/kegg_result"
	GO_OUT   = "/data/output/go_result"
	//expression
	EXPRESSION_OUT = "/data/output/expression_result"

	BSA_GENOME_PREFIX = "/data/input/index/Pepper.chr"
)

type Plugin interface {
	Build(ctx context.Context) error
	Name() string
}

const (
	MaxIdleConns        int = 100
	MaxIdleConnsPerHost int = 100
	IdleConnTimeout     int = 90
)
const (
	R1Sample          = "_R1.fastq.gz"
	R1SampleEx        = "_R1.fq.gz"
	R1SampleFq        = "_1.fastq.gz"
	R1SampleFqEx      = "_1.fq.gz"
	R1SampleClean     = "_R1.clean.fastq.gz"
	R1SampleCleanEx   = "_R1.clean.fq.gz"
	R1SampleCleanFq   = "_1.clean.fastq.gz"
	R1SampleCleanFqEx = "_1.clean.fq.gz"
)

type SampleType uint8

func (s SampleType) Type() string {
	switch s {
	case 1:
		return R1Sample
	case 2:
		return R1SampleEx
	case 3:
		return R1SampleFq
	case 4:
		return R1SampleFqEx
	}
	return ""
}
func (s SampleType) CleanType() string {
	switch s {
	case 1:
		return R1SampleClean
	case 2:
		return R1SampleCleanEx
	case 3:
		return R1SampleCleanFq
	case 4:
		return R1SampleCleanFqEx
	}
	return ""
}
