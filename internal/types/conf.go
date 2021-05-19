package types

type SummaryMeta struct {
	TotalReads      float64 `json:"total_reads"`
	TotalBases      float64 `json:"total_bases"`
	Q20Bases        float64 `json:"q20_bases"`
	Q30Bases        float64 `json:"q30_bases"`
	Q20Rate         float64 `json:"q20_rate"`
	Q30Rate         float64 `json:"q30_rate"`
	Read1MeanLength int     `json:"read1_mean_length"`
	Read2MeanLength int     `json:"read2_mean_length"`
	GCContent       float64 `json:"gc_content"`
}
type Summary struct {
	BeforeFiltering *SummaryMeta `json:"before_filtering"`
	AfterFiltering  *SummaryMeta `json:"after_filtering"`
}
type QualityContentCurves struct {
	QualityCurves *QualityCurves `json:"quality_curves"`
	ContentCurves *ContentCurves `json:"content_curves"`
}
type SummaryDocs struct {
	Summary     *Summary              `json:"summary"`
	Read1Before *QualityContentCurves `json:"read1_before_filtering"`
	Read2Before *QualityContentCurves `json:"read2_before_filtering"`

	Read1After *QualityContentCurves `json:"read1_after_filtering"`
	Read2After *QualityContentCurves `json:"read2_after_filtering"`
}

type QualityCurves struct {
	A    []float64 `json:"A"`
	T    []float64 `json:"T"`
	C    []float64 `json:"C"`
	G    []float64 `json:"G"`
	Mean []float64 `json:"mean"`
}
type ContentCurvesMeta struct {
	ContentCurves *ContentCurves `json:"content_curves"`
}
type ContentCurves struct {
	A  []float64 `json:"A"`
	T  []float64 `json:"T"`
	C  []float64 `json:"C"`
	G  []float64 `json:"G"`
	N  []float64 `json:"N"`
	GC []float64 `json:"GC"`
}

type Project struct {
	Name         string `json:"name"`
	Number       string `json:"number"`
	Organization string `json:"organization"`
	Time         string `json:"time"`
}
type ProjectDesc struct {
	Number       string `json:"number"`
	Proposal     string `json:"proposal"`
	Type         string `json:"type"`
	Species      string `json:"species"`
	SampleForm   string `json:"sampleForm"`
	SampleNumber string `json:"sampleNumber"`
	Platform     string `json:"platform"`
	DataSize     string `json:"data_size"`
	Analyse      string `json:"analyse"`
	Finished     string `json:"finished"`
}
type Text struct {
	Project     Project     `json:"project"`
	ProjectDesc ProjectDesc `json:"project_desc"`
}
type SampleItem struct {
	Type string `json:"type"`
}
type Tables struct {
}
type Images struct {
}
