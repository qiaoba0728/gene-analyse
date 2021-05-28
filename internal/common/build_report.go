package common

import (
	"fmt"
	"github.com/jszwec/csvutil"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"io/ioutil"
	"log"
	"strings"
)

const (
	SummaryTemplate = "|%s|%s|%s|%s|%s|%s|%s|"
)
const (
	reportBeforeContentTemplate = `<details><summary>%s</summary><center><img src="../asset/%s_content_before.png">注：横坐标是Reads中碱基位置（5'->3'），横坐标是该位点某碱基含量</center></details>`
	reportAfterContentTemplate  = `<details><summary>%s</summary><center><img src="../asset/%s_content_after.png">注：横坐标是Reads中碱基位置（5'->3'），横坐标是该位点某碱基含量</center></details>`

	reportBeforeQualityTemplate = `<details><summary>%s</summary><center><img src="../asset/%s_quality_before.png">注：横坐标是Reads中碱基位置（5'->3'），横坐标是该位点某碱基质量值</center></details>`
	reportAfterQualityTemplate  = `<details><summary>%s</summary><center><img src="../asset/%s_quality_after.png">注：横坐标是Reads中碱基位置（5'->3'），横坐标是该位点某碱基质量值</center></details>`
)

type ReportSummaryMeta struct {
	TotalReads      float64 `csv:"total_reads"`
	TotalBases      float64 `csv:"total_bases"`
	Q20Bases        float64 `csv:"q20_bases"`
	Q30Bases        float64 `csv:"q30_bases"`
	Q20Rate         float64 `csv:"q20_rate"`
	Q30Rate         float64 `csv:"q30_rate"`
	Read1MeanLength int     `csv:"read1_mean_length"`
	Read2MeanLength int     `csv:"read2_mean_length"`
	GCContent       float64 `csv:"gc_content"`
	Name            string  `csv:"name"`
}
type ReportQualityCurves struct {
	Index int     `csv:"index"`
	A     float64 `csv:"A"`
	T     float64 `csv:"T"`
	C     float64 `csv:"C"`
	G     float64 `csv:"G"`
	Mean  float64 `csv:"mean"`
}
type ReportContentCurves struct {
	Index int     `csv:"index"`
	A     float64 `csv:"A"`
	T     float64 `csv:"T"`
	C     float64 `csv:"C"`
	G     float64 `csv:"G"`
	N     float64 `csv:"N"`
	GC    float64 `csv:"GC"`
}

func BuildReport() error {
	beforeSummary := make([]*ReportSummaryMeta, 0)
	afterSummary := make([]*ReportSummaryMeta, 0)
	read1BeforeQualityCurves := make([]*ReportQualityCurves, 0)
	read2BeforeQualityCurves := make([]*ReportQualityCurves, 0)
	read1AfterQualityCurves := make([]*ReportQualityCurves, 0)
	read2AfterQualityCurves := make([]*ReportQualityCurves, 0)

	read1BeforeContentCurves := make([]*ReportContentCurves, 0)
	read2BeforeContentCurves := make([]*ReportContentCurves, 0)
	read1AfterContentCurves := make([]*ReportContentCurves, 0)
	read2AfterContentCurves := make([]*ReportContentCurves, 0)

	reportBeforeContentTemplates := make([]string, 0)
	reportAfterContentTemplates := make([]string, 0)
	reportBeforeQualityTemplates := make([]string, 0)
	reportAfterQualityTemplates := make([]string, 0)
	files, _ := ioutil.ReadDir(types.FASTP_OUT)
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".json") {
			docs, err := conf.GetSummary(fmt.Sprintf("%s/%s", types.FASTP_OUT, v.Name()))
			if err != nil {
				return err
			}
			name := strings.TrimSuffix(v.Name(), ".json")
			before := &ReportSummaryMeta{
				TotalReads:      docs.Summary.BeforeFiltering.TotalReads,
				TotalBases:      docs.Summary.BeforeFiltering.TotalBases,
				Q20Bases:        docs.Summary.BeforeFiltering.Q20Bases,
				Q30Bases:        docs.Summary.BeforeFiltering.Q30Bases,
				Q20Rate:         docs.Summary.BeforeFiltering.Q20Rate,
				Q30Rate:         docs.Summary.BeforeFiltering.Q30Rate,
				Read1MeanLength: docs.Summary.BeforeFiltering.Read1MeanLength,
				Read2MeanLength: docs.Summary.BeforeFiltering.Read2MeanLength,
				GCContent:       docs.Summary.BeforeFiltering.GCContent,
				Name:            name,
			}
			beforeSummary = append(beforeSummary, before)
			after := &ReportSummaryMeta{
				TotalReads:      docs.Summary.AfterFiltering.TotalReads,
				TotalBases:      docs.Summary.AfterFiltering.TotalBases,
				Q20Bases:        docs.Summary.AfterFiltering.Q20Bases,
				Q30Bases:        docs.Summary.AfterFiltering.Q30Bases,
				Q20Rate:         docs.Summary.AfterFiltering.Q20Rate,
				Q30Rate:         docs.Summary.AfterFiltering.Q30Rate,
				Read1MeanLength: docs.Summary.AfterFiltering.Read1MeanLength,
				Read2MeanLength: docs.Summary.AfterFiltering.Read2MeanLength,
				GCContent:       docs.Summary.AfterFiltering.GCContent,
				Name:            name,
			}
			afterSummary = append(afterSummary, after)
			for i := 0; i < len(docs.Read1Before.QualityCurves.A); i++ {
				read1BeforeQuality := &ReportQualityCurves{
					Index: i,
					A:     docs.Read1Before.QualityCurves.A[i],
					T:     docs.Read1Before.QualityCurves.T[i],
					C:     docs.Read1Before.QualityCurves.C[i],
					G:     docs.Read1Before.QualityCurves.G[i],
					Mean:  docs.Read1Before.QualityCurves.Mean[i],
				}
				read2BeforeQuality := &ReportQualityCurves{
					Index: i,
					A:     docs.Read2Before.QualityCurves.A[i],
					T:     docs.Read2Before.QualityCurves.T[i],
					C:     docs.Read2Before.QualityCurves.C[i],
					G:     docs.Read2Before.QualityCurves.G[i],
					Mean:  docs.Read2Before.QualityCurves.Mean[i],
				}
				read1AfterQuality := &ReportQualityCurves{
					Index: i,
					A:     docs.Read2After.QualityCurves.A[i],
					T:     docs.Read2After.QualityCurves.T[i],
					C:     docs.Read2After.QualityCurves.C[i],
					G:     docs.Read2After.QualityCurves.G[i],
					Mean:  docs.Read2After.QualityCurves.Mean[i],
				}
				read2AftereQuality := &ReportQualityCurves{
					Index: i,
					A:     docs.Read2After.QualityCurves.A[i],
					T:     docs.Read2After.QualityCurves.T[i],
					C:     docs.Read2After.QualityCurves.C[i],
					G:     docs.Read2After.QualityCurves.G[i],
					Mean:  docs.Read2After.QualityCurves.Mean[i],
				}
				read1BeforeQualityCurves = append(read1BeforeQualityCurves, read1BeforeQuality)
				read2BeforeQualityCurves = append(read2BeforeQualityCurves, read2BeforeQuality)
				read1AfterQualityCurves = append(read1AfterQualityCurves, read1AfterQuality)
				read2AfterQualityCurves = append(read2AfterQualityCurves, read2AftereQuality)
				read1BeforeContent := &ReportContentCurves{
					Index: i,
					A:     docs.Read1Before.ContentCurves.A[i],
					T:     docs.Read1Before.ContentCurves.T[i],
					C:     docs.Read1Before.ContentCurves.C[i],
					G:     docs.Read1Before.ContentCurves.G[i],
					GC:    docs.Read1Before.ContentCurves.GC[i],
				}
				read2BeforeContent := &ReportContentCurves{
					Index: i,
					A:     docs.Read2Before.ContentCurves.A[i],
					T:     docs.Read2Before.ContentCurves.T[i],
					C:     docs.Read2Before.ContentCurves.C[i],
					G:     docs.Read2Before.ContentCurves.G[i],
					GC:    docs.Read2Before.ContentCurves.GC[i],
				}
				read1AfterContent := &ReportContentCurves{
					Index: i,
					A:     docs.Read1After.ContentCurves.A[i],
					T:     docs.Read1After.ContentCurves.T[i],
					C:     docs.Read1After.ContentCurves.C[i],
					G:     docs.Read1After.ContentCurves.G[i],
					GC:    docs.Read1After.ContentCurves.GC[i],
				}
				read2AfterContent := &ReportContentCurves{
					Index: i,
					A:     docs.Read2After.ContentCurves.A[i],
					T:     docs.Read2After.ContentCurves.T[i],
					C:     docs.Read2After.ContentCurves.C[i],
					G:     docs.Read2After.ContentCurves.G[i],
					GC:    docs.Read2After.ContentCurves.GC[i],
				}
				read1BeforeContentCurves = append(read1BeforeContentCurves, read1BeforeContent)
				read2BeforeContentCurves = append(read2BeforeContentCurves, read2BeforeContent)
				read1AfterContentCurves = append(read1AfterContentCurves, read1AfterContent)
				read2AfterContentCurves = append(read2AfterContentCurves, read2AfterContent)
			}
			b, err := csvutil.Marshal(read1BeforeQualityCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read1BeforeQualityCurves.csv"), b, 0644)
			b, err = csvutil.Marshal(read2BeforeQualityCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read2BeforeQualityCurves.csv"), b, 0644)
			b, err = csvutil.Marshal(read1AfterQualityCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read1AfterQualityCurves.csv"), b, 0644)
			b, err = csvutil.Marshal(read2AfterQualityCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read2AfterQualityCurves.csv"), b, 0644)

			b, err = csvutil.Marshal(read1BeforeContentCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read1BeforeContentCurves.csv"), b, 0644)
			b, err = csvutil.Marshal(read2BeforeContentCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read2BeforeContentCurves.csv"), b, 0644)
			b, err = csvutil.Marshal(read1AfterContentCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read1AfterContentCurves.csv"), b, 0644)
			b, err = csvutil.Marshal(read2AfterContentCurves)
			if err != nil {
				log.Println("csv marshal error")
				return err
			}
			_ = ioutil.WriteFile(fmt.Sprintf("%s/%s-%s", types.FASTP_OUT, name, "read2AfterContentCurves.csv"), b, 0644)

			//build html]
			temp := fmt.Sprintf(reportBeforeContentTemplate, name, name)
			reportBeforeContentTemplates = append(reportBeforeContentTemplates, temp)
			temp = fmt.Sprintf(reportAfterContentTemplate, name, name)
			reportAfterContentTemplates = append(reportAfterContentTemplates, temp)
			temp = fmt.Sprintf(reportBeforeQualityTemplate, name, name)
			reportBeforeQualityTemplates = append(reportBeforeQualityTemplates, temp)
			temp = fmt.Sprintf(reportAfterQualityTemplate, name, name)
			reportAfterQualityTemplates = append(reportAfterQualityTemplates, temp)
		}
	}
	result := make([]string, 0)
	for _, v := range beforeSummary {
		temp := fmt.Sprintf(SummaryTemplate, v.Name, utils.ToString(v.TotalReads), utils.ToString(v.TotalBases), utils.ToString(v.Q20Bases), utils.ToString(v.Q30Bases), fmt.Sprintf("%.4f%%", v.Q20Rate*100.0), fmt.Sprintf("%.4f%%", v.Q30Rate*100.0))
		result = append(result, temp)
	}
	b := strings.Join(result, "\n")
	_ = ioutil.WriteFile(fmt.Sprintf("%s/%s", types.FASTP_OUT, "beforeSummary.template"), []byte(b), 0644)
	result = make([]string, 0)
	for _, v := range afterSummary {
		temp := fmt.Sprintf(SummaryTemplate, v.Name, utils.ToString(v.TotalReads), utils.ToString(v.TotalBases), utils.ToString(v.Q20Bases), utils.ToString(v.Q30Bases), fmt.Sprintf("%.4f%%", v.Q20Rate*100.0), fmt.Sprintf("%.4f%%", v.Q30Rate*100.0))
		result = append(result, temp)
	}
	b = strings.Join(result, "\n")
	_ = ioutil.WriteFile(fmt.Sprintf("%s/%s", types.FASTP_OUT, "afterSummary.template"), []byte(b), 0644)
	result = make([]string, 0)
	result = append(result, reportBeforeContentTemplates...)
	result = append(result, reportAfterContentTemplates...)
	result = append(result, reportBeforeQualityTemplates...)
	result = append(result, reportAfterQualityTemplates...)
	b = strings.Join(result, "\n")
	_ = ioutil.WriteFile(fmt.Sprintf("%s/%s", types.FASTP_OUT, "reportTemplate.template"), []byte(b), 0644)
	return nil
}
