package conf

import (
	"github.com/stretchr/testify/assert"
	"testing"
)

func TestGetSummary(t *testing.T) {
	summary, err := GetSummary("I://code//go//mod//gene-analyse//fastp.json")
	assert.Nil(t, err)
	t.Log(summary.Summary.AfterFiltering)
	t.Log(summary.Summary.BeforeFiltering)
	t.Log(summary.Read1After.ContentCurves)
	t.Log(len(summary.Read1After.ContentCurves.C))
	t.Log(summary.Read1After.QualityCurves)
	t.Log(summary.Read1Before.ContentCurves)
	t.Log(summary.Read1Before.QualityCurves)
	t.Log(summary.Read2Before.ContentCurves)
	t.Log(summary.Read2Before.QualityCurves)
	t.Log(summary.Read2After.ContentCurves)
	t.Log(summary.Read2After.QualityCurves)
}
