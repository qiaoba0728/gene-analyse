package conf

import (
	"encoding/json"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"io/ioutil"
)

func GetSummary(path string) (*types.SummaryDocs, error) {
	summary := &types.SummaryDocs{}
	buf, err := ioutil.ReadFile(path)
	if err != nil {
		return summary, err
	}
	err = json.Unmarshal(buf, summary)
	return summary, err
}
