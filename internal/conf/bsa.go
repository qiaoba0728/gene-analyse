package conf

import (
	"encoding/json"
	"io/ioutil"
	"sync"
)

type BSAConfig struct {
	Fastp    int      `json:"fastp"`
	Hisat2   int      `json:"hisat2"`
	Samtools int      `json:"samtools"`
	Subject  string   `json:"subject"`
	Body     string   `json:"body"`
	Sends    []string `json:"sends"`
}

var (
	bsaConfigPath string
	bsaConfig     BuildConfig
	bsaOnce       sync.Once
)

func InitBSAConfig(path string) {
	bsaConfigPath = path
}
func GetBSAConfig() *BuildConfig {
	bsaOnce.Do(func() {
		jsonFile, err := ioutil.ReadFile(bsaConfigPath)
		if err != nil {
			panic(err)
		}
		err = json.Unmarshal(jsonFile, &bsaConfig)
		if err != nil {
			panic(err)
		}
	})
	return &bsaConfig
}
