package conf

import (
	"encoding/json"
	"io/ioutil"
	"sync"
)

type BuildConfig struct {
	Fastp    int      `json:"fastp"`
	Hisat2   int      `json:"hisat2"`
	Samtools int      `json:"samtools"`
	Subject  string   `json:"subject"`
	Body     string   `json:"body"`
	Sends    []string `json:"sends"`
}

var (
	buildConfigPath string
	buildConfig     BuildConfig
	buildOnce       sync.Once
)

func InitBuildConfig(path string) {
	buildConfigPath = path
}
func GetBuildConfig() *BuildConfig {
	buildOnce.Do(func() {
		jsonFile, err := ioutil.ReadFile(buildConfigPath)
		if err != nil {
			panic(err)
		}
		err = json.Unmarshal(jsonFile, &buildConfig)
		if err != nil {
			panic(err)
		}
	})
	return &buildConfig
}
