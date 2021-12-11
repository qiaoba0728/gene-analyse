package conf

import (
	"encoding/json"
	"io/ioutil"
	"sync"
)

type SampleConfig struct {
}

var (
	sampleConfigPath string
	sampleConfig     BuildConfig
	sampleOnce       sync.Once
)

func InitSampleConfig(path string) {
	sampleConfigPath = path
}
func GetSampleConfig() *BuildConfig {
	sampleOnce.Do(func() {
		jsonFile, err := ioutil.ReadFile(sampleConfigPath)
		if err != nil {
			panic(err)
		}
		err = json.Unmarshal(jsonFile, &sampleConfig)
		if err != nil {
			panic(err)
		}
	})
	return &sampleConfig
}
