package conf

import (
	"encoding/json"
	"io/ioutil"
	"sync"
)

//config
type Group struct {
	Start         string `json:"start"`         //起始列
	End           string `json:"end"`           //结束列
	Name          string `json:"name"`          //生成图片名字
	StartRepeated string `json:"startRepeated"` //重复数目
	EndRepeated   string `json:"endRepeated"`   //重复数目
	Output        string `json:"output"`        //生成的文件目录
	Richer        bool   `json:"richer"`        //go kegg分析
}
type DiffGroup struct {
	Name   string   `json:"name"` //生成图片名字
	Groups []string `json:"groups"`
	Output string   `json:"output"` //生成的文件目录
	Richer bool     `json:"richer"` //go kegg分析
}
type Config struct {
	GeneDB string  `json:"geneDB"` //['Rsativus']
	Factor float64 // 0.05
	//File string `json:"file"`			//基因文件位置
	Group     []*Group   `json:"group"`
	DiffGroup *DiffGroup `json:"diffGroup"`
}

var (
	configPath string
	conf       Config
	once       sync.Once
)

func InitConfig(path string) {
	configPath = path
}
func GetConfig() *Config {
	once.Do(func() {
		jsonFile, err := ioutil.ReadFile(configPath)
		if err != nil {
			panic(err)
		}
		err = json.Unmarshal(jsonFile, &conf)
		if err != nil {
			panic(err)
		}
	})
	return &conf
}
