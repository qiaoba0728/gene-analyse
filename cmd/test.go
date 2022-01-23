package cmd

import (
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/spf13/cobra"
)

var testCmd = &cobra.Command{
	Use:   "test",
	Short: "test gene-analyse plugin test",
	Long:  "run gene-analyse plugin build test data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
		//conf.InitBuildConfig("I://code//go//src//gene-analyse//build.json")
	},
	Run: func(cmd *cobra.Command, args []string) {
		//conf.GetBuildConfig()
		//l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		//random := build.NewRandomPlugin(l)
		//random.Build(context.Background())
		err := common.BuildReport()
		if err != nil {
			panic(err)
		}
	},
}
