package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
)

var (
	gatkConfigFile string
)

func init() {
	gatkCmd.Flags().StringVar(&gatkConfigFile, "conf", "/work/build.json", "config file (default is ./config.yml)")
}

var gatkCmd = &cobra.Command{
	Use:   "gatk",
	Short: "gene-analyse plugin",
	Long:  "run gene-analyse plugin build upstream data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//if utils.IsExist("/data/build.json") {
		//	conf.InitBuildConfig("/data/build.json")
		//} else {
		//	conf.InitBuildConfig(gatkConfigFile)
		//}
	},
	Run: func(cmd *cobra.Command, args []string) {
		var err error
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		if len(args) == 1 && args[0] == "result" {
			plugin := build.NewGATKResultPlugin(l.Named("result"))
			if err = plugin.Build(context.Background()); err != nil {
				l.Error("gatk result data build ", zap.Error(err))
				return
			}
		} else if len(args) == 1 && args[0] == "sample" {
			plugin := build.NewGATKSingleResultPlugin(l.Named("result_sample"))
			if err = plugin.Build(context.Background()); err != nil {
				l.Error("gatk bsa result data build ", zap.Error(err))
				return
			}
		} else {
			plugin := build.NewGATKPlugin(l.Named("build"))
			if err = plugin.Build(context.Background()); err != nil {
				l.Error("gatk data build ", zap.Error(err))
				return
			}
		}
	},
}
