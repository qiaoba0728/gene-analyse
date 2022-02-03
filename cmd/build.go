package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
)

var (
	buildConfigFile string
)

func init() {
	mainCmd.Flags().StringVar(&buildConfigFile, "conf", "/work/build.json", "config file (default is ./config.yml)")
}

var mainCmd = &cobra.Command{
	Use:   "build",
	Short: "gene-analyse plugin",
	Long:  "run gene-analyse plugin build upstream data",
	PreRun: func(cmd *cobra.Command, args []string) {
		if utils.IsExist("/data/build.json") {
			conf.InitBuildConfig("/data/build.json")
		} else {
			conf.InitBuildConfig(buildConfigFile)
		}
	},
	Run: func(cmd *cobra.Command, args []string) {
		var err error
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		plugin := build.NewUpStreamPlugin(l.Named("build"))
		if err = plugin.Build(context.Background()); err != nil {
			l.Error("data build ", zap.Error(err))
			return
		} else {
			l.Info("data build success")
			plugin = build.NewExpressionPlugin(l.Named("expression"))
			if err = plugin.Build(context.Background()); err != nil {
				l.Error("data build expression", zap.Error(err))
				return
			} else {
				l.Info("data build expression success")
			}
		}
	},
}
