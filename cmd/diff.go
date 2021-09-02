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
	cfgFile string
)

func init() {
	diffCmd.Flags().StringVar(&cfgFile, "conf", "/data/config.json", "config file (default is ./config.yml)")
}

var diffCmd = &cobra.Command{
	Use:   "diff",
	Short: "gene-analyse plugin",
	Long:  "run gene-analyse plugin build picture",
	PreRun: func(cmd *cobra.Command, args []string) {
		conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		plugin := build.NewGenePlugin(conf.GetConfig(), l.Named("diff"))
		err := plugin.Build(context.Background())
		if err != nil {
			l.Error("build fail", zap.Error(err))
		}
	},
}
