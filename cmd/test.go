package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
)

var testCmd = &cobra.Command{
	Use:   "test",
	Short: "test gene-analyse plugin",
	Long:  "run gene-analyse plugin build test data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		plugin := build.NewReportPlugin(l)
		if err := plugin.Build(context.Background()); err != nil {
			l.Error("report  fail", zap.Error(err))
		} else {
			l.Info("report finished")
		}
	},
}
