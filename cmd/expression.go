package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
)

//获取表达矩阵相关的数据
var expressionCmd = &cobra.Command{
	Use:   "expression",
	Short: "gene-analyse plugin",
	Long:  "run gene-analyse plugin build tpm and fkpm and counts",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		plugin := build.NewExpressionPlugin(l)
		if err := plugin.Build(context.Background()); err != nil {
			l.Error("expression  fail", zap.Error(err))
		} else {
			l.Info("expression finished")
		}
	},
}
