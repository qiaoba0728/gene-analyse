package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
)

var randomCmd = &cobra.Command{
	Use:   "random",
	Short: "gene-analyse plugin",
	Long:  "run gene-analyse plugin build random data",
	PreRun: func(cmd *cobra.Command, args []string) {

	},
	Run: func(cmd *cobra.Command, args []string) {
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		random := build.NewRandomPlugin(l)
		err := random.Build(context.Background())
		if err != nil {
			l.Error("data build ", zap.Error(err))
		}
	},
}