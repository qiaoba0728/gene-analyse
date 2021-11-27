package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
	"os"
	"os/exec"
	"path"
)

var reportCmd = &cobra.Command{
	Use:   "report",
	Short: "test gene-analyse plugin",
	Long:  "run gene-analyse plugin build report data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		if len(args) == 1 && args[0] == "stat" {
			err := utils.WriteFile("build_stat.R", scripts.BUILD_STAT)
			if err != nil {
				l.Named("build_stat").Error("cmd run fail", zap.Error(err))
				return
			}
			cmd := exec.Command("/bin/sh", "-x", "/work/build_stat.sh")
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			l.Named("build_stat").Info("cmd run ", zap.String("cmd", cmd.String()))
			if err := cmd.Run(); err != nil {
				l.Named("build_stat").Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
				return
			}
			wd, _ := os.Getwd()
			cmd = exec.Command("Rscript", path.Join(wd, "script", "build_stat.R"))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			l.Named("build_stat").Info("cmd run ", zap.String("cmd", cmd.String()))
			if err := cmd.Run(); err != nil {
				l.Named("build_stat").Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
				return
			}
		} else {
			plugin := build.NewReportPlugin(l)
			if err := plugin.Build(context.Background()); err != nil {
				l.Error("report  fail", zap.Error(err))
			} else {
				l.Info("report finished")
			}
		}
	},
}
