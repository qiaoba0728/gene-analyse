package cmd

import (
	"context"
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
	"strings"
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
		} else if len(args) == 1 && args[0] == "read" {
			err := utils.WriteFile("read.R", scripts.READ_REPORT)
			if err != nil {
				l.Named("read").Error("cmd run fail", zap.Error(err))
				return
			}
			files, err := ioutil.ReadDir(types.REPORT_OUT)
			if err != nil {
				l.Named("read").Error("cmd run fail", zap.Error(err))
				return
			}
			for _, v := range files {
				if strings.HasSuffix(v.Name(), "_read_distribution.log") {
					temp := strings.TrimSuffix(v.Name(), "_read_distribution.log")
					cmd := exec.Command("/bin/sh", "-c", fmt.Sprintf("cat %s/%s_read_distribution.log | tail -n +5 | head -n 11 > %s/%s_read_distribution_plot.txt", types.REPORT_OUT, temp, types.REPORT_OUT, temp))
					l.Named("read").Info("run read_quality", zap.String("cmd", cmd.String()))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						l.Named("read").Error("read_quality bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					wd, _ := os.Getwd()
					cmd = exec.Command("Rscript", path.Join(wd, "script", "read.R"), fmt.Sprintf("%s/%s_read_distribution_plot.txt", types.REPORT_OUT, temp), fmt.Sprintf("%s/%s_read_distribution_plot.png", types.REPORT_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					l.Named("build_stat").Info("cmd run ", zap.String("cmd", cmd.String()))
					if err := cmd.Run(); err != nil {
						l.Named("read").Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
				}
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
