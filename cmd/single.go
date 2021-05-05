package cmd

import (
	"context"
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
	"io/ioutil"
)

var (
	singleConfigFile string
)

func init() {
	singleCmd.Flags().StringVar(&singleConfigFile, "conf", "/work/build.json", "config file (default is ./config.yml)")
}

var singleCmd = &cobra.Command{
	Use:   "single",
	Short: "gene-analyse plugin",
	Long:  "run gene-analyse single plugin build",
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
		plugin := build.NewSinglePlugin(l.Named("single"))
		if err = plugin.Build(context.Background()); err != nil {
			l.Error("data single", zap.Error(err))
			return
		} else {
			l.Info("data build single success")
		}
		email := common.NewQQEmail(l.Named("email"))
		if err != nil {
			err = email.SendEmail("494340090@qq.com", conf.GetBuildConfig().Sends, "ngybgeetypogbibc", "smtp.qq.com", "hx", conf.GetBuildConfig().Subject, err.Error(), 25, false)
		} else {
			inputFiles := make([]string, 0)
			files, _ := ioutil.ReadDir(types.SINGLE_OUT)
			for _, v := range files {
				temp := fmt.Sprintf("%s/%s", types.SINGLE_OUT, v.Name())
				inputFiles = append(inputFiles, temp)
			}
			err := utils.ZipFiles("result.zip", inputFiles, "/data/output", ".")
			if err != nil {
				l.Error("zip error", zap.Error(err))
			}
			err = email.SendEmailWithPoolAndFile(conf.GetBuildConfig().Sends, "494340090@qq.com", "ngybgeetypogbibc", "smtp.qq.com", "hx", conf.GetBuildConfig().Subject, conf.GetBuildConfig().Body, 25, "/bsa/result.zip")
		}
		if err != nil {
			l.Error("send email", zap.Error(err))
		} else {
			l.Info("send email success")
		}
	},
}
