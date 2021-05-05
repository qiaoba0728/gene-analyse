package cmd

import (
	"context"
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
	"io/ioutil"
	"strings"
)

var (
	bsaConfigFile string
)

func init() {
	bsaCmd.Flags().StringVar(&bsaConfigFile, "conf", "/bsa/bsa.json", "config file (default is ./config.yml)")
}

var bsaCmd = &cobra.Command{
	Use:   "bsa",
	Short: "gene-analyse plugin bsa",
	Long:  "run gene-analyse plugin build bsa data",
	PreRun: func(cmd *cobra.Command, args []string) {
		if utils.IsExist("/data/bsa.json") {
			conf.InitBSAConfig("/data/bsa.json")
		} else {
			conf.InitBSAConfig(bsaConfigFile)
		}
	},
	Run: func(cmd *cobra.Command, args []string) {
		var err error
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		if len(args) == 1 && args[0] == "rna" {
			plugin := build.NewBsaRNAPlugin(l.Named("RNA"))
			if err = plugin.Build(context.Background()); err != nil {
				l.Error("data bsa rna build ", zap.Error(err))
			} else {
				l.Info("wait bsa rna finished......")
			}
		} else {
			plugin := build.NewBsaPlugin(l.Named("DNA"))
			if err = plugin.Build(context.Background()); err != nil {
				l.Error("data bsa build ", zap.Error(err))
			} else {
				l.Info("wait bsa dna finished......")
			}
		}
		email := common.NewQQEmail(l.Named("email"))
		if err != nil {
			err = email.SendEmail("494340090@qq.com", conf.GetBSAConfig().Sends, "ngybgeetypogbibc", "smtp.qq.com", "hx", conf.GetBSAConfig().Subject, err.Error(), 25, false)
		} else {
			inputFiles := make([]string, 0)
			files, _ := ioutil.ReadDir("/data/output/")
			for _, v := range files {
				if strings.HasSuffix(v.Name(), ".pdf") {
					temp := fmt.Sprintf("/data/output/%s", v.Name())
					inputFiles = append(inputFiles, temp)
				}
			}
			err := utils.ZipFiles("result.zip", inputFiles, "/data/output", ".")
			if err != nil {
				l.Error("zip error", zap.Error(err))
			}
			//err = email.SendEmail("494340090@qq.com",conf.GetBSAConfig().Sends, "ngybgeetypogbibc", "smtp.qq.com","hx", conf.GetBSAConfig().Subject, conf.GetBSAConfig().Body, 25,false)
			err = email.SendEmailWithPoolAndFile(conf.GetBSAConfig().Sends, "494340090@qq.com", "ngybgeetypogbibc", "smtp.qq.com", "hx", conf.GetBSAConfig().Subject, conf.GetBSAConfig().Body, 25, "/bsa/result.zip")
		}
		if err != nil {
			l.Error("send email", zap.Error(err))
		} else {
			l.Info("send email success")
		}
	},
}
