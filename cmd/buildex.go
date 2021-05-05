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
	buildExConfigFile string
)

func init() {
	buildExCmd.Flags().StringVar(&buildExConfigFile, "conf", "/work/build.json", "config file (default is ./config.yml)")
}

var buildExCmd = &cobra.Command{
	Use:   "buildex",
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
		plugin := build.NewBuildExPlugin(l.Named("build"))
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
		//email := common.NewQQEmail(l.Named("email"))
		//if err != nil {
		//	err = email.SendEmail("494340090@qq.com",conf.GetBuildConfig().Sends, "ngybgeetypogbibc", "smtp.qq.com","hx", conf.GetBuildConfig().Subject, err.Error(), 25,false)
		//}else {
		//	inputFiles := make([]string,0)
		//	files,_ := ioutil.ReadDir("/data/output/clean")
		//	for _,v := range  files {
		//		if strings.HasSuffix(v.Name(),".csv") {
		//			temp := fmt.Sprintf("/data/output/clean/%s",v.Name())
		//			inputFiles = append(inputFiles,temp)
		//		}
		//	}
		//	files,_ = ioutil.ReadDir("/data/output/expression_result")
		//	for _,v := range  files {
		//		if strings.HasSuffix(v.Name(),".csv") {
		//			temp := fmt.Sprintf("/data/output/expression_result/%s",v.Name())
		//			inputFiles = append(inputFiles,temp)
		//		}
		//	}
		//	files,_ = ioutil.ReadDir("/data/output/report_result")
		//	for _,v := range  files {
		//		temp := fmt.Sprintf("/data/output/report_result/%s",v.Name())
		//		inputFiles = append(inputFiles,temp)
		//	}
		//	err := utils.ZipFiles("result.zip",inputFiles,"/data/output",".")
		//	if err != nil {
		//		l.Error("zip error",zap.Error(err))
		//	}
		//	//err = email.SendEmail("494340090@qq.com",conf.GetBSAConfig().Sends, "ngybgeetypogbibc", "smtp.qq.com","hx", conf.GetBSAConfig().Subject, conf.GetBSAConfig().Body, 25,false)
		//	err = email.SendEmailWithPoolAndFile(conf.GetBuildConfig().Sends,"494340090@qq.com", "ngybgeetypogbibc", "smtp.qq.com","hx", conf.GetBuildConfig().Subject, conf.GetBuildConfig().Body, 25,"/bsa/result.zip")
		//}
		//if err !=  nil {
		//	l.Error("send email",zap.Error(err))
		//}else {
		//	l.Info("send email success")
		//}
	},
}
