package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
	sync2 "sync"
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
		//if utils.IsExist(types.KEGG_OUT) {
		//	cmd := exec.Command("rm","-f",types.KEGG_OUT)
		//	cmd.Stdout = os.Stdout
		//	cmd.Stderr = os.Stderr
		//	if err :=  cmd.Run();err != nil {
		//		l.Error("delete kegg out fail",zap.Error(err))
		//	}
		//}
		//exeCmd := exec.Command("mkdir","-p",types.KEGG_OUT)
		//exeCmd.Stdout = os.Stdout
		//exeCmd.Stderr = os.Stderr
		//if err :=  exeCmd.Run();err != nil {
		//	l.Error("create kegg out fail",zap.Error(err))
		//}
		//if utils.IsExist(types.KEGG_OUT) {
		//	cmd := exec.Command("rm","-f",types.GO_OUT)
		//	cmd.Stdout = os.Stdout
		//	cmd.Stderr = os.Stderr
		//	if err :=  cmd.Run();err != nil {
		//		l.Error("delete go out fail",zap.Error(err))
		//	}
		//}
		//exeCmd = exec.Command("mkdir","-p",types.GO_OUT)
		//exeCmd.Stdout = os.Stdout
		//exeCmd.Stderr = os.Stderr
		//if err :=  exeCmd.Run();err != nil {
		//	l.Error("create go out fail",zap.Error(err))
		//}
		var wg sync2.WaitGroup
		for _, v := range conf.GetConfig().Group {
			plugin := build.NewGenePlugin(&conf.Group{
				Name:          v.Name,
				Start:         v.Start,
				End:           v.End,
				StartRepeated: v.StartRepeated,
				EndRepeated:   v.EndRepeated,
				Output:        v.Output,
			}, l.Named("diff"))
			wg.Add(1)
			go func() {
				defer wg.Done()
				err := plugin.Build(context.Background())
				if err != nil {
					l.Error("build fail", zap.Error(err))
				}
			}()
		}
		l.Info("wait diff finished")
		wg.Wait()
	},
}
