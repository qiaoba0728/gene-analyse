package cmd

import (
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/spf13/cobra"
	"os"
	"os/exec"
)

var testCmd = &cobra.Command{
	Use:   "test",
	Short: "test gene-analyse plugin",
	Long:  "run gene-analyse plugin build test data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
		//l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		//plugin := build.NewReportPlugin(l)
		//if err := plugin.Build(context.Background()); err != nil {
		//	l.Error("report  fail", zap.Error(err))
		//} else {
		//	l.Info("report finished")
		//}
		command := exec.Command("gatk", "VariantFiltration",
			//"--java-options", `"-Xmx15G -Djava.io.tmpdir=./"`,
			"-V", fmt.Sprintf("%s/%s.indel.vcf.gz", types.GATK_OUT, "10Z"),
			"--filter-expression", `'QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'`,
			"--filter-name", "PASS",
			"-O", fmt.Sprintf("%s/%s.indel.filter.vcf.gz", types.GATK_OUT, "10Z"))
		command.Stdout = os.Stdout
		command.Stderr = os.Stderr
		//g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
		if err := command.Run(); err != nil {
			//g.logger.Error("run gatk VariantFiltration bam", zap.Error(err), zap.String("cmd", cmd.String()))
			fmt.Println("cmd", command.String())
			return
		}
	},
}
