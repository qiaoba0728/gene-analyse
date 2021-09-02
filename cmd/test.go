package cmd

import (
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
	"io/ioutil"
	"os"
	"os/exec"
	"strings"
	"time"
)

var testCmd = &cobra.Command{
	Use:   "test",
	Short: "test gene-analyse plugin",
	Long:  "run gene-analyse plugin build test data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
		//g, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		//plugin := build.NewReportPlugin(l)
		//if err := plugin.Build(context.Background()); err != nil {
		//	l.Error("report  fail", zap.Error(err))
		//} else {
		//	l.Info("report finished")
		//}
		//command := exec.Command("gatk", "VariantFiltration",
		//	//"--java-options", `"-Xmx15G -Djava.io.tmpdir=./"`,
		//	"-V", fmt.Sprintf("%s/%s.indel.vcf.gz", types.GATK_OUT, "W"),
		//	"--filter-expression", "'QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'",
		//	"--filter-name", "PASS",
		//	"-O", fmt.Sprintf("%s/%s.indel.filter.vcf.gz", "/work", "10Z"))
		//command := exec.Command("/bin/sh", "-c",
		//	fmt.Sprintf(`gatk VariantFiltration -V %s/%s.indel.vcf.gz --filter-expression 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name PASS -O %s/%s.indel.filter.vcf.gz`,
		//		types.GATK_OUT, "W",
		//		"/work", "10Z"))
		//command.Stdout = os.Stdout
		//command.Stderr = os.Stderr
		////g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
		//if err := command.Run(); err != nil {
		//	//g.logger.Error("run gatk VariantFiltration bam", zap.Error(err), zap.String("cmd", cmd.String()))
		//	fmt.Println("cmd", command.String())
		//	return
		//}
		//vcfs := ""
		//files, err := ioutil.ReadDir(types.GATK_OUT)
		//if err != nil {
		//	return
		//}
		//for _, v := range files {
		//	if strings.HasSuffix(v.Name(), ".g.vcf") {
		//		//vcfs = append(vcfs,fmt.Sprintf("%s/%s",types.GATK_OUT,v.Name()))
		//		vcfs = vcfs + fmt.Sprintf("-V %s/%s ", types.GATK_OUT, v.Name())
		//	}
		//}
		//if vcfs != "" {
		//	start := time.Now()
		//	temp := "merge"
		//	start = time.Now()
		//	vcfs = "gatk CombineGVCFs -R gene.fa " + vcfs + " -O " + fmt.Sprintf("%s.g.vcf", temp)
		//	cmd := exec.Command("/bin/sh", "-c", vcfs)
		//	cmd.Stdout = os.Stdout
		//	cmd.Stderr = os.Stderr
		//	if err = cmd.Run(); err != nil {
		//		g.Error("run gatk CombineGVCFs bam", zap.Error(err), zap.String("cmd", cmd.String()))
		//		return
		//	}
		//	g.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
		//	start = time.Now()
		//	cmd = exec.Command("gatk", "GenotypeGVCFs", "-R", "gene.fa",
		//		"-V", fmt.Sprintf("%s/%s.g.vcf", types.GATK_OUT, temp),
		//		"-O", fmt.Sprintf("%s/%s.vcf", types.GATK_OUT, temp))
		//	cmd.Stdout = os.Stdout
		//	cmd.Stderr = os.Stderr
		//	if err = cmd.Run(); err != nil {
		//		g.Error("run gatk GenotypeGVCFs bam", zap.Error(err), zap.String("cmd", cmd.String()))
		//		return
		//	}
		//
		//	cmd = exec.Command("bgzip", "-f", fmt.Sprintf("%s/%s.vcf", types.GATK_OUT, temp))
		//	cmd.Stdout = os.Stdout
		//	cmd.Stderr = os.Stderr
		//	if err = cmd.Run(); err != nil {
		//		g.Error("run bgzip bam", zap.Error(err), zap.String("cmd", cmd.String()))
		//		return
		//	}
		//	g.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
		//	//tabix -p vcf ${file}.vcf.gz
		//	start = time.Now()
		//	cmd = exec.Command("tabix", "-p", "vcf", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_OUT, temp))
		//	cmd.Stdout = os.Stdout
		//	cmd.Stderr = os.Stderr
		//	if err = cmd.Run(); err != nil {
		//		g.Error("run tabix bam", zap.Error(err), zap.String("cmd", cmd.String()))
		//		return
		//	}
		//	g.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
		//}
	},
}
