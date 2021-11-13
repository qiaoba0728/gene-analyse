package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qianlnk/pgbar"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
	"strings"
	"sync"
	"time"
)

type gatkResultPlugin struct {
	bar *pgbar.Pgbar
	tp  types.SampleType
	//samples []string
	logger *zap.Logger
}

func NewGATKResultPlugin(logger *zap.Logger) types.Plugin {
	return &gatkResultPlugin{
		logger: logger,
		bar:    pgbar.New("gatk"),
	}
}
func (g *gatkResultPlugin) check() {
	if b := utils.IsExist(types.LOG); !b {
		cmd := exec.Command("mkdir", "-p", types.LOG)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create log dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.LOG)
		if err != nil {
			g.logger.Error("existed log dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("log", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.FASTP_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.FASTP_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create clean dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.FASTP_OUT)
		if err != nil {
			g.logger.Error("existed clean dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("fastp out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.REPORT_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.REPORT_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create report dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.REPORT_OUT)
		if err != nil {
			g.logger.Error("existed report dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("report out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.HISAT2_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.HISAT2_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create hisat2 dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.HISAT2_OUT)
		if err != nil {
			g.logger.Error("existed hisat2 dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("hisat2 out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.SORTED_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.SORTED_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create hisat2 sorted dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.SORTED_OUT)
		if err != nil {
			g.logger.Error("existed hisat2 sorted dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("hisat2 sorted out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.GATK_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.GATK_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create gatk dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.GATK_OUT)
		if err != nil {
			g.logger.Error("existed gatk dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("gatk out", zap.Strings("files", names))
		}
	}
}
func (g *gatkResultPlugin) Name() string {
	return "gatkResultPlugin"
}
func (g *gatkResultPlugin) Build(ctx context.Context) error {
	if err := g.index(); err != nil {
		return err
	}
	if err := g.merge(); err != nil {
		return err
	}
	if err := g.result(); err != nil {
		return err
	}
	g.logger.Info("gatk build finished")
	return nil
}
func (g *gatkResultPlugin) merge() error {
	vcfs := ""
	files, err := ioutil.ReadDir(types.GATK_OUT)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".g.vcf") {
			//vcfs = append(vcfs,fmt.Sprintf("%s/%s",types.GATK_OUT,v.Name()))
			vcfs = vcfs + fmt.Sprintf(" -V %s/%s ", types.GATK_OUT, v.Name())
		}
	}
	if vcfs != "" {
		temp := "merge"
		start := time.Now()
		vcfs = "gatk --java-options '-Xmx30G' CombineGVCFs -R gene.fa " + vcfs + " -O " + fmt.Sprintf("%s/%s.g.vcf", types.GATK_OUT, temp)
		cmd := exec.Command("/bin/sh", "-c", vcfs)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
		if err = cmd.Run(); err != nil {
			g.logger.Error("run gatk CombineGVCFs bam", zap.Error(err), zap.String("cmd", cmd.String()))
			return err
		}
		start = time.Now()
		cmd = exec.Command("gatk", "GenotypeGVCFs", "-R", "gene.fa",
			"-V", fmt.Sprintf("%s/%s.g.vcf", types.GATK_OUT, temp),
			"-O", fmt.Sprintf("%s/%s.vcf", types.GATK_OUT, temp))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
		if err = cmd.Run(); err != nil {
			g.logger.Error("run gatk GenotypeGVCFs bam", zap.Error(err), zap.String("cmd", cmd.String()))
			return err
		}

		cmd = exec.Command("bgzip", "-f", fmt.Sprintf("%s/%s.vcf", types.GATK_OUT, temp))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
		if err = cmd.Run(); err != nil {
			g.logger.Error("run bgzip bam", zap.Error(err), zap.String("cmd", cmd.String()))
			return err
		}
		g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
		//tabix -p vcf ${file}.vcf.gz
		start = time.Now()
		cmd = exec.Command("tabix", "-p", "vcf", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_OUT, temp))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
		if err = cmd.Run(); err != nil {
			g.logger.Error("run tabix bam", zap.Error(err), zap.String("cmd", cmd.String()))
			return err
		}
		g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
	}
	return nil
}
func (g *gatkResultPlugin) getfa() (string, error) {
	files, err := ioutil.ReadDir(types.REFERENCES)
	if err != nil {
		return "", err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".fasta") || strings.HasSuffix(v.Name(), ".fa") {
			return path.Join(types.REFERENCES, v.Name()), nil
		}
	}
	return "", errors.New("fa not find")
}
func (g *gatkResultPlugin) index() error {
	fa, err := g.getfa()
	if err != nil {
		return err
	}
	cmd := exec.Command("cp", "-r", fa, "./gene.fa")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	g.logger.Info("index run ", zap.String("cmd", cmd.String()))
	if err := cmd.Run(); err != nil {
		g.logger.Error("create index dir", zap.Error(err))
		return err
	}
	cmd = exec.Command("bwa", "index", "gene.fa")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	g.logger.Info("index run ", zap.String("cmd", cmd.String()))
	if err := cmd.Run(); err != nil {
		g.logger.Error("create index dir", zap.Error(err))
		return err
	}
	cmd = exec.Command("samtools", "faidx", "gene.fa")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	g.logger.Info("index run ", zap.String("cmd", cmd.String()))
	if err := cmd.Run(); err != nil {
		g.logger.Error("create index dir", zap.Error(err))
		return err
	}
	cmd = exec.Command("gatk", "CreateSequenceDictionary", "-R",
		"gene.fa", "-O", "gene.dict")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	g.logger.Info("index run ", zap.String("cmd", cmd.String()))
	if err := cmd.Run(); err != nil {
		g.logger.Error("create index dir", zap.Error(err))
		return err
	}
	return nil
}
func (g *gatkResultPlugin) result() error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
		//thread int
	)
	//fa,err := g.getfa()
	//if err != nil {
	//	return err
	//}
	//wd, _ = os.Getwd()
	//env := os.Getenv("THREAD")
	//if env != "" {
	//	n, _ := strconv.Atoi(env)
	//	thread = n
	//} else {
	//	thread = 4
	//}
	pool, err = ants.NewPool(8)
	if err != nil {
		return err
	}
	//bar := e.bar.NewBar("featurecounts",len(files))
	files, err = ioutil.ReadDir(types.GATK_OUT)
	if err != nil {
		return err
	}
	bar := g.bar.NewBar("vcf -> result", len(files))
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".vcf.gz") {
			wg.Add(1)
			name := v.Name()
			//input := path.Join(types.GATK_OUT, name)
			if err = pool.Submit(func() {
				var wgVCF sync.WaitGroup
				wgVCF.Add(2)
				temp := strings.TrimSuffix(name, ".vcf.gz")
				defer func() {
					bar.Add(1)
					wg.Done()
					g.logger.Info("build gatk file success", zap.String("names", name))
				}()
				go func() {
					defer wgVCF.Done()
					start := time.Now()
					//gatk SelectVariants -select-type SNP -V ${file}.vcf.gz -O ${file}.snp.vcf.gz
					cmd := exec.Command("gatk", "SelectVariants", "-select-type",
						"SNP", "-V", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_OUT, temp), "-O", fmt.Sprintf("%s/%s.snp.vcf.gz", types.GATK_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						g.logger.Error("run gatk SelectVariants bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
					start = time.Now()
					cmd = exec.Command("/bin/sh", "-c",
						fmt.Sprintf(`gatk VariantFiltration -V %s/%s.snp.vcf.gz --filter-expression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name PASS -O %s/%s.snp.filter.vcf.gz`,
							types.GATK_OUT, temp,
							types.GATK_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						g.logger.Error("run gatk VariantFiltration bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
				}()
				go func() {
					defer wgVCF.Done()
					//gatk SelectVariants -select-type INDEL -V ${file}.vcf.gz -O ${file}.indel.vcf.gz
					start := time.Now()
					cmd := exec.Command("gatk", "SelectVariants",
						//"--java-options", `"-Xmx15G -Djava.io.tmpdir=./"`,
						"-select-type", "INDEL", "-V", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_OUT, temp),
						"-O", fmt.Sprintf("%s/%s.indel.vcf.gz", types.GATK_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						g.logger.Error("run gatk SelectVariants bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
					start = time.Now()
					cmd = exec.Command("/bin/sh", "-c",
						fmt.Sprintf(`gatk VariantFiltration -V %s/%s.indel.vcf.gz --filter-expression 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name PASS -O %s/%s.indel.filter.vcf.gz`,
							types.GATK_OUT, temp,
							types.GATK_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						g.logger.Error("run gatk VariantFiltration bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
				}()
				g.logger.Error("gatk waiting snp and indel finished")
				wgVCF.Wait()
				//gatk MergeVcfs -I ${file}.snp.filter.vcf.gz -I ${file}.indel.filter.vcf.gz -O ${file}.filter.vcf.gz
				start := time.Now()
				cmd := exec.Command("gatk", "MergeVcfs",
					//"--java-options", `"-Xmx15G -Djava.io.tmpdir=./"`,
					"-I", fmt.Sprintf("%s/%s.snp.filter.vcf.gz", types.GATK_OUT, temp),
					"-I", fmt.Sprintf("%s/%s.indel.filter.vcf.gz", types.GATK_OUT, temp),
					"-O", fmt.Sprintf("%s/%s.filter.vcf.gz", types.GATK_OUT, temp))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					g.logger.Error("run gatk MergeVcfs bam", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
				g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
			}); err != nil {
				g.logger.Error("pool run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
