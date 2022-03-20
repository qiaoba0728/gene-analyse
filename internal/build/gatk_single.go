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

type gatkSingleResultPlugin struct {
	bar *pgbar.Pgbar
	tp  types.SampleType
	//samples []string
	logger *zap.Logger
}

func NewGATKSingleResultPlugin(logger *zap.Logger) types.Plugin {
	return &gatkSingleResultPlugin{
		logger: logger,
		bar:    pgbar.New("gatk_single"),
	}
}
func (g *gatkSingleResultPlugin) check() {
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
	if b := utils.IsExist(types.GATK_G_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.GATK_G_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create gatk dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.GATK_G_OUT)
		if err != nil {
			g.logger.Error("existed gatk dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("gatk out vcf", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.GATK_SINGLE_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.GATK_SINGLE_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create gatk single sorted dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.GATK_SINGLE_OUT)
		if err != nil {
			g.logger.Error("existed gatk single sorted dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("gatk single out", zap.Strings("files", names))
		}
	}
}
func (g *gatkSingleResultPlugin) getGTF() (string, error) {
	files, err := ioutil.ReadDir(types.REFERENCES)
	if err != nil {
		return "", err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".gtf") {
			return path.Join(types.REFERENCES, v.Name()), nil
		}
	}
	return "", errors.New("gtf not find")
}
func (g *gatkSingleResultPlugin) buildBed12() error {
	var (
		gtf string
		err error
	)
	if utils.IsExist(fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12")) {
		g.logger.Info("bed file has existed")
		return nil
	}
	if gtf, err = g.getGTF(); err != nil {
		return err
	}
	f, err := os.Create(fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"))
	if err != nil {
		g.logger.Error("read_distribution create fail", zap.Error(err))
		return err
	}
	defer f.Close()
	cmd := exec.Command("/work/gtf2bed12.perl", gtf)
	cmd.Stdout = f
	cmd.Stderr = os.Stderr
	g.logger.Info("run gtf to bed", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		g.logger.Error("run build bed12", zap.Error(err), zap.String("cmd", cmd.String()))
	}
	return err
}
func (g *gatkSingleResultPlugin) Name() string {
	return "gatkResultPlugin"
}
func (g *gatkSingleResultPlugin) Build(ctx context.Context) error {
	g.check()
	go func() {
		err := g.buildBed12()
		if err != nil {
			g.logger.Error("build bed fail", zap.Error(err))
		}
	}()
	if err := g.index(); err != nil {
		return err
	}
	if err := g.merge(); err != nil {
		return err
	}
	if err := g.result(); err != nil {
		return err
	}
	//if err := g.geneDepthCoverage(types.SORTED_OUT); err != nil {
	//	return err
	//}
	g.logger.Info("gatk build finished")
	return nil
}
func (g *gatkSingleResultPlugin) merge() error {
	files, err := ioutil.ReadDir(types.GATK_G_OUT)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".g.vcf") && v.Name() != "merge.g.vcf" {
			start := time.Now()
			temp := strings.TrimSuffix(v.Name(), ".g.vcf")
			cmd := exec.Command("gatk", "GenotypeGVCFs", "-R", "gene.fa",
				"-V", fmt.Sprintf("%s/%s.g.vcf", types.GATK_G_OUT, temp),
				"-O", fmt.Sprintf("%s/%s.vcf", types.GATK_SINGLE_OUT, temp))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
			if err = cmd.Run(); err != nil {
				g.logger.Error("run gatk GenotypeGVCFs bam", zap.Error(err), zap.String("cmd", cmd.String()))
				return err
			}

			cmd = exec.Command("bgzip", "-f", fmt.Sprintf("%s/%s.vcf", types.GATK_SINGLE_OUT, temp))
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
			cmd = exec.Command("tabix", "-p", "vcf", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_SINGLE_OUT, temp))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
			if err = cmd.Run(); err != nil {
				g.logger.Error("run tabix bam", zap.Error(err), zap.String("cmd", cmd.String()))
				return err
			}
			g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
		}
	}
	return nil
}
func (g *gatkSingleResultPlugin) getfa() (string, error) {
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
func (g *gatkSingleResultPlugin) index() error {
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
func (g *gatkSingleResultPlugin) result() error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
		//thread int
	)
	pool, err = ants.NewPool(8)
	if err != nil {
		return err
	}
	//bar := e.bar.NewBar("featurecounts",len(files))
	files, err = ioutil.ReadDir(types.GATK_SINGLE_OUT)
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
						"SNP", "-V", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_SINGLE_OUT, temp), "-O", fmt.Sprintf("%s/%s.snp.vcf.gz", types.GATK_SINGLE_OUT, temp))
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
							types.GATK_SINGLE_OUT, temp,
							types.GATK_SINGLE_OUT, temp))
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
						"-select-type", "INDEL", "-V", fmt.Sprintf("%s/%s.vcf.gz", types.GATK_SINGLE_OUT, temp),
						"-O", fmt.Sprintf("%s/%s.indel.vcf.gz", types.GATK_SINGLE_OUT, temp))
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
							types.GATK_SINGLE_OUT, temp,
							types.GATK_SINGLE_OUT, temp))
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
					"-I", fmt.Sprintf("%s/%s.snp.filter.vcf.gz", types.GATK_SINGLE_OUT, temp),
					"-I", fmt.Sprintf("%s/%s.indel.filter.vcf.gz", types.GATK_SINGLE_OUT, temp),
					"-O", fmt.Sprintf("%s/%s.filter.vcf.gz", types.GATK_SINGLE_OUT, temp))
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
func (g *gatkResultPlugin) singleReport() error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
	)
	pool, err = ants.NewPool(8)
	if err != nil {
		return err
	}
	//bar := e.bar.NewBar("featurecounts",len(files))
	files, err = ioutil.ReadDir(types.SINGLE_OUT)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".indel.filter.vcf.gz") {
			wg.Add(1)
			source := v.Name()
			temp := strings.TrimSuffix(v.Name(), ".indel.filter.vcf.gz")
			if err := pool.Submit(func() {
				defer func() {
					wg.Done()
					g.logger.Info("build report success")
				}()
				cmd := exec.Command("bin/bash", "-c", fmt.Sprintf("cp %s/%s .", types.SINGLE_OUT, source))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
				cmd = exec.Command("bin/bash", "-c", fmt.Sprintf("gunzip %s/%s", types.SINGLE_OUT, source))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
				cmd = exec.Command("python", "gatk_indel.py", fmt.Sprintf("%s/%s", types.SINGLE_OUT, source), fmt.Sprintf("%s/%s_indel.report", types.REPORT_OUT, temp))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
				}
				cmd = exec.Command("bin/bash", "-c", "rm", "-rf", source)
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
			}); err != nil {
				g.logger.Error("run submmit fail", zap.Error(err))
			}
		}
		if strings.HasSuffix(v.Name(), ".snp.filter.vcf.gz") {
			wg.Add(1)
			source := v.Name()
			temp := strings.TrimSuffix(v.Name(), ".snp.filter.vcf.gz")
			if err := pool.Submit(func() {
				defer func() {
					wg.Done()
					g.logger.Info("build report success")
				}()
				cmd := exec.Command("bin/bash", "-c", fmt.Sprintf("cp %s/%s .", types.SINGLE_OUT, source))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
				cmd = exec.Command("bin/bash", "-c", fmt.Sprintf("gunzip %s/%s", types.SINGLE_OUT, source))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
				cmd = exec.Command("python", "gatk_snp.py", fmt.Sprintf("%s/%s", types.SINGLE_OUT, source), fmt.Sprintf("%s/%s_snp.report", types.REPORT_OUT, temp))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
				}
				cmd = exec.Command("bin/bash", "-c", "rm", "-rf", source)
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("run cmd", zap.String("cmd", cmd.String()))
				if err := cmd.Run(); err != nil {
					g.logger.Error("run cp", zap.Error(err), zap.String("cmd", cmd.String()))
					return
				}
			}); err != nil {
				g.logger.Error("run submmit fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
