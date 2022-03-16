package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qianlnk/pgbar"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"
	"sync"
	"time"
)

const (
	faFile = "/work/gene.fa"
)

type gatkPlugin struct {
	bar *pgbar.Pgbar
	tp  types.SampleType
	//samples []string
	logger *zap.Logger
}

func NewGATKPlugin(logger *zap.Logger) types.Plugin {
	return &gatkPlugin{
		logger: logger,
		bar:    pgbar.New("gatk"),
	}
}
func (g *gatkPlugin) check() {
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
	if b := utils.IsExist(types.FASTP_QC_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.FASTP_QC_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create clean dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.FASTP_QC_OUT)
		if err != nil {
			g.logger.Error("existed clean dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("fastqc out", zap.Strings("files", names))
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
}
func (g *gatkPlugin) Name() string {
	return "gatkPlugin"
}
func (g *gatkPlugin) Build(ctx context.Context) error {
	g.check()
	done := utils.Wait(func() error {
		var wg sync.WaitGroup
		wg.Add(2)
		go func() {
			defer wg.Done()
			if err := g.index(); err != nil {
				g.logger.Error("gene index build fail", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if b := utils.IsExist(types.INPUT); b {
				if b := utils.IsExist(types.FASTP_OUT); !b || (b && utils.Files(types.FASTP_OUT) == 0) {
					if err := g.fastq(types.INPUT); err != nil {
						g.logger.Error("create clean data fail", zap.Error(err))
					} else {
						err = common.BuildReport()
						if err != nil {
							g.logger.Error("build report data fail", zap.Error(err))
						} else {

						}
					}
				} else {
					g.logger.Info("fastq data has build")
					files, err := ioutil.ReadDir(types.FASTP_OUT)
					if err != nil {
						g.logger.Error("read clean data fail", zap.Error(err))
					} else {
						names := make([]string, 0)
						for _, v := range files {
							names = append(names, v.Name())
						}
						g.logger.Warn("clean data ", zap.Strings("files", names))
					}
				}
			} else {
				g.logger.Error("fastq input not existed")
			}
		}()
		wg.Wait()
		if b := utils.IsExist(types.FASTP_OUT); b {
			if b := utils.IsExist(types.HISAT2_OUT); !b || (b && utils.Files(types.HISAT2_OUT) == 0) {
				if err := g.fastqc(types.FASTP_OUT); err != nil {
					g.logger.Error("clean -> fastqc data fail", zap.Error(err))
					return err
				}
			} else {
				g.logger.Info("fastqc data has build")
				files, err := ioutil.ReadDir(types.FASTP_OUT)
				if err != nil {
					g.logger.Error("read fastqc data fail", zap.Error(err))
				} else {
					names := make([]string, 0)
					for _, v := range files {
						names = append(names, v.Name())
					}
					g.logger.Warn("fastqc data ", zap.Strings("files", names))
				}
			}
		} else {
			g.logger.Error("fastqc input not existed")
		}

		//clean -> sam
		if b := utils.IsExist(types.FASTP_OUT); b {
			if b := utils.IsExist(types.HISAT2_OUT); !b || (b && utils.Files(types.HISAT2_OUT) == 0) {
				if err := g.bwa(types.FASTP_OUT); err != nil {
					g.logger.Error("clean -> hisat data fail", zap.Error(err))
					return err
				}
			} else {
				g.logger.Info("hisat2 data has build")
				files, err := ioutil.ReadDir(types.FASTP_OUT)
				if err != nil {
					g.logger.Error("read hisat2 data fail", zap.Error(err))
				} else {
					names := make([]string, 0)
					for _, v := range files {
						names = append(names, v.Name())
					}
					g.logger.Warn("hisat2 data ", zap.Strings("files", names))
				}
			}
		} else {
			g.logger.Error("hisat2 input not existed")
		}
		//sam -> bam
		if b := utils.IsExist(types.HISAT2_OUT); b {
			if b := utils.IsExist(types.SORTED_OUT); !b || (b && utils.Files(types.SORTED_OUT) == 0) {
				if err := g.sort(types.HISAT2_OUT); err != nil {
					g.logger.Error("create hisat2 sorted data fail", zap.Error(err))
					return err
				} else {
					g.logger.Info("create hisat2 sorted data success")
				}
			}
		} else {
			g.logger.Error("hisat2 dir input not existed")
		}
		//bam -> count
		if b := utils.IsExist(types.GATK_OUT); b {
			if err := g.buildVCF(); err != nil {
				g.logger.Error("create vcf  fail", zap.Error(err))
				return err
			} else {
				g.logger.Info("create vcf success")
			}
			//if err := g.result(); err != nil {
			//	g.logger.Error("select vcf fail", zap.Error(err))
			//	return err
			//} else {
			//	g.logger.Info("select vcf success")
			//}
		} else {
			g.logger.Error("hisat2 sorted input not existed")
		}
		return nil
	})
	select {
	case err := <-done:
		return err
	case <-ctx.Done():
		return errors.New("time out")
	}
}
func (g *gatkPlugin) getfa() (string, error) {
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
func (g *gatkPlugin) index() error {
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
func (g *gatkPlugin) fastq(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	//fastp -i "${temp}".R1.fastq.gz -I "${temp}".R2.fastq.gz -o ${fastp_out}/"${temp}".R1.clean.fastq.gz -O ${fastp_out}/"${temp}".R2.clean.fastq.gz > ${log}/fastp.log
	//bar := u.bar.NewBar("fastq --> clean.fastq",len(files) / 2)
	pool, err := ants.NewPool(8)
	if err != nil {
		return err
	}
	var (
		tp types.SampleType
		wg sync.WaitGroup
	)
	for _, v := range files {
		switch {
		case strings.HasSuffix(v.Name(), types.R1Sample):
			tp = types.SampleType(1)
		case strings.HasSuffix(v.Name(), types.R1SampleEx):
			tp = types.SampleType(2)
		case strings.HasSuffix(v.Name(), types.R1SampleFq):
			tp = types.SampleType(3)
		case strings.HasSuffix(v.Name(), types.R1SampleFqEx):
			tp = types.SampleType(4)
		default:
			g.logger.Info("file name", zap.String("name", v.Name()))
			continue
		}
		name := v.Name()
		//sample := strings.TrimSuffix(name, tp.Type())
		//g.samples = append(g.samples, sample)
		wg.Add(1)
		log.Println("clean source file", name)
		if err = pool.Submit(func() {
			defer func() {
				//bar.Add(1)
				wg.Done()
				g.logger.Info("build clean file success", zap.String("name", name))
			}()
			temp := strings.TrimSuffix(name, tp.Type())
			g.logger.Info("fastq -> clean", zap.String("source", name), zap.String("target", fmt.Sprintf("%s_R2.clean.fastq.gz", temp)))
			cmd := exec.Command("fastp", "-i",
				fmt.Sprintf("%s/%s%s", types.INPUT, temp, tp.Type()),
				"-I", fmt.Sprintf("%s/%s%s", types.INPUT, temp, strings.Replace(tp.Type(), "1", "2", -1)),
				"-o", fmt.Sprintf("%s/%s%s", types.FASTP_OUT, temp, tp.CleanType()),
				"-O", fmt.Sprintf("%s/%s%s", types.FASTP_OUT, temp, strings.Replace(tp.CleanType(), "1", "2", -1)),
				"-h", fmt.Sprintf("%s/%s.html", types.FASTP_OUT, temp),
				"-j", fmt.Sprintf("%s/%s.json", types.FASTP_OUT, temp))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			g.logger.Info("fastp run ", zap.String("cmd", cmd.String()))
			if err = cmd.Run(); err != nil {
				g.logger.Error("create clean file", zap.Error(err), zap.String("cmd", cmd.String()))
			}
		}); err != nil {
			g.logger.Error("task submit fail", zap.Error(err))
		}
	}
	g.logger.Info("fastq -> clean waiting")
	wg.Wait()
	g.logger.Info("fastq -> clean finished")
	g.tp = tp
	return nil
}
func (g *gatkPlugin) bwa(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	pool, err := ants.NewPool(4)
	if err != nil {
		return err
	}
	var (
		tp     types.SampleType
		wg     sync.WaitGroup
		thread string
	)
	if thread = os.Getenv("THREAD"); thread == "" {
		thread = "8"
	}
	for _, v := range files {
		switch {
		case strings.HasSuffix(v.Name(), types.R1SampleClean):
			tp = types.SampleType(1)
		case strings.HasSuffix(v.Name(), types.R1SampleCleanEx):
			tp = types.SampleType(2)
		case strings.HasSuffix(v.Name(), types.R1SampleCleanFq):
			tp = types.SampleType(3)
		case strings.HasSuffix(v.Name(), types.R1SampleCleanFqEx):
			tp = types.SampleType(4)
		default:
			g.logger.Info("file name", zap.String("name", v.Name()))
			continue
		}
		name := v.Name()
		//sample := strings.TrimSuffix(name, tp.Type())
		//g.samples = append(g.samples, sample)
		wg.Add(1)
		if err = pool.Submit(func() {
			defer func() {
				//bar.Add(1)
				wg.Done()
				g.logger.Info("build bwa file success", zap.String("name", name))
			}()
			temp := strings.TrimSuffix(name, tp.CleanType())
			g.logger.Info("clean -> bwa", zap.String("source", name),
				zap.String("target1", fmt.Sprintf("%s%s", temp, tp.Type())),
				zap.String("target2", fmt.Sprintf("%s%s", temp, strings.Replace(tp.Type(), "1", "2", -1))))
			cmd := exec.Command("bwa", "mem", "-t", thread, "-R",
				fmt.Sprintf(`@RG\tID:group_%s\tLB:library_%s\tPL:illumina\tSM:sample_%s`, temp, temp, temp),
				faFile, fmt.Sprintf("%s/%s%s", types.FASTP_OUT, temp, tp.CleanType()),
				fmt.Sprintf("%s/%s%s", types.FASTP_OUT, temp, strings.Replace(tp.CleanType(), "1", "2", -1)),
				"-o", fmt.Sprintf("%s/%s.sam", types.HISAT2_OUT, temp))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			g.logger.Info("bwa run ", zap.String("cmd", cmd.String()))
			if err = cmd.Run(); err != nil {
				g.logger.Error("create bwa file", zap.Error(err), zap.String("cmd", cmd.String()))
			}
		}); err != nil {
			g.logger.Error("task submit fail", zap.Error(err))
		}
	}
	g.logger.Info("clean -> bwa waiting")
	wg.Wait()
	g.logger.Info("clean -> bwa finished")
	g.tp = tp
	return nil
}
func (g *gatkPlugin) sort(dir string) error {
	var (
		//bar *pgbar.Bar
		wg sync.WaitGroup
	)
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	pool, err := ants.NewPool(8)
	if err != nil {
		return err
	}
	if !utils.SufFiles(types.SORTED_OUT, ".bam") {
		//bar := u.bar.NewBar("sam --> bam",len(files))
		for _, v := range files {
			if strings.HasSuffix(v.Name(), ".sam") {
				name := v.Name()
				wg.Add(1)
				if err = pool.Submit(func() {
					defer func() {
						//bar.Add(1)
						wg.Done()
						g.logger.Info("build sam file success", zap.String("name", name))
					}()
					temp := strings.TrimSuffix(name, ".sam")
					g.logger.Info("sam -> bam", zap.String("source", fmt.Sprintf("%s", name)), zap.String("target", fmt.Sprintf("%s.bam", temp)))
					cmd := exec.Command("samtools", "view",
						"-b", fmt.Sprintf("%s/%s", types.HISAT2_OUT, name),
						"-o", fmt.Sprintf("%s/%s.bam", types.SORTED_OUT, temp))
					//cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					g.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
					if err = cmd.Run(); err != nil {
						g.logger.Error("samtools run fail", zap.Error(err), zap.String("cmd", cmd.String()))
					}
					cmd = exec.Command("samtools", "sort", "-@",
						"4", fmt.Sprintf("%s/%s.bam", types.SORTED_OUT, temp),
						"-o", fmt.Sprintf("%s/%s.sorted.bam", types.SORTED_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					g.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
					if err = cmd.Run(); err != nil {
						g.logger.Error("samtools sorted run fail", zap.Error(err))
					}
					cmd = exec.Command("samtools", "index", "-@",
						"4", fmt.Sprintf("%s/%s.sorted.bam", types.SORTED_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					g.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
					if err = cmd.Run(); err != nil {
						g.logger.Error("samtools index run fail", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					f, err := os.Create(fmt.Sprintf("%s/%s.report", types.REPORT_OUT, temp))
					if err != nil {
						g.logger.Error("samtools create fail", zap.Error(err))
						return
					}
					defer f.Close()
					cmd = exec.Command("samtools", "flagstat",
						fmt.Sprintf("%s/%s.sorted.bam", types.SORTED_OUT, temp))
					cmd.Stdout = f
					cmd.Stderr = os.Stderr
					g.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
					if err = cmd.Run(); err != nil {
						g.logger.Error("samtools index run fail", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
				}); err != nil {
					g.logger.Error("samtools submit task run fail", zap.Error(err))
				}
			}
		}
		g.logger.Info("sam -> bam waiting")
		wg.Wait()
	} else {
		g.logger.Info("create bam success")
	}
	return nil
}
func (g *gatkPlugin) buildVCF() error {
	var (
		err    error
		pool   *ants.Pool
		files  []os.FileInfo
		wg     sync.WaitGroup
		thread int
	)
	//fa,err := g.getfa()
	//if err != nil {
	//	return err
	//}
	//wd, _ = os.Getwd()
	env := os.Getenv("THREAD")
	if env != "" {
		n, _ := strconv.Atoi(env)
		thread = n
	} else {
		thread = 4
	}
	pool, err = ants.NewPool(thread)
	if err != nil {
		return err
	}
	if !utils.IsExist(fmt.Sprintf("%s/merge.g.vcf", types.GATK_OUT)) {
		files, err = ioutil.ReadDir(types.SORTED_OUT)
		if err != nil {
			return err
		}
		bar := g.bar.NewBar("sort -> gatk", len(files)/3)
		for _, v := range files {
			if strings.HasSuffix(v.Name(), ".sorted.bam") {
				wg.Add(1)
				name := v.Name()
				input := path.Join(types.SORTED_OUT, name)
				if err = pool.Submit(func() {
					temp := strings.TrimSuffix(name, ".sorted.bam")
					defer func() {
						bar.Add(1)
						wg.Done()
						g.logger.Info("build vcf file success", zap.String("names", name))
					}()
					//标记重复序列
					start := time.Now()
					cmd := exec.Command("gatk", "MarkDuplicates",
						"-I", input, "-O", fmt.Sprintf("%s/%s.markdup.bam", types.GATK_OUT, temp),
						"-M", fmt.Sprintf("%s/%s.markdup_metrics.txt", types.GATK_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						g.logger.Error("run gatk MarkDuplicates", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
					cmd = exec.Command("samtools", "index",
						fmt.Sprintf("%s/%s.markdup.bam", types.GATK_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					g.logger.Info("run cmd ", zap.String("cmd", cmd.String()))
					if err = cmd.Run(); err != nil {
						g.logger.Error("run samtools bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
					start = time.Now()
					cmd = exec.Command("gatk", "HaplotypeCaller", "-R", "gene.fa",
						//"--java-options", `"-Xmx15G -Djava.io.tmpdir=./"`,
						"--emit-ref-confidence", "GVCF", "-I", fmt.Sprintf("%s/%s.markdup.bam", types.GATK_OUT, temp),
						"-O", fmt.Sprintf("%s/%s.g.vcf", types.GATK_G_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						g.logger.Error("run gatk HaplotypeCaller bam", zap.Error(err), zap.String("cmd", cmd.String()))
						return
					}
					g.logger.Info("run cmd finished", zap.String("cmd", cmd.String()), zap.Duration("lost", time.Since(start)))
				}); err != nil {
					g.logger.Error("pool run fail", zap.Error(err))
				}
			}
		}
		wg.Wait()
	}
	return nil
}
func (g *gatkPlugin) fastqc(dir string) error {
	//fastqc -t 12 -o out_path sample1_1.fq sample1_2.fq
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	//bar := u.bar.NewBar("clean.fastq --> sam",len(files) / 2 )
	//todo 内存占用较大（单个5G）
	pool, err := ants.NewPool(4)
	if err != nil {
		return err
	}
	var (
		wg sync.WaitGroup
		//flag bool
		tp types.SampleType
	)
	for _, v := range files {
		if g.tp == 0 {
			switch {
			case strings.HasSuffix(v.Name(), types.R1SampleClean):
				tp = types.SampleType(1)
			case strings.HasSuffix(v.Name(), types.R1SampleCleanEx):
				tp = types.SampleType(2)
			case strings.HasSuffix(v.Name(), types.R1SampleCleanFq):
				tp = types.SampleType(3)
			case strings.HasSuffix(v.Name(), types.R1SampleCleanFqEx):
				tp = types.SampleType(4)
			default:
				g.logger.Info("file name", zap.String("name", v.Name()))
				continue
			}
			g.tp = tp
		} else {
			tp = g.tp
		}
		name := v.Name()
		g.logger.Info("start fastqc", zap.String("name", name), zap.String("tp", g.tp.CleanType()))
		wg.Add(1)
		if err = pool.Submit(func() {
			defer func() {
				//bar.Add(1)
				wg.Done()
				//log.Println(name, "build map file success", types.HISAT2_OUT, name)
			}()
			if strings.HasSuffix(name, g.tp.CleanType()) {
				temp := strings.TrimSuffix(name, g.tp.CleanType())
				cmd := exec.Command("fastqc", "-t", "4",
					"-o", types.FASTP_QC_OUT,
					fmt.Sprintf("%s/%s%s", types.FASTP_OUT, temp, g.tp.CleanType()),
					fmt.Sprintf("%s/%s%s", types.FASTP_OUT, temp, strings.Replace(g.tp.CleanType(), "1", "2", -1)))
				//cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				g.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					g.logger.Error("fastqc run fail", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}
		}); err != nil {
			g.logger.Error("hisat2 submit task run fail", zap.Error(err))
		}
	}
	if tp == 0 {
		g.logger.Error("please check raw data(fastqc)")
		return errors.New("please check clean data")
	}
	g.logger.Info("clean -> hisat2 waiting")
	wg.Wait()
	g.logger.Info("clean -> hisat2 finished")
	//解压
	files, err = ioutil.ReadDir(types.FASTP_QC_OUT)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".zip") {
			err = utils.UnZip(fmt.Sprintf("%s/%s", types.FASTP_QC_OUT, v.Name()), types.FASTP_QC_OUT)
			if err != nil {
				return err
			}
		}
	}
	return nil
}
