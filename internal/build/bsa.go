package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"
	"sync"
)

const (
	bam1 = "/data/output/smsnpMapper_out/smhisat2_out1/acc.sorted.bam"
	bam2 = "/data/output/smsnpMapper_out/smhisat2_out2/acc.sorted.bam"
)

type bsaPlugin struct {
	tp      types.SampleType
	samples []string
	logger  *zap.Logger
}

func NewBsaPlugin(logger *zap.Logger) types.Plugin {
	return &bsaPlugin{
		logger: logger,
	}
}
func (g *bsaPlugin) check() error {
	//check log
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
			g.logger.Warn("log dir", zap.Strings("files", names))
		}
	}
	//check index
	if b := utils.IsExist(types.GENOME_INDEX); !b {
		cmd := exec.Command("mkdir", "-p", types.GENOME_INDEX)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create index dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.GENOME_INDEX)
		if err != nil {
			g.logger.Error("index dir existed", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("index dir", zap.Strings("files", names))
		}
	}
	return nil
}
func (g *bsaPlugin) Name() string {
	return "bsaPlugin"
}
func (g *bsaPlugin) Build(ctx context.Context) error {
	if err := g.check(); err != nil {
		return err
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
	if b := utils.IsExist(types.BSA_DNA_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.BSA_DNA_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create bsa dna dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.BSA_DNA_OUT)
		if err != nil {
			g.logger.Error("existed bsa dna dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("bsa dna out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.FASTP_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.FASTP_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("create clean dir existed", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.FASTP_OUT)
		if err != nil {
			g.logger.Error("create clean dir existed", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			g.logger.Warn("create clean dir", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.INPUT); b {
		if b := utils.IsExist(types.FASTP_OUT); !b || (b && utils.Files(types.FASTP_OUT) == 0) {
			if err := g.fastq(types.INPUT); err != nil {
				g.logger.Error("create clean existed", zap.Error(err))
			} else {
				if err = common.BuildReport(); err != nil {
					g.logger.Error("fastp -> report", zap.Error(err))
				}
			}
		} else {
			g.logger.Info("clean file existed")
			files, err := ioutil.ReadDir(types.FASTP_OUT)
			if err != nil {
				g.logger.Error("create clean fail", zap.Error(err))
			} else {
				names := make([]string, 0)
				for _, v := range files {
					names = append(names, v.Name())
				}
				g.logger.Warn("clean dir", zap.Strings("files", names))
			}
		}
	} else {
		if b := utils.IsExist(types.FASTP_OUT); !b {
			g.logger.Warn("fastp input not existed")
			return errors.New("please check raw data")
		}
	}
	return g.pipeline()
}

func (g *bsaPlugin) fastq(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	//fastp -i "${temp}".R1.fastq.gz -I "${temp}".R2.fastq.gz -o ${fastp_out}/"${temp}".R1.clean.fastq.gz -O ${fastp_out}/"${temp}".R2.clean.fastq.gz > ${log}/fastp.log
	//bar := u.bar.NewBar("fastq --> clean.fastq",len(files) / 2)
	pool, err := ants.NewPool(4)
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
		sample := strings.TrimSuffix(name, tp.Type())
		g.samples = append(g.samples, sample)
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
				"-j", fmt.Sprintf("%s/%s.json", types.FASTP_OUT, temp),
				">", types.LOG)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err = cmd.Run(); err != nil {
				g.logger.Error("create clean file", zap.Error(err), zap.String("cmd", cmd.String()))
			}
		}); err != nil {
			g.logger.Error("task submit fail", zap.Error(err))
		}
	}
	g.logger.Info("fastq -> clean waiting")
	wg.Wait()
	if len(g.samples) != 2 {
		g.logger.Error("sample num is error", zap.Strings("samples", g.samples))
		return errors.New("sample num is error")
	}
	g.logger.Info("fastq -> clean finished")
	g.tp = tp
	return nil
}
func (g *bsaPlugin) getfa() (string, error) {
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
func (g *bsaPlugin) fastqc(dir string) error {
	//fastqc -t 12 -o out_path sample1_1.fq sample1_2.fq
	cmd := exec.Command("fastqc", "-t", "4",
		"-o", types.REPORT_OUT,
		fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[0], g.tp.CleanType()),
		fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[0], strings.Replace(g.tp.CleanType(), "1", "2", -1)),
		fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[1], g.tp.CleanType()),
		fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[1], strings.Replace(g.tp.CleanType(), "1", "2", -1)))
	g.logger.Info("dna fastqc", zap.String("cmd", cmd.String()))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		g.logger.Error("fastqc", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	//解压
	files, err := ioutil.ReadDir(types.REPORT_OUT)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".zip") {
			err = utils.UnZip(fmt.Sprintf("%s/%s", types.REPORT_OUT, v.Name()), types.REPORT_OUT)
			if err != nil {
				return err
			}
		}
	}
	return nil
}
func (g *bsaPlugin) pipeline() error {
	fa, err := g.getfa()
	if err != nil {
		return err
	}
	var tp types.SampleType
	if len(g.samples) != 2 {
		g.logger.Warn("sample num is error", zap.Strings("samples", g.samples))
		files, err := ioutil.ReadDir(types.FASTP_OUT)
		if err != nil {
			return err
		}
		for _, v := range files {
			if strings.HasSuffix(v.Name(), types.R1SampleClean) {
				sample := strings.TrimSuffix(v.Name(), types.R1SampleClean)
				g.samples = append(g.samples, sample)
				tp = types.SampleType(1)
			} else if strings.HasSuffix(v.Name(), types.R1SampleCleanFqEx) {
				sample := strings.TrimSuffix(v.Name(), types.R1SampleCleanFqEx)
				g.samples = append(g.samples, sample)
				tp = types.SampleType(4)
			} else if strings.HasSuffix(v.Name(), types.R1SampleCleanEx) {
				sample := strings.TrimSuffix(v.Name(), types.R1SampleCleanEx)
				g.samples = append(g.samples, sample)
				tp = types.SampleType(2)
			} else if strings.HasSuffix(v.Name(), types.R1SampleCleanFq) {
				sample := strings.TrimSuffix(v.Name(), types.R1SampleCleanFq)
				g.samples = append(g.samples, sample)
				tp = types.SampleType(3)
			}
		}
		if len(g.samples) != 2 {
			return errors.New("sample num is error")
		}
		g.tp = tp
	}
	//g.tp = tp
	//build clean qc
	err = g.fastqc(types.FASTP_OUT)
	if err != nil {
		return err
	}
	var (
		step   string
		window string
		thread string
	)
	if step = os.Getenv("STEP"); step == "" {
		step = "500000"
	}
	if window = os.Getenv("WINDOW"); window == "" {
		window = "1000000"
	}
	if thread = os.Getenv("THREAD"); thread == "" {
		thread = "10"
	}
	cmd := exec.Command("DNA_BSA_pipeline.pl", "-f", fa,
		"-b", types.BSA_GENOME_PREFIX,
		"-3", fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[0], g.tp.CleanType()),
		"-4", fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[0], strings.Replace(g.tp.CleanType(), "1", "2", -1)),
		"-5", fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[1], g.tp.CleanType()),
		"-6", fmt.Sprintf("%s/%s%s", types.FASTP_OUT, g.samples[1], strings.Replace(g.tp.CleanType(), "1", "2", -1)),
		"-a", thread, "-s", step, "-r", window, "-1", g.samples[0], "-2", g.samples[1])
	g.logger.Info("dna pipeline", zap.String("cmd", cmd.String()))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	g.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		g.logger.Error("dna pipeline fail", zap.Error(err))
		return err
	}
	return nil
}

func (g *bsaPlugin) report() {
	fa, err := g.getfa()
	if err != nil {
		g.logger.Error("samtools mpileup fail", zap.Error(err))
		return
	}
	var wg sync.WaitGroup
	go func() {
		defer wg.Done()
		cmd := exec.Command("/bin/sh", "-c",
			fmt.Sprintf(`samtools mpileup -f %s %s |perl -alne '{$pos=int($F[1]/1000); $key="$F[0]\t$pos";$GC{$key}++ if $F[2]=~/[GC]/;$counts_sum{$key}+=$F[3];$number{$key}++;}END{print "$_\t$number{$_}\t$GC{$_}\t$counts_sum{$_}" foreach sort{$a<=>$b} keys %number}' > %s/%s.txt`,
				fa, bam1,
				types.REPORT_OUT, "bam1"))
		g.logger.Info("samtools mpileup", zap.String("cmd", cmd.String()))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			g.logger.Error("samtools mpileup fail", zap.Error(err))
			return
		}
	}()
	go func() {
		defer wg.Done()
		cmd := exec.Command("/bin/sh", "-c",
			fmt.Sprintf(`samtools mpileup -f %s %s |perl -alne '{$pos=int($F[1]/1000); $key="$F[0]\t$pos";$GC{$key}++ if $F[2]=~/[GC]/;$counts_sum{$key}+=$F[3];$number{$key}++;}END{print "$_\t$number{$_}\t$GC{$_}\t$counts_sum{$_}" foreach sort{$a<=>$b} keys %number}' > %s/%s.txt`,
				fa, bam2,
				types.REPORT_OUT, "bam2"))
		g.logger.Info("samtools mpileup", zap.String("cmd", cmd.String()))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			g.logger.Error("samtools mpileup fail", zap.Error(err))
			return
		}
	}()
	wg.Wait()
}
