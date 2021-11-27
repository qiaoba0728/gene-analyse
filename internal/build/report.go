package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
	"strings"
	"sync"
)

const (
//bam1 = "/data/output/smsnpMapper_out/smhisat2_out1/acc.sorted.bam"
//bam2 = "/data/output/smsnpMapper_out/smhisat2_out2/acc.sorted.bam"
)

type reportPlugin struct {
	logger *zap.Logger
}

func NewReportPlugin(logger *zap.Logger) types.Plugin {
	return &reportPlugin{
		logger: logger,
	}
}

func (r *reportPlugin) Build(ctx context.Context) error {
	var (
		err error
	)
	if err = r.buildBed12(); err != nil {
		return err
	}
	var wg sync.WaitGroup
	if utils.IsExist(types.SORTED_OUT) {
		wg.Add(5)
		go func() {
			defer wg.Done()
			if err = r.readDistribution(types.SORTED_OUT); err != nil {
				r.logger.Error("readDistribution", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.geneBodyCoverage(types.SORTED_OUT); err != nil {
				r.logger.Error("geneBodyCoverage", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.rpkmSaturation(types.SORTED_OUT); err != nil {
				r.logger.Error("rpkmSaturation", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.geneDepthCoverage(types.SORTED_OUT); err != nil {
				r.logger.Error("geneDepthCoverage", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.stringTieHandle(types.SORTED_OUT); err != nil {
				r.logger.Error("stringTieHandle", zap.Error(err))
			}
		}()
	} else {
		wg.Add(4)
		go func() {
			defer wg.Done()
			if err = r.readDistributionEx(); err != nil {
				r.logger.Error("readDistribution", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.geneBodyCoverageEx(); err != nil {
				r.logger.Error("geneBodyCoverage", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.rpkmSaturationEx(); err != nil {
				r.logger.Error("rpkmSaturation", zap.Error(err))
			}
		}()
		go func() {
			defer wg.Done()
			if err = r.geneDepthCoverageEx(types.SORTED_OUT); err != nil {
				r.logger.Error("geneDepthCoverageEx", zap.Error(err))
			}
		}()
	}
	wg.Wait()
	return nil
}
func (r *reportPlugin) readDistribution(dir string) error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
	)
	pool, err = ants.NewPool(2)
	if err != nil {
		return err
	}
	files, err = ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			//r := path.Join(wd,"script","featurecounts.R")
			input := path.Join(dir, name)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				f, err := os.Create(fmt.Sprintf("%s/%s_read_distribution.log", types.REPORT_OUT, temp))
				if err != nil {
					r.logger.Error("read_distribution create fail", zap.Error(err))
					return
				}
				defer f.Close()
				cmd := exec.Command("read_distribution.py", "-i", input,
					"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"))
				defer func() {
					//bar.Add(1)
					wg.Done()
					r.logger.Info("read_distribution file success", zap.String("names", name))
				}()
				r.logger.Info("run readDistribution", zap.String("cmd", cmd.String()))
				cmd.Stdout = f
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					r.logger.Error("read_distribution bam", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				r.logger.Error("pool run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
func (r *reportPlugin) readDistributionEx() error {
	var (
		err  error
		pool *ants.Pool
		wg   sync.WaitGroup
	)
	pool, err = ants.NewPool(2)
	if err != nil {
		return err
	}
	wg.Add(2)
	err = pool.Submit(func() {
		defer func() {
			wg.Done()
		}()
		f, err := os.Create(fmt.Sprintf("%s/%s_read_distribution.log", types.REPORT_OUT, "bam1"))
		if err != nil {
			r.logger.Error("read_distribution create fail", zap.Error(err))
			return
		}
		defer f.Close()
		cmd := exec.Command("read_distribution.py", "-i", bam1,
			"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"))
		r.logger.Info("run read_quality", zap.String("cmd", cmd.String()))
		cmd.Stdout = f
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("read_quality bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	})
	if err != nil {
		return err
	}
	err = pool.Submit(func() {
		defer func() {
			wg.Done()
		}()
		f, err := os.Create(fmt.Sprintf("%s/%s_read_distribution.log", types.REPORT_OUT, "bam2"))
		if err != nil {
			r.logger.Error("read_distribution create fail", zap.Error(err))
			return
		}
		defer f.Close()
		cmd := exec.Command("read_distribution.py", "-i", bam2,
			"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"))
		r.logger.Info("run read_quality", zap.String("cmd", cmd.String()))
		cmd.Stdout = f
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("read_quality bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	})
	wg.Wait()
	return nil
}
func (r *reportPlugin) geneBodyCoverage(dir string) error {
	var (
		err    error
		pool   *ants.Pool
		files  []os.FileInfo
		wg     sync.WaitGroup
		inputs []string
	)
	pool, err = ants.NewPool(4)
	if err != nil {
		return err
	}
	//bar := e.bar.NewBar("featurecounts",len(files))
	files, err = ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			//r := path.Join(wd,"script","featurecounts.R")
			input := path.Join(dir, name)
			inputs = append(inputs, input)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				cmd := exec.Command("geneBody_coverage.py", "-i", input,
					"-f", "png",
					"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
					"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, temp))
				defer func() {
					//bar.Add(1)
					wg.Done()
					r.logger.Info("geneBody_coverage file success", zap.String("name", name))
				}()
				r.logger.Info("run geneBodyCoverage", zap.String("cmd", cmd.String()))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					r.logger.Error("geneBody_coverage bam", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				r.logger.Error("pool run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	cmd := exec.Command("geneBody_coverage.py", "-i", strings.Join(inputs, ","),
		"-f", "png",
		"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
		"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "all"))
	r.logger.Info("run geneBodyCoverage", zap.String("cmd", cmd.String()))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err = cmd.Run(); err != nil {
		r.logger.Error("build all geneBody_coverage bam", zap.Error(err), zap.String("cmd", cmd.String()))
	}
	return nil
}
func (r *reportPlugin) geneBodyCoverageEx() error {
	var (
		err    error
		inputs []string
	)
	inputs = []string{bam1, bam2}
	cmd := exec.Command("geneBody_coverage.py", "-i", strings.Join(inputs, ","),
		"-f", "png",
		"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
		"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "all"))
	r.logger.Info("run geneBodyCoverage", zap.String("cmd", cmd.String()))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err = cmd.Run(); err != nil {
		r.logger.Error("build all geneBody_coverage bam", zap.Error(err), zap.String("cmd", cmd.String()))
	}
	return nil
}
func (r *reportPlugin) geneDepthCoverage(dir string) error {
	var (
		err    error
		pool   *ants.Pool
		files  []os.FileInfo
		wg     sync.WaitGroup
		inputs []string
	)
	pool, err = ants.NewPool(4)
	if err != nil {
		return err
	}
	//bar := e.bar.NewBar("featurecounts",len(files))
	files, err = ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			//r := path.Join(wd,"script","featurecounts.R")
			input := path.Join(dir, name)
			inputs = append(inputs, input)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				if !utils.IsExist(fmt.Sprintf("%s/%s", types.REPORT_OUT, temp)) {
					cmd := exec.Command("mkdir", "-p", fmt.Sprintf("%s/%s", types.REPORT_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						r.logger.Error("mkdir bam", zap.Error(err), zap.String("cmd", cmd.String()))
					}
				}
				cmd := exec.Command("/work/bamdst", "-p",
					fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
					"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, temp),
					input)
				defer func() {
					//bar.Add(1)
					wg.Done()
					r.logger.Info("geneDepthCoverage file success", zap.String("name", name))
				}()
				r.logger.Info("run geneDepthCoverage", zap.String("cmd", cmd.String()))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					r.logger.Error("geneDepthCoverage bam", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				r.logger.Error("pool run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
func (r *reportPlugin) geneDepthCoverageEx(dir string) error {
	var (
		err error
	)
	var wg sync.WaitGroup
	wg.Add(2)
	go func() {
		defer wg.Done()
		if !utils.IsExist(fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam1")) {
			cmd := exec.Command("mkdir", "-p", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam1"))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err = cmd.Run(); err != nil {
				r.logger.Error("mkdir bam", zap.Error(err), zap.String("cmd", cmd.String()))
			}
		}
		cmd := exec.Command("/work/bamdst", "-p",
			fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
			"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam1"),
			bam1)
		r.logger.Info("run geneDepthCoverageEx", zap.String("cmd", cmd.String()))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("build all geneDepthCoverageEx bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	}()
	go func() {
		defer wg.Done()
		if !utils.IsExist(fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam2")) {
			cmd := exec.Command("mkdir", "-p", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam2"))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err = cmd.Run(); err != nil {
				r.logger.Error("mkdir bam", zap.Error(err), zap.String("cmd", cmd.String()))
			}
		}
		cmd := exec.Command("/work/bamdst", "-p",
			fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
			"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam2"),
			bam2)
		r.logger.Info("run geneDepthCoverageEx", zap.String("cmd", cmd.String()))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("build all geneDepthCoverageEx bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	}()
	wg.Wait()
	return nil
}
func (r *reportPlugin) getGTF() (string, error) {
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
func (r *reportPlugin) buildBed12() error {
	var (
		gtf string
		err error
	)
	if utils.IsExist(fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12")) {
		r.logger.Info("bed file has existed")
		return nil
	}
	if gtf, err = r.getGTF(); err != nil {
		return err
	}
	f, err := os.Create(fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"))
	if err != nil {
		r.logger.Error("read_distribution create fail", zap.Error(err))
		return err
	}
	defer f.Close()
	cmd := exec.Command("/work/gtf2bed12.perl", gtf)
	cmd.Stdout = f
	cmd.Stderr = os.Stderr
	r.logger.Info("run gtf to bed", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		r.logger.Error("run build bed12", zap.Error(err), zap.String("cmd", cmd.String()))
	}
	return err
}
func (r *reportPlugin) readQuality(dir string) error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
	)
	pool, err = ants.NewPool(1)
	if err != nil {
		return err
	}
	files, err = ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			//r := path.Join(wd,"script","featurecounts.R")
			input := path.Join(dir, name)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				cmd := exec.Command("read_quality.py", "-i", input,
					"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, temp))
				defer func() {
					//bar.Add(1)
					wg.Done()
					r.logger.Info("read_quality file success", zap.String("name", name))
				}()
				r.logger.Info("run read_quality", zap.String("cmd", cmd.String()))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					r.logger.Error("geneBody_coverage bam", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				r.logger.Error("pool run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
func (r *reportPlugin) readQualityEx() error {
	var (
		err  error
		pool *ants.Pool
		wg   sync.WaitGroup
	)
	pool, err = ants.NewPool(2)
	if err != nil {
		return err
	}
	wg.Add(2)
	err = pool.Submit(func() {
		defer func() {
			wg.Done()
		}()
		cmd := exec.Command("read_quality.py", "-i", bam1,
			"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam1"))
		r.logger.Info("run read_quality", zap.String("cmd", cmd.String()))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("read_quality bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	})
	if err != nil {
		return err
	}
	err = pool.Submit(func() {
		defer func() {
			wg.Done()
		}()
		cmd := exec.Command("read_quality.py", "-i", bam2,
			"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam2"))
		r.logger.Info("run read_quality", zap.String("cmd", cmd.String()))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("read_quality bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	})
	wg.Wait()
	return err
}
func (r *reportPlugin) rpkmSaturation(dir string) error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
	)
	pool, err = ants.NewPool(2)
	if err != nil {
		return err
	}
	files, err = ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			//r := path.Join(wd,"script","featurecounts.R")
			input := path.Join(dir, name)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				cmd := exec.Command("RPKM_saturation.py", "-i", input,
					//"-f", "png",
					"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
					"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, temp))
				defer func() {
					//bar.Add(1)
					wg.Done()
					r.logger.Info("RPKM_saturation file success", zap.String("name", name))
				}()
				r.logger.Info("run RPKM_saturation", zap.String("cmd", cmd.String()))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					r.logger.Error("RPKM_saturation bam", zap.Error(err), zap.String("cmd", cmd.String()))
				} else {
					//出png的图片
					cmd = exec.Command("sed", "-i", "s/pdf/png/g",
						fmt.Sprintf("%s/%s.saturation.r", types.REPORT_OUT, temp))
					cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					if err = cmd.Run(); err != nil {
						r.logger.Error("RPKM_saturation bam", zap.Error(err), zap.String("cmd", cmd.String()))
					} else {
						cmd = exec.Command("Rscript", fmt.Sprintf("%s/%s.saturation.r", types.REPORT_OUT, temp))
						cmd.Stdout = os.Stdout
						cmd.Stderr = os.Stderr
						if err = cmd.Run(); err != nil {
							r.logger.Error("RPKM_saturation bam", zap.Error(err), zap.String("cmd", cmd.String()))
						}
					}
				}
			}); err != nil {
				r.logger.Error("pool RPKM_saturation run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
func (r *reportPlugin) rpkmSaturationEx() error {
	var (
		err  error
		pool *ants.Pool
		wg   sync.WaitGroup
	)
	pool, err = ants.NewPool(2)
	if err != nil {
		return err
	}
	wg.Add(2)
	err = pool.Submit(func() {
		defer func() {
			wg.Done()
		}()
		cmd := exec.Command("RPKM_saturation.py", "-i", bam1,
			"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
			"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam1"))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("RPKM_saturation bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	})
	if err != nil {
		return err
	}
	err = pool.Submit(func() {
		defer func() {
			wg.Done()
		}()
		cmd := exec.Command("RPKM_saturation.py", "-i", bam2,
			"-r", fmt.Sprintf("%s/%s", types.REFERENCES, "gtf.bed12"),
			"-o", fmt.Sprintf("%s/%s", types.REPORT_OUT, "bam2"))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err = cmd.Run(); err != nil {
			r.logger.Error("RPKM_saturation bam", zap.Error(err), zap.String("cmd", cmd.String()))
		}
	})
	wg.Wait()
	return nil
}

func (r *reportPlugin) stringTieHandle(dir string) error {
	var (
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
		gtf   string
	)
	pool, err = ants.NewPool(2)
	if err != nil {
		return err
	}
	files, err = ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			//r := path.Join(wd,"script","featurecounts.R")
			input := path.Join(dir, name)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				cmd := exec.Command("stringtie", input,
					//"-f", "png",
					"-l", temp,
					"-o", fmt.Sprintf("%s/%s.transcripts.gtf", types.REPORT_OUT, temp))
				defer func() {
					//bar.Add(1)
					wg.Done()
					r.logger.Info("build transcripts file success", zap.String("name", name))
				}()
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				r.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					r.logger.Error("build transcripts bam", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				r.logger.Error("pool transcripts run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	cmd := exec.Command("/bin/sh", "-c", fmt.Sprintf("find %s -name *.transcripts.gtf > %s/merglist.txt", types.REPORT_OUT, types.REPORT_OUT))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	r.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		r.logger.Error("build transcripts bam", zap.Error(err), zap.String("cmd", cmd.String()))
	}
	if gtf, err = r.getGTF(); err != nil {
		return err
	}
	cmd = exec.Command("/bin/sh", "-c", fmt.Sprintf("stringtie --merge -p 8 -G %s -o %s/stringtie_merged.gtf %s/merglist.txt", gtf, types.REPORT_OUT, types.REPORT_OUT))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	r.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		r.logger.Error("build merge bam", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	cmd = exec.Command("/bin/sh", "-c", fmt.Sprintf("gffcompare -r %s -G %s/stringtie_merged.gtf", gtf, types.REPORT_OUT))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	r.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		r.logger.Error("gffcompare gtf", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	cmd = exec.Command("/bin/sh", "-c", fmt.Sprintf("cp -r /work/gffcmp.* %s", types.REPORT_OUT))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	r.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
	if err = cmd.Run(); err != nil {
		r.logger.Error("cp gtf", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	return nil
}
func (r *reportPlugin) Name() string {
	return "reportPlugin"
}
