package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qianlnk/pgbar"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
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

type expressionPlugin struct {
	bar    *pgbar.Pgbar
	logger *zap.Logger
}

func NewExpressionPlugin(logger *zap.Logger) types.Plugin {
	return &expressionPlugin{
		bar:    pgbar.New("expression\n"),
		logger: logger,
	}
}

//检查注释文件是否存在
func (e *expressionPlugin) check() error {
	wd, _ := os.Getwd()
	scriptDir := path.Join(wd, "script")
	if !utils.IsExist(scriptDir) {
		cmd := exec.Command("mkdir", "-p", scriptDir)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			e.logger.Error("create script fail", zap.Error(err))
			return err
		}
	}
	err := utils.WriteFile("featurecounts.R", scripts.FEATURE_SCRIPT)
	if err != nil {
		return err
	}
	err = utils.WriteFile("matrix_count.R", scripts.MATRIX_COUNT)
	if err != nil {
		return err
	}
	err = utils.WriteFile("matrix_fpkm.R", scripts.MATRIX_FPKM)
	if err != nil {
		return err
	}
	err = utils.WriteFile("fast_report.R", scripts.FAST_REPORT)
	if err != nil {
		return err
	}
	err = utils.WriteFile("fpkm_report.R", scripts.FPKM_REPORT)
	if err != nil {
		return err
	}
	err = utils.WriteFile("matrix_tpm.R", scripts.MATRIX_TPM)
	if err != nil {
		return err
	}
	err = utils.WriteFile("pca_count.R", scripts.PCA_REPORT)
	if err != nil {
		return err
	}
	err = utils.WriteFile("freq_count.R", scripts.FREQ_REPORT)
	if err != nil {
		return err
	}
	err = utils.WriteFile("heatmap_count.R", scripts.HEATMAP_REPORT)
	if err != nil {
		return err
	}
	if b := utils.IsExist(types.SORTED_OUT); !b {
		e.logger.Error("sorted file not existed")
		return errors.New("sorted file not existed")
	} else {
		files, err := ioutil.ReadDir(types.SORTED_OUT)
		if err != nil {
			e.logger.Error("read sorted file fail", zap.Error(err))
			return err
		} else {
			count := 0
			for _, v := range files {
				if strings.HasSuffix(v.Name(), ".sorted.bam.bai") {
					count++
				}
			}
			length := len(files)
			if length%3 != 0 || count != length/3 {
				e.logger.Error("sorted file num error", zap.Int("files", len(files)), zap.Int("count", count))
				return errors.New("sorted file num error")
			}
		}
	}
	if b := utils.IsExist(types.EXPRESSION_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.EXPRESSION_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			e.logger.Error("create expression fail", zap.Error(err))
			return err
		}
	} else {
		files, err := ioutil.ReadDir(types.EXPRESSION_OUT)
		if err != nil {
			e.logger.Error("read expression file fail", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			e.logger.Warn("expression files", zap.Strings("names", names))
		}
		cmd := exec.Command("rm", "-rf", fmt.Sprintf("%s/%s", types.EXPRESSION_OUT, "*"))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			e.logger.Error("delete expression dir fail", zap.Error(err))
			return err
		}
	}
	return nil
}
func (e *expressionPlugin) Name() string {
	return "expressionPlugin"
}
func (e *expressionPlugin) Build(ctx context.Context) error {
	if err := e.check(); err != nil {
		return err
	}
	if utils.IsExist(fmt.Sprintf("%s/%s", types.EXPRESSION_OUT, "merge.count")) {
		return errors.New("merge is existed")
	}
	if utils.Files(types.EXPRESSION_OUT) == 0 {
		if err := e.feature(); err != nil {
			return err
		}
	} else {
		files, err := ioutil.ReadDir(types.EXPRESSION_OUT)
		if err != nil {
			e.logger.Error("read expression file fail", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			e.logger.Warn("merge expression files", zap.Strings("names", names))
		}
	}
	if err := e.merge(); err != nil {
		return err
	}
	if err := e.matrix(); err != nil {
		return err
	}
	if err := e.buildReport(); err != nil {
		return err
	}
	return nil
}
func (e *expressionPlugin) getGTF() (string, error) {
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
func (e *expressionPlugin) feature() error {
	var (
		gtf   string
		wd    string
		err   error
		pool  *ants.Pool
		files []os.FileInfo
		wg    sync.WaitGroup
	)
	if gtf, err = e.getGTF(); err != nil {
		return err
	}
	wd, _ = os.Getwd()
	pool, err = ants.NewPool(20)
	if err != nil {
		return err
	}
	//bar := e.bar.NewBar("featurecounts",len(files))
	files, err = ioutil.ReadDir(types.SORTED_OUT)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			wg.Add(1)
			name := v.Name()
			r := path.Join(wd, "script", "featurecounts.R")
			input := path.Join(types.SORTED_OUT, name)
			if err = pool.Submit(func() {
				temp := strings.TrimSuffix(name, ".sorted.bam")
				cmd := exec.Command("Rscript", r, input,
					gtf,
					"8",
					fmt.Sprintf("%s/%s", types.EXPRESSION_OUT, temp),
					">", types.LOG)
				defer func() {
					//bar.Add(1)
					wg.Done()
					e.logger.Info("build featurecounts file success", zap.String("names", name))
				}()
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				if err = cmd.Run(); err != nil {
					e.logger.Error("run sorted bam", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				e.logger.Error("pool run fail", zap.Error(err))
			}
		}
	}
	wg.Wait()
	return nil
}
func (e *expressionPlugin) merge() error {
	//mergePath := path.Join(types.EXPRESSION_OUT,"merge.count")
	command := fmt.Sprintf(`cd %s;perl -lne 'if ($ARGV=~/(.*).count/){print "$1\t$_"}' *.count > %s`, types.EXPRESSION_OUT, "merge.count")
	cmd := exec.Command("/bin/sh", "-c", command)
	//cmd := exec.Command("perl","-lne",`if ($ARGV=~/(.*).count/){print "$1\t$_"}`,
	//	path.Join(types.EXPRESSION_OUT,"*.count"),
	//	">",mergePath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	command = fmt.Sprintf(`cd %s;awk -F"\t" '{print $1"\t"$2"\t"$3}' %s > %s`, types.EXPRESSION_OUT, "merge.count", "gene.count")
	cmd = exec.Command("/bin/sh", `-c`, command)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	command = fmt.Sprintf(`cd %s;awk -F"\t" '{print $1"\t"$2"\t"$5}' %s > %s`, types.EXPRESSION_OUT, "merge.count", "gene.tpm")
	cmd = exec.Command("/bin/sh", `-c`, command)
	//cmd = exec.Command("awk",`-F"\t"`,`{print $1"\t"$2"\t"$5}`,
	//	mergePath,
	//	">",path.Join(types.EXPRESSION_OUT,"gene.tpm"))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	command = fmt.Sprintf(`cd %s;awk -F"\t" '{print $1"\t"$2"\t"$4}' %s > %s`, types.EXPRESSION_OUT, "merge.count", "gene.fpkm")
	cmd = exec.Command("/bin/sh", `-c`, command)
	//cmd = exec.Command("awk",`-F"\t"`,`{print $1"\t"$2"\t"$4}`,
	//	mergePath,
	//	">",path.Join(types.EXPRESSION_OUT,"gene.fpkm"))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	return nil
}
func (e *expressionPlugin) matrix() error {
	wd, _ := os.Getwd()
	input := path.Join(types.EXPRESSION_OUT, "gene.count")
	output := path.Join(types.EXPRESSION_OUT, "gene_count.csv")
	outputNumber := path.Join(types.EXPRESSION_OUT, "gene_count_number.csv")
	outputXls := path.Join(types.EXPRESSION_OUT, "gene_count.xls")
	cmd := exec.Command("Rscript", path.Join(wd, "script", "matrix_count.R"),
		input, output, outputNumber, outputXls)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	input = path.Join(types.EXPRESSION_OUT, "gene.tpm")
	output = path.Join(types.EXPRESSION_OUT, "gene_tpm.csv")
	outputXls = path.Join(types.EXPRESSION_OUT, "gene_tpm.xls")
	cmd = exec.Command("Rscript", path.Join(wd, "script", "matrix_tpm.R"),
		input, output, outputXls)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	input = path.Join(types.EXPRESSION_OUT, "gene.fpkm")
	output = path.Join(types.EXPRESSION_OUT, "gene_fpkm.csv")
	outputXls = path.Join(types.EXPRESSION_OUT, "gene_fpkm.xls")
	cmd = exec.Command("Rscript", path.Join(wd, "script", "matrix_fpkm.R"),
		input, output, outputXls)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	return utils.BuildConfig(path.Join(types.EXPRESSION_OUT, "gene_count.csv"), 3, "/data/build_config.json")
}
func (e *expressionPlugin) buildReport() error {
	wd, _ := os.Getwd()
	files, _ := ioutil.ReadDir(types.FASTP_OUT)
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".json") {
			name := strings.TrimSuffix(v.Name(), ".json")
			cmd := exec.Command("Rscript", path.Join(wd, "script", "fast_report.R"), name)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
				return err
			}
		}
	}
	cmd := exec.Command("Rscript", path.Join(wd, "script", "fpkm_report.R"))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	cmd = exec.Command("Rscript", path.Join(wd, "script", "pca_count.R"))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	cmd = exec.Command("Rscript", path.Join(wd, "script", "freq_count.R"))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	cmd = exec.Command("Rscript", path.Join(wd, "script", "heatmap_count.R"))
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		e.logger.Error("bash run ", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	return nil
}
