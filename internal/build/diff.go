package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"
)

type genePlugin struct {
	config *conf.Config `json:"conf"` //配置文件
	logger *zap.Logger
}

func NewGenePlugin(config *conf.Config, logger *zap.Logger) types.Plugin {
	return &genePlugin{
		config: config,
		logger: logger,
	}
}
func (g *genePlugin) check(group *conf.Group) error {
	if !utils.IsExist(group.Output) {
		cmd := exec.Command("mkdir", "-p", group.Output)
		if err := cmd.Run(); err != nil {
			g.logger.Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
			return err
		}
	}
	err := utils.WriteFile("diff_matrix.R", scripts.DIFF_MATRIX)
	if err != nil {
		g.logger.Error("cmd run fail", zap.Error(err))
		return err
	}
	err = utils.WriteFile("nomode_kegg.R", scripts.NOMODO_KEGG)
	if err != nil {
		g.logger.Error("cmd run fail", zap.Error(err))
		return err
	}
	err = utils.WriteFile("nomode_go.R", scripts.NOMODE_GO)
	if err != nil {
		g.logger.Error("cmd run fail", zap.Error(err))
		return err
	}
	err = utils.WriteFile("insertsect_go.R", scripts.INSERTSECT)
	if err != nil {
		g.logger.Error("cmd run fail", zap.Error(err))
		return err
	}
	return nil
}
func (g *genePlugin) Name() string {
	return "genePlugin"
}
func (g *genePlugin) getGTF() (string, error) {
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
func (g *genePlugin) Build(ctx context.Context) error {
	wd, _ := os.Getwd()
	inputs := make([]string, 0)
	for _, v := range g.config.Group {
		if err := g.check(v); err != nil {
			return err
		}
		if v.StartRepeated != "1" && v.EndRepeated != "1" {
			cmd := exec.Command("Rscript", path.Join(wd, "script", "diff_matrix.R"), path.Join(types.EXPRESSION_OUT, "gene_count.csv"),
				v.Start, v.End,
				v.StartRepeated,
				v.EndRepeated,
				v.Output,
				v.Name)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				log.Println("diff error", err.Error())
				return err
			}
			inputs = append(inputs, path.Join(v.Output, fmt.Sprintf("diffexpr-%s-0.05.txt", v.Name)))
			cmd = exec.Command("Rscript", path.Join(wd, "script", "nomode_go.R"),
				path.Join(v.Output, fmt.Sprintf("diffexpr-%s-0.05.txt", v.Name)),
				fmt.Sprintf("%s/%s", v.Output, v.Name))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				log.Println("diff error", err.Error())
				return err
			}
			cmd = exec.Command("Rscript", path.Join(wd, "script", "nomode_kegg.R"),
				path.Join(v.Output, fmt.Sprintf("diffexpr-%s-0.05.txt", v.Name)),
				fmt.Sprintf("%s/%s", v.Output, v.Name))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				g.logger.Error("diff error", zap.String("cmd", cmd.String()), zap.Error(err))
				return err
			}
		} else {
			//gfold diff -s1 sample1 -s2 sample2 -suf .read_cnt -o result.diff
			cmd := exec.Command("gfold", "diff",
				"-s1", fmt.Sprintf("%s/%s", types.SINGLE_OUT, v.Start),
				"-s2", fmt.Sprintf("%s/%s", types.SINGLE_OUT, v.End),
				"-o", fmt.Sprintf("%s/%s", v.Output, v.Name))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				g.logger.Error("run fail", zap.String("cmd", cmd.String()), zap.Error(err))
			}
		}
		g.logger.Info("diff finished......")
	}
	if g.config.DiffGroup != nil && len(inputs) > 0 {
		cmd := exec.Command("Rscript", path.Join(wd, "script", "insertsect_go.R"), strings.Join(inputs, ","), fmt.Sprintf("%s/all", g.config.DiffGroup.Output))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("run fail", zap.String("cmd", cmd.String()), zap.Error(err))
			return err
		}
		cmd = exec.Command("Rscript", path.Join(wd, "script", "nomode_go.R"),
			fmt.Sprintf("%s/all_merge.txt", g.config.DiffGroup.Output),
			fmt.Sprintf("%s/%s", g.config.DiffGroup.Output, g.config.DiffGroup.Name))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			log.Println("diff error", err.Error())
			return err
		}
		cmd = exec.Command("Rscript", path.Join(wd, "script", "nomode_kegg.R"),
			fmt.Sprintf("%s/all_merge.txt", g.config.DiffGroup.Output),
			fmt.Sprintf("%s/%s", g.config.DiffGroup.Output, g.config.DiffGroup.Name))
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			g.logger.Error("diff error", zap.String("cmd", cmd.String()), zap.Error(err))
			return err
		}
	}
	return nil
}
