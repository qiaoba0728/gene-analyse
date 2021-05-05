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
	group  *conf.Group `json:"conf"` //配置文件
	logger *zap.Logger
}

func NewGenePlugin(group *conf.Group, logger *zap.Logger) types.Plugin {
	return &genePlugin{
		group:  group,
		logger: logger,
	}
}
func (g *genePlugin) check() error {
	cmd := exec.Command("rm", "-rf", g.group.Output)
	if err := cmd.Run(); err != nil {
		log.Println("cmd run rm err", err)
		g.logger.Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	cmd = exec.Command("mkdir", "-p", g.group.Output)
	if err := cmd.Run(); err != nil {
		g.logger.Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	err := utils.WriteFile("diff_matrix.R", scripts.DIFF_MATRIX)
	if err != nil {
		g.logger.Error("cmd run fail", zap.Error(err), zap.String("cmd", cmd.String()))
		return err
	}
	//err = utils.WriteFile("nomode_kegg.R", scripts.NOMODE_KEGG)
	//if err != nil {
	//	g.logger.Error("cmd run fail",zap.Error(err),zap.String("cmd",cmd.String()))
	//	return err
	//}
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
	if err := g.check(); err != nil {
		return err
	}
	wd, _ := os.Getwd()
	done := utils.Wait(func() error {
		if g.group.StartRepeated != "1" && g.group.EndRepeated != "1" {
			cmd := exec.Command("Rscript", path.Join(wd, "script", "diff_matrix.R"), path.Join(types.EXPRESSION_OUT, "gene_count.csv"),
				g.group.Start, g.group.End,
				g.group.StartRepeated,
				g.group.EndRepeated,
				g.group.Output,
				g.group.Name)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				log.Println("diff error", err.Error())
				return err
			}
		} else {
			//gfold diff -s1 sample1 -s2 sample2 -suf .read_cnt -o result.diff
			cmd := exec.Command("gfold", "diff",
				"-s1", fmt.Sprintf("%s/%s", types.SINGLE_OUT, g.group.Start),
				"-s2", fmt.Sprintf("%s/%s", types.SINGLE_OUT, g.group.End),
				"-o", fmt.Sprintf("%s/%s", g.group.Output, g.group.Name))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				g.logger.Error("run fail", zap.String("cmd", cmd.String()), zap.Error(err))
			}
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
