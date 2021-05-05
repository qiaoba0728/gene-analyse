package build

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
)

type goPlugin struct {
	list []string
}

func NewGoPlugin(group *conf.Group) types.Plugin {
	return &genePlugin{
		group: group,
	}
}
func (g *goPlugin) check() error {
	err := utils.WriteFile("diff_matrix.R", scripts.DIFF_MATRIX)
	if err != nil {
		return err
	}
	return nil
}
func (g *goPlugin) Name() string {
	return "goPlugin"
}
func (g *goPlugin) Build(ctx context.Context) error {
	return nil
}
