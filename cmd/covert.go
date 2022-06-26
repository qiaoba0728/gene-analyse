package cmd

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"go.uber.org/zap"
)

var (
	input    string
	contract string //合同名称
	name     string //客户姓名
	order    string //合同编号
	address  string //客户单位
)
var covertCmd = &cobra.Command{
	Use:   "covert",
	Short: "gene-analyse plugin",
	Long:  "covert html to owner",
	PreRun: func(cmd *cobra.Command, args []string) {
	},
	Run: func(cmd *cobra.Command, args []string) {
		l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
		plugin := build.NewCovertPlugin(input, name, order, contract, address, l.Named("covert"))
		err := plugin.Build(context.Background())
		if err != nil {
			l.Error("build fail", zap.Error(err))
		}
	},
}

func init() {
	covertCmd.Flags().StringVarP(&input, "input", "i", "./test.zip", "use zip file")
	covertCmd.Flags().StringVarP(&name, "name", "n", "刘静", "use zip file")
	covertCmd.Flags().StringVarP(&contract, "contract", "c", "测序分析样品检测报告", "合同名称")
	covertCmd.Flags().StringVarP(&order, "order", "o", "Q202201016", "合同编号")
	covertCmd.Flags().StringVarP(&address, "address", "a", "华中科技大学同济医院", "客户单位")
}
