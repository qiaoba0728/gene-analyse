package cmd

import (
	"github.com/spf13/cobra"
)

var (
	list string
)

func init() {
	mainCmd.Flags().StringVar(&list, "list", "./list.txt", "config file (default is ./list.txt)")
}

var keggCmd = &cobra.Command{
	Use:   "kegg",
	Short: "kegg gene-analyse plugin",
	Long:  "run gene-analyse plugin build kegg data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {

	},
}
