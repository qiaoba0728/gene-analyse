package cmd

import "github.com/spf13/cobra"

var goCmd = &cobra.Command{
	Use:   "go",
	Short: "go gene-analyse plugin",
	Long:  "run gene-analyse plugin build go data",
	PreRun: func(cmd *cobra.Command, args []string) {
		//conf.InitConfig(cfgFile)
	},
	Run: func(cmd *cobra.Command, args []string) {
	},
}
