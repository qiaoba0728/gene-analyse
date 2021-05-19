package cmd

import (
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"log"
)

var (
	file string
)

func init() {
	syncCmd.Flags().StringVar(&file, "db", "./db.txt", "config file (default is ./db.txt)")
}

var syncCmd = &cobra.Command{
	Use:   "sync",
	Short: "gene sync plugin",
	Long:  "run sync data plugin",
	PreRun: func(cmd *cobra.Command, args []string) {
		//err := common.GetMysqlDB().CreateTables(&kegg.GeneKO{},&kegg.KO{})
		//if err != nil {
		//	panic(err)
		//}
	},
	Run: func(cmd *cobra.Command, args []string) {
		if result, err := utils.ReadKEGG(file); err != nil {
			log.Println("read fail", err.Error())
			return
		} else {
			koInput := make(map[string]struct{})
			for _, v := range result {
				if _, ok := koInput[v]; !ok {
					koInput[v] = struct{}{}
				}
			}
			log.Println("result:", len(result), "koInput:", len(koInput))
			//batchData := make([]*kegg.GeneKO,0)
		}
	},
}
