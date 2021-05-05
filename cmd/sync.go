package cmd

import (
	//"github.com/qiaoba0728/gene-analyse/internal/common"
	//"github.com/qiaoba0728/gene-analyse/internal/dao/kegg"
	"github.com/jszwec/csvutil"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"github.com/spf13/cobra"
	"io/ioutil"
	"log"
	"os"
	"path"
	"sync/atomic"
	"time"
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
			//result gene -> ko
			err := sync(result, koInput)
			if err != nil {
				log.Println("sync fail", err.Error())
				return
			}
			//batchData := make([]*kegg.GeneKO,0)
		}
	},
}

//ko -- name
type term2name struct {
	Term string `csv:"term"`
	Name string `csv:"name"`
}

//ko -- gene
type k2ko struct {
	Term string `csv:"term"`
	Gene string `csv:"gene"`
}

//func sync(input map[string]string,kosInput map[string]struct{}) error {
//	var(
//		//resp []*types.Kegg
//		err error
//		count uint32
//		wg sync2.WaitGroup
//		desc map[string]string
//		kom map[string][]string    // ko -> []koid
//		mu sync2.Mutex
//	)
//	desc = make(map[string]string)
//	kom = make(map[string][]string)
//	log.Println("start get kegg db",len(kosInput))
//	pool, _ := ants.NewPool(10)
//	client := utils.NewClient()
//	for k,_ := range kosInput {
//		wg.Add(1)
//		kId := k
//		if err = pool.Submit(func() {
//			var (
//				results []*types.Kegg
//				err error
//			)
//			defer func() {
//				wg.Done()
//				atomic.AddUint32(&count,1)
//				log.Println("get kegg count:",kId,atomic.LoadUint32(&count),len(kosInput),len(kom))
//			}()
//			for {
//				results,err = utils.GetKeggInfo(kId,client)
//				if err != nil {
//					time.Sleep(1 * time.Second)
//					log.Println("@hx restart get kegg:",kId,err.Error())
//					continue
//				}
//				break
//			}
//			if len(results) == 0 {
//				log.Println("get kegg empty:",kId)
//				return
//			}
//			mu.Lock()
//			kos := make([]string,len(results))
//			// k -- [ko1,ko2]
//			for key,val := range results {
//				if _,ok := desc[val.KOId];!ok {
//					desc[val.KOId] = val.Description
//				}
//				kos[key] = val.KOId
//				log.Println("key:",key,"val:",val)
//			}
//			if _,ok := kom[k];!ok {
//				kom[k] = kos
//			}
//			mu.Unlock()
//		});err != nil {
//			log.Println("pool submit error")
//		}
//	}
//	wg.Wait()
//	for pool.Running() != 0 {
//		log.Println("waiting sync finished",pool.Running())
//		time.Sleep(3 * time.Second)
//	}
//	log.Println("sync finished")
//	log.Println("term2name:",len(desc),len(kosInput),len(kom))
//	data := make([]*term2name,0)
//	for k,v := range desc{
//		temp := &term2name{
//			Term: k,
//			Name: v,
//		}
//		data = append(data,temp)
//	}
//	b, err := csvutil.Marshal(data)
//	if err != nil {
//		log.Println("csv marshal error")
//		return err
//	}
//	wd,_ := os.Getwd()
//	_ = ioutil.WriteFile(path.Join(wd,"term2name.csv"),b,0644)
//	log.Println("write term2name success")
//	geneko := make([]*k2ko,0)
//	for k,v := range input {
//		if vals,ok := kom[v];ok {
//			for _,val := range vals {
//				temp := &k2ko{
//					Term:   val,
//					Gene: k,
//				}
//				geneko = append(geneko,temp)
//			}
//		}else {
//			log.Println("not existed",v,k)
//		}
//	}
//	log.Println("term2gene:",len(geneko),len(input))
//	b, err = csvutil.Marshal(geneko)
//	if err != nil {
//		log.Println("csv marshal error")
//		return err
//	}
//	_ = ioutil.WriteFile(path.Join(wd,"term2gene.csv"),b,0644)
//	log.Println("write term2gene success")
//	return nil
//}
func sync(input map[string]string, kosInput map[string]struct{}) error {
	var (
		//resp []*types.Kegg
		err   error
		count uint32
		//wg sync2.WaitGroup
		desc map[string]string
		kom  map[string][]string // ko -> []koid
		//mu sync2.Mutex
	)
	desc = make(map[string]string)
	kom = make(map[string][]string)
	log.Println("start get kegg db", len(kosInput))
	//pool, _ := ants.NewPool(100,ants.WithMaxBlockingTasks(100))
	client := utils.NewClient()
	for k, _ := range kosInput {
		//wg.Add(1)
		kId := k
		var (
			results []*types.Kegg
			err     error
		)
		for {
			results, err = utils.GetKeggInfo(kId, client)
			if err != nil {
				time.Sleep(1 * time.Second)
				log.Println("restart get kegg:", kId)
				continue
			}
			break
		}
		if len(results) == 0 {
			log.Println("get kegg empty:", kId)
			continue
		}
		kos := make([]string, len(results))
		// k -- [ko1,ko2]
		for key, val := range results {
			desc[val.KOId] = val.Description
			kos[key] = val.KOId
		}
		kom[k] = kos
		atomic.AddUint32(&count, 1)
		log.Println("get kegg count:", atomic.LoadUint32(&count), len(kosInput))
	}
	//wg.Wait()
	log.Println("sync finished")
	log.Println("term2name:", len(desc), len(kosInput))
	data := make([]*term2name, 0)
	for k, v := range desc {
		temp := &term2name{
			Term: k,
			Name: v,
		}
		data = append(data, temp)
	}
	b, err := csvutil.Marshal(data)
	if err != nil {
		log.Println("csv marshal error")
		return err
	}
	wd, _ := os.Getwd()
	_ = ioutil.WriteFile(path.Join(wd, "term2name.csv"), b, 0644)
	log.Println("write term2name success")
	geneko := make([]*k2ko, 0)
	for k, v := range input {
		if vals, ok := kom[v]; ok {
			for _, val := range vals {
				temp := &k2ko{
					Term: val,
					Gene: k,
				}
				geneko = append(geneko, temp)
			}
		} else {
			log.Println("not existed", v, k)
		}
	}
	log.Println("term2gene:", len(geneko), len(input))
	b, err = csvutil.Marshal(geneko)
	if err != nil {
		log.Println("csv marshal error")
		return err
	}
	_ = ioutil.WriteFile(path.Join(wd, "term2gene.csv"), b, 0644)
	log.Println("write term2gene success")
	return nil
}
