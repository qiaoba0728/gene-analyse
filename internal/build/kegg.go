package build

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
)

type KeggResult struct {
	KO          string
	Description string
}
type keggPlugin struct {
}

func NewKeggPlugin() types.Plugin {
	return &keggPlugin{}
}
func (g *keggPlugin) check() error {
	err := utils.WriteFile("nomode_kegg.R", scripts.NOMODE_KEGG)
	if err != nil {
		return err
	}
	return nil
}
func (g *keggPlugin) Name() string {
	return "keggPlugin"
}
func (g *keggPlugin) Build(ctx context.Context) error {
	//list := g.getList()
	//log.Println("Build:",list,len(list))
	return nil
}

//func (g *keggPlugin)getList() map[string][]string{
//	//kos := make([]string,0)
//	ks := make([]string,0)
//	resp := make(map[string][]string)
//	for _,v := range g.list {
//		temp := kegg.GeneKO{
//			GeneId:  v,
//		}
//		if err := temp.FindById(context.Background());err != nil {
//			log.Println("find error",err.Error(),v)
//			continue
//		}
//		ks = append(ks,temp.KId)
//	}
//	ks = utils.DeleteRepeat(ks)
//	log.Println(ks,len(ks))
//	for _,v := range ks {
//		temp := kegg.KO{
//			KId:  v,
//		}
//		if results,err := temp.FindByKId(context.Background());err != nil {
//			log.Println("find by kid error",err.Error(),v)
//			continue
//		}else {
//			for _,val := range results {
//				if _,ok := resp[val.KOId];ok {
//					resp[val.KOId] = append(resp[val.KOId],val.KId)
//				}else {
//					resp[val.KOId] = []string{val.KId}
//				}
//			}
//		}
//	}
//	for k,_ := range resp {
//		resp[k] = utils.DeleteRepeat(resp[k])
//	}
//	return resp
//}
