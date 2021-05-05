package kegg

import (
	"context"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/stretchr/testify/assert"
	"testing"
	"time"
)

func TestGetMysqlDB(t *testing.T) {
	err := common.GetMysqlDB().CreateTables(&GeneKO{}, &KO{})
	assert.Nil(t, err)
}
func TestKO_BatchCreate(t *testing.T) {
	k := &GeneKO{
		GeneId:  "xxx",
		KId:     "yyy",
		Created: time.Time{},
		Updated: time.Time{},
	}
	result := make([]*GeneKO, 0)
	result = append(result, k)
	err := k.BatchCreate(context.Background(), result)
	assert.Nil(t, err)
}
