package kegg

import (
	"context"
	"errors"
	"github.com/go-xorm/xorm"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"time"
)

type GeneKO struct {
	Id      int64     `xorm:"pk autoincr"`
	GeneId  string    `xorm:"gene_id"`
	KId     string    `xorm:"kid"`
	Created time.Time `xorm:"created"`
	Updated time.Time `xorm:"updated"`
}

func (g *GeneKO) Create(ctx context.Context) error {
	_, err := common.GetMysqlDB().InsertOne(g)
	return err
}
func (g *GeneKO) BatchCreate(ctx context.Context, beans ...interface{}) error {
	_, err := common.GetMysqlDB().Insert(beans...)
	return err
}
func (g *GeneKO) FindById(ctx context.Context) error {
	var (
		err   error
		db    *xorm.Engine
		query *xorm.Session
		b     bool
	)
	db = common.GetMysqlDB()
	query = db.Where("gene_id=?", g.GeneId)
	b, err = query.Get(g)
	if !b {
		return errors.New("record is not existed")
	}
	return err
}

type KO struct {
	Id          int64     `xorm:"pk autoincr"`
	KOId        string    `xorm:"ko_id"`
	KId         string    `xorm:"kid"`
	Description string    `xorm:"description"`
	Created     time.Time `xorm:"created"`
	Updated     time.Time `xorm:"updated"`
}

func (k *KO) Create(ctx context.Context) error {
	_, err := common.GetMysqlDB().Insert(k)
	return err
}
func (k *KO) BatchCreate(ctx context.Context, beans ...interface{}) error {
	_, err := common.GetMysqlDB().Insert(beans...)
	return err
}
func (k *KO) FindByKId(ctx context.Context) ([]*KO, error) {
	var (
		err    error
		db     *xorm.Engine
		query  *xorm.Session
		result []*KO
	)
	db = common.GetMysqlDB()
	query = db.Where("kid=?", k.KId)
	err = query.Find(&result)
	return result, err
}
func (k *KO) FindByKO(ctx context.Context) ([]*KO, error) {
	var (
		err    error
		db     *xorm.Engine
		query  *xorm.Session
		result []*KO
	)
	db = common.GetMysqlDB()
	query = db.Where("ko_id=?", k.KOId)
	err = query.Find(&result)
	return result, err
}

type KODesc struct {
	KOId        string    `xorm:"pk 'ko_id''"`
	Description string    `xorm:"description"`
	Created     time.Time `xorm:"created"`
	Updated     time.Time `xorm:"updated"`
}

func (k *KODesc) Create(ctx context.Context) error {
	_, err := common.GetMysqlDB().Insert(k)
	return err
}
