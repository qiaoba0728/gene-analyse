package common

import (
	_ "github.com/go-sql-driver/mysql"
	"github.com/go-xorm/xorm"
	"sync"
)

var con *xorm.Engine
var once sync.Once

func GetMysqlDB() *xorm.Engine {
	once.Do(func() {
		orm, err := xorm.NewEngine("mysql", "root:123456@tcp(192.168.150.200:3306)/kegg?charset=utf8")
		if err != nil {
			panic(err)
		}
		//orm.ShowSQL(true)
		con = orm
	})
	return con
}
