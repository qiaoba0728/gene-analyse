package utils

import (
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qianlnk/pgbar"
	"github.com/qiaoba0728/gene-analyse/internal/scripts"
	"github.com/stretchr/testify/assert"
	"go.uber.org/zap"
	"sync"
	"testing"
	"time"
)

func TestWait(t *testing.T) {
	var wg sync.WaitGroup
	wg.Add(4)
	bar := pgbar.New("多线程进度条")
	u := bar.NewBar("test", 4)
	pool, err := ants.NewPool(2)
	assert.Nil(t, err)
	err = pool.Submit(func() {
		fmt.Println("xxx")
		time.Sleep(3 * time.Second)
		u.Add(1)
		wg.Done()
	})
	assert.Nil(t, err)
	err = pool.Submit(func() {
		fmt.Println("xxx")
		u.Add(1)
		time.Sleep(3 * time.Second)
		wg.Done()
	})
	assert.Nil(t, err)
	err = pool.Submit(func() {
		fmt.Println("xxx")
		time.Sleep(3 * time.Second)
		u.Add(1)
		wg.Done()
	})
	assert.Nil(t, err)
	err = pool.Submit(func() {
		fmt.Println("xxx")
		time.Sleep(3 * time.Second)
		u.Add(1)
		wg.Done()
	})
	assert.Nil(t, err)
	fmt.Println("wait")
	wg.Wait()
}
func TestWriteFile(t *testing.T) {
	err := WriteFile("test.R", scripts.FEATURE_SCRIPT)
	assert.Nil(t, err)
}

func TestReadKEGG(t *testing.T) {
	result, err := ReadKEGG("E://data//feature//query.ko")
	assert.Nil(t, err)
	t.Log(len(result))
	for k, v := range result {
		t.Logf("k:%s,v:%s", k, v)
	}
}
func TestRemoveBlank(t *testing.T) {
	res := RemoveFirstBlank("      xxxx xxx")
	t.Log(res)
}
func TestGetKeggInfo(t *testing.T) {
	//client := NewClient()
	//result, err := GetKeggInfo("K12843", client)
	//assert.Nil(t, err)
	//t.Log(len(result))
	//for _, v := range result {
	//	t.Logf("%s,%s", v.KOId, v.Description)
	//}
}
func TestCompressZip(t *testing.T) {
	l, _ := NewZapLogger(&Opt{LogOutput: CONSOLE})
	err := CompressZip("../gene-analyse", "result.zip", l.Named("zip"))
	if err != nil {
		l.Error("zip error", zap.Error(err))
	}
}

func TestZipFiles(t *testing.T) {
	err := ZipFiles("result.zip", []string{"I://code/go/mod/github.com/qiaoba0728/gene-analyse/internal/types/kegg.go", "I://code/go/mod/github.com/qiaoba0728/gene-analyse/internal/types/types.go"}, "I://code/go/mod/github.com/qiaoba0728/gene-analyse/internal", ".")
	assert.Nil(t, err)
}
func TestUnZip(t *testing.T) {
	err := UnZip("D://data//P_R1.clean_fastqc.zip", "D://data//bsa//")
	if err != nil {
		fmt.Println(err)
	}
}
func TestDelSuf(t *testing.T) {
	err := DelSuf("D://data", ".txt")
	assert.Nil(t, err)
}

func TestRandFloats(t *testing.T) {
	f := RandFloats(1000, 2000, 10)
	for _, v := range f {
		t.Log(v)
	}
}

func TestBuildConfig(t *testing.T) {
	err := BuildConfig("D://gene_count.csv", 3, "D://build.json")
	assert.Nil(t, err)
}
