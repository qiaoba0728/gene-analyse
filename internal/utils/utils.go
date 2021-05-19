package utils

import (
	"archive/zip"
	"bufio"
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"go.uber.org/zap"
	"io"
	"io/ioutil"
	"log"
	"net"
	"net/http"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"strings"
	"time"
)

const (
	http_content_type = "application/json;charset=utf-8"
)

func Wait(f func() error) <-chan error {
	done := make(chan error)
	go func() {
		done <- f()
	}()
	return done
}

// IsExist checks whether a file or directory exists.
// It returns false when the file or directory does not exist.
func IsExist(f string) bool {
	_, err := os.Stat(f)
	return err == nil || os.IsExist(err)
}

func Files(f string) int {
	files, err := ioutil.ReadDir(f)
	if err != nil {
		panic(err)
	}
	return len(files)
}
func SufFiles(f, suf string) bool {
	files, err := ioutil.ReadDir(f)
	if err != nil {
		panic(err)
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), suf) {
			return true
		}
	}
	return false
}
func WriteFile(name, data string) error {
	dir, _ := os.Getwd()
	d := path.Join(dir, "script")
	if b := IsExist(d); !b {
		cmd := exec.Command("mkdir", "-p", d)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			log.Println("=================> mkdir err:", err.Error())
		}
	}
	return ioutil.WriteFile(path.Join(d, name), []byte(data), 0644)
}

func ReadKEGG(file string) (map[string]string, error) {
	result := make(map[string]string)
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	rd := bufio.NewReader(f)
	for {
		line, err := rd.ReadString('\n') //以'\n'为结束符读入一行
		if err != nil || io.EOF == err {
			break
		}
		split := strings.Split(line, "\t")
		if len(split) != 2 {
			continue
		}
		key := strings.Replace(split[0], " ", "", -1)
		key = strings.Replace(key, "\n", "", -1)
		value := strings.Replace(split[1], " ", "", -1)
		value = strings.Replace(value, "\n", "", -1)
		result[key] = value
	}
	return result, nil
}
func RemoveFirstBlank(param string) string {
	index := 0
	for _, v := range param {
		if v == ' ' {
			index++
			continue
		} else {
			break
		}
	}
	return strings.Replace(param[index:], "\n", "", -1)
}

//func GetKeggInfo(id string, client *http.Client) ([]*types.Kegg, error) {
//	result := make([]*types.Kegg, 0)
//	url := fmt.Sprintf("http://rest.kegg.jp/get/%s", id)
//	log.Println("post:", url)
//	resp, err := client.Post(url, http_content_type, nil)
//	if err != nil {
//		log.Println("http post error:", err.Error())
//		return result, err
//	}
//	defer resp.Body.Close()
//	rd := bufio.NewReader(resp.Body)
//	for {
//		line, err := rd.ReadString('\n') //以'\n'为结束符读入一行
//		if err != nil || io.EOF == err {
//			break
//		}
//		reg := regexp.MustCompile(`\sko[0-9]{5}\s`)
//		if reg == nil {
//			continue
//		}
//		match := reg.FindAllStringSubmatch(line, -1)
//		if len(match) != 0 {
//			index := match[0][0]
//			n := strings.Index(line, index)
//			res := line[n+len(index):]
//			//log.Println(n,index,res)
//			temp := &types.Kegg{
//				KOId:        strings.Replace(index, " ", "", -1),
//				Description: RemoveFirstBlank(res),
//			}
//			result = append(result, temp)
//		}
//	}
//	return result, nil
//}
func ReadCSV(file string) []string {
	result := make([]string, 0)
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	rd := bufio.NewReader(f)
	for {
		line, err := rd.ReadString('\n') //以'\n'为结束符读入一行
		if err != nil || io.EOF == err {
			break
		}
		result = append(result, line)
	}
	return result
}

func NewClient() *http.Client {
	producer := &http.Client{
		Transport: &http.Transport{
			Proxy: http.ProxyFromEnvironment,
			DialContext: (&net.Dialer{
				Timeout:   3 * time.Second,
				KeepAlive: 30 * time.Second,
			}).DialContext,
			MaxIdleConns:        types.MaxIdleConns,
			MaxIdleConnsPerHost: types.MaxIdleConnsPerHost,
			IdleConnTimeout:     time.Duration(types.IdleConnTimeout) * time.Second,
		},
	}
	return producer
}
func DeleteRepeat(list []string) []string {
	mapdata := make(map[string]interface{})
	if len(list) <= 0 {
		return nil
	}
	for _, v := range list {
		mapdata[v] = struct{}{}
	}
	var datas []string
	for k, _ := range mapdata {
		if k == "" {
			continue
		}
		datas = append(datas, k)
	}
	return datas
}

func CompressZip(src, zipName string, logger *zap.Logger) (err error) {
	dir, err := ioutil.ReadDir(src)
	if err != nil {
		logger.Error("ioutil.ReadDir err:", zap.Error(err))
		return err
	}
	if len(dir) == 0 {
		logger.Info(src + " is empty dir!")
		return nil
	}
	// 预防：旧文件无法覆盖
	os.RemoveAll(zipName)

	// 创建：zip文件
	zipfile, _ := os.Create(zipName)
	defer zipfile.Close()

	// 打开：zip文件
	archive := zip.NewWriter(zipfile)
	defer archive.Close()

	// 遍历路径信息
	filepath.Walk(src, func(path string, info os.FileInfo, _ error) error {
		// 如果是源路径，提前进行下一个遍历
		if path == src {
			return nil
		}
		// 获取：文件头信息
		header, _ := zip.FileInfoHeader(info)
		header.Name = strings.TrimPrefix(path, src+`\`)
		// 判断：文件是不是文件夹
		if info.IsDir() {
			header.Name += `/`
		} else {
			// 设置：zip的文件压缩算法
			header.Method = zip.Deflate
		}
		// 创建：压缩包头部信息
		writer, _ := archive.CreateHeader(header)
		if !info.IsDir() {
			file, _ := os.Open(path)
			defer file.Close()
			io.Copy(writer, file)
		}
		return nil
	})
	return nil
}
func ZipFiles(filename string, files []string, oldform, newform string) error {
	newZipFile, err := os.Create(filename)
	if err != nil {
		return err
	}

	defer newZipFile.Close()
	zipWriter := zip.NewWriter(newZipFile)
	defer zipWriter.Close()

	// 把files添加到zip中
	for _, file := range files {
		zipfile, err := os.Open(file)
		if err != nil {
			return err
		}
		//defer zipfile.Close()
		info, err := zipfile.Stat()
		if err != nil {
			zipfile.Close()
			return err
		}
		header, err := zip.FileInfoHeader(info)
		if err != nil {
			return err
		}
		header.Name = strings.Replace(file, oldform, newform, -1)
		header.Method = zip.Deflate
		writer, err := zipWriter.CreateHeader(header)
		if err != nil {
			zipfile.Close()
			return err
		}

		if _, err = io.Copy(writer, zipfile); err != nil {
			zipfile.Close()
			return err
		}
		zipfile.Close()
	}
	return nil
}
func UnZip(file string, target string) error {
	r, err := zip.OpenReader(file)
	if err != nil {
		return err
	}
	for _, k := range r.Reader.File {
		if k.FileInfo().IsDir() {
			file := fmt.Sprintf("%s/%s", target, k.Name)
			err := os.MkdirAll(file, 0644)
			if err != nil {
				fmt.Println(err)
			}
			continue
		}
		r, err := k.Open()
		if err != nil {
			continue
		}
		file := fmt.Sprintf("%s/%s", target, k.Name)
		NewFile, err := os.Create(file)
		if err != nil {
			fmt.Println(err)
			continue
		}
		io.Copy(NewFile, r)
		NewFile.Close()
		r.Close()
	}
	return nil
}
