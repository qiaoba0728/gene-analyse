package build

import (
	"context"
	"errors"
	"fmt"
	"github.com/panjf2000/ants"
	"github.com/qianlnk/pgbar"
	"github.com/qiaoba0728/gene-analyse/internal/common"
	"github.com/qiaoba0728/gene-analyse/internal/conf"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"
	"sync"
)

type buildExPlugin struct {
	//bar *pgbar.Pgbar
	logger *zap.Logger
}

func NewBuildExPlugin(logger *zap.Logger) types.Plugin {
	return &buildExPlugin{
		//bar:pgbar.New("RNA-seq"),
		logger: logger,
	}
}
func (u *buildExPlugin) Name() string {
	return "buildExPlugin"
}
func (u *buildExPlugin) check() {
	if b := utils.IsExist(types.GENOME_INDEX); !b {
		cmd := exec.Command("mkdir", "-p", types.GENOME_INDEX)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create index dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.GENOME_INDEX)
		if err != nil {
			u.logger.Error("existed index", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("gene index", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.LOG); !b {
		cmd := exec.Command("mkdir", "-p", types.LOG)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create log dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.LOG)
		if err != nil {
			u.logger.Error("existed log dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("log", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.FASTP_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.FASTP_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create clean dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.FASTP_OUT)
		if err != nil {
			u.logger.Error("existed clean dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("fastp out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.REPORT_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.REPORT_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create report dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.REPORT_OUT)
		if err != nil {
			u.logger.Error("existed report dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("report out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.HISAT2_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.HISAT2_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create hisat2 dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.HISAT2_OUT)
		if err != nil {
			u.logger.Error("existed hisat2 dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("hisat2 out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.SORTED_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.SORTED_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create hisat2 sorted dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.SORTED_OUT)
		if err != nil {
			u.logger.Error("existed hisat2 sorted dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("hisat2 sorted out", zap.Strings("files", names))
		}
	}
	if b := utils.IsExist(types.EXPRESSION_OUT); !b {
		cmd := exec.Command("mkdir", "-p", types.EXPRESSION_OUT)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			u.logger.Error("create expression dir", zap.Error(err))
		}
	} else {
		files, err := ioutil.ReadDir(types.EXPRESSION_OUT)
		if err != nil {
			u.logger.Error("existed expression dir", zap.Error(err))
		} else {
			names := make([]string, 0)
			for _, v := range files {
				names = append(names, v.Name())
			}
			u.logger.Warn("expression out", zap.Strings("files", names))
		}
	}
}
func (u *buildExPlugin) Build(ctx context.Context) error {
	u.check()
	done := utils.Wait(func() error {
		var wg sync.WaitGroup
		wg.Add(2)
		go func() {
			defer wg.Done()
			if b := utils.IsExist(types.GENOME_INDEX); !b || (b && utils.Files(types.GENOME_INDEX) == 0) {
				if b := utils.IsExist(types.REFERENCES); b {
					if err := u.index(types.REFERENCES); err != nil {
						u.logger.Error("gene index build fail", zap.Error(err))
					}
				} else {
					u.logger.Error("gene index has existed")
				}
			} else {
				u.logger.Warn("gene source index has existed")
			}
		}()
		go func() {
			defer wg.Done()
			if b := utils.IsExist(types.INPUT); b {
				if b := utils.IsExist(types.FASTP_OUT); !b || (b && utils.Files(types.FASTP_OUT) == 0) {
					if err := u.fastq(types.INPUT); err != nil {
						u.logger.Error("create clean data fail", zap.Error(err))
					} else {
						err = common.BuildReport()
						if err != nil {
							u.logger.Error("build report data fail", zap.Error(err))
						} else {

						}
					}
				} else {
					u.logger.Info("fastq data has build")
					files, err := ioutil.ReadDir(types.FASTP_OUT)
					if err != nil {
						u.logger.Error("read clean data fail", zap.Error(err))
					} else {
						names := make([]string, 0)
						for _, v := range files {
							names = append(names, v.Name())
						}
						u.logger.Warn("clean data ", zap.Strings("files", names))
					}
				}
			} else {
				u.logger.Error("fastq input not existed")
			}
		}()
		wg.Wait()
		//clean -> sam
		if b := utils.IsExist(types.FASTP_OUT); b {
			if b := utils.IsExist(types.HISAT2_OUT); !b || (b && utils.Files(types.HISAT2_OUT) == 0) {
				if err := u.bowtie2(types.FASTP_OUT); err != nil {
					u.logger.Error("clean -> hisat data fail", zap.Error(err))
					return err
				}
			} else {
				u.logger.Info("hisat2 data has build")
				files, err := ioutil.ReadDir(types.FASTP_OUT)
				if err != nil {
					u.logger.Error("read hisat2 data fail", zap.Error(err))
				} else {
					names := make([]string, 0)
					for _, v := range files {
						names = append(names, v.Name())
					}
					u.logger.Warn("hisat2 data ", zap.Strings("files", names))
				}
			}
		} else {
			u.logger.Error("hisat2 input not existed")
		}
		//sam -> bam
		if b := utils.IsExist(types.HISAT2_OUT); b {
			if err := u.sort(types.HISAT2_OUT); err != nil {
				u.logger.Error("create hisat2 sorted data fail", zap.Error(err))
				return err
			} else {
				u.logger.Info("create hisat2 sorted data success")
			}
		} else {
			u.logger.Error("hisat2 sorted input not existed")
		}
		//bam -> count
		//if b := utils.IsExist(types.SORTED_OUT);b {
		//	if b := utils.IsExist(types.EXPRESSION_OUT);!b  ||(b && utils.Files(types.EXPRESSION_OUT) == 0){
		//		if err := u.feature(types.SORTED_OUT);err != nil {
		//			log.Println("=================> feature count err:",err.Error())
		//			return err
		//		}else {
		//			log.Println("=================> feature count file existed")
		//		}
		//	}else {
		//		log.Println("=================> feature count has build")
		//		files, err := ioutil.ReadDir(types.EXPRESSION_OUT)
		//		if err != nil {
		//			log.Println("=================> existed feature count:",err.Error())
		//		}else {
		//			names := make([]string,0)
		//			for _,v := range files {
		//				names = append(names,v.Name())
		//			}
		//			log.Println(types.EXPRESSION_OUT,names)
		//		}
		//	}
		//}else {
		//	log.Println("=================> feature count file existed")
		//}
		return nil
	})
	select {
	case err := <-done:
		return err
	case <-ctx.Done():
		return errors.New("time out")
	}
}
func (u *buildExPlugin) fastq(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	//fastp -i "${temp}".R1.fastq.gz -I "${temp}".R2.fastq.gz -o ${fastp_out}/"${temp}".R1.clean.fastq.gz -O ${fastp_out}/"${temp}".R2.clean.fastq.gz > ${log}/fastp.log
	//bar := u.bar.NewBar("fastq --> clean.fastq",len(files) / 2)
	pool, err := ants.NewPool(conf.GetBuildConfig().Fastp)
	if err != nil {
		return err
	}
	var wg sync.WaitGroup
	var flag bool
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".R1.fastq.gz") {
			flag = true
			name := v.Name()
			wg.Add(1)
			if err = pool.Submit(func() {
				defer func() {
					//bar.Add(1)
					wg.Done()
					u.logger.Info("build clean file success", zap.String("name", name))
				}()
				temp := strings.TrimSuffix(name, ".R1.fastq.gz")
				u.logger.Info("fastq -> clean", zap.String("source1", name), zap.String("source2", fmt.Sprintf("%s.R2.clean.fastq.gz", temp)))
				cmd := exec.Command("fastp", "-i",
					fmt.Sprintf("%s/%s.R1.fastq.gz", types.INPUT, temp),
					"-I", fmt.Sprintf("%s/%s.R2.fastq.gz", types.INPUT, temp),
					"-o", fmt.Sprintf("%s/%s.R1.clean.fastq.gz", types.FASTP_OUT, temp),
					"-O", fmt.Sprintf("%s/%s.R2.clean.fastq.gz", types.FASTP_OUT, temp),
					"-h", fmt.Sprintf("%s/%s.html", types.FASTP_OUT, temp),
					"-j", fmt.Sprintf("%s/%s.json", types.FASTP_OUT, temp),
					">", types.LOG)
				//cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					u.logger.Error("fastq run fail", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				u.logger.Error("submit task run fail", zap.Error(err))
			}
		}
	}
	if !flag {
		u.logger.Error("please check raw data(.R1.fastq.gz)")
		return errors.New("please check raw data")
	}
	u.logger.Info("fastq -> clean waiting")
	wg.Wait()
	u.logger.Info("fastq -> clean finished")
	return nil
}

// build fa -> Species
func (u *buildExPlugin) index(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".fa") || strings.HasSuffix(v.Name(), ".fasta") {
			name := v.Name()
			cmd := exec.Command("bowtie2-build", fmt.Sprintf("%s/%s", types.REFERENCES, name),
				types.GENOME_PREFIX,
				"1",
				">", fmt.Sprintf("%s/%s", types.LOG, "bowtie2-build.log"), "2",
				">", fmt.Sprintf("%s/%s", types.LOG, "bowtie2-build.err"))
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
			if err = cmd.Run(); err != nil {
				u.logger.Error("bowtie2-build task run fail", zap.Error(err), zap.String("name", name), zap.String("cmd", cmd.String()))
			} else {
				u.logger.Info("build gene index file success")
			}
		}
	}
	u.logger.Info("waiting index file success")
	return nil
}

//clean -> sam
func (u *buildExPlugin) bowtie2(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	//bar := u.bar.NewBar("clean.fastq --> sam",len(files) / 2 )
	//todo 内存占用较大（单个5G）
	pool, err := ants.NewPool(conf.GetBuildConfig().Hisat2)
	if err != nil {
		return err
	}
	var wg sync.WaitGroup
	var flag bool
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".R1.clean.fastq.gz") {
			flag = true
			name := v.Name()
			wg.Add(1)
			if err = pool.Submit(func() {
				defer func() {
					//bar.Add(1)
					wg.Done()
					log.Println(name, "build map file success", types.HISAT2_OUT, name)
				}()
				temp := strings.TrimSuffix(name, ".R1.clean.fastq.gz")
				cmd := exec.Command("bowtie2", "-p", "8",
					"--no-unal", "-x", types.GENOME_PREFIX,
					"-1", fmt.Sprintf("%s/%s.%s", types.FASTP_OUT, temp, "R1.clean.fastq.gz"),
					"-2", fmt.Sprintf("%s/%s.%s", types.FASTP_OUT, temp, "R2.clean.fastq.gz"),
					"-S", fmt.Sprintf("%s/%s.sam", types.HISAT2_OUT, temp),
					">", fmt.Sprintf("%s/%s.err", types.LOG, temp))
				//cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					u.logger.Error("hisat2 run fail", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				u.logger.Error("hisat2 submit task run fail", zap.Error(err))
			}
		}
	}
	if !flag {
		u.logger.Error("please check raw data(*R1.clean.fastq.gz)")
		return errors.New("please check clean data")
	}
	wg.Wait()
	u.logger.Info("clean -> hisat2 waiting")
	wg.Wait()
	u.logger.Info("clean -> hisat2 finished")
	return nil
}

// sam -> sort
func (u *buildExPlugin) sort(dir string) error {
	var (
		bar *pgbar.Bar
		wg  sync.WaitGroup
	)
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	pool, err := ants.NewPool(conf.GetBuildConfig().Samtools)
	if err != nil {
		return err
	}
	if !utils.SufFiles(types.SORTED_OUT, ".bam") {
		//bar := u.bar.NewBar("sam --> bam",len(files))
		for _, v := range files {
			if strings.HasSuffix(v.Name(), ".sam") {
				name := v.Name()
				wg.Add(1)
				if err = pool.Submit(func() {
					defer func() {
						//bar.Add(1)
						wg.Done()
						u.logger.Info("build sam file success", zap.String("name", name))
					}()
					temp := strings.TrimSuffix(name, ".sam")
					u.logger.Info("sam -> bam", zap.String("source", fmt.Sprintf("%s", name)), zap.String("target", fmt.Sprintf("%s.bam", temp)))
					cmd := exec.Command("samtools", "view",
						"-bS", fmt.Sprintf("%s/%s", types.HISAT2_OUT, name),
						"-o", fmt.Sprintf("%s/%s.bam", types.SORTED_OUT, temp))
					//cmd.Stdout = os.Stdout
					cmd.Stderr = os.Stderr
					u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
					if err = cmd.Run(); err != nil {
						u.logger.Error("samtools run fail", zap.Error(err), zap.String("cmd", cmd.String()))
					}
				}); err != nil {
					u.logger.Error("samtools submit task run fail", zap.Error(err))
				}
			}
		}
		u.logger.Info("sam -> bam waiting")
		wg.Wait()
	} else {
		u.logger.Info("create bam success")
	}
	pool, err = ants.NewPool(conf.GetBuildConfig().Samtools)
	if err != nil {
		return err
	}
	files, err = ioutil.ReadDir(types.SORTED_OUT)
	if err != nil {
		return err
	}
	//bar = u.bar.NewBar("bam --> bam sort",len(files))
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".bam") && !strings.HasSuffix(v.Name(), ".sorted.bam") {
			name := v.Name()
			wg.Add(1)
			if err = pool.Submit(func() {
				defer func() {
					bar.Add(1)
					wg.Done()
					u.logger.Info("build bam sorted file success", zap.String("name", name))
				}()
				temp := strings.TrimSuffix(name, ".bam")
				cmd := exec.Command("samtools", "sort", "-@",
					"4", fmt.Sprintf("%s/%s", types.SORTED_OUT, name),
					"-o", fmt.Sprintf("%s/%s.sorted.bam", types.SORTED_OUT, temp))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					u.logger.Error("samtools sorted run fail", zap.Error(err))
				}
				cmd = exec.Command("samtools", "index", "-@",
					"4", fmt.Sprintf("%s/%s.sorted.bam", types.SORTED_OUT, temp))
				//cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					u.logger.Error("samtools index run fail", zap.Error(err), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				u.logger.Error("samtools submit task run fail", zap.Error(err))
			}
		}
	}
	u.logger.Info("sam -> sort bam and index bam waiting")
	wg.Wait()
	u.logger.Info("sam -> sort bam and index bam finished")
	return nil
}
func (u *buildExPlugin) getGTF() (string, error) {
	files, err := ioutil.ReadDir(types.REFERENCES)
	if err != nil {
		return "", err
	}
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".gtf") {
			return path.Join(types.REFERENCES, v.Name()), nil
		}
	}
	return "", errors.New("gtf not find")
}
func (u *buildExPlugin) feature(dir string) error {
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return err
	}
	pool, err := ants.NewPool(10)
	if err != nil {
		return err
	}
	//bar := u.bar.NewBar("bam --> count",len(files) / 3)
	gtf, err := u.getGTF()
	if err != nil {
		return err
	}
	//bams := make([]string,0)
	var wg sync.WaitGroup
	for _, v := range files {
		if strings.HasSuffix(v.Name(), ".sorted.bam") {
			//bams = append(bams,fmt.Sprintf("%s/%s",types.SORTED_OUT,v.Name()))
			name := v.Name()
			wg.Add(1)
			if err = pool.Submit(func() {
				defer func() {
					//bar.Add(1)
					wg.Done()
					log.Println(name, "build count file finished")
				}()
				//featureCounts -p -Q 10 -s 0 -T "$cpu" -a "$gtf" -o $feature_out/"${file}".counts.txt "${file}"
				temp := strings.TrimSuffix(name, ".sorted.bam")
				cmd := exec.Command("featureCounts", "-p", "-Q",
					"10", "-s", "0", "-T", "4", "-a", gtf, "-o",
					fmt.Sprintf("%s/%s.counts.txt", types.EXPRESSION_OUT, temp),
					fmt.Sprintf("%s/%s", types.SORTED_OUT, name))
				cmd.Stdout = os.Stdout
				cmd.Stderr = os.Stderr
				u.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
				if err = cmd.Run(); err != nil {
					log.Println(name, err.Error(), zap.String("cmd", cmd.String()))
				}
			}); err != nil {
				u.logger.Error("pool run err", zap.Error(err))
			}
		}
	}
	wg.Wait()
	//featureCounts -p -Q 10 -s 0 -T 4 -a gtf -o test.counts.txt {file}
	//cmd := exec.Command("featureCounts","-p","-Q",
	//	"10","-s","0","-T","4","-a",gtf,"-o",
	//	fmt.Sprintf("%s/%s.all.counts.txt",types.EXPRESSION_OUT,"feature"),
	//	fmt.Sprintf("%s/%s.sorted.bam",types.SORTED_OUT,"*"))
	//cmd.Stdout = os.Stdout
	//cmd.Stderr = os.Stderr
	//if err =  cmd.Run();err != nil {
	//	return err
	//}
	return nil
}
