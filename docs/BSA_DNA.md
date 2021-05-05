## BSA-DNA/RNA 下机数据处理（不包括数据处理）
测序公司测序完成后会返回我们需要的下机数据

说明:基于docker构建整体的工具链（镜像仓库）
### 目录设置
- 新建工作目录test(可根据测序样品设置自己设置),例如 /home/hexing/test
- 在工作目录下新建input目录（用于存放测序数据(/home/hexing/test/input/raw_data/*.fast.q.gz),基因注释文件(home/hexing/test/input/references/*.fasta)）


### 下机数据处理

目前镜像版本
- DNA处理镜像地址 dockerhub.qingcloud.com/iot_demo/r-base:v0.5
- RNA处理镜像地址 dockerhub.qingcloud.com/iot_demo/r-base:v0.7
- 图形参数处理镜像地址 dockerhub.qingcloud.com/iot_demo/r-base:plot

查看镜像

```
sudo docker images
```

运行镜像生成容器任务

```
sudo docker run -it -d --name test -v /home/hexing/test:/data dockerhub.qingcloud.com/iot_demo/r-base:v0.5
```

或者

```
sudo docker run -it -d --name test -v /home/hexing/test:/data {imageId}
```

可以查看任务进度：

```
sudo docker logs -f {containerId}
```

如果想修改图形生成的参数
```
sudo docker run -it -d --name test -v /home/hexing/test:/data --env type=rna --env step=500000 --env window=1000000 --env thread=16 {imageId}
```

### 结果生成

目录结构
```bash
├── input
│   ├── index
│   ├── raw_data
│   │   ├── P-1.R1.fastq.gz
│   │   ├── P-1.R2.fastq.gz
│   │   ├── P-2.R1.fastq.gz
│   │   ├── P-2.R2.fastq.gz
│   │   ├── P-3.R1.fastq.gz
│   │   ├── P-3.R2.fastq.gz
│   │   ├── W-1.R1.fastq.gz
│   │   ├── W-1.R2.fastq.gz
│   │   ├── W-2.R1.fastq.gz
│   │   ├── W-2.R2.fastq.gz
│   │   ├── W-3.R1.fastq.gz
│   │   └── W-3.R2.fastq.gz
│   └── references
│       ├── Qing.chr.v1.fa
├── log
└── output
```

#### input（输入文件）
- index（参考基因组生成的索引）
- raw_data（输入下机测序数据fastq）
- references（参考基因组fasta文件）

#### output（输出文件）
包括生成的图片和相关的数据文件
#### log（运行时日志）
目前生成一些运行时日志
