## RNA-seq 下机数据处理（不包括数据处理）
测序公司测序完成后会返回我们需要的下机数据
- fastp(质控和过滤)
- hisat2-build(参考基因组建立索引)
- hisat2(mapping到参考基因组)
- samtools(比对结果压缩排序以及构建索引)

说明:基于docker构建整体的工具链（镜像仓库）
### 目录设置
- 新建工作目录test(可根据测序样品设置自己设置),例如 /home/hexing/test
- 在工作目录下新建input目录（用于存放测序数据(/home/hexing/test/input/raw_data/*.fast.q.gz),基因注释文件(/home/hexing/test/input/references/*.gtf,/home/hexing/test/input/references/*.fa)）

示例：
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
│       ├── Qing.gene.update.gtf
├── log
└── output
    ├── clean
    ├── sorted_result
    ├── hisat2_result
    └── expression_result
```
+ input 任务输入相关的文件
    + raw_data 存放下机数据（fastq.gz）
    + references 物种的参考基因组（fa）和基因注释文件（gtf）
    + index 参开基因组的索引文件（参考基因组生成）
+ log 任务中间日志文件
+ output
    + clean 质控处理之后的文件（clean.fastq.gz）
    + hisat2_result 序列跟参考基因组比对的结果文件（sam）
    + sorted_result 比对结果文件排序建立索引（sorted.bam）
    + expression_result 表达量矩阵和fpkm、tpm数据

### 下机数据处理

目前镜像版本（dockerhub.qingcloud.com/iot_demo/r-base:v0.6）
查看镜像

```
sudo docker images
```

运行镜像生成容器任务

```
sudo docker run -it -d --name test -v /home/hexing/test:/data dockerhub.qingcloud.com/iot_demo/r-base:v0.6
```

或者

```
sudo docker run -it -d --name test -v /home/hexing/test:/data {imageId}
```

可以查看任务进度：

```
sudo docker logs -f {containerId}
```


### 结果生成

目录结构

#### input（输入文件）
- index（参考基因组生成的索引）
- raw_data（输入下机测序数据fastq）
- references（参考基因组gtf文件和fa文件）

#### output（输出文件）
- clean（质控和过滤后的文件）
- hisat2_result（比对到参考基因组生成的文件）
- sorted_result（压缩排序建立索引后的文件）
- expression_result(表达矩阵数据)
#### log（运行时日志）
目前生成一些运行时日志


## 数据分析作图

### 配置文件
由于差异分析会涉及到样本之间的比较，目前会使用配置文件的方式设计比较组
```json
{
  "group": [
    {
      "start": "1",
      "end": "6",
      "name": "G14H-G21H",
      "repeated": "3",
      "output": "/data/output/diff"
    },
    {
      "start": "6",
      "end": "9",
      "name": "G21H-G21L",
      "repeated": "3",
      "output": "/data/output/diff"
    }
  ]
}
```
参数说明（A vs B）：
- group 表示当前分组设置（数组中每一个元素代表一个分组比较）
- start 表示根据表达矩阵设置的组A的起始列
- end 表示根据表达矩阵设置的组B的起始列
- repeated 表示当前表达矩阵的单个样本的重复数
- output 表示生成的相关图片的位置（默认在工作目录下的）

配置文件需要在根目录下的创建，命名为config.json

### 运行
由于出图会很快我们不需要后台运行直接进入容器虚拟环境运行相关命令即可

目前镜像版本（dockerhub.qingcloud.com/iot_demo/base:v0.3）
查看镜像

```
sudo docker images
```

运行镜像生成容器任务

```
sudo docker run -it --rm --name test -v /home/hexing/test:/data dockerhub.qingcloud.com/iot_demo/base:v0.3 /bin/bash
```

进入到容器之后 执行
```
gene-analyse diff
```

执行完成之后退出环境，查看结果
```
exit
```