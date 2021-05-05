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

### 下机数据处理

目前镜像版本（dockerhub.qingcloud.com/iot_demo/conda:v0.2）
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



for var in *_R2.fq.gz; do mv "$var" "${var%_R2.fq.gz}.R2.fastq.gz"; done