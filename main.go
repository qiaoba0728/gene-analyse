package main

import (
	"flag"
	"github.com/qiaoba0728/gene-analyse/cmd"
)

func main() {
	flag.Parse()
	//defer log.Flush()
	cmd.Execute()
}
