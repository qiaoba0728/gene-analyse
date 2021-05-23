package build

import (
	"fmt"
	"os/exec"
	"testing"
)

func TestGenePlugin_Build(t *testing.T) {
	//res := strings.Replace(types.R1Sample, "1", "2", -1)
	//fmt.Println(res)
	cmd := exec.Command("/bin/bash", "-c", "echo", "xxxx")
	fmt.Println(cmd.String())
	cmd.Run()
}
