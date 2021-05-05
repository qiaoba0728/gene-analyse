package build

import (
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"strings"
	"testing"
)

func TestGenePlugin_Build(t *testing.T) {
	res := strings.Replace(types.R1Sample, "1", "2", -1)
	fmt.Println(res)
}
