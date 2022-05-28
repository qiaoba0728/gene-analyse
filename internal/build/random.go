package build

import (
	"context"
	"encoding/json"
	"fmt"
	"github.com/360EntSecGroup-Skylar/excelize"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"os"
)

const (
	max    = 1000
	six    = "6G"
	eight  = "8G"
	ten    = "10G"
	twel   = "12G"
	fift   = "15G"
	sixten = "16G"
	twenty = "20G"
)

type config struct {
	Data [][]string `json:"data"`
}
type randomPlugin struct {
	tp     string
	logger *zap.Logger
}

func NewRandomPlugin(tp string, logger *zap.Logger) types.Plugin {
	return &randomPlugin{
		logger: logger,
		tp:     tp,
	}
}
func (r *randomPlugin) Build(ctx context.Context) error {
	c := config{Data: make([][]string, max)}
	for i := 0; i < max; i++ {
		temp := r.buildItem()
		c.Data[i] = append(c.Data[i], temp...)
	}
	return r.writeExcel(c.Data)
}
func (r *randomPlugin) Name() string {
	return "randomPlugin"
}

func (r *randomPlugin) buildItem() []string {
	var (
		readNo         float64
		basesBp        float64
		results        []string
		readN          float64
		readQ20        float64
		readQ30        float64
		cleanReadsRate float64
		cleanDataRate  float64
		cleanReadsNo   float64
		cleanDataBp    float64
		items          []float64
	)
	// I
	items = utils.RandFloats(91.00, 94.50, 1)
	cleanReadsRate = items[0]
	// J
	cleanDataRate = cleanReadsRate
	// G
	switch r.tp {
	case six:
		items = utils.RandFloats(41000000, 53000000, 1)
	case eight:
		items = utils.RandFloats(54000000, 70000000, 1)
	case ten:
		items = utils.RandFloats(68000000, 82000000, 1)
	case twel:
		items = utils.RandFloats(80000000, 100000000, 1)
	case fift:
		items = utils.RandFloats(100000000, 120000000, 1)
	case sixten:
		items = utils.RandFloats(108000000, 125000000, 1)
	case twenty:
		items = utils.RandFloats(135000000, 160000000, 1)
	default:
		items = utils.RandFloats(2160000000, 2500000000, 1)
	}
	temp := int(items[0])
	if temp%2 != 0 {
		temp = temp + 1
	}
	cleanReadsNo = float64(temp)
	// H
	cleanDataBp = cleanReadsNo * 150

	// E
	items = utils.RandFloats(92, 95, 1)
	readQ30 = items[0]

	// C
	items = utils.RandFloats(0.007500, 0.009500, 1)
	readN = items[0]
	// A
	temp = int(cleanReadsNo * 100 / cleanReadsRate)
	if temp%2 != 0 {
		temp = temp + 1
	}
	readNo = float64(temp)

	basesBp = readNo * 150

	items = utils.RandFloats(96.00, 98.00, 1)
	// D
	readQ20 = items[0]

	results = append(results, fmt.Sprintf("%.f", readNo), fmt.Sprintf("%.f", basesBp),
		fmt.Sprintf("%.6f", readN),
		fmt.Sprintf("%.2f", readQ20),
		fmt.Sprintf("%.2f", readQ30),
		fmt.Sprintf("%.f", cleanReadsNo),
		fmt.Sprintf("%.f", cleanDataBp),
		fmt.Sprintf("%.2f", cleanReadsRate),
		fmt.Sprintf("%.2f", cleanDataRate))
	return results
}

func (r *randomPlugin) writeJson(data [][]string) error {
	resp, _ := json.Marshal(data)
	fp, err := os.OpenFile("random.json", os.O_RDWR|os.O_CREATE, 0755)
	if err != nil {
		r.logger.Error("open fail ", zap.Error(err))
		return err
	}
	defer fp.Close()
	_, err = fp.Write(resp)
	return err
}
func (r *randomPlugin) writeExcel(data [][]string) error {
	xlsx := excelize.NewFile()
	sheet := "Sheet1"
	index := xlsx.NewSheet(sheet)

	// build header
	xlsx.SetCellValue(sheet, fmt.Sprintf("A%d", 1), "readNo")
	xlsx.SetCellValue(sheet, fmt.Sprintf("B%d", 1), "basesBp")
	xlsx.SetCellValue(sheet, fmt.Sprintf("C%d", 1), "readN")
	xlsx.SetCellValue(sheet, fmt.Sprintf("D%d", 1), "readQ20")
	xlsx.SetCellValue(sheet, fmt.Sprintf("E%d", 1), "readQ30")

	xlsx.SetCellValue(sheet, fmt.Sprintf("G%d", 1), "cleanReadsNo")
	xlsx.SetCellValue(sheet, fmt.Sprintf("H%d", 1), "cleanDataBp")
	xlsx.SetCellValue(sheet, fmt.Sprintf("I%d", 1), "cleanReadsRate")
	xlsx.SetCellValue(sheet, fmt.Sprintf("J%d", 1), "cleanDataRate")
	for k, v := range data {
		xlsx.SetCellValue(sheet, fmt.Sprintf("A%d", k+2), v[0])
		xlsx.SetCellValue(sheet, fmt.Sprintf("B%d", k+2), v[1])
		xlsx.SetCellValue(sheet, fmt.Sprintf("C%d", k+2), v[2])
		xlsx.SetCellValue(sheet, fmt.Sprintf("D%d", k+2), v[3])
		xlsx.SetCellValue(sheet, fmt.Sprintf("E%d", k+2), v[4])

		xlsx.SetCellValue(sheet, fmt.Sprintf("G%d", k+2), v[5])
		xlsx.SetCellValue(sheet, fmt.Sprintf("H%d", k+2), v[6])
		xlsx.SetCellValue(sheet, fmt.Sprintf("I%d", k+2), v[7])
		xlsx.SetCellValue(sheet, fmt.Sprintf("J%d", k+2), v[8])
	}
	xlsx.SetActiveSheet(index)
	return xlsx.SaveAs(fmt.Sprintf("./%s_random.xlsx", r.tp))
}
