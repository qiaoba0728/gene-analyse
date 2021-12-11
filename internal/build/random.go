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

const max = 1000

type config struct {
	Data [][]string `json:"data"`
}
type randomPlugin struct {
	logger *zap.Logger
}

func NewRandomPlugin(logger *zap.Logger) types.Plugin {
	return &randomPlugin{
		logger: logger,
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
	)
	items := utils.RandFloats(55047245, 65047245, 1)
	readNo = items[0]
	rates := utils.RandFloats(149, 151, 1)
	basesBp = readNo * rates[0]

	items = utils.RandFloats(0.000090, 0.000100, 1)
	readN = items[0]
	items = utils.RandFloats(0.9600, 0.9800, 1)
	readQ20 = items[0]
	items = utils.RandFloats(0.9200, 0.9500, 1)
	readQ30 = items[0]

	items = utils.RandFloats(0.9600, 0.9800, 1)
	cleanReadsRate = items[0]
	items = utils.RandFloats(0.9200, 0.9400, 1)
	cleanDataRate = items[0]
	cleanReadsNo = cleanReadsRate * readNo
	cleanDataBp = cleanDataRate * basesBp

	results = append(results, fmt.Sprintf("%.f", readNo), fmt.Sprintf("%.f", basesBp),
		fmt.Sprintf("%.6f", readN),
		fmt.Sprintf("%.4f", readQ20),
		fmt.Sprintf("%.4f", readQ30),
		fmt.Sprintf("%.f", cleanReadsNo),
		fmt.Sprintf("%.f", cleanDataBp),
		fmt.Sprintf("%.4f", cleanReadsRate),
		fmt.Sprintf("%.4f", cleanDataRate))
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
	return xlsx.SaveAs("./random.xlsx")
}
