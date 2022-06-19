package main

import (
	"context"
	"fmt"
	"fyne.io/fyne"
	"fyne.io/fyne/app"
	"fyne.io/fyne/container"
	"fyne.io/fyne/dialog"
	"fyne.io/fyne/widget"
	"github.com/flopp/go-findfont"
	"github.com/qiaoba0728/gene-analyse/internal/build"
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"go.uber.org/zap"
	"net/url"
	"os"
	"path"
	"strings"
)

func main() {
	logger, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
	a := app.New()
	w := a.NewWindow("zqz-tool")
	w.Resize(fyne.NewSize(800, 300))
	w.CenterOnScreen()
	output := widget.NewLabel("日志输出:")
	//合同编号
	orderEntry := widget.NewEntry()
	orderEntry.SetPlaceHolder("Q202206029")
	orderEntry.OnChanged = func(content string) {
		logger.Info("orderEntry", zap.String("order", orderEntry.Text))
		output.SetText(fmt.Sprintf("合同编号:%s", orderEntry.Text))
	}
	//项目名称
	projecNameEntry := widget.NewEntry()
	projecNameEntry.SetPlaceHolder("2个家蚕BSA测序分析")
	projecNameEntry.OnChanged = func(content string) {
		logger.Info("projecNameEntry", zap.String("projectName", projecNameEntry.Text))
		output.SetText(fmt.Sprintf("项目名称:%s", projecNameEntry.Text))
	}

	//客户单位
	accountEntry := widget.NewEntry()
	accountEntry.SetPlaceHolder("湖北省农业科学院经济作物研究所")
	accountEntry.OnChanged = func(content string) {
		logger.Info("accountEntry", zap.String("account", accountEntry.Text))
		output.SetText(fmt.Sprintf("客户单位:%s", projecNameEntry.Text))
	}

	userEntry := widget.NewEntry()
	userEntry.SetPlaceHolder("xxx(客户姓名)")
	userEntry.OnChanged = func(content string) {
		logger.Info("userEntry", zap.String("user", userEntry.Text))
		output.SetText(fmt.Sprintf("客户姓名:%s", projecNameEntry.Text))
	}
	inputFile := widget.NewLabel("输入压缩文件(zip):")
	input := widget.NewEntry()
	openButton := widget.NewButton("open", func() {
		dialog.ShowFileOpen(func(closer fyne.URIReadCloser, err error) {
			if err != nil {
				dialog.ShowError(err, w)
				return
			}
			if closer == nil {
				return
			}
			u, err := url.Parse(closer.URI().String())
			if err != nil {
				dialog.ShowError(err, w)
				return
			}
			logger.Info("get path", zap.String("path", u.String()))
			input.SetText(fmt.Sprintf("%s%s", u.Host, u.RawPath))
		}, w)
	})
	buildButton := widget.NewButton("重新生成", func() {
		output.SetText(fmt.Sprintf("%s \n%s %s %s %s",
			input.Text,
			orderEntry.Text, projecNameEntry.Text, accountEntry.Text, userEntry.Text))
		plugin := build.NewCovertPlugin(input.Text, projecNameEntry.Text, orderEntry.Text, userEntry.Text, accountEntry.Text, logger.Named("covert"))
		err := plugin.Build(context.Background())
		if err != nil {
			output.SetText(err.Error())
		} else {
			wd, _ := os.Getwd()
			output.SetText(fmt.Sprintf("输出文件目录: %s", path.Join(wd, orderEntry.Text)))
		}
	})

	head := container.NewBorder(nil, nil, inputFile, openButton, input)
	ctnt := container.NewVBox(head, orderEntry, projecNameEntry, accountEntry, userEntry, buildButton, output)
	w.SetContent(ctnt)

	w.ShowAndRun()
}
func init() {
	fontPaths := findfont.List()
	for _, path := range fontPaths {
		if strings.Contains(path, "simkai.ttf") {
			os.Setenv("FYNE_FONT", path)
			break
		}
	}
}
