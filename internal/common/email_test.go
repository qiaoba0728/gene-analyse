package common

import (
	"github.com/qiaoba0728/gene-analyse/internal/utils"
	"testing"
)

func TestNewQQEmail(t *testing.T) {
	l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
	q := NewQQEmail(l)
	var to = []string{"771473657@qq.com"}
	from := "494340090@qq.com"
	nickname := "贺兴"
	secret := "ngybgeetypogbibc"
	host := "smtp.qq.com"
	port := 25
	subject := "转录组测序分析"
	body := "人类基因测序分析结果"
	if err := q.SendEmailWithPoolAndFile(to, from, secret, host, subject, body, nickname, port, "./config/build.json"); err != nil {
		t.Log("发送失败: ", err)
	} else {
		t.Log("发送成功")
	}
}
func TestQqEmail_SendEmail(t *testing.T) {
	l, _ := utils.NewZapLogger(&utils.Opt{LogOutput: utils.CONSOLE})
	q := NewQQEmail(l)
	var to = []string{"771473657@qq.com"}
	from := "494340090@qq.com"
	nickname := "贺兴"
	secret := "ngybgeetypogbibc"
	host := "smtp.qq.com"
	port := 25
	subject := "转录组测序分析"
	body := "人类基因测序分析结果"
	if err := q.SendEmail(from, to, secret, host, nickname, subject, body, port, false); err != nil {
		t.Log("发送失败: ", err)
	} else {
		t.Log("发送成功")
	}
}
