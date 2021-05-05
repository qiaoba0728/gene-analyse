package common

import (
	"crypto/tls"
	"fmt"
	"github.com/jordan-wright/email"
	"go.uber.org/zap"
	"net/smtp"
	"time"
)

var (
	pool      *email.Pool
	maxClient int = 10
)

type Email interface {
	SendEmail(from string, to []string, secret string, host string, nickname string, subject string, body string, port int, ssl bool) error
	SendEmailWithFile(from string, to []string, secret string, host string, nickname string, subject string, body string, port int, ssl bool, attach string) error
	SendEmailWithPool(to []string, from, secret, host, subject, body, nickname string, port int) error
	SendEmailWithPoolAndFile(to []string, from, secret, host, subject, body, nickname string, port int, attach string) error
}

type qqEmail struct {
	logger *zap.Logger
}

func NewQQEmail(logger *zap.Logger) Email {
	return &qqEmail{
		logger: logger,
	}
}

func (q *qqEmail) SendEmail(from string, to []string, secret string, host string, nickname string, subject string, body string, port int, ssl bool) error {
	auth := smtp.PlainAuth("", from, secret, host)
	e := email.NewEmail()
	e.From = fmt.Sprintf("%s<%s>", nickname, from)
	e.To = to
	e.Subject = subject
	e.HTML = []byte(body)
	hostAddr := fmt.Sprintf("%s:%d", host, port)
	if ssl {
		return e.SendWithTLS(hostAddr, auth, &tls.Config{ServerName: host})
	}
	return e.Send(hostAddr, auth)
}

func (q *qqEmail) SendEmailWithFile(from string, to []string, secret string, host string, nickname string, subject string, body string, port int, ssl bool, attach string) error {
	auth := smtp.PlainAuth("", from, secret, host)
	e := email.NewEmail()
	e.From = fmt.Sprintf("%s<%s>", nickname, from)
	e.To = to
	e.Subject = subject
	e.HTML = []byte(body)
	if attach != "" {
		_, _ = e.AttachFile(attach)
	}
	hostAddr := fmt.Sprintf("%s:%d", host, port)
	if ssl {
		return e.SendWithTLS(hostAddr, auth, &tls.Config{ServerName: host})
	}
	return e.Send(hostAddr, auth)
}

func (q *qqEmail) SendEmailWithPool(to []string, from, secret, host, subject, body, nickname string, port int) (err error) {
	hostAddr := fmt.Sprintf("%s:%d", host, port)
	auth := smtp.PlainAuth("", from, secret, host)
	if pool == nil {
		pool, err = email.NewPool(hostAddr, maxClient, auth)
		if err != nil {
			q.logger.Error("SendEmailWithPool fail", zap.Error(err))
		}
	}
	e := &email.Email{
		From:    fmt.Sprintf("%s<%s>", nickname, from),
		To:      to,
		Subject: subject,
		Text:    []byte(body),
	}
	return pool.Send(e, 5*time.Second)
}

func (q *qqEmail) SendEmailWithPoolAndFile(to []string, from, secret, host, subject, body, nickname string, port int, attach string) (err error) {
	hostAddr := fmt.Sprintf("%s:%d", host, port)
	auth := smtp.PlainAuth("", from, secret, host)
	if pool == nil {
		pool, err = email.NewPool(hostAddr, maxClient, auth)
		if err != nil {
			q.logger.Error("SendEmailWithPoolAndFile fail", zap.Error(err))
		}
	}
	e := &email.Email{
		From:    fmt.Sprintf("%s<%s>", nickname, from),
		To:      to,
		Subject: subject,
		Text:    []byte(body),
	}
	if attach != "" {
		_, _ = e.AttachFile(attach)
	}
	return pool.Send(e, 5*time.Second)
}

//163.com:
//
//POP3服务器地址:pop.163.com（端口：110）
//
//SMTP服务器地址:smtp.163.com（端口：25）
//
//
//
//126邮箱：
//
//POP3服务器地址:pop.126.com（端口：110）
//
//SMTP服务器地址:smtp.126.com（端口：25）
//
//
//
//139邮箱：
//
//POP3服务器地址：POP.139.com（端口：110）
//
//SMTP服务器地址：SMTP.139.com(端口：25)
//
//
//
//QQ邮箱：
//
//POP3服务器地址：pop.qq.com（端口：110）
//
//SMTP服务器地址：smtp.qq.com （端口：25）
//
//
//
//QQ企业邮箱 ：
//
//POP3服务器地址：pop.exmail.qq.com （SSL启用 端口：995）
//
//SMTP服务器地址：smtp.exmail.qq.com（SSL启用 端口：587/465）
//
//
//
//gmail(google.com) ：
//
//POP3服务器地址:pop.gmail.com（SSL启用 端口：995）
//
//SMTP服务器地址:smtp.gmail.com（SSL启用 端口：587）
//
//
//
//Foxmail：
//
//POP3服务器地址:POP.foxmail.com（端口：110）
//
//SMTP服务器地址:SMTP.foxmail.com（端口：25）
//
//
//
//sina.com:
//
//POP3服务器地址:pop3.sina.com.cn（端口：110）
//
//SMTP服务器地址:smtp.sina.com.cn（端口：25）
//
//
//
//sinaVIP：
//
//POP3服务器:pop3.vip.sina.com （端口：110）
//
//SMTP服务器:smtp.vip.sina.com （端口：25）
//
//
//
//sohu.com:
//
//POP3服务器地址:pop3.sohu.com（端口：110）
//
//SMTP服务器地址:smtp.sohu.com（端口：25）
//
//
//
//yahoo.com:
//
//POP3服务器地址:pop.mail.yahoo.com
//
//SMTP服务器地址:smtp.mail.yahoo.com
//
//
//
//yahoo.com.cn:
//
//POP3服务器地址:pop.mail.yahoo.com.cn（端口：995）
//
//SMTP服务器地址:smtp.mail.yahoo.com.cn（端口：587  ）
//
//
//
//HotMail ：
//
//POP3服务器地址：pop3.live.com （端口：995）
//
//SMTP服务器地址：smtp.live.com （端口：587）
//
//
//
//263.net:
//
//POP3服务器地址:pop3.263.net（端口：110）
//
//SMTP服务器地址:smtp.263.net（端口：25）
//
//
//
//263.net.cn:
//
//POP3服务器地址:pop.263.net.cn（端口：110）
//
//SMTP服务器地址:smtp.263.net.cn（端口：25）
//
//
//
//x263.net:
//
//POP3服务器地址:pop.x263.net（端口：110）
//
//SMTP服务器地址:smtp.x263.net（端口：25）
//
//
//
//21cn.com:
//
//POP3服务器地址:pop.21cn.com（端口：110）
//
//SMTP服务器地址:smtp.21cn.com（端口：25）
//
//
//china.com:
//
//POP3服务器地址:pop.china.com（端口：110）
//
//SMTP服务器地址:smtp.china.com（端口：25）
//
//
//tom.com:
//
//POP3服务器地址:pop.tom.com（端口：110）
//
//SMTP服务器地址:smtp.tom.com（端口：25）
//
//
//etang.com:
//
//POP3服务器地址:pop.etang.com
//
//SMTP服务器地址:smtp.etang.com
//
//常用邮箱的服务器地址与端口
