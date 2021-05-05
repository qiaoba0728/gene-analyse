/**
 * @Author: hexing
 * @Description:
 * @File:  log_test
 * @Version: 1.0.0
 * @Date: 20-3-24 下午6:15
 */

package utils

import (
	"go.uber.org/zap"
	"testing"
)

func TestNewZapLogger(t *testing.T) {
	l, _ := NewZapLogger(&Opt{LogOutput: CONSOLE})
	l.Named("test").Info("hello world", zap.String("err", "error"))
}
