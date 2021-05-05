package utils

import (
	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
	"gopkg.in/natefinch/lumberjack.v2"
	"io"
	"os"
	"path/filepath"
)

type Output int

const (
	FILE Output = iota
	CONSOLE
	ALL
)

type Level int

const (
	DEBUG Level = iota
	INFO
	WARN
	ERROR
	FATAL
)

func (l *Level) Unmarshal(text string) {
	switch text {
	case "debug", "DEBUG":
		*l = DEBUG
	case "info", "INFO", "": // make the zero value useful
		*l = INFO
	case "warn", "WARN":
		*l = WARN
	case "error", "ERROR":
		*l = ERROR
	case "fatal", "FATAL":
		*l = FATAL
	default:
		*l = INFO
	}
}

type Opt struct {
	LogPath    string
	LogName    string
	LogLevel   string
	MaxSize    int
	MaxBackup  int
	MaxAge     int
	LogOutput  Output
	JsonFormat bool
}

func NewZapLogger(opt *Opt) (*zap.Logger, io.Writer) {
	var (
		writer   io.Writer
		core     zapcore.Core
		zapLevel zapcore.Level
	)
	switch opt.LogOutput {
	case FILE:
		writer = newFileWriter(opt)
	case CONSOLE:
		writer = os.Stdout
	case ALL:
		writer = io.MultiWriter(newFileWriter(opt), os.Stdout)
	}
	encoderConfig := zap.NewProductionEncoderConfig()
	encoderConfig.EncodeTime = zapcore.ISO8601TimeEncoder
	var l Level
	l.Unmarshal(opt.LogLevel)
	zapLevel = asZapLevel(l)
	if !opt.JsonFormat {
		core = zapcore.NewCore(zapcore.NewConsoleEncoder(encoderConfig), zapcore.AddSync(writer), zapLevel)
	} else {
		core = zapcore.NewCore(zapcore.NewJSONEncoder(encoderConfig), zapcore.AddSync(writer), zapLevel)
	}

	if zapLevel == zapcore.DebugLevel {
		//return zap.New(core, zap.AddCaller()), writer
		return zap.New(core, zap.Fields(zap.String("source", "gene"),
			zap.String("appid", "analyse")), zap.AddCaller()), writer
	}
	return zap.New(core, zap.Fields(zap.String("source", "gene"),
		zap.String("appid", "analyse"))), writer
}
func newFileWriter(opt *Opt) io.Writer {
	if opt.LogPath == "" {
		opt.LogPath = os.TempDir()
	}
	if opt.MaxSize <= 0 {
		opt.MaxSize = 100
	}
	if opt.MaxBackup <= 0 {
		opt.MaxBackup = 10
	}
	if opt.MaxAge <= 0 {
		opt.MaxAge = 28
	}
	return &lumberjack.Logger{
		Filename:   filepath.Join(opt.LogPath, opt.LogName+".log"),
		MaxSize:    opt.MaxSize,
		MaxBackups: opt.MaxBackup,
		MaxAge:     opt.MaxAge,
		LocalTime:  true,
		Compress:   true,
	}
}

func asZapLevel(level Level) zapcore.Level {
	zapLevel := zap.InfoLevel
	switch level {
	case DEBUG:
		zapLevel = zap.DebugLevel
	case INFO:
		zapLevel = zap.InfoLevel
	case WARN:
		zapLevel = zap.WarnLevel
	case ERROR:
		zapLevel = zap.ErrorLevel
	case FATAL:
		zapLevel = zap.FatalLevel
	}
	return zapLevel
}
