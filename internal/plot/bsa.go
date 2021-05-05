package plot

import (
	"context"
	"errors"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"go.uber.org/zap"
	"os"
	"os/exec"
)

type bsapPlugin struct {
	input  string
	logger *zap.Logger
}

func NewBsapPlugin(input string, logger *zap.Logger) types.Plugin {
	return &bsapPlugin{
		input:  input,
		logger: logger,
	}
}
func (b *bsapPlugin) Build(ctx context.Context) error {
	step := os.Getenv("step")
	window := os.Getenv("window")
	thread := os.Getenv("thread")
	if step == "" || window == "" || thread == "" {
		return errors.New("env is not set")
	}
	cmd := exec.Command("Rscript", "/bsa/BSA_permutation_parallel.R", b.input, window, step, thread)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	b.logger.Info("cmd run ", zap.String("cmd", cmd.String()))
	if err := cmd.Run(); err != nil {
		b.logger.Error("Rscript sorted run fail", zap.Error(err))
		return err
	}
	b.logger.Info("plot finished", zap.String("step", step), zap.String("window", window), zap.String("thread", thread))
	return nil
}
func (b *bsapPlugin) Name() string {
	return "bsap"
}
