package plot

import (
	"context"
	"errors"
	"fmt"
	"github.com/qiaoba0728/gene-analyse/internal/types"
	"go.uber.org/zap"
	"os"
	"os/exec"
	"time"
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
	// default
	step := os.Getenv("STEP")
	if step == "" {
		step = "10000"
	}
	window := os.Getenv("WINDOW")
	if window == "" {
		window = "50000"
	}
	thread := os.Getenv("THREAD")
	if step == "" || window == "" || thread == "" {
		return errors.New("env is not set")
	}
	dir := fmt.Sprintf("/data/output/%s/", fmt.Sprintf("%s_%s_%d", window, step, time.Now().Unix()))
	cmd := exec.Command("mkdir", "-p", dir)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		b.logger.Error("create bsa dna dir", zap.Error(err))
	}
	cmd = exec.Command("Rscript", "/work/BSA_permutation_parallel.R", b.input, window, step, thread, dir)
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
