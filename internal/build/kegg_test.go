package build

import (
	"context"
	"github.com/stretchr/testify/assert"
	"testing"
)

func TestNewKeggPlugin(t *testing.T) {
	p := NewKeggPlugin()
	err := p.Build(context.Background())
	assert.Nil(t, err)
}
