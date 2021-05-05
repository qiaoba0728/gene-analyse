package common

import (
	"github.com/stretchr/testify/assert"
	"testing"
)

func TestBuildReport(t *testing.T) {
	err := BuildReport()
	assert.Nil(t, err)
}
