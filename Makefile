# gene Makefile

GOCMD = GO111MODULE=on go
VERSION = v0.0.1
GIT_DIR = "gene-analyse"
GIT_COMMIT=$(shell git rev-parse HEAD)
GIT_BRANCH=$(shell git name-rev --name-only HEAD)
BUILD_DATE=$(shell date '+%Y-%m-%d-%H:%M:%S')
GORUN = $(GOCMD) run -ldflags "-X $(GIT_DIR)/pkg/version.GitCommit=${GIT_COMMIT} -X $(GIT_DIR)/pkg/version.BuildDate=${BUILD_DATE}"
LDFLAGS=-X $(GIT_DIR)/pkg/version.Version=$(VERSION) -X $(GIT_DIR)/pkg/version.BuildDate=$(BUILD_DATE) -X $(GIT_DIR)/pkg/version.GitCommit=$(GIT_COMMIT) -X $(GIT_DIR)/pkg/version.GitBranch=$(GIT_BRANCH)

GOBUILD = $(GOCMD) build -gcflags "all=-N -l" -ldflags "$(LDFLAGS)"
BINNAME = gene-analyse

## help: Help for this project
help: Makefile
	@echo "Usage:\n  make [command]"
	@echo
	@echo "Available Commands:"
	@sed -n 's/^##//p' $< | column -t -s ':' |  sed -e 's/^/ /'
run:
	@echo "---------------------------"
	@echo "-		 Run			 -"
	@echo "---------------------------"
	@$(GORUN) . serve  --conf ./config.yml
dev:
	@rm -rf bin/
	@mkdir bin/
	@echo "---------------------------"
	@echo "-		build...		 -"
	@$(GOBUILD)	-o bin/linux/$(BINNAME)
	@echo "-	 build(linux)...	 -"
	@CGO_ENABLED=0 GOOS=linux GOARCH=amd64  $(GOBUILD) -o bin/linux/$(BINNAME)
	@echo "-	builds completed!	-"
	@echo "---------------------------"
	@TEST=`./build.sh`
	@echo $(TEST)
build:
	@rm -rf bin/
	@mkdir bin/
	@echo "---------------------------"
	@echo "-		build...		 -"
	@$(GOBUILD)	-o bin/linux/$(BINNAME)
	@echo "-	 build(linux)...	 -"
	@CGO_ENABLED=0 GOOS=linux GOARCH=amd64  $(GOBUILD) -o bin/linux/$(BINNAME)
	@echo "-	builds completed!	-"
	@echo "---------------------------"
	@bin/linux/$(BINNAME) version
clean:
	@@rm -rf bin/*
	@echo "-	clean completed!	-"
docker: build
	docker build -t github.com/qiaoba0728/gene-analyse:v0.0.1 .


.PHONY: install generate

-include .dev/*.makefile