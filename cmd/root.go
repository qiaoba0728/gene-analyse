/*
 * Copyright (C) 2019 Yunify, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this work except in compliance with the License.
 * You may obtain a copy of the License in the LICENSE file, or at:
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package cmd

import (
	"fmt"
	"github.com/spf13/cobra"
	"os"
)

var rootCmd = &cobra.Command{
	Use:   "gene-analyse",
	Short: "This is gene-analyse service.",
	Long:  `This is gene-analyse service.`,
}

//Execute execute the root command
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

func init() {
	//rootCmd.PersistentFlags().StringVarP(&metricsBackend, "metrics", "m", "expvar", "Metrics backend (expvar|prometheus)")
	//rand.Seed(int64(time.Now().Nanosecond()))
	//logger, _ = log.NewProduction(
	//	logf.AddStacktrace(log.FatalLevel),
	//)
	rootCmd.AddCommand(diffCmd)
	rootCmd.AddCommand(mainCmd)
	rootCmd.AddCommand(expressionCmd)
	rootCmd.AddCommand(syncCmd)
	rootCmd.AddCommand(testCmd)
	rootCmd.AddCommand(bsaCmd)
	rootCmd.AddCommand(singleCmd)
	rootCmd.AddCommand(buildExCmd)
	rootCmd.AddCommand(bsapCmd)
	rootCmd.AddCommand(gatkCmd)
	rootCmd.AddCommand(reportCmd)
	rootCmd.AddCommand(randomCmd)
	cobra.OnInitialize(onInitialize)
	//flag.CommandLine.AddGoFlagSet(goflag.CommandLine)
} //

// onInitialize is called before the command is executed.
func onInitialize() {
}
