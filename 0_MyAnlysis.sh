#! bin/bash

mkdir -p MyAnalysisDir
cd MyAnalysisDir
root -l -q ../MyAnalysis.C
