# !bin/bash

#FILE="/Users/hanseopark/alice/work/Main_Analysis/TriggerCombination_PFR.dat"
FILE="/Users/hanseopark/alice/work/Main_Analysis/TriggerCombination_PFR_bla.dat"
#FILE="/Users/hanseopark/alice/work/Main_Analysis/TriggerCombination_PFR_Only.dat"
MODE=4
TRIGGERNUMBER=1
SUFFIX="eps"
ISDATA="data"
ISMC="MC"
ENERGY="13TeV_Jets"
PERIOD="LHC16,17,18x"
PILEUPAPPLIED=kTRUE
MAXPTPI0=20
#AVERAGEDPI0=kTRUE
AVERAGEDPI0=kFALSE
ENABLEETA=kTRUE
MAXPTETA=40
AVERAGEDETA=kFALSE
V2CLUSTER=KFALSE
NAMEFILE=""
CLUSTEROUTPUT=kTRUE
FILEINPUT=""

function ProduceFinalResultsPatchedTriggers()
{
	root -l -q TaskV1/ProduceFinalResultsPatchedTriggers.C\+\($1\)
}

#OPTION=\"$FILE\"\,\"$MODE\"\,\"$TRIGGERNUMBER\"\,\"$SUFFIX\"\,\"$ISDATA\"\,\"$ENERGY\"\,\"$PERIOD\"\,\"$PILEUPAPPLIED\"\,\"$MAXPTPI0\"\,\"$AVERAGEDPI0\"\,\"$ENABLEETA\"\,\"$MAXPTETA\"\,\"$AVERAGEDETA\"\,\"$V2CLUSTER\"\,\"$NAMEFILE\"\,\"$CLUSTEROUTPUT\"\,\"$FILEINPUT\"
OPTION=\"$FILE\"\,$MODE\,$TRIGGERNUMBER\,\"$SUFFIX\"\,\"$ISDATA\"\,\"$ENERGY\"\,\"$PERIOD\"\,kTRUE\,$MAXPTPI0\,kFALSE\,kTRUE\,$MAXPTETA\,kFALSE\,kFALSE\,\"$NAMEFILE\"\,kTRUE\,\"$FILEINPUT\"
#root -l -b -x -q 'TaskV1/ProduceFinalResultsPatchedTriggers.C++("/Users/hanseopark/alice/work/Main_Analysis/200430_0_pp_13TeV_EMCal/TriggerCombination_PFR.dat",4,3,"eps","data","13TeV","LHC16,17,18x", kTRUE, 26, kTRUE, kTRUE, 50, kTRUE, kFALSE, "", kTRUE, "")'

echo -e "  To perform ProduceFinalResultsPatchedTriggers  "
ProduceFinalResultsPatchedTriggers $OPTION
