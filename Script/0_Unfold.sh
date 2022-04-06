#! bin/bash

#############################################################333

ENERGY="13TeV"
CUTSELECTION="00010113_411791107l032230000_3l631031000000d0"
SUFFIX="eps"
UNFOLDINGFILE="/home/alidock/alice/work/Data/pp_13TeV/MCMerged/GammaCalo_913.root"
crystal=Gaussian
DIRECTPHOTON="No"
ADVMESONQA="AdvancedMesonQA"
BINSPTPI0=18 #Jet is limited by pt <10GeV
BINSPTETA=6 #Jet is limited by pt <10GeV
MODE=4
USETHNSPARSE=0
CORRFSETTING=""
function ExtractSignal()
{
       root -x -q -l -b  TaskV1/ExtractSignalV2.C\+\($1\,$MODE\,$USETHNSPARSE\,-1\,\"$CORRFSETTING\"\)
}
UnfoldingEnergy=$ENERGY
UnfoldingEnergy+="_Unfolding_AsData"
mkdir -p $CUTSELECTION/$UnfoldingEnergy/$SUFFIX
OPTIONSPI0MC=\"Pi0\"\,\"$UNFOLDINGFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$UnfoldingEnergy\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
ExtractSignal $OPTIONSPI0MC

OPTIONSETAMC=\"Eta\"\,\"$UNFOLDINGFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$UnfoldingEnergy\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
ExtractSignal $OPTIONSETAMC

UnfoldingEnergy=$ENERGY
UnfoldingEnergy+="_Unfolding_Missed"
mkdir -p $CUTSELECTION/$UnfoldingEnergy/$SUFFIX
OPTIONSPI0MC=\"Pi0\"\,\"$UNFOLDINGFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$UnfoldingEnergy\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
ExtractSignal $OPTIONSPI0MC

OPTIONSETAMC=\"Eta\"\,\"$UNFOLDINGFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$UnfoldingEnergy\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
ExtractSignal $OPTIONSETAMC

UnfoldingEnergy=$ENERGY
UnfoldingEnergy+="_Unfolding_Reject"
mkdir -p $CUTSELECTION/$UnfoldingEnergy/$SUFFIX
OPTIONSPI0MC=\"Pi0\"\,\"$UNFOLDINGFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$UnfoldingEnergy\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
ExtractSignal $OPTIONSPI0MC

OPTIONSETAMC=\"Eta\"\,\"$UNFOLDINGFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$UnfoldingEnergy\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
ExtractSignal $OPTIONSETAMC

root -b -x -l -q Jet_Unfolding_Macro.C\+\(\"$CUTSELECTION\"\,\"$ENERGY\"\,$MODE\,$BINSPTPI0\,$BINSPTETA\)



