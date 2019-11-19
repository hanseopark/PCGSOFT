#! bin/bash

############################################################

ENERGY="13TeV"
#FILE="./CutStudies/13TeV/Pi0_data_SystematicErrorCuts.root"

FILEPI0="/home/alidock/alice/work/Main_Analysis/191118_0_pp_13TeV_EMCal/CutStudies/13TeV/Pi0_data_SystematicErrorCuts.root"
FILEETA="/home/alidock/alice/work/Main_Analysis/191118_0_pp_13TeV_EMCal/CutStudies/13TeV/Eta_data_SystematicErrorCuts.root"
FILEPI0ETA="/home/alidock/alice/work/Main_Analysis/191118_0_pp_13TeV_EMCal/CutStudies/13TeV/Pi0EtaBinning_data_SystematicErrorCuts.root"
MESONPI0="Pi0"
MESONETA="ETA"
MESONPI0ETA="Pi0EtaBinning"
NUMBEROFPTBIN=19
NUMBERCUTSTUDIES=8
STARTPTSYS=0
NAME="pp"
NAMEOUTPUT="INT7"
SUFFIX="eps"
MODE=4

function SystematicError()
{
	root -l -x -b -q FinaliseSystematicErrorsCalo_ppV2.C\+\($1\)
}

OPTIONPI0=\"$FILEPI0\"\,\"$ENERGY\"\,\"$MESONPI0\"\,$NUMBEROFPTBIN\,$NUMBERCUTSTUDIES\,$STARTPTSYS\,\"$NAME\"\,\"$NAMEOUTPUT\",\"$SUFFIX\"\,$MODE
#OPTIONETA=\"$FILEETA\"\,\"$ENERGY\"\,\"$MESONETA\"\,$NUMBEROFPTBIN\,\,$NUMBERCUTSTUDIES\,$STARTPTSYS\,\"$NAME\"\,\"$NAMEOUTPUT\"\,\"$SUFFIX\"\,$MODE
#OPTIONPI0ETA=\"$FILEPI0ETA\"\,\"$ENERGY\"\,\"$MESONPI0ETA\"\,$NUMBEROFPTBIN\,\,$NUMBERCUTSTUDIES\,$STARTPTSYS\,\"$NAME\"\,\"$NAMEOUTPUT\"\,\"$SUFFIX\"\,$MODE

echo -e"PI0"
SystematicError $OPTIONPI0

#echo -e"Eta"
#SystematicError $OPTIONETA

#echo -e"Pi0EtaBinning"
#SystematicError $OPTIONPI0ETA
