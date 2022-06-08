#! bin/bash

#######################################################

FILEPI0="${PCGDIR}/Result/pp/EMC/SystematicErrors/SystematicErrorAveragedSingleEMCEMC_Pi0_13TeVINT7_2021_09_30.dat"
FILEETA="${PCGDIR}/Result/pp/EMC/SystematicErrors/SystematicErrorAveragedSingleEMCEMC_Pi0_13TeVINT7_2021_09_30.dat"

FILEEMC="${PCGDIR}/Result/pp/EMC/data_EMCAL-EMCALResultsFullCorrection_PP.root"
FILEPHOS="${PCGDIR}/Result/pp/PHOS/data_PHOS-PHOSResultsFullCorrection_PP.root"
FILEPCM="${PCGDIR}/Result/pp/PCM/data_PCMResultsFullCorrection_PP.root"
FILEEMCPCM="${PCGDIR}/Result/pp/EMCPCM/data_PCM-EMCALResultsFullCorrection_PP.root"

COMBMODE="systems"
MESONPI0="Pi0"
MESONETA="Eta"
ENERGY="13TeV"
MODE=4
SUFFIX="eps"
ISSTATCORR=kFALSE
CENTRALITY=""
EVENTCUT=""
NOMPi0=3
NOMEta=2
MinPtPi0=1.4
MaxPtPi0=25.0
MinPtEta=2.0
MaxPtEta=40.0

function CorrelationFactors()
{
root -x -l -q -b ComputeCorrelationFactors.C\+\($1\)
}

function CombineMeasurement()
{
root -x -l -q -b CombineMesonMeasurements13TeV_Jets.C\+\($1\)
}
OPTIONPI0=\"$FILEPI0\"\,\"$COMBMODE\"\,\"$MESONPI0\"\,\"$ENERGY\"\,$MODE\,\"$SUFFIX\"\,$ISSTATCORR\,\"$CENTRALITY\"\,\"$EVENTCUT\"
OPTIONETA=\"$FILEETA\"\,\"$COMBMODE\"\,\"$MESONETA\"\,\"$ENERGY\"\,\"$MODE\"\,\"$SUFFIX\"\,\"$ISSTATCORR\"\,\"$CENTRALITY\"\,\"$EVENTCUT\"
COMBINEOPTION=\"$FILEPCM\"\,\"$FILEPHOS\"\,\"$FILEEMC\"\,\"\"\,\"\"\,\"$FILEEMCPCM\"\,\"\"\,$NOMPi0\,$NOMEta\,$MinPtPi0\,$MaxPtPi0\,$MinPtEta\,$MaxPtEta\,\"$SUFFIX\"\,\"\"

#echo -e "Pi0"
#CorrelationFactors $OPTIONPI0

#echo -e "Eta"
#CorrelationFactors $OPTIONETA

CombineMeasurement $COMBINEOPTION

