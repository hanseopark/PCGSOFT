#echo -e "00010113_411790007l032230000_2l631031000000d0" > CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_912.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_912.root eps  < answers_EMCal.txt
#echo -e "00010113_411790607l032230000_2l631031000000d0" > CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_912.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_912.root eps  < answers_EMCal.txt

# It is needs before main analysis
MESON="Pi0"
ENERGY="13TeV"
DataMainFile="/home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaConvV1_500.root"
DataAddFile=""
CUTNUMBER="00010113_00200009327000008250400000_2152103500000000"
NPTBINS=50
SUFFIX="eps"

root -b -x -q -l TaskV1/AnalyseDCADist.C+\(\"$MESON\"\,\"$DataMainFile\"\,\"$DataAddFile\"\,\"$CUTNUMBER\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"\"\,kFALSE\,$NPTBINS\,0\)

# main anaylsis
echo -e "00010113_00200009327000008250400000_2152103500000000" > CutSelection.log #I have not changed the MC data yet
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaConvV1_500.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaConvV1_500.root eps  < answers_PCM.txt
