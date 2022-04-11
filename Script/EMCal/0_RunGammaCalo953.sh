# Cluster shape cut
echo -e "00010113_411790607l032220000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal.txt

echo -e "00010113_411790607l032250000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal.txt

echo -e "00010113_411790607l032220000_2l631031000000d0\n00010113_411790607l032250000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal_953_shape.txt

# TM Cut
echo -e "00010113_411790607e032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal.txt

echo -e "00010113_411790607g032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal.txt

echo -e "00010113_411790607h032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal.txt

echo -e "00010113_411790607e032230000_2l631031000000d0\n00010113_411790607g032230000_2l631031000000d0\n00010113_411790607h032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_953.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_953.root eps  < answers_EMCal_953_TM.txt
