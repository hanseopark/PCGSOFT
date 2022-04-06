# Opening angle cut
echo -e "0008d113_411790607l032230000_2l631031000000b0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_961.txt

echo -e "0008d113_411790607l032230000_2l631031000000g0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_961.txt

echo -e "0008d113_411790607l032230000_2l631031000000a0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_961.txt

echo -e "0008d113_411790607l032230000_2l631031000000b0\n0008d113_411790607l032230000_2l631031000000g0\n0008d113_411790607l032230000_2l631031000000a0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_964_opening.txt

# Alpha cut 
echo -e "0008d113_411790607l032230000_2l631041000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_961.txt

echo -e "0008d113_411790607l032230000_2l631051000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_961.txt

echo -e "0008d113_411790607l032230000_2l631061000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_961.txt

echo -e "0008d113_411790607l032230000_2l631041000000d0\n0008d113_411790607l032230000_2l631051000000d0\n0008d113_411790607l032230000_2l631061000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_964.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_964.root eps  < answers_EMCal_964_alpha.txt
