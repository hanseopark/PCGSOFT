# Opening angle cut
echo -e "00010113_411790607l032230000_2l631031000000b0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal.txt

echo -e "00010113_411790607l032230000_2l631031000000g0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal.txt

echo -e "00010113_411790607l032230000_2l631031000000a0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal.txt

echo -e "00010113_411790607l032230000_2l631031000000b0\n00010113_411790607l032230000_2l631031000000g0\n00010113_411790607l032230000_2l631031000000a0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal_954_opening.txt

# Alpha cut 
#echo -e "00010113_411790607l032230000_2l631041000000d0" > CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal.txt
#
#echo -e "00010113_411790607l032230000_2l631051000000d0" > CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal.txt
#
#echo -e "00010113_411790607l032230000_2l631061000000d0" > CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal.txt

echo -e "00010113_411790607l032230000_2l631041000000d0\n00010113_411790607l032230000_2l631051000000d0\n00010113_411790607l032230000_2l631061000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_954.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_954.root eps  < answers_EMCal_954_alpha.txt
