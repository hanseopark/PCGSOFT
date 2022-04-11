#echo -e "00010113_411790007l032230000_2l631031000000d0" > CutSelection.log
echo -e "00010113_411790607l032230000_2l631031000000d0" > CutSelection.log
#echo -e "00010113_411791107l032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/DataMerged/GammaCalo_912.root /home/alidock/alice/work/Data/pp_13TeV/MCMerged/GammaCalo_912.root eps  < answers_EMCal.txt
