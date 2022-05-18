# Energy Calibration
echo -e "00010113_411790607l022230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_951.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_951.root eps  < Answers/answers_EMCal.txt

echo -e "00010113_411790607l042230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_951.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_951.root eps  < Answers/answers_EMCal.txt

echo -e "00010113_411790607l052230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_951.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_951.root eps  < Answers/answers_EMCal.txt

# Energy cluster cut
echo -e "00010113_411790607l022230000_2l631031000000d0\n00010113_411790607l042230000_2l631031000000d0\n00010113_411790607l052230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_951.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_951.root eps  < Answers/answers_EMCal_951_energy.txt

# NCell cut
echo -e "00010113_411790607l031230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_951.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_951.root eps  < Answers/answers_EMCal.txt

echo -e "00010113_411790607l033230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_951.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_951.root eps  < Answers/answers_EMCal.txt

echo -e "00010113_411790607l031230000_2l631031000000d0\n00010113_411790607l033230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff bla eps  < Answers/answers_EMCal_951_ncells.txt


# # Alpha
#echo -e "80000113_00200009327000008250400000_2444451041013200000_0163103100000010\n80000113_00200009327000008250400000_2444451041013200000_0163105100000010\n80000113_00200009327000008250400000_2444451041013200000_0163106100000010" > CutSelectionAlpha.log
#    echo -e "Alpha\nLHC13bc\n3\nY\npPb5\nY\n/home/mike/2_pPb_EMC/0_analysis/170803_final_EMC/CocktailEMC_4Mio.root\n0.80\nN\nY\nY" > Answers/answersAlpha.txt
#    cp CutSelectionAlpha.log CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff bla eps < Answers/answersAlpha.txt
