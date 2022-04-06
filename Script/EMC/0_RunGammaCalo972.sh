# Timing cut
echo -e "0008e113_411790605l032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_972.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_972.root eps  < answers_EMCal.txt

echo -e "0008e113_411790609l032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_972.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_972.root eps  < answers_EMCal.txt

echo -e "0008e113_411790608l032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_972.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_972.root eps  < answers_EMCal.txt

echo -e "0008e113_411790605l032230000_2l631031000000d0\n0008e113_411790609l032230000_2l631031000000d0\n0008e113_411790608l032230000_2l631031000000d0" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_972.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_972.root eps  < answers_EMCal_972_timing.txt
# # Alpha
#echo -e "80000113_00200009327000008250400000_2444451041013200000_0163103100000010\n80000113_00200009327000008250400000_2444451041013200000_0163105100000010\n80000113_00200009327000008250400000_2444451041013200000_0163106100000010" > CutSelectionAlpha.log
#    echo -e "Alpha\nLHC13bc\n3\nY\npPb5\nY\n/home/mike/2_pPb_EMC/0_analysis/170803_final_EMC/CocktailEMC_4Mio.root\n0.80\nN\nY\nY" > answersAlpha.txt
#    cp CutSelectionAlpha.log CutSelection.log
#    bash start_FullMesonAnalysis_TaskV3.sh -dgammaOff bla eps < answersAlpha.txt
