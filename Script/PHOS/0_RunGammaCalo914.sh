#MB
echo -e "00010113_2446651044012300000_2163103100000010" > CutSelection.log
    bash start_FullMesonAnalysis_TaskV3.sh -gammaOff /home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_914.root /home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_914.root eps  < answers_PHOS.txt
