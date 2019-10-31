//////////////////////////////////////////////////////////////
///////// My Analysis ///////////////////////////////////////
/////////////////////////////////////////////////////////////


void MyAnalysis(){

TString EventCuts = "00010113";
TString CaloCuts  = "411790007l032230000";
TString MesonCuts = "2l631031000000d0";
TStirng Sufix = "13TeV";

TString CutSelection = Form("%s_%s_%s",EvnetCuts.Data,CaloCuts.Data,MesonCuts.Data);



TFile *file1 = new TFile(Form("%s/%s/Pi0_MC_GammaConvV1CorrectionHistos_00010113_411790007l032230000_2l631031000000d0.root",CutSelection.Data,Sufix.Data),"read");

TH1D* histoMesonEffi;
TH1D* histoTrueMesonEffi;
TH1D* histoAcceptance;

histoMesonEffi = (TH1D*)file1->FindObject("MesonEffiPt");
histoTrueMesonEffi = (TH1D*)file1->FindObject("TrueMesonEffiPt");
histoAcceptance = (TH1D*)file1->FindObject("fMCMesonAccepPt");

histoMesonEffi->Draw();





}
