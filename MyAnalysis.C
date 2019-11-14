//////////////////////////////////////////////////////////////
///////// My Analysis ///////////////////////////////////////
/////////////////////////////////////////////////////////////


void MyAnalysis(){

	TString EventCuts = "00010113";
	TString CaloCuts  = "411790607l032230000";
	TString MesonCuts = "2l631031000000d0";
	TString Sufix = "13TeV";

	TString CutSelection = Form("%s_%s_%s",EventCuts.Data(),CaloCuts.Data(),MesonCuts.Data());

	TCanvas *c1 = new TCanvas("ratio Pi0/Eta","ratio Pi0/Eta",1200,800);


	TFile *Pi0EtaBinning_MC 	= new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Pi0EtaBinning_Data 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_MC 			= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Data 		= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	
	TFile *Pi0EtaBinning_MC_InclusiveJet 	= new TFile(Form("../%s/%s_InclusiveJet/Pi0EtaBinning_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Pi0EtaBinning_Data_InclusiveJet 	= new TFile(Form("../%s/%s_InclusiveJet/Pi0EtaBinning_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_MC_InclusiveJet 			= new TFile(Form("../%s/%s_InclusiveJet/Eta_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Data_InclusiveJet 		= new TFile(Form("../%s/%s_InclusiveJet/Eta_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");

	TFile *Pi0EtaBinning_Factor 		 = new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConv_OnlyCorrectionFactor_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Factor	    		 = new TFile(Form("../%s/%s/Eta_MC_GammaConv_OnlyCorrectionFactor_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	
	TFile *Pi0EtaBinning_Factor_InclusiveJet = new TFile(Form("../%s/%s_InclusiveJet/Pi0EtaBinning_MC_GammaConv_OnlyCorrectionFactor_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Factor_InclusiveJet	    	 = new TFile(Form("../%s/%s_InclusiveJet/Eta_MC_GammaConv_OnlyCorrectionFactor_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	

	TH1D* histoTrueMesonEffii;

	//	TH1D* histoMesonEffi;
	//	TH1D* histoAcceptance;
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_Data;
	//	TH1D* histoCorrectedYieldJetNormalEff;
	//	TH1D* histoMesonEffi;
	//	TH1D* histoAcceptance;
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_MC;
	//	TH1D* histoCorrectedYieldJetNormalEff;
	//	TH1D* histoMesonEffi;
	//	TH1D* histoAcceptance;
	TH1D* histoCorrectedYieldNormalEff_Eta_Data;
	//	TH1D* histoCorrectedYieldJetNormalEff;
	//	TH1D* histoMesonEffi;
	//	TH1D* histoAcceptance;
	TH1D* histoCorrectedYieldNormalEff_Eta_MC;
	//	TH1D* histoCorrectedYieldJetNormalEff;

	//	Int_t BinNumber = histoCorrectedYieldNormalEff_Eta_Data->GetNBinsX();
	TH1D* histoAcceptance_Pi0EtaBinning;
	TH1D* histoAcceptance_Eta;
	TH1D* histoAcceptance_Pi0EtaBinning_InclusiveJet;
	TH1D* histoAcceptance_Eta_InclusiveJet;

	histoCorrectedYieldNormalEff_Pi0EtaBinning_Data 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0EtaBinning_MC		= (TH1D*) Pi0EtaBinning_MC->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_Data 			= (TH1D*) Eta_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_MC 			= (TH1D*) Eta_MC->Get("CorrectedYieldNormEff");
	
	histoAcceptance_Pi0EtaBinning = (TH1D*) Pi0EtaBinning_Factor -> Get("fMCMesonAccepPt");
	histoAcceptance_Eta		= (TH1D*) Eta_Factor -> Get("fMCMesonAccepPt");
	histoAcceptance_Pi0EtaBinning_InclusiveJet = (TH1D*) Pi0EtaBinning_Factor_InclusiveJet->Get("fMCMesonAccepPt");
	histoAcceptance_Eta_InclusiveJet = (TH1D*) Eta_Factor_InclusiveJet -> Get("fMCMesonAccepPt");

	///////// ratio Pi0/Eta Caculation ///////////////////	

	TFile *MyAnalysis = new TFile("MyAnalysis.root","recreate");
	
	TH1D* historatioPi0Eta_Data;	
	TH1D* historatioPi0Eta_MC;

	historatioPi0Eta_Data = (TH1D*) histoCorrectedYieldNormalEff_Eta_Data->Clone();
	historatioPi0Eta_Data->Divide(historatioPi0Eta_Data ,histoCorrectedYieldNormalEff_Pi0EtaBinning_Data,1.,1.,"");

	historatioPi0Eta_MC = (TH1D*) histoCorrectedYieldNormalEff_Eta_MC->Clone();
	historatioPi0Eta_MC->Divide(historatioPi0Eta_MC ,histoCorrectedYieldNormalEff_Pi0EtaBinning_MC,1.,1.,"");

	//	histoCorrectedYieldNormalEff_Pi0EtaBinning_Data->Draw();
	c1->cd(1);
	historatioPi0Eta_Data->SetXTitle("p_T");
	historatioPi0Eta_Data->SetYTitle("Pi0/Eta");
	TLegend *l1 = new TLegend(0.1,0.8,0.2,0.9);
	l1 -> AddEntry(historatioPi0Eta_Data,"Data","f");
	l1 -> AddEntry(historatioPi0Eta_MC,"MC","f");
	historatioPi0Eta_Data->Draw();
	historatioPi0Eta_MC -> Draw("same");
	l1 -> Draw("same");
	
	c1->Write();	
	
	////////// Acceptance in jets/ Inclusive Jet //////////////////

	TH1D* histoAcceptance_ratio_Pi0EtaBinning;
	TH1D* histoAcceptance_ratio_Eta;


	histoAcceptance_ratio_Pi0EtaBinning = (TH1D*) histoAcceptance_Pi0EtaBinning -> Clone();
	histoAcceptance_ratio_Pi0EtaBinning->Divide(histoAcceptance_ratio_Pi0EtaBinning, histoAcceptance_Pi0EtaBinning_InclusiveJet, 1., 1., "");

	histoAcceptance_ratio_Eta = (TH1D*) histoAcceptance_Eta -> Clone();
	histoAcceptance_ratio_Eta->Divide(histoAcceptance_ratio_Eta, histoAcceptance_Eta_InclusiveJet, 1., 1., "");

	histoAcceptance_ratio_Pi0EtaBinning->SetTitle("Pi0 ratio Acceptance in jets / Inclusive Jet");
	histoAcceptance_ratio_Eta->SetTitle("Eta ratio Acceptance in jet / InclusiveJet");
	histoAcceptance_ratio_Pi0EtaBinning->SetYTitle("Ratio A_Pi0");
	histoAcceptance_ratio_Eta->SetYTitle("Ratio A_Eta");
	histoAcceptance_ratio_Pi0EtaBinning->SetAxisRange(0.,2.,"Y");	
	histoAcceptance_ratio_Eta->SetAxisRange(0.,2.,"Y");	
	
	histoAcceptance_ratio_Pi0EtaBinning->Write();
	histoAcceptance_ratio_Eta->Write();

	
	///////////////////////////////////////////////////////////////
	//


	MyAnalysis->Close();	

}
