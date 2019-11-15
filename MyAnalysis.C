//////////////////////////////////////////////////////////////
///////// My Analysis ///////////////////////////////////////
/////////////////////////////////////////////////////////////

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"




void MyAnalysis(){

	////////// Basic Setting //////////////

	StyleSettingsThesis();
 	SetPlotStyle();

//	TString dateForOutput   = ReturnDateStringForOutput(); //date in format: YYYY_MM_DD
// 	TString collisionSystem = ReturnFullCollisionsSystem(); //i.e. "pp, #sqrt{s} = 7 TeV"

	TString nameDet = "ENCal";

	Color_t  colorDet = GetDefaultColorDiffDetectors(nameDet.Data(), kFALSE, kFALSE, kTRUE);
 	Style_t  markerStyleDet = GetDefaultMarkerStyleDiffDetectors(nameDet.Data(), kFALSE);
 	Size_t   markerSizeDet = GetDefaultMarkerSizeDiffDetectors(nameDet.Data(), kFALSE);

	Color_t colorEnergy       = GetColorDefaultColor("13TeV", "", "");
 	Style_t markerStyleEnergy = GetDefaultMarkerStyle("13TeV", "", "");
 	Size_t markerSizeEnergy   = GetDefaultMarkerSize("13TeV", "", "");

// 	DrawGammaSetMarkerTGraphAsym(dummgyTGrAsym , 2, styleLineJETPHOX, colorJETPHOX, colorJETPHOX, widthLinesBoxes, kTRUE, colorJETPHOX,kTRUE);	

	TCanvas* canvas_ratio_Pi0Eta = new TCanvas("canvas_ratio_Pi0Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Pi0Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_ratio_Acceptance_Pi0 = new TCanvas("canvas_ratio_Acceptance_Pi0","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Acceptance_Pi0, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Acceptance_Eta = new TCanvas("canvas_ratio_Acceptance_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Acceptance_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Efficiency_Pi0EtaBinning = new TCanvas("canvas_ratio_Efficiecny_Pi0EtaBinning","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Efficiency_Pi0EtaBinning, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Efficiency_Eta = new TCanvas("canvas_ratio_Efficiency_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Efficiency_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TF1 * f1 = new TF1 ("f1", "1", 0, 50); // constant fuction y = 1	
	
	Int_t textSizeLabelsPixel = 900*0.04;
	Int_t expectedLinesInLegend = 1;
//	TLegend* legendExample   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
// 	legendExample->AddEntry(historatioPi0Eta_Data,"Data","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...

	/////////////////////////////////////////


	TString EventCuts = "00010113";
	TString CaloCuts  = "411790607l032230000";
	TString MesonCuts = "2l631031000000d0";
	TString Sufix = "13TeV";

	TString CutSelection = Form("%s_%s_%s",EventCuts.Data(),CaloCuts.Data(),MesonCuts.Data());

	// Read File for ratio Eta / Pi0
	TFile *Pi0EtaBinning_MC 	= new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Pi0EtaBinning_Data 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_MC 			= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Data 		= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	
	TFile *Pi0EtaBinning_MC_InclusiveJet 	= new TFile(Form("../%s/%s_InclusiveJet/Pi0EtaBinning_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Pi0EtaBinning_Data_InclusiveJet 	= new TFile(Form("../%s/%s_InclusiveJet/Pi0EtaBinning_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_MC_InclusiveJet 		= new TFile(Form("../%s/%s_InclusiveJet/Eta_MC_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Data_InclusiveJet 		= new TFile(Form("../%s/%s_InclusiveJet/Eta_data_GammaConvV1Correction_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");

	// Read File for Acceptance and Efficiency in jet / InclusiveJet
	TFile *Pi0EtaBinning_Histos 		 = new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Histos	    		 = new TFile(Form("../%s/%s/Eta_MC_GammaConvV1CorrectionHistos_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	
	TFile *Pi0EtaBinning_Histos_InclusiveJet = new TFile(Form("../%s/%s_InclusiveJet/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");
	TFile *Eta_Histos_InclusiveJet	    	 = new TFile(Form("../%s/%s_InclusiveJet/Eta_MC_GammaConvV1CorrectionHistos_00010113_411790607l032230000_2l631031000000d0.root",CutSelection.Data(),Sufix.Data()),"read");

	
	// histo def
	TH1D* histoTrueMesonEffii;
	
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_Data;
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_MC;
	TH1D* histoCorrectedYieldNormalEff_Eta_Data;
	TH1D* histoCorrectedYieldNormalEff_Eta_MC;
	
	TH1D* histoAcceptance_Pi0EtaBinning;
	TH1D* histoAcceptance_Eta;
	TH1D* histoAcceptance_Pi0EtaBinning_InclusiveJet;
	TH1D* histoAcceptance_Eta_InclusiveJet;

	TH1D* histoEfficiency_Pi0EtaBinning;
	TH1D* histoEfficiency_Pi0EtaBinning_True;
	TH1D* histoEfficiency_Eta;
	TH1D* histoEfficiency_Eta_True;
	TH1D* histoEfficiency_Pi0EtaBinning_InclusiveJet;
	TH1D* histoEfficiency_Pi0EtaBinning_True_InclusiveJet;
	TH1D* histoEfficiency_Eta_InclusiveJet;
	TH1D* histoEfficiency_Eta_True_InclusiveJet;
	

	// Read histo for ratio Eta / Pi0
	histoCorrectedYieldNormalEff_Pi0EtaBinning_Data 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0EtaBinning_MC		= (TH1D*) Pi0EtaBinning_MC->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_Data 			= (TH1D*) Eta_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_MC 			= (TH1D*) Eta_MC->Get("CorrectedYieldNormEff");
	
	// Read histo for Acceptance in jet / Inclusive Jet
	histoAcceptance_Pi0EtaBinning 				= (TH1D*) Pi0EtaBinning_Histos->Get("fMCMesonAccepPt");
	histoAcceptance_Eta					= (TH1D*) Eta_Histos->Get("fMCMesonAccepPt");
	histoAcceptance_Pi0EtaBinning_InclusiveJet 		= (TH1D*) Pi0EtaBinning_Histos_InclusiveJet->Get("fMCMesonAccepPt");
	histoAcceptance_Eta_InclusiveJet 			= (TH1D*) Eta_Histos_InclusiveJet->Get("fMCMesonAccepPt");

	// Read histo for Efficiency in jet / Inclusive Jet
	histoEfficiency_Pi0EtaBinning				= (TH1D*) Pi0EtaBinning_Histos->Get("MesonEffiPt");
	histoEfficiency_Pi0EtaBinning_True			= (TH1D*) Pi0EtaBinning_Histos->Get("TrueMesonEffiPt");
	histoEfficiency_Eta					= (TH1D*) Eta_Histos->Get("MesonEffiPt");
	histoEfficiency_Eta_True				= (TH1D*) Eta_Histos->Get("TrueMesonEffiPt");
	histoEfficiency_Pi0EtaBinning_InclusiveJet		= (TH1D*) Pi0EtaBinning_Histos_InclusiveJet->Get("MesonEffiPt");
	histoEfficiency_Pi0EtaBinning_True_InclusiveJet		= (TH1D*) Pi0EtaBinning_Histos_InclusiveJet->Get("TrueMesonEffiPt");
	histoEfficiency_Eta_InclusiveJet			= (TH1D*) Eta_Histos_InclusiveJet->Get("MesonEffiPt");
	histoEfficiency_Eta_True_InclusiveJet			= (TH1D*) Eta_Histos_InclusiveJet->Get("TrueMesonEffiPt");


	///////// ratio Pi0/Eta Caculation ///////////////////	

	TFile *MyAnalysis = new TFile("MyAnalysis.root","recreate");
	
	TH1D* historatioPi0Eta_Data;	
	TH1D* historatioPi0Eta_MC;

	historatioPi0Eta_Data = (TH1D*) histoCorrectedYieldNormalEff_Eta_Data->Clone();
	historatioPi0Eta_Data->Divide(historatioPi0Eta_Data ,histoCorrectedYieldNormalEff_Pi0EtaBinning_Data,1.,1.,"");

	historatioPi0Eta_MC = (TH1D*) histoCorrectedYieldNormalEff_Eta_MC->Clone();
	historatioPi0Eta_MC->Divide(historatioPi0Eta_MC ,histoCorrectedYieldNormalEff_Pi0EtaBinning_MC,1.,1.,"");

	//	histoCorrectedYieldNormalEff_Pi0EtaBinning_Data->Draw();
	canvas_ratio_Pi0Eta->cd();
	
	DrawGammaSetMarker(historatioPi0Eta_Data, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(historatioPi0Eta_MC, 20, 1., kRed, kRed);
	historatioPi0Eta_Data->SetXTitle("p_{T}");
	historatioPi0Eta_Data->SetYTitle("#frac{#eta}{#pi_{0}}");
 	
	TLegend* legendratioPi0Eta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendratioPi0Eta->AddEntry(historatioPi0Eta_Data,"Data","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta->AddEntry(historatioPi0Eta_MC,"MC","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...

	historatioPi0Eta_Data->Draw();
	historatioPi0Eta_MC  ->Draw("same");
 	legendratioPi0Eta    ->Draw("same");
	
	canvas_ratio_Pi0Eta->Write();	
	
	////////// Acceptance in jets/ Inclusive Jet //////////////////

	TH1D* histoAcceptance_ratio_Pi0EtaBinning;
	TH1D* histoAcceptance_ratio_Eta;

	histoAcceptance_ratio_Pi0EtaBinning = (TH1D*) histoAcceptance_Pi0EtaBinning -> Clone();
	histoAcceptance_ratio_Pi0EtaBinning->Divide(histoAcceptance_ratio_Pi0EtaBinning, histoAcceptance_Pi0EtaBinning_InclusiveJet, 1., 1., "");

	histoAcceptance_ratio_Eta = (TH1D*) histoAcceptance_Eta -> Clone();
	histoAcceptance_ratio_Eta->Divide(histoAcceptance_ratio_Eta, histoAcceptance_Eta_InclusiveJet, 1., 1., "");
	
	histoAcceptance_ratio_Pi0EtaBinning->SetTitle("Pi0 ratio Acceptance in jets / Inclusive Jet");
	histoAcceptance_ratio_Eta->SetTitle("Eta ratio Acceptance in jet / InclusiveJet");
	histoAcceptance_ratio_Pi0EtaBinning->SetYTitle("Ratio A_{#pi_{0}}");
	histoAcceptance_ratio_Eta->SetYTitle("Ratio A_{#eta}");
	histoAcceptance_ratio_Pi0EtaBinning->SetAxisRange(0.,2.,"Y");	
	histoAcceptance_ratio_Eta->SetAxisRange(0.,2.,"Y");	
	DrawGammaSetMarker(histoAcceptance_ratio_Pi0EtaBinning, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoAcceptance_ratio_Eta, 20, 1., kBlack, kBlack);
	
	canvas_ratio_Acceptance_Pi0->cd();
	
	histoAcceptance_ratio_Pi0EtaBinning->Draw();
	f1->Draw("same");	

	canvas_ratio_Acceptance_Pi0->Write();

	canvas_ratio_Acceptance_Eta->cd();

	histoAcceptance_ratio_Eta->Draw();
	f1->Draw("same");

	canvas_ratio_Acceptance_Eta->Write();

	histoAcceptance_ratio_Pi0EtaBinning->Write();
	histoAcceptance_ratio_Eta->Write();

	/////////// Efficiency in jets / Inclusive Jet ///////////////
	

	TH1D* histoEfficiency_ratio_Pi0EtaBinning;
	TH1D* histoEfficiency_ratio_Pi0EtaBinning_True;
	TH1D* histoEfficiency_ratio_Eta;
	TH1D* histoEfficiency_ratio_Eta_True;
	
	canvas_ratio_Efficiency_Pi0EtaBinning->cd(); // Draw the Efficiency ratio_Pi0 in jet / Inclusive jet in canvas

	histoEfficiency_ratio_Pi0EtaBinning = (TH1D*) histoEfficiency_Pi0EtaBinning -> Clone();
	histoEfficiency_ratio_Pi0EtaBinning->Divide(histoEfficiency_ratio_Pi0EtaBinning, histoEfficiency_Pi0EtaBinning_InclusiveJet, 1., 1., "");

	histoEfficiency_ratio_Pi0EtaBinning_True = (TH1D*) histoEfficiency_Pi0EtaBinning_True -> Clone();
	histoEfficiency_ratio_Pi0EtaBinning_True->Divide(histoEfficiency_ratio_Pi0EtaBinning_True, histoEfficiency_Pi0EtaBinning_True_InclusiveJet, 1., 1., "");

	DrawGammaSetMarker(histoEfficiency_ratio_Pi0EtaBinning, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoEfficiency_ratio_Pi0EtaBinning_True, 20, 1., kGreen, kGreen);
	histoEfficiency_ratio_Pi0EtaBinning->SetAxisRange(0.,2.,"Y");
	histoEfficiency_ratio_Pi0EtaBinning->SetYTitle("Ratio #epsilon_{#pi_{0}}");
	histoEfficiency_ratio_Pi0EtaBinning->SetXTitle("p_{T}");
	
	TLegend* legendratioEfficiency_Pi0EtaBinning   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendratioEfficiency_Pi0EtaBinning->AddEntry(histoEfficiency_ratio_Pi0EtaBinning, "#frac{Rec #epsilon_{jet}}{Rec #epsilon_{MB}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioEfficiency_Pi0EtaBinning->AddEntry(histoEfficiency_ratio_Pi0EtaBinning_True, "#frac{True #epsilon_{Jet}}{True #epsilon_{MB}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	
	histoEfficiency_ratio_Pi0EtaBinning->Draw();
	histoEfficiency_ratio_Pi0EtaBinning_True->Draw("same");	
	legendratioEfficiency_Pi0EtaBinning->Draw("same");
	f1->Draw("same");
	canvas_ratio_Efficiency_Pi0EtaBinning->Write();

	canvas_ratio_Efficiency_Eta->cd(); // Draw the Efficiency ratio_Eta in jet / Inclusive Jet in canvas
	
	histoEfficiency_ratio_Eta = (TH1D*) histoEfficiency_Eta -> Clone();
	histoEfficiency_ratio_Eta->Divide(histoEfficiency_ratio_Eta, histoEfficiency_Eta_InclusiveJet, 1., 1., "");

	histoEfficiency_ratio_Eta_True = (TH1D*) histoEfficiency_Eta_True -> Clone();
	histoEfficiency_ratio_Eta_True->Divide(histoEfficiency_ratio_Eta_True, histoEfficiency_Eta_True_InclusiveJet, 1., 1., "");
	
	DrawGammaSetMarker(histoEfficiency_ratio_Eta, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoEfficiency_ratio_Eta_True, 20, 1., kGreen, kGreen);
	histoEfficiency_ratio_Eta->SetAxisRange(0.,2.,"Y");
	histoEfficiency_ratio_Eta->SetYTitle("Ratio #epsilon_{#eta}");
	histoEfficiency_ratio_Eta->SetXTitle("p_{T}");

	TLegend* legendratioEfficiency_Eta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendratioEfficiency_Eta->AddEntry(histoEfficiency_ratio_Eta, "#frac{Rec #epsilon_{jet}}{Rec #epsilon_{MB}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioEfficiency_Eta->AddEntry(histoEfficiency_ratio_Eta_True, "#frac{Rec #epsilon_{jet}}{Rec #epsilon_{MB}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	histoEfficiency_ratio_Eta->Draw();
	histoEfficiency_ratio_Eta_True->Draw("same");
	legendratioEfficiency_Eta->Draw("same");
	f1->Draw("same");
	canvas_ratio_Efficiency_Eta->Write();

	histoEfficiency_ratio_Pi0EtaBinning->Write();
	histoEfficiency_ratio_Pi0EtaBinning_True->Write();
	histoEfficiency_ratio_Eta->Write();
	histoEfficiency_ratio_Eta_True->Write();

	
	/////// Save the Canvas file /////////////////////////////////
	//
	canvas_ratio_Pi0Eta->SaveAs("RatioPi0Eta.eps");
	canvas_ratio_Acceptance_Pi0->SaveAs("RatioAcceptancePi0.eps");
	canvas_ratio_Acceptance_Eta->SaveAs("RatioAcceptanceEta.eps");
	canvas_ratio_Efficiency_Pi0EtaBinning->SaveAs("RatioEfficiencyPi0.eps");
	canvas_ratio_Efficiency_Eta->SaveAs("RatioEfficiencyEta.eps");

	//////////////////////////////////////////////////////////////
	MyAnalysis->Close();	

}
