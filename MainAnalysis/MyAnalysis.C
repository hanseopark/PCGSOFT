//////////////////////////////////////////////////////////////
///////// My Analysis ///////////////////////////////////////
/////////////////////////////////////////////////////////////

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"
//#include "MyCommonHeaders/Filipad.h"



void MyAnalysis(){

    const Int_t MaxNumberOfFiles    = 13;
	////////// Basic Setting //////////////
	TLatex *t = new TLatex();
	t->SetTextAlign(12);
	t->SetTextSize(0.05);
	
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

	TCanvas* ctest = new TCanvas("ctest","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(ctest, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Pi0Eta = new TCanvas("canvas_ratio_Pi0Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Pi0Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_CorrectedYield_Pi0 = new TCanvas("canvas_CorrectedYield_Pi0","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_CorrectedYield_Pi0, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_CorrectedYield_Pi0_MC = new TCanvas("canvas_CorrectedYield_Pi0_MC","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_CorrectedYield_Pi0_MC, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_CorrectedYield_Eta = new TCanvas("canvas_CorrectedYield_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_CorrectedYield_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Pi0Eta_trigger = new TCanvas("canvas_ratio_Pi0Eta_trigger","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Pi0Eta_trigger, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_ratio_diffright = new TCanvas("canvas_ratio_diffrigt","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_diffright, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)


//	TCanvas* canvas_ratio_DiffDefault_Pi0[6] = nullptr;
//
	TCanvas* canvas_ratio_DiffDefault_Pi0 = new TCanvas("canvas_ratio_DiffDefault_Pi0","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_DiffDefault_Pi0, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Acceptance_Pi0 = new TCanvas("canvas_ratio_Acceptance_Pi0","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Acceptance_Pi0, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Acceptance_Eta = new TCanvas("canvas_ratio_Acceptance_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Acceptance_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Efficiency_Pi0EtaBinning = new TCanvas("canvas_ratio_Efficiecny_Pi0EtaBinning","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Efficiency_Pi0EtaBinning, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_ratio_Efficiency_Eta = new TCanvas("canvas_ratio_Efficiency_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_ratio_Efficiency_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_Acceptance = new TCanvas("canvas_Acceptance","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_Acceptance, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_Efficiency = new TCanvas("canvas_Efficiency","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_Efficiency, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_trigger_Pi0Yield = new TCanvas("canvas_trigger_Pi0Yield","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_trigger_Pi0Yield, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_trigger_Pi0Yield_Inclusive = new TCanvas("canvas_trigger_Pi0Yield_Inclusive","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_trigger_Pi0Yield_Inclusive, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_trigger_EtaYield = new TCanvas("canvas_trigger_EtaYield","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_trigger_EtaYield, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_triggerClusterData = new TCanvas("canvas_triggerClusterData","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_triggerClusterData, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_triggerClusterData_Inclusive = new TCanvas("canvas_triggerClusterData_Inclusive","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_triggerClusterData_Inclusive, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_triggerClusterMC = new TCanvas("canvas_triggerClusterMC","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_triggerClusterMC, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_trigger_Eta = new TCanvas("canvas_trigger_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_trigger_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	
	TCanvas* canvas_triggerPi0YieldDataRatio = new TCanvas("canvas_triggerPi0YieldDataRatio","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_triggerPi0YieldDataRatio, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_triggerPi0YieldDataRatio_Inclusive = new TCanvas("canvas_triggerPi0YieldDataRatio_Inclusive","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_triggerPi0YieldDataRatio_Inclusive, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TCanvas* canvas_triggerPi0YieldMCRatio = new TCanvas("canvas_triggerPi0YieldMCRatio","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_triggerPi0YieldMCRatio, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TF1 * f1 = new TF1 ("f1", "1", 0, 50); // constant fuction y = 1	
	
	Int_t textSizeLabelsPixel = 900*0.04;
	Int_t expectedLinesInLegend = 1;
//	TLegend* legendExample   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
// 	legendExample->AddEntry(historatioPi0Eta_Data,"Data","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...

	/////////////////////////////////////////


	TString EventCuts = "00010113";
	TString EventCuts_EG1 = "0008d113";
	TString EventCuts_EG2 = "0008e113";
	TString CaloCuts  = "411790607l032230000";
	TString MesonCuts = "2l631031000000d0";
	TString CaloCuts_Inclusive = "411790607l032230000";
	TString CaloCuts_Pi0_Inclusive = "411791206f032230000";
	//TString CaloCuts_Eta_Inclusive = "411791106f032230000";
	TString MesonCuts_Inclusive = "0l631031000000d0";
	TString MesonCuts_Inclusive_trigger = "01631031000000d0";
	TString Energy = "13TeV";

	TString CutSelection = Form("%s_%s_%s",EventCuts.Data(),CaloCuts.Data(),MesonCuts.Data());
	TString CutSelection_EG1 = Form("%s_%s_%s",EventCuts_EG1.Data(),CaloCuts.Data(),MesonCuts.Data());
	TString CutSelection_EG2 = Form("%s_%s_%s",EventCuts_EG2.Data(),CaloCuts.Data(),MesonCuts.Data());
	TString CutSelection_Inclusive = Form("%s_%s_%s",EventCuts.Data(),CaloCuts_Inclusive.Data(),MesonCuts_Inclusive.Data());
	TString CutSelection_INT7_Inclusive = Form("%s_%s_%s",EventCuts.Data(),CaloCuts_Pi0_Inclusive.Data(),MesonCuts_Inclusive_trigger.Data());
	TString CutSelection_EG1_Inclusive = Form("%s_%s_%s",EventCuts_EG1.Data(),CaloCuts_Pi0_Inclusive.Data(),MesonCuts_Inclusive_trigger.Data());
	TString CutSelection_EG2_Inclusive = Form("%s_%s_%s",EventCuts_EG2.Data(),CaloCuts_Pi0_Inclusive.Data(),MesonCuts_Inclusive_trigger.Data());
	//TString CutSelection_Eta_Inclusive = Form("%s_%s_%s",EventCuts.Data(),CaloCuts_Eta_Inclusive.Data(),MesonCuts_Inclusive.Data());
	
	// Read File for compare cluster energy with each trigger.
	TFile *GammaCalo_912			= new TFile("/home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_912.root");
	TFile *GammaCalo_950			= new TFile("/home/alidock/alice/work/Data/pp_13TeV/Injet/DataMerged/GammaCalo_950.root");
	TFile *GammaCalo_950_MC			= new TFile("/home/alidock/alice/work/Data/pp_13TeV/Injet/MCMerged/GammaCalo_950.root");

	// Read File for ratio Eta / pi0 included preliminary datas which are in jet and Inclusive at 5TeV.
	TFile *PreliminaryData		= new TFile("/home/alidock/alice/work/Data/CombinedResultsPaperPP5TeV_2019_12_11.root");
	
	// Read File for ratio Eta / Pi0.
	TFile *Pi0EtaBinning_MC 	= new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1Correction_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Pi0EtaBinning_Data 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1Correction_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Eta_MC 			= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1Correction_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Eta_Data 		= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	
	TFile *Pi0EtaBinning_MC_Inclusive 	= new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1Correction_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");
	TFile *Pi0EtaBinning_Data_Inclusive 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1Correction_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");
	TFile *Eta_MC_Inclusive 		= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1Correction_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");
	TFile *Eta_Data_Inclusive 		= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");

	// Read File for Acceptance and Efficiency in jet / Inclusive
	TFile *Pi0_Histos 		 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1CorrectionHistos_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Pi0_Histos_Inclusive 			= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1CorrectionHistos_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");
	TFile *Pi0EtaBinning_Histos             = new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Eta_Histos	    		= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1CorrectionHistos_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	
	TFile *Pi0EtaBinning_Histos_Inclusive = new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");
	TFile *Eta_Histos_Inclusive	    	 = new TFile(Form("../%s/%s/Eta_MC_GammaConvV1CorrectionHistos_%s.root",CutSelection_Inclusive.Data(),Energy.Data(),CutSelection_Inclusive.Data()),"read");

	// Read File to compare each trigger
	TFile *Pi0_Data 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1Correction_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Pi0_Data_EG1 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1Correction_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0_Data_EG2 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1Correction_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Eta_Data_EG1 	= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Eta_Data_EG2 	= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Pi0_MC 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Pi0_MC_EG1 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0_MC_EG2 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Eta_MC_EG1 	= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1Correction_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Eta_MC_EG2 	= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1Correction_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Pi0EtaBinning_Data_EG1 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1Correction_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0EtaBinning_Data_EG2 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1Correction_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");

	TFile *Pi0_Data_Inclusive		= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1Correction_%s.root",CutSelection_INT7_Inclusive.Data(),Energy.Data(),CutSelection_INT7_Inclusive.Data()),"read");
	TFile *Pi0_Data_EG1_Inclusive 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1Correction_%s.root",CutSelection_EG1_Inclusive.Data(),Energy.Data(),CutSelection_EG1_Inclusive.Data()),"read");
	TFile *Pi0_Data_EG2_Inclusive 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1Correction_%s.root",CutSelection_EG2_Inclusive.Data(),Energy.Data(),CutSelection_EG2_Inclusive.Data()),"read");
//	TFile *Eta_Data_EG1_Inclusive 	= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_%s.root",CutSelection_EG1_Inclusive.Data(),Energy.Data(),CutSelection_EG1_Inclusive.Data()),"read");
//	TFile *Eta_Data_EG2_Inclusive 	= new TFile(Form("../%s/%s/Eta_data_GammaConvV1Correction_%s.root",CutSelection_EG2_Inclusive.Data(),Energy.Data(),CutSelection_EG2_Inclusive.Data()),"read");
	// Read Uncorrected File to compare each trigger
	TFile *Pi0_UncorrectedData 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Pi0_UncorrectedData_EG1 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0_UncorrectedData_EG2 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Eta_UncorrectedData_EG1 	= new TFile(Form("../%s/%s/Eta_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Eta_UncorrectedData_EG2 	= new TFile(Form("../%s/%s/Eta_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Pi0EtaBinning_UncorrectedData_EG1 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0EtaBinning_UncorrectedData_EG2 	= new TFile(Form("../%s/%s/Pi0EtaBinning_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	
	TFile *Pi0_UncorrectedData_Inclusive		= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_INT7_Inclusive.Data(),Energy.Data(),CutSelection_INT7_Inclusive.Data()),"read");
	TFile *Pi0_UncorrectedData_EG1_Inclusive 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1_Inclusive.Data(),Energy.Data(),CutSelection_EG1_Inclusive.Data()),"read");
	TFile *Pi0_UncorrectedData_EG2_Inclusive 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2_Inclusive.Data(),Energy.Data(),CutSelection_EG2_Inclusive.Data()),"read");

	TFile *Pi0_UncorrectedMC 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection.Data(),Energy.Data(),CutSelection.Data()),"read");
	TFile *Pi0_UncorrectedMC_EG1 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0_UncorrectedMC_EG2 	= new TFile(Form("../%s/%s/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Eta_UncorrectedMC_EG1 	= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Eta_UncorrectedMC_EG2 	= new TFile(Form("../%s/%s/Eta_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");
	TFile *Pi0EtaBinning_UncorrectedMC_EG1 	= new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG1.Data(),Energy.Data(),CutSelection_EG1.Data()),"read");
	TFile *Pi0EtaBinning_UncorrectedMC_EG2 	= new TFile(Form("../%s/%s/Pi0EtaBinning_MC_GammaConvV1WithoutCorrection_%s.root",CutSelection_EG2.Data(),Energy.Data(),CutSelection_EG2.Data()),"read");

	// histo def
	TH1D* histoCorrectedYieldNormalEff_Right_Normal_Pi0EtaBinning;
	TH1D* histoCorrectedYieldNormalEff_Right_Narrow_Pi0EtaBinning;
	TH1D* histoCorrectedYieldNormalEff_Right_Wide_Pi0EtaBinning;
	TH1D* histoCorrectedYieldNormalEff_Left_Normal_Pi0EtaBinning;
	TH1D* histoCorrectedYieldNormalEff_Left_Narrow_Pi0EtaBinning;
	TH1D* histoCorrectedYieldNormalEff_Left_Wide_Pi0EtaBinning;

	TH1D* histoTrueMesonEffii;
	
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_Data;
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_MC;
	TH1D* histoCorrectedYieldNormalEff_Eta_INT7_Data;
	TH1D* histoCorrectedYieldNormalEff_Eta_INT7_MC;
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_EG1_Data;
	TH1D* histoCorrectedYieldNormalEff_Eta_EG1_Data;
	TH1D* histoCorrectedYieldNormalEff_Pi0EtaBinning_EG2_Data;
	TH1D* histoCorrectedYieldNormalEff_Eta_EG2_Data;
	TH1D* histoCorrectedYieldNormalEff_Pi0_INT7_Data;
	TH1D* histoCorrectedYieldNormalEff_Pi0_EG1_Data;
	TH1D* histoCorrectedYieldNormalEff_Pi0_EG2_Data;
	
	TH1D *histoCorrectedYieldNormalEff_Narrow_Pi0EtaBinning_Data;
	TH1D *histoCorrectedYieldNormalEff_Wide_Pi0EtaBinning_Data;
	TH1D *histoCorrectedYieldNormalEff_Narrow_Eta_Data;
	TH1D *histoCorrectedYieldNormalEff_Wide_Eta_Data;
	TH1D *histoCorrectedYieldNormalEff_Narrow_Pi0EtaBinning_MC;
	TH1D *histoCorrectedYieldNormalEff_Wide_Pi0EtaBinning_MC;
	TH1D *histoCorrectedYieldNormalEff_Narrow_Eta_MC;
	TH1D *histoCorrectedYieldNormalEff_Wide_Eta_MC;

	TH1D* histoAcceptance_Pi0EtaBinning;
	TH1D* histoAcceptance_Eta;
	TH1D* histoAcceptance_Pi0EtaBinning_Inclusive;
	TH1D* histoAcceptance_Eta_Inclusive;

	TH1D* histoEfficiency_Pi0;
	TH1D* histoEfficiency_Pi0_Inclusive;
	TH1D* histoEfficiency_Pi0EtaBinning;
	TH1D* histoEfficiency_Pi0EtaBinning_True;
	TH1D* histoEfficiency_Eta;
	TH1D* histoEfficiency_Eta_True;
	TH1D* histoEfficiency_Pi0EtaBinning_Inclusive;
	TH1D* histoEfficiency_Pi0EtaBinning_True_Inclusive;
	TH1D* histoEfficiency_Eta_Inclusive;
	TH1D* histoEfficiency_Eta_True_Inclusive;
		
	TH1D* histoYieldMeson_Pi0_Data_INT7;
	TH1D* histoYieldMeson_Pi0_Data_EG1;
	TH1D* histoYieldMeson_Pi0_Data_EG2;
	
		
	TH1D* histoYieldMeson_Pi0_Data_INT7_Inclusive;
	TH1D* histoYieldMeson_Pi0_Data_EG1_Inclusive;
	TH1D* histoYieldMeson_Pi0_Data_EG2_Inclusive;

	TH1D* histoYieldMeson_Eta_Data_INT7;
	TH1D* histoYieldMeson_Eta_Data_EG1;
	TH1D* histoYieldMeson_Eta_Data_EG2;
	
	TH1D* histoYieldMeson_Pi0_MC_INT7;
	TH1D* histoYieldMeson_Pi0_MC_EG1;
	TH1D* histoYieldMeson_Pi0_MC_EG2;
	
	TH1D* histoYieldMeson_Eta_MC_INT7;
	TH1D* histoYieldMeson_Eta_MC_EG1;
	TH1D* histoYieldMeson_Eta_MC_EG2;


	TH1D* histoYieldMeson_Pi0EtaBinning_Data_EG1;
	TH1D* histoYieldMeson_Pi0EtaBinning;

	TH1D* histoClusterEnergyINT7;
	TH1D* histoClusterEnergyEG1;
	TH1D* histoClusterEnergyEG2;

	TH1D* histoClusterEnergyINT7_MC;
	TH1D* histoClusterEnergyEG1_MC;
	TH1D* histoClusterEnergyEG2_MC;

	TH1D* histoNEventsINT7_Data;
	TH1D* histoNEventsEG1_Data;
	TH1D* histoNEventsEG2_Data;
	
	TH1D* histoNEventsINT7_MC;
	TH1D* histoNEventsEG1_MC;
	TH1D* histoNEventsEG2_MC;

	// Read histo for ratio Eta / Pi0
	histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_Data 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_MC	= (TH1D*) Pi0EtaBinning_MC->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_INT7_Data 		= (TH1D*) Eta_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_INT7_MC 		= (TH1D*) Eta_MC->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0EtaBinning_EG1_Data 	= (TH1D*) Pi0EtaBinning_Data_EG1->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_EG1_Data 		= (TH1D*) Eta_Data_EG1->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0EtaBinning_EG2_Data 	= (TH1D*) Pi0EtaBinning_Data_EG2->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Eta_EG2_Data 		= (TH1D*) Eta_Data_EG2->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0_INT7_Data 	= (TH1D*) Pi0_Data->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0_EG1_Data 	= (TH1D*) Pi0_Data_EG1->Get("CorrectedYieldNormEff");
	histoCorrectedYieldNormalEff_Pi0_EG2_Data 	= (TH1D*) Pi0_Data_EG2->Get("CorrectedYieldNormEff");
	
	// Read histo for narrow/default wide/default
	
	histoCorrectedYieldNormalEff_Narrow_Pi0EtaBinning_Data =  (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffNarrow");
	histoCorrectedYieldNormalEff_Wide_Pi0EtaBinning_Data   =  (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffWide");
	histoCorrectedYieldNormalEff_Narrow_Eta_Data	          =  (TH1D*) Eta_Data->Get("CorrectedYieldNormEffNarrow");
	histoCorrectedYieldNormalEff_Wide_Eta_Data	          =  (TH1D*) Eta_Data->Get("CorrectedYieldNormEffWide");
	histoCorrectedYieldNormalEff_Narrow_Pi0EtaBinning_MC   =  (TH1D*) Pi0EtaBinning_MC->Get("CorrectedYieldNormEffNarrow");
	histoCorrectedYieldNormalEff_Wide_Pi0EtaBinning_MC     =  (TH1D*) Pi0EtaBinning_MC->Get("CorrectedYieldNormEffWide");
	histoCorrectedYieldNormalEff_Narrow_Eta_MC	          =  (TH1D*) Eta_MC->Get("CorrectedYieldNormEffNarrow");
	histoCorrectedYieldNormalEff_Wide_Eta_MC	          =  (TH1D*) Eta_MC->Get("CorrectedYieldNormEffWide");

	// Read histo for ratio of normal, narrow, wide / right and normal,narrow, wide / left
	
	histoCorrectedYieldNormalEff_Right_Normal_Pi0EtaBinning = (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEff");	
	histoCorrectedYieldNormalEff_Right_Narrow_Pi0EtaBinning = (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffNarrow");	
	histoCorrectedYieldNormalEff_Right_Wide_Pi0EtaBinning 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffWide");	
	histoCorrectedYieldNormalEff_Left_Normal_Pi0EtaBinning 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffLeft");	
	histoCorrectedYieldNormalEff_Left_Narrow_Pi0EtaBinning 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffLeftNarrow");	
	histoCorrectedYieldNormalEff_Left_Wide_Pi0EtaBinning 	= (TH1D*) Pi0EtaBinning_Data->Get("CorrectedYieldNormEffLeftWide");	

	// Read histo for Acceptance in jet / Inclusive Jet
	histoAcceptance_Pi0EtaBinning 			= (TH1D*) Pi0EtaBinning_Histos->Get("fMCMesonAccepPt");
	histoAcceptance_Eta				= (TH1D*) Eta_Histos->Get("fMCMesonAccepPt");
	histoAcceptance_Pi0EtaBinning_Inclusive 		= (TH1D*) Pi0EtaBinning_Histos_Inclusive->Get("fMCMesonAccepPt");
	histoAcceptance_Eta_Inclusive 				= (TH1D*) Eta_Histos_Inclusive->Get("fMCMesonAccepPt");

	// Read histo for Efficiency in jet / Inclusive Jet
	histoEfficiency_Pi0				= (TH1D*) Pi0_Histos->Get("MesonEffiPt");
	histoEfficiency_Pi0_Inclusive				= (TH1D*) Pi0_Histos_Inclusive->Get("MesonEffiPt");
	histoEfficiency_Pi0EtaBinning			= (TH1D*) Pi0EtaBinning_Histos->Get("MesonEffiPt");
	histoEfficiency_Pi0EtaBinning_True		= (TH1D*) Pi0EtaBinning_Histos->Get("TrueMesonEffiPt");
	histoEfficiency_Eta				= (TH1D*) Eta_Histos->Get("MesonEffiPt");
	histoEfficiency_Eta_True			= (TH1D*) Eta_Histos->Get("TrueMesonEffiPt");
	histoEfficiency_Pi0EtaBinning_Inclusive		= (TH1D*) Pi0EtaBinning_Histos_Inclusive->Get("MesonEffiPt");
	histoEfficiency_Pi0EtaBinning_True_Inclusive		= (TH1D*) Pi0EtaBinning_Histos_Inclusive->Get("TrueMesonEffiPt");
	histoEfficiency_Eta_Inclusive				= (TH1D*) Eta_Histos_Inclusive->Get("MesonEffiPt");
	histoEfficiency_Eta_True_Inclusive			= (TH1D*) Eta_Histos_Inclusive->Get("TrueMesonEffiPt");
	
	// Read histo to compare each trigger
	histoYieldMeson_Pi0_Data_INT7 	= (TH1D*) Pi0_Data->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0_Data_EG1 	= (TH1D*) Pi0_Data_EG1->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0_Data_EG2 	= (TH1D*) Pi0_Data_EG2->Get("histoYieldMesonPerEvent");

	histoYieldMeson_Pi0_Data_INT7_Inclusive 	= (TH1D*) Pi0_Data_Inclusive->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0_Data_EG1_Inclusive		= (TH1D*) Pi0_Data_EG1_Inclusive->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0_Data_EG2_Inclusive	 	= (TH1D*) Pi0_Data_EG2_Inclusive->Get("histoYieldMesonPerEvent");

	histoYieldMeson_Eta_Data_INT7 	= (TH1D*) Eta_Data->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Eta_Data_EG1 	= (TH1D*) Eta_Data_EG1->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Eta_Data_EG2 	= (TH1D*) Eta_Data_EG2->Get("histoYieldMesonPerEvent");

//	histoYieldMeson_Eta_Data_INT7_Inclusive 	= (TH1D*) Eta_Data_Inclusive->Get("histoYieldMesonPerEvent");
//	histoYieldMeson_Eta_Data_EG1_Inclusive		= (TH1D*) Eta_Data_EG1_Inclusive->Get("histoYieldMesonPerEvent");
//	histoYieldMeson_Eta_Data_EG2_Inclusive	 	= (TH1D*) Eta_Data_EG2_Inclusive->Get("histoYieldMesonPerEvent");

	histoYieldMeson_Pi0_MC_INT7 	= (TH1D*) Pi0_MC->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0_MC_EG1 	= (TH1D*) Pi0_MC_EG1->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0_MC_EG2 	= (TH1D*) Pi0_MC_EG2->Get("histoYieldMesonPerEvent");

	histoYieldMeson_Eta_MC_INT7 	= (TH1D*) Eta_MC->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Eta_MC_EG1 	= (TH1D*) Eta_MC_EG1->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Eta_MC_EG2 	= (TH1D*) Eta_MC_EG2->Get("histoYieldMesonPerEvent");

	histoYieldMeson_Pi0EtaBinning_Data_EG1 	= (TH1D*) Pi0EtaBinning_Data_EG1->Get("histoYieldMesonPerEvent");
	histoYieldMeson_Pi0EtaBinning 	= (TH1D*) Pi0EtaBinning_Data_EG2->Get("histoYieldMesonPerEvent");
	
	TList *TopDir_950                   = (TList*) GammaCalo_950->Get("GammaCalo_950");
        TList *HistosGammaConversion_INT7   = (TList*)TopDir_950->FindObject(Form("Cut Number %s",CutSelection.Data()));
    	TList *ESDContainer_INT7                 = (TList*) HistosGammaConversion_INT7->FindObject(Form("%s ESD histograms",CutSelection.Data()));
	
	TList *TopDir_912                  = (TList*) GammaCalo_950->Get("GammaCalo_950");
        TList *HistosGammaConversion_EG1   = (TList*)TopDir_912->FindObject(Form("Cut Number %s",CutSelection_EG1.Data()));
    	TList *ESDContainer_EG1                 = (TList*) HistosGammaConversion_EG1->FindObject(Form("%s ESD histograms",CutSelection_EG1.Data()));
	
        TList *HistosGammaConversion_EG2   = (TList*)TopDir_912->FindObject(Form("Cut Number %s",CutSelection_EG2.Data()));
    	TList *ESDContainer_EG2                 = (TList*) HistosGammaConversion_EG2->FindObject(Form("%s ESD histograms",CutSelection_EG2.Data()));
	histoClusterEnergyINT7				= (TH1D*) ESDContainer_INT7->FindObject("ClusGamma_E");
	histoClusterEnergyEG1				= (TH1D*) ESDContainer_EG1->FindObject("ClusGamma_E");
	histoClusterEnergyEG2				= (TH1D*) ESDContainer_EG2->FindObject("ClusGamma_E");
	histoNEventsINT7_Data				= (TH1D*) ESDContainer_INT7->FindObject("NEvents");
	histoNEventsEG1_Data				= (TH1D*) ESDContainer_EG1->FindObject("NEvents");
	histoNEventsEG2_Data				= (TH1D*) ESDContainer_EG2->FindObject("NEvents");
		
	TList *TopDir_950_MC                   = (TList*) GammaCalo_950_MC->Get("GammaCalo_950");
        TList *HistosGammaConversion_INT7_MC   = (TList*)TopDir_950_MC->FindObject(Form("Cut Number %s",CutSelection.Data()));
        TList *HistosGammaConversion_EG1_MC   = (TList*)TopDir_950_MC->FindObject(Form("Cut Number %s",CutSelection_EG1.Data()));
        TList *HistosGammaConversion_EG2_MC   = (TList*)TopDir_950_MC->FindObject(Form("Cut Number %s",CutSelection_EG2.Data()));
			TList *ESDContainer_INT7_MC                = (TList*) HistosGammaConversion_INT7_MC->FindObject(Form("%s ESD histograms",CutSelection.Data()));
			TList *ESDContainer_EG1_MC                 = (TList*) HistosGammaConversion_EG1_MC->FindObject(Form("%s ESD histograms",CutSelection_EG1.Data()));
			TList *ESDContainer_EG2_MC                 = (TList*) HistosGammaConversion_EG2_MC->FindObject(Form("%s ESD histograms",CutSelection_EG2.Data()));
	histoClusterEnergyINT7_MC				= (TH1D*) ESDContainer_INT7_MC->FindObject("ClusGamma_E");
	histoClusterEnergyEG1_MC				= (TH1D*) ESDContainer_EG1_MC->FindObject("ClusGamma_E");
	histoClusterEnergyEG2_MC				= (TH1D*) ESDContainer_EG2_MC->FindObject("ClusGamma_E");
	histoNEventsINT7_MC				= (TH1D*) ESDContainer_INT7_MC->FindObject("NEvents");
	histoNEventsEG1_MC				= (TH1D*) ESDContainer_EG1_MC->FindObject("NEvents");
	histoNEventsEG2_MC				= (TH1D*) ESDContainer_EG2_MC->FindObject("NEvents");
	
	Int_t NEvtDataINT7 = histoNEventsINT7_Data->GetBinContent(1);	
	Int_t NEvtDataEG1 = histoNEventsEG1_Data->GetBinContent(1);	
	Int_t NEvtDataEG2 = histoNEventsEG2_Data->GetBinContent(1);	
	TH1D* histoClusterEnergyINT7perEventData; 
	TH1D* histoClusterEnergyEG1perEventData;
	TH1D* histoClusterEnergyEG2perEventData;
	histoClusterEnergyINT7perEventData = (TH1D*) Pi0_UncorrectedData->Get("ClusterEPerEvent");;
	histoClusterEnergyEG1perEventData = (TH1D*) Pi0_UncorrectedData_EG1->Get("ClusterEPerEvent");;
	histoClusterEnergyEG2perEventData = (TH1D*) Pi0_UncorrectedData_EG2->Get("ClusterEPerEvent");;

	TH1D* histoClusterEnergyINT7perEventDataInclusive;
	TH1D* histoClusterEnergyEG1perEventDataInclusive;
	TH1D* histoClusterEnergyEG2perEventDataInclusive;
	histoClusterEnergyINT7perEventDataInclusive = (TH1D*) Pi0_UncorrectedData_Inclusive->Get("ClusterEPerEvent");
	histoClusterEnergyEG1perEventDataInclusive = (TH1D*) Pi0_UncorrectedData_EG1_Inclusive->Get("ClusterEPerEvent");
	histoClusterEnergyEG2perEventDataInclusive = (TH1D*) Pi0_UncorrectedData_EG2_Inclusive->Get("ClusterEPerEvent");

	TH1D* histoClusterEnergyINT7perEventMC; 
	TH1D* histoClusterEnergyEG1perEventMC;
	TH1D* histoClusterEnergyEG2perEventMC;
	histoClusterEnergyINT7perEventMC = (TH1D*) Pi0_UncorrectedMC->Get("ClusterEPerEvent");;
	histoClusterEnergyEG1perEventMC = (TH1D*) Pi0_UncorrectedMC_EG1->Get("ClusterEPerEvent");;
	histoClusterEnergyEG2perEventMC = (TH1D*) Pi0_UncorrectedMC_EG2->Get("ClusterEPerEvent");;
	histoNEventsINT7_MC				= (TH1D*) Pi0_UncorrectedMC->Get("NEvents");
	histoNEventsEG1_MC				= (TH1D*) Pi0_UncorrectedMC_EG1->Get("NEvents");
	histoNEventsEG2_MC				= (TH1D*) Pi0_UncorrectedMC_EG2->Get("NEvents");

	Int_t NEvtMCINT7 = histoNEventsINT7_MC->GetBinContent(1);	
	Int_t NEvtMCEG1 = histoNEventsEG1_MC->GetBinContent(1);	
	Int_t NEvtMCEG2 = histoNEventsEG2_MC->GetBinContent(1);	

				///////// read prelininary ratio Eta / Pi0 /////////////////////
	
	TGraph* graphratioPi0Eta_preliminary;
	TGraph* graphratioPi0Eta_Inclusive_preliminary;
	
	graphratioPi0Eta_preliminary	 = (TGraph*) PreliminaryData->Get("Eta5TeVJets/graphRatioEtaToPi0Comb5TeVTotErr"); 
	graphratioPi0Eta_Inclusive_preliminary	 = (TGraph*) PreliminaryData->Get("Eta5TeVMB/graphRatioEtaToPi0Comb5TeVTotErr"); 

	///////// ratio Eta/Pi0 Calculation ///////////////////	

	TFile *MyAnalysis = new TFile("MyAnalysis.root","recreate");

	////////////ratio Eta/Pi0 /////////////////
	canvas_ratio_Pi0Eta->cd();
	TH1D* historatioPi0Eta_Data;	
	TH1D* historatioPi0Eta_MC;
	TH1D* historatioPi0Eta_Data_EG1;
	TH1D* historatioPi0Eta_Data_EG2;	

	historatioPi0Eta_Data = (TH1D*) histoCorrectedYieldNormalEff_Eta_INT7_Data->Clone();
	historatioPi0Eta_Data->Divide(historatioPi0Eta_Data ,histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_Data,1.,1.,"");
	historatioPi0Eta_MC = (TH1D*) histoCorrectedYieldNormalEff_Eta_INT7_MC->Clone();
	historatioPi0Eta_MC->Divide(historatioPi0Eta_MC ,histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_MC,1.,1.,"");
	
	historatioPi0Eta_Data_EG1 = (TH1D*) histoCorrectedYieldNormalEff_Eta_EG1_Data->Clone();
	historatioPi0Eta_Data_EG1->Divide(historatioPi0Eta_Data_EG1 ,histoCorrectedYieldNormalEff_Pi0EtaBinning_EG1_Data,1.,1.,"");
	historatioPi0Eta_Data_EG2 = (TH1D*) histoCorrectedYieldNormalEff_Eta_EG2_Data->Clone();
	historatioPi0Eta_Data_EG2->Divide(historatioPi0Eta_Data_EG2 ,histoCorrectedYieldNormalEff_Pi0EtaBinning_EG2_Data,1.,1.,"");
	
	Int_t n = 9;
	Double_t x[n];
	Double_t y[n];
	Double_t xerror[n];	
	Double_t yerror[n];
	Double_t x_EG1[n];
	Double_t y_EG1[n];
	Double_t xerror_EG1[n];	
	Double_t yerror_EG1[n];
	Double_t x_EG2[n];
	Double_t y_EG2[n];
	Double_t xerror_EG2[n];	
	Double_t yerror_EG2[n];


	for (Int_t i=0; i <n; i ++)
	{
//	x[i] = historatioPi0Eta_Data->GetBinCenter(i);
//	y[i] = historatioPi0Eta_Data->GetBinContent(i);
//	xerror[i] = historatioPi0Eta_Data->GetBinWidth(i)/2;
//	yerror[i] = historatioPi0Eta_Data->GetBinError(i);
	x[i] = historatioPi0Eta_Data->GetBinCenter(i+1);
	y[i] = historatioPi0Eta_Data->GetBinContent(i+1);
	xerror[i] = historatioPi0Eta_Data->GetBinWidth(i+1)/2;
	yerror[i] = historatioPi0Eta_Data->GetBinError(i+1);
	
	x_EG1[i] = historatioPi0Eta_Data_EG1->GetBinCenter(i+1);
	y_EG1[i] = historatioPi0Eta_Data_EG1->GetBinContent(i+1);
	xerror_EG1[i] = historatioPi0Eta_Data_EG1->GetBinWidth(i+1)/2;
	yerror_EG1[i] = historatioPi0Eta_Data_EG1->GetBinError(i+1);
	
	x_EG2[i] = historatioPi0Eta_Data_EG2->GetBinCenter(i+1);
	y_EG2[i] = historatioPi0Eta_Data_EG2->GetBinContent(i+1);
	xerror_EG2[i] = historatioPi0Eta_Data_EG2->GetBinWidth(i+1)/2;
	yerror_EG2[i] = historatioPi0Eta_Data_EG2->GetBinError(i+1);

	cout << x[i] << " " << y[i] << " " << xerror[i] << " " << yerror[i]<< endl;
	}
	TGraphErrors* graphratioPi0Eta = new TGraphErrors(n,x,y,xerror,yerror);
	TGraphErrors* graphratioPi0Eta_EG1 = new TGraphErrors(n,x_EG1,y_EG1,xerror_EG1,yerror_EG1);
	TGraphErrors* graphratioPi0Eta_EG2 = new TGraphErrors(n,x_EG2,y_EG2,xerror_EG2,yerror_EG2);
	
	graphratioPi0Eta_preliminary->SetTitle("#sqrt{s} = 5.02 TeV inside jets preliminary");
	graphratioPi0Eta_preliminary->SetMarkerStyle(24);
	graphratioPi0Eta_preliminary->SetMarkerSize(2);
	graphratioPi0Eta_preliminary->SetMarkerColor(kGray+2);
	graphratioPi0Eta_preliminary->SetDrawOption("AP");
	graphratioPi0Eta_preliminary->SetLineWidth(1);
	graphratioPi0Eta_preliminary->SetLineColor(kGray+1);
	graphratioPi0Eta_preliminary->SetFillStyle(0);


	graphratioPi0Eta_Inclusive_preliminary->SetTitle("#sqrt{s} = 5.02 TeV inclussive preliminary");
	graphratioPi0Eta_Inclusive_preliminary->SetMarkerStyle(25);
	graphratioPi0Eta_Inclusive_preliminary->SetMarkerSize(2);
	graphratioPi0Eta_Inclusive_preliminary->SetMarkerColor(kGray+3);
	graphratioPi0Eta_Inclusive_preliminary->SetDrawOption("AP");
	graphratioPi0Eta_Inclusive_preliminary->SetLineWidth(1);
	graphratioPi0Eta_Inclusive_preliminary->SetLineColor(kGray+3);
	graphratioPi0Eta_Inclusive_preliminary->SetFillStyle(0);

	graphratioPi0Eta->SetTitle("#sqrt{s} = 13 TeV inside jets(My Analysis)");
	graphratioPi0Eta->SetMarkerStyle(23);
	graphratioPi0Eta->SetMarkerSize(2);
	graphratioPi0Eta->SetMarkerColor(4);
	graphratioPi0Eta->SetDrawOption("A");
	graphratioPi0Eta->SetLineWidth(1);
	graphratioPi0Eta->SetLineColor(4);
	graphratioPi0Eta->SetFillStyle(0);

	graphratioPi0Eta_EG1->SetTitle("#sqrt{s} = 13 TeV inside jets(EG1)");
	graphratioPi0Eta_EG1->SetMarkerStyle(22);
	graphratioPi0Eta_EG1->SetMarkerSize(2);
	graphratioPi0Eta_EG1->SetMarkerColor(2);
	graphratioPi0Eta_EG1->SetDrawOption("A");
	graphratioPi0Eta_EG1->SetLineWidth(1);
	graphratioPi0Eta_EG1->SetLineColor(2);
	graphratioPi0Eta_EG1->SetFillStyle(0);

	graphratioPi0Eta_EG2->SetTitle("#sqrt{s} = 13 TeV inside jets(EG2)");
	graphratioPi0Eta_EG2->SetMarkerStyle(21);
	graphratioPi0Eta_EG2->SetMarkerSize(2);
	graphratioPi0Eta_EG2->SetMarkerColor(1);
	graphratioPi0Eta_EG2->SetDrawOption("A");
	graphratioPi0Eta_EG2->SetLineWidth(1);
	graphratioPi0Eta_EG2->SetLineColor(1);
	graphratioPi0Eta_EG2->SetFillStyle(0);
	
	graphratioPi0Eta_preliminary->Write();
	graphratioPi0Eta_Inclusive_preliminary->Write();
	graphratioPi0Eta->Write();
	graphratioPi0Eta_EG1->Write();
	graphratioPi0Eta_EG2->Write();
	
	TMultiGraph *mgratio = new TMultiGraph();
	mgratio->SetTitle(";#it{p}_{T} (GeV/#it{c});#eta/#pi^{0}");
	mgratio->Add(graphratioPi0Eta);
	mgratio->Add(graphratioPi0Eta_preliminary);
	mgratio->Add(graphratioPi0Eta_Inclusive_preliminary);
	mgratio->Add(graphratioPi0Eta_EG1);
	mgratio->Add(graphratioPi0Eta_EG2);
	mgratio->GetHistogram()->GetXaxis()->SetRangeUser(0.,20.);
	mgratio->GetHistogram()->GetYaxis()->SetRangeUser(0.,1.);
	TLegend* legendratioPi0Eta   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
 	legendratioPi0Eta->AddEntry(graphratioPi0Eta,"#sqrt{s} = 13 TeV inside jets (My Analysis, INT7)","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta->AddEntry(graphratioPi0Eta_EG1,"#sqrt{s} = 13 TeV inside jets (My Analysis, EG1)","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta->AddEntry(graphratioPi0Eta_EG2,"#sqrt{s} = 13 TeV inside jets (My Analysis, EG2)","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta->AddEntry(graphratioPi0Eta_preliminary,"#sqrt{s} = 5.02 TeV inside jets preliminary","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta->AddEntry(graphratioPi0Eta_Inclusive_preliminary,"#sqrt{s} = 5.02 TeV inclussive preliminary","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	mgratio->Write();
	TLine *l0 = new TLine(0,0.4,20,0.4);	
	l0->SetLineStyle(2);
	TLine *l1 = new TLine(0,0.6,20,0.6);	
	l1->SetLineStyle(2);
	mgratio->Draw("ap");	
	l0->Draw();
	l1->Draw();
	legendratioPi0Eta->Draw("same");
	canvas_ratio_Pi0Eta->Write();	

	TMultiGraph *mgratio2 = new TMultiGraph();
	mgratio2->SetTitle(";#it{p}_{T} (GeV/#it{c});#eta/#pi^{0}");
	mgratio2->Add(graphratioPi0Eta);
	mgratio2->Add(graphratioPi0Eta_EG1);
	mgratio2->Add(graphratioPi0Eta_EG2);
	mgratio2->GetHistogram()->GetXaxis()->SetRangeUser(0.,20.);
	mgratio2->GetHistogram()->GetYaxis()->SetRangeUser(0.,1.);
	TLegend* legendratioPi0Eta_trigger   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
 	legendratioPi0Eta_trigger->AddEntry(graphratioPi0Eta,"#sqrt{s} = 13 TeV inside jets (INT7)","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta_trigger->AddEntry(graphratioPi0Eta_EG1,"EG1","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioPi0Eta_trigger->AddEntry(graphratioPi0Eta_EG2,"EG2","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	mgratio2->Write();
	canvas_ratio_Pi0Eta_trigger->cd();
	mgratio2->Draw("ap");
	legendratioPi0Eta_trigger->Draw("same");
	canvas_ratio_Pi0Eta_trigger->Write();
 
	///////// To compare each trigger (INT7, EG1, EG2) ////

	// Data
	canvas_triggerClusterData->cd();	
	canvas_triggerClusterData->SetLogy();
	TH1D* histoCompareEG1INT7;
	TH1D* histoCompareEG2INT7;
	TH1D* histoCompareEG1EG2;
	TH2D *href7 = new TH2D("href7", "", 100, 0, 100, 100, 1e-1, 1e+6);
	TF1 *fEnergyEG1INT7 = new TF1 ("fEnergyEG1INT7", "[0]", 0, 100); 
	TF1 *fEnergyEG2INT7 = new TF1 ("fEnergyEG2INT7", "[0]", 0, 100); 
	TF1 *fEnergyEG1EG2  = new TF1 ("fEnergyEG1EG2", "[0]", 0, 100); 

	histoCompareEG1INT7 = (TH1D*) histoClusterEnergyEG1perEventData->Clone();
	histoCompareEG2INT7 = (TH1D*) histoClusterEnergyEG2perEventData->Clone();
	histoCompareEG1EG2 = (TH1D*) histoClusterEnergyEG1perEventData->Clone();
	histoCompareEG1INT7->Divide(histoCompareEG1INT7,histoClusterEnergyINT7perEventData, 1, 1, "");
	histoCompareEG2INT7->Divide(histoCompareEG2INT7,histoClusterEnergyINT7perEventData, 1, 1, "");
	histoCompareEG1EG2->Divide(histoCompareEG1EG2,histoClusterEnergyEG2perEventData, 1, 1, "");
	cout<< "EG1/INT7" <<endl;
	histoCompareEG1INT7->Fit(fEnergyEG1INT7, "", "", 11.9, 40);
	cout<< "parameter: " << fEnergyEG1INT7->GetParameter(0) << endl;
	//Double_t RFEG1_Data = fEnergyEG1INT7->GetParameter(0)*NEvtDataINT7/NEvtDataEG1;
	Double_t RFEG1_Data = 439.23*11.820;
	cout<< "EG1 RF: " << RFEG1_Data << endl;
	cout<< "EG2/INT7" <<endl;
	histoCompareEG2INT7->Fit(fEnergyEG2INT7, "", "", 8.1, 16);
	cout<< "EG2/INT7" <<endl;
	//Double_t RFEG2_Data = fEnergyEG2INT7->GetParameter(0)*NEvtDataINT7/NEvtDataEG2;
	Double_t RFEG2_Data = 439.23;
	cout<< "EG2 RF: " << RFEG2_Data << endl;
	cout<< "parameter: " << fEnergyEG2INT7->GetParameter(0) << endl;
	cout<< "EG1/EG2" <<endl;
	histoCompareEG1EG2->Fit(fEnergyEG1EG2, "", "", 11.9, 40);
	cout<< "parameter: " << fEnergyEG1EG2->GetParameter(0) << endl;
	cout << "Nevent of INT7: " << NEvtDataINT7 <<endl; 
	cout << "Nevent of EG1 : " << NEvtDataEG1 <<endl; 
	cout << "Nevent of EG2 : " << NEvtDataEG2 <<endl; 

	DrawGammaSetMarker(histoCompareEG1INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2, 26, 1., 1, 1);
	href7->SetXTitle("#it{E} (GeV)"); 
	href7->SetYTitle("Ratio of cluster per event for data"); 
	//histoCompareEG1INT7->SetAxisRange(0.,10000,"Y");	
	TLegend* legendComparetrigger   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger->AddEntry(histoCompareEG1INT7, "#frac{EG1}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger->AddEntry(histoCompareEG2INT7, "#frac{EG2}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger->AddEntry(histoCompareEG1EG2, "#frac{EG1}{EG2}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	TLatex *latexComparetrigger = new TLatex();
	latexComparetrigger->SetTextAlign(12);
	latexComparetrigger->SetTextSize(0.05);

	href7->Draw();
	histoCompareEG1INT7->Draw("same");
	histoCompareEG2INT7->Draw("same");
	histoCompareEG1EG2->Draw("same");
	legendComparetrigger->Draw("same");
	latexComparetrigger->DrawLatex(60, 0.5, Form("EG1/INT7: %.3e", fEnergyEG1INT7->GetParameter(0) ));
	latexComparetrigger->DrawLatex(60, 1.2, Form("EG2/INT7: %.3e", fEnergyEG2INT7->GetParameter(0) ));
	histoCompareEG1INT7->Write();	
	canvas_triggerClusterData->Write();

//// Data and EG1 trigger
	
	canvas_triggerClusterData_Inclusive->cd();	
	canvas_triggerClusterData_Inclusive->SetLogy();
	TH1D* histoCompareEG1INT7_Inclusive;
	TH1D* histoCompareEG2INT7_Inclusive;
	TH1D* histoCompareEG1EG2_Inclusive;
	TH2D *hhref7 = new TH2D("href7", "", 100, 0, 100, 100, 1e-1, 1e+6);
	TF1 *fEnergyEG1INT7_Inclusive = new TF1 ("fEnergyEG1INT7_Inclusive", "[0]", 0, 100); 
	TF1 *fEnergyEG2INT7_Inclusive = new TF1 ("fEnergyEG2INT7_Inclusive", "[0]", 0, 100); 
	TF1 *fEnergyEG1EG2_Inclusive  = new TF1 ("fEnergyEG1EG2_Inclusive", "[0]", 0, 100); 

	histoCompareEG1INT7_Inclusive = (TH1D*) histoClusterEnergyEG1perEventDataInclusive->Clone();
	histoCompareEG2INT7_Inclusive = (TH1D*) histoClusterEnergyEG2perEventDataInclusive->Clone();
	histoCompareEG1EG2_Inclusive = (TH1D*) histoClusterEnergyEG1perEventDataInclusive->Clone();
	histoCompareEG1INT7_Inclusive->Divide(histoCompareEG1INT7_Inclusive,histoClusterEnergyINT7perEventDataInclusive, 1, 1, "");
	histoCompareEG2INT7_Inclusive->Divide(histoCompareEG2INT7_Inclusive,histoClusterEnergyINT7perEventDataInclusive, 1, 1, "");
	histoCompareEG1EG2_Inclusive->Divide(histoCompareEG1EG2_Inclusive,histoClusterEnergyEG2perEventDataInclusive, 1, 1, "");
	cout<< "EG1/INT7" <<endl;
	histoCompareEG1INT7_Inclusive->Fit(fEnergyEG1INT7_Inclusive, "", "", 11.9, 40);
	cout<< "parameter: " << fEnergyEG1INT7_Inclusive->GetParameter(0) << endl;
	//Double_t RFEG1_Data = fEnergyEG1INT7->GetParameter(0)*NEvtDataINT7/NEvtDataEG1;
	//Double_t RFEG1_Data = 439.23*11.820;
	Double_t RFEG1_Data_Inclusive = fEnergyEG1INT7_Inclusive->GetParameter(0);
	Double_t RFEG1_DataError_Inclusive = fEnergyEG1INT7_Inclusive->GetParError(0);
	cout<< "EG1 RF: " << RFEG1_Data << endl;
	cout<< "EG2/INT7" <<endl;
	histoCompareEG2INT7_Inclusive->Fit(fEnergyEG2INT7_Inclusive, "", "", 8.1, 16);
	cout<< "EG2/INT7" <<endl;
	//Double_t RFEG2_Data = fEnergyEG2INT7->GetParameter(0)*NEvtDataINT7/NEvtDataEG2;
	//Double_t RFEG2_Data = 439.23;
	Double_t RFEG2_Data_Inclusive = fEnergyEG2INT7_Inclusive->GetParameter(0);
	Double_t RFEG2_DataError_Inclusive = fEnergyEG2INT7_Inclusive->GetParError(0);
	cout<< "EG2 RF: " << RFEG2_Data << endl;
	cout<< "parameter: " << fEnergyEG2INT7_Inclusive->GetParameter(0) << endl;
	cout<< "EG1/EG2" <<endl;
	histoCompareEG1EG2_Inclusive->Fit(fEnergyEG1EG2_Inclusive, "", "", 11.9, 40);
	cout<< "parameter: " << fEnergyEG1EG2->GetParameter(0) << endl;

	DrawGammaSetMarker(histoCompareEG1INT7_Inclusive, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7_Inclusive, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2_Inclusive, 26, 1., 1, 1);
	hhref7->SetXTitle("#it{E} (GeV)"); 
	hhref7->SetYTitle("Ratio of cluster per event for data about EG1 trigger which is inclusive measurement"); 
	//histoCompareEG1INT7->SetAxisRange(0.,10000,"Y");	
	TLegend* legendComparetrigger_Inclusive   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_Inclusive->AddEntry(histoCompareEG1INT7_Inclusive, "#frac{EG1}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Inclusive->AddEntry(histoCompareEG2INT7_Inclusive, "#frac{EG2}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Inclusive->AddEntry(histoCompareEG1EG2_Inclusive, "#frac{EG1}{EG2}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
//	TLatex *latexComparetrigger = new TLatex();
//	latexComparetrigger->SetTextAlign(12);
//	latexComparetrigger->SetTextSize(0.05);

	hhref7->Draw();
	histoCompareEG1INT7_Inclusive->Draw("same");
	histoCompareEG2INT7_Inclusive->Draw("same");
	histoCompareEG1EG2_Inclusive->Draw("same");
	legendComparetrigger_Inclusive->Draw("same");
	latexComparetrigger->DrawLatex(60, 0.5, Form("EG1/INT7: %.3e", fEnergyEG1INT7_Inclusive->GetParameter(0) ));
	latexComparetrigger->DrawLatex(60, 1.2, Form("EG2/INT7: %.3e", fEnergyEG2INT7_Inclusive->GetParameter(0) ));
	histoCompareEG1INT7_Inclusive->Write();	
	canvas_triggerClusterData_Inclusive->Write();



////MC

	canvas_triggerClusterMC->cd();	
	canvas_triggerClusterMC->SetLogy();
	TH1D* histoCompareEG1INT7_MC;
	TH1D* histoCompareEG2INT7_MC;
	TH1D* histoCompareEG1EG2_MC;
	TH2D *href8 = new TH2D("href8", "", 100, 0, 100, 100, 1e-1, 1e+6);
	TF1 *fEnergyEG1INT7_MC = new TF1 ("fEnergyEG1INT7_MC", "[0]", 0, 100); 
	TF1 *fEnergyEG2INT7_MC = new TF1 ("fEnergyEG2INT7_MC", "[0]", 0, 100); 
	TF1 *fEnergyEG1EG2_MC  = new TF1 ("fEnergyEG1EG2_MC", "[0]", 0, 100); 

	histoCompareEG1INT7_MC = (TH1D*) histoClusterEnergyEG1perEventMC->Clone();
	histoCompareEG2INT7_MC = (TH1D*) histoClusterEnergyEG2perEventMC->Clone();
	histoCompareEG1EG2_MC = (TH1D*) histoClusterEnergyEG1perEventMC->Clone();
	histoCompareEG1INT7_MC->Divide(histoCompareEG1INT7_MC,histoClusterEnergyINT7perEventMC, 1, 1, "");
	histoCompareEG2INT7_MC->Divide(histoCompareEG2INT7_MC,histoClusterEnergyINT7perEventMC, 1, 1, "");
	histoCompareEG1EG2_MC->Divide(histoCompareEG1EG2_MC,histoClusterEnergyEG2perEventMC, 1, 1, "");
	cout<< "MC: EG1/INT7" <<endl;
	histoCompareEG1INT7_MC->Fit(fEnergyEG1INT7_MC, "", "", 10, 100);
	Double_t RFEG1Pi0ClusterMC = fEnergyEG1INT7_MC->GetParameter(0);
	cout<< "parameter: " << fEnergyEG1INT7_MC->GetParameter(0) << endl;
	cout<< "MC: EG2/INT7" <<endl;
	histoCompareEG2INT7_MC->Fit(fEnergyEG2INT7_MC, "", "", 5, 100);
	Double_t RFEG2Pi0ClusterMC = fEnergyEG2INT7_MC->GetParameter(0);
	cout<< "parameter: " << fEnergyEG2INT7_MC->GetParameter(0) << endl;
	cout<< "MC: EG1/EG2" <<endl;
	histoCompareEG1EG2_MC->Fit(fEnergyEG1EG2_MC, "", "", 10, 100);
	cout<< "parameter: " << fEnergyEG1EG2_MC->GetParameter(0) << endl;

	DrawGammaSetMarker(histoCompareEG1INT7_MC, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7_MC, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2_MC, 26, 1., 1, 1);
	href8->SetXTitle("#it{E} (GeV)"); 
	href8->SetYTitle("Ratio of Cluster per event for MC"); 
	TLegend* legendComparetrigger_MC   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_MC->AddEntry(histoCompareEG1INT7_MC, "#frac{EG1}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_MC->AddEntry(histoCompareEG2INT7_MC, "#frac{EG2}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_MC->AddEntry(histoCompareEG1EG2_MC, "#frac{EG1}{EG2}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	TLatex *latexComparetrigger_MC = new TLatex();
	latexComparetrigger_MC->SetTextAlign(12);
	latexComparetrigger_MC->SetTextSize(0.05);

	href8->Draw();
	histoCompareEG1INT7_MC->Draw("same");
	histoCompareEG2INT7_MC->Draw("same");
	histoCompareEG1EG2_MC->Draw("same");
	legendComparetrigger_MC->Draw("same");
	latexComparetrigger_MC->DrawLatex(60, 0.5, Form("EG1/INT7: %.3e", fEnergyEG1INT7_MC->GetParameter(0) ));
	latexComparetrigger_MC->DrawLatex(60, 1.2, Form("EG2/INT7: %.3e", fEnergyEG2INT7_MC->GetParameter(0) ));
	histoCompareEG1INT7_MC->Write();	
	canvas_triggerClusterMC->Write();

	
	canvas_triggerPi0YieldDataRatio->cd();	
	canvas_triggerPi0YieldDataRatio->SetLogy();
	TH1D* histoCompareEG1INT7Pi0YieldData;
	TH1D* histoCompareEG2INT7Pi0YieldData;
	TH1D* histoCompareEG1EG2Pi0YieldData;
	TH2D *href5 = new TH2D("href5", "", 100, 0, 20, 100, 0, 1e+2);	
	href5->SetXTitle("#it{p}_{T} (GeV/#it{c})"); 
	href5->SetYTitle("Ratio of raw yield per event"); 
	histoCompareEG1INT7Pi0YieldData = (TH1D*) histoYieldMeson_Pi0_Data_EG1->Clone();
	histoCompareEG2INT7Pi0YieldData = (TH1D*) histoYieldMeson_Pi0_Data_EG2->Clone();
	histoCompareEG1EG2Pi0YieldData = (TH1D*) histoYieldMeson_Pi0_Data_EG1->Clone();
 	histoCompareEG1INT7Pi0YieldData->Divide(histoCompareEG1INT7Pi0YieldData,histoYieldMeson_Pi0_Data_INT7, 1, 1, "");
	histoCompareEG2INT7Pi0YieldData->Divide(histoCompareEG2INT7Pi0YieldData,histoYieldMeson_Pi0_Data_INT7, 1, 1, "");
	histoCompareEG1EG2Pi0YieldData->Divide(histoCompareEG1EG2Pi0YieldData,histoYieldMeson_Pi0_Data_EG2, 1, 1, "");
	TF1 *fYieldEG1INT7Pi0Data = new TF1("fYieldEG1INT7","[0]", 0, 100);
	TF1 *fYieldEG2INT7Pi0Data = new TF1("fYieldEG2INT7","[0]", 0, 100);
	TF1 *fYieldEG1EG2Pi0Data  = new TF1("fYieldEG1EG2","[0]", 0, 100);
	histoCompareEG1INT7Pi0YieldData->Fit("fYieldEG1INT7","","",10,20);	
	histoCompareEG2INT7Pi0YieldData->Fit("fYieldEG2INT7","","",8,20);	
	histoCompareEG1EG2Pi0YieldData ->Fit("fYieldEG1EG2","","",10,20);	
	Double_t RFEG1Pi0YieldData = fYieldEG1INT7Pi0Data->GetParameter(0);
	Double_t RFEG2Pi0YieldData = fYieldEG2INT7Pi0Data->GetParameter(0);
	Double_t RFEGX2Pi0YieldData = fYieldEG1EG2Pi0Data->GetParameter(0);
	Double_t RFEG1Pi0YieldDataError = fYieldEG1INT7Pi0Data->GetParError(0);
	Double_t RFEG2Pi0YieldDataError = fYieldEG2INT7Pi0Data->GetParError(0);
	fYieldEG1INT7Pi0Data->Write("EG1INT7RejectionFactorOfYield");		
	fYieldEG2INT7Pi0Data->Write("EG2INT7RejectionFactorOfYield");		
	fYieldEG1EG2Pi0Data->Write("EG1EG2RejectionFactorOfYield");		

	DrawGammaSetMarker(histoCompareEG1INT7Pi0YieldData, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7Pi0YieldData, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2Pi0YieldData, 26, 1., 1, 1);
	TLegend* legendComparetriggerPi0YieldData   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetriggerPi0YieldData->AddEntry(histoCompareEG1INT7Pi0YieldData, "#frac{EG1}{INT7}", "p"); 
 	legendComparetriggerPi0YieldData->AddEntry(histoCompareEG2INT7Pi0YieldData, "#frac{EG2}{INT7}", "p"); 
 	legendComparetriggerPi0YieldData->AddEntry(histoCompareEG1EG2Pi0YieldData, "#frac{EG1}{EG2}", "p"); 

	TLatex *latexComparetriggerPi0YieldData = new TLatex();
	latexComparetriggerPi0YieldData->SetTextAlign(12);
	latexComparetriggerPi0YieldData->SetTextSize(0.05);

	TLatex *latexComparetriggerEtaYieldData = new TLatex();
	latexComparetriggerEtaYieldData->SetTextAlign(12);
	latexComparetriggerEtaYieldData->SetTextSize(0.05);

	href5->Draw();
	histoCompareEG1INT7Pi0YieldData->Draw("samep");
	histoCompareEG2INT7Pi0YieldData->Draw("samep");
	histoCompareEG1EG2Pi0YieldData->Draw("samep");
	legendComparetriggerPi0YieldData->Draw("same");
	latexComparetriggerPi0YieldData->DrawLatex(10, 1.5, Form("EG1/INT7: %.3e", RFEG1Pi0YieldData ));
	latexComparetriggerPi0YieldData->DrawLatex(10, 2.0, Form("EG2/INT7: %.3e", RFEG2Pi0YieldData ));
	histoCompareEG1INT7Pi0YieldData->Write();	
	canvas_triggerPi0YieldDataRatio->Write();
	
	canvas_triggerPi0YieldDataRatio_Inclusive->cd();	
	canvas_triggerPi0YieldDataRatio_Inclusive->SetLogy();
	TH1D* histoCompareEG1INT7Pi0YieldData_Inclusive;
	TH1D* histoCompareEG2INT7Pi0YieldData_Inclusive;
	TH1D* histoCompareEG1EG2Pi0YieldData_Inclusive;
	TH2D *href5_Inclusive = new TH2D("href5_Inclusive", "", 100, 0, 20, 100, 7e-2, 1e+4);	
	href5_Inclusive->SetXTitle("#it{p}_{T} (GeV/#it{c})"); 
	href5_Inclusive->SetYTitle("Ratio of raw yield per event in inclusive measurement"); 
	histoCompareEG1INT7Pi0YieldData_Inclusive = (TH1D*) histoYieldMeson_Pi0_Data_EG1_Inclusive->Clone();
	histoCompareEG2INT7Pi0YieldData_Inclusive = (TH1D*) histoYieldMeson_Pi0_Data_EG2_Inclusive->Clone();
	histoCompareEG1EG2Pi0YieldData_Inclusive = (TH1D*) histoYieldMeson_Pi0_Data_EG1_Inclusive->Clone();
 	histoCompareEG1INT7Pi0YieldData_Inclusive->Divide(histoCompareEG1INT7Pi0YieldData_Inclusive,histoYieldMeson_Pi0_Data_INT7_Inclusive, 1, 1, "");
	histoCompareEG2INT7Pi0YieldData_Inclusive->Divide(histoCompareEG2INT7Pi0YieldData_Inclusive,histoYieldMeson_Pi0_Data_INT7_Inclusive, 1, 1, "");
	histoCompareEG1EG2Pi0YieldData_Inclusive->Divide(histoCompareEG1EG2Pi0YieldData_Inclusive,histoYieldMeson_Pi0_Data_EG2_Inclusive, 1, 1, "");
	TF1 *fYieldEG1INT7Pi0Data_Inclusive = new TF1("fYieldEG1INT7_Inclusive","[0]", 0, 100);
	TF1 *fYieldEG2INT7Pi0Data_Inclusive = new TF1("fYieldEG2INT7_Inclusive","[0]", 0, 100);
	TF1 *fYieldEG1EG2Pi0Data_Inclusive  = new TF1("fYieldEG1EG2_Inclusive","[0]", 0, 100);
	histoCompareEG1INT7Pi0YieldData_Inclusive->Fit("fYieldEG1INT7_Inclusive","","",10,20);	
	histoCompareEG2INT7Pi0YieldData_Inclusive->Fit("fYieldEG2INT7_Inclusive","","",5,20);	
	histoCompareEG1EG2Pi0YieldData_Inclusive ->Fit("fYieldEG1EG2_Inclusive","","",10,20);	
	Double_t RFEG1Pi0YieldData_Inclusive = fYieldEG1INT7Pi0Data_Inclusive->GetParameter(0);
	Double_t RFEG2Pi0YieldData_Inclusive = fYieldEG2INT7Pi0Data_Inclusive->GetParameter(0);

	DrawGammaSetMarker(histoCompareEG1INT7Pi0YieldData_Inclusive, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7Pi0YieldData_Inclusive, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2Pi0YieldData_Inclusive, 26, 1., 1, 1);
	TLegend* legendComparetriggerPi0YieldData_Inclusive   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetriggerPi0YieldData_Inclusive->AddEntry(histoCompareEG1INT7Pi0YieldData_Inclusive, "#frac{EG1}{INT7}", "p"); 
 	legendComparetriggerPi0YieldData_Inclusive->AddEntry(histoCompareEG2INT7Pi0YieldData_Inclusive, "#frac{EG2}{INT7}", "p"); 
 	legendComparetriggerPi0YieldData_Inclusive->AddEntry(histoCompareEG1EG2Pi0YieldData_Inclusive, "#frac{EG1}{EG2}", "p"); 

	TLatex *latexComparetriggerPi0YieldData_Inclusive = new TLatex();
	latexComparetriggerPi0YieldData_Inclusive->SetTextAlign(12);
	latexComparetriggerPi0YieldData_Inclusive->SetTextSize(0.05);

	href5_Inclusive->Draw();
	histoCompareEG1INT7Pi0YieldData_Inclusive->Draw("samep");
	histoCompareEG2INT7Pi0YieldData_Inclusive->Draw("samep");
	histoCompareEG1EG2Pi0YieldData_Inclusive->Draw("samep");
	legendComparetriggerPi0YieldData_Inclusive->Draw("same");
	latexComparetriggerPi0YieldData_Inclusive->DrawLatex(10, 1.5, Form("EG1/INT7: %.3e", RFEG1Pi0YieldData_Inclusive ));
	latexComparetriggerPi0YieldData_Inclusive->DrawLatex(10, 2.0, Form("EG2/INT7: %.3e", RFEG2Pi0YieldData_Inclusive ));
	histoCompareEG1INT7Pi0YieldData_Inclusive->Write();	
	canvas_triggerPi0YieldDataRatio_Inclusive->Write();


	
	canvas_triggerPi0YieldMCRatio->cd();	
	canvas_triggerPi0YieldMCRatio->SetLogy();
	TH1D* histoCompareEG1INT7Pi0YieldMC;
	TH1D* histoCompareEG2INT7Pi0YieldMC;
	TH1D* histoCompareEG1EG2Pi0YieldMC;
	TH2D *hhref5 = new TH2D("hhref5", "", 100, 0, 20, 100, 1e-1, 1e+6);	
	hhref5->SetXTitle("#it{p}_{T} (GeV/#it{c})"); 
	hhref5->SetYTitle("Ratio of raw yield per event"); 
	histoCompareEG1INT7Pi0YieldMC = (TH1D*) histoYieldMeson_Pi0_MC_EG1->Clone();
	histoCompareEG2INT7Pi0YieldMC = (TH1D*) histoYieldMeson_Pi0_MC_EG2->Clone();
	histoCompareEG1EG2Pi0YieldMC = (TH1D*) histoYieldMeson_Pi0_MC_EG1->Clone();
 	histoCompareEG1INT7Pi0YieldMC->Divide(histoCompareEG1INT7Pi0YieldMC,histoYieldMeson_Pi0_MC_INT7, 1, 1, "");
	histoCompareEG2INT7Pi0YieldMC->Divide(histoCompareEG2INT7Pi0YieldMC,histoYieldMeson_Pi0_MC_INT7, 1, 1, "");
	histoCompareEG1EG2Pi0YieldMC->Divide(histoCompareEG1EG2Pi0YieldMC,histoYieldMeson_Pi0_MC_EG2, 1, 1, "");
	TF1 *fYieldEG1INT7Pi0MC = new TF1("fYieldEG1INT7","[0]", 0, 100);
	TF1 *fYieldEG2INT7Pi0MC = new TF1("fYieldEG2INT7","[0]", 0, 100);
	TF1 *fYieldEG1EG2Pi0MC  = new TF1("fYieldEG1EG2","[0]", 0, 100);
	histoCompareEG1INT7Pi0YieldMC->Fit("fYieldEG1INT7","","",12,20);	
	histoCompareEG2INT7Pi0YieldMC->Fit("fYieldEG2INT7","","",6,20);	
	histoCompareEG1EG2Pi0YieldMC ->Fit("fYieldEG1EG2","","",12,20);	

	Double_t RFEG1Pi0YieldMC = fYieldEG1INT7Pi0MC->GetParameter(0);
	Double_t RFEG2Pi0YieldMC = fYieldEG2INT7Pi0MC->GetParameter(0);
	DrawGammaSetMarker(histoCompareEG1INT7Pi0YieldMC, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7Pi0YieldMC, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2Pi0YieldMC, 26, 1., 1, 1);
	TLegend* legendComparetriggerPi0YieldMC   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetriggerPi0YieldMC->AddEntry(histoCompareEG1INT7Pi0YieldMC, "#frac{EG1}{INT7}", "p"); 
 	legendComparetriggerPi0YieldMC->AddEntry(histoCompareEG2INT7Pi0YieldMC, "#frac{EG2}{INT7}", "p"); 
 	legendComparetriggerPi0YieldMC->AddEntry(histoCompareEG1EG2Pi0YieldMC, "#frac{EG1}{EG2}", "p"); 

	TLatex *latexComparetriggerPi0YieldMC = new TLatex();
	latexComparetriggerPi0YieldMC->SetTextAlign(12);
	latexComparetriggerPi0YieldMC->SetTextSize(0.05);

	hhref5->Draw();
	histoCompareEG1INT7Pi0YieldMC->Draw("samep");
	histoCompareEG2INT7Pi0YieldMC->Draw("samep");
	histoCompareEG1EG2Pi0YieldMC->Draw("samep");
	legendComparetriggerPi0YieldMC->Draw("same");
	latexComparetriggerPi0YieldMC->DrawLatex(10, 0.5, Form("EG1/INT7: %.3e", RFEG1Pi0YieldMC ));
	latexComparetriggerPi0YieldMC->DrawLatex(10, 1, Form("EG2/INT7: %.3e", RFEG2Pi0YieldMC ));
	histoCompareEG1INT7Pi0YieldMC->Write();	
	canvas_triggerPi0YieldMCRatio->Write();


	// Pi0 Correctedyield with each trigger
	canvas_CorrectedYield_Pi0->cd();
	canvas_CorrectedYield_Pi0->SetLogy();
	TH2D* href = new TH2D("href", "", 100, 0 ,20, 10, 1e-8, 1e-1);
	href->GetYaxis()->SetTitle("Corrected Yield");
	href->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	
	TH1D *histoCorrectedYieldNormalEff_Pi0_INT7_Data_Clone;
	TH1D *histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone;
	TH1D *histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone;
	
	histoCorrectedYieldNormalEff_Pi0_INT7_Data_Clone = (TH1D*) histoCorrectedYieldNormalEff_Pi0_INT7_Data->Clone();
	histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone = (TH1D*) histoCorrectedYieldNormalEff_Pi0_EG1_Data->Clone();
	histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone = (TH1D*) histoCorrectedYieldNormalEff_Pi0_EG2_Data->Clone();

	DrawGammaSetMarker(histoCorrectedYieldNormalEff_Pi0_INT7_Data_Clone, 24, 1., 4,4 );
	DrawGammaSetMarker(histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone, 25, 1.,2 ,2 );
	DrawGammaSetMarker(histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone, 26, 1., 1, 1);

	TLegend* legendCorrectedYieldPi0   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
 	legendCorrectedYieldPi0->AddEntry(histoCorrectedYieldNormalEff_Pi0_INT7_Data_Clone,"INT7","p"); 
 	legendCorrectedYieldPi0->AddEntry(histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone,"EG1","p"); 
 	legendCorrectedYieldPi0->AddEntry(histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone,"EG2","p");
	histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone->Scale(1/RFEG1Pi0YieldData); // N_{Event}/N_{EG1}
	histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone->Scale(1/RFEG2Pi0YieldData); // N_{Event}/N_{EG2}
	//histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone->Divide(fEnergyEG1INT7,1);
//	histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone->Scale(1/RFEG1_Data); // N_{Event}/N_{EG1}
//	histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone->Scale(NEvtDataEG1/NEvtDataINT7); // N_{Event}/N_{EG1}
//	histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone->Scale(1/RFEG2_Data);
//	histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone->Scale(NEvtDataEG2/NEvtDataINT7); // N_{Event}/N_{EG1}
	href->Draw();
	histoCorrectedYieldNormalEff_Pi0_INT7_Data_Clone->Draw("samep");
	histoCorrectedYieldNormalEff_Pi0_EG1_Data_Clone->Draw("samep");
	histoCorrectedYieldNormalEff_Pi0_EG2_Data_Clone->Draw("samep");
	legendCorrectedYieldPi0->Draw("samep");

	// Eta Correctedyield with each trigger
	canvas_CorrectedYield_Eta->cd();
	canvas_CorrectedYield_Eta->SetLogy();
	TH2D* href2 = new TH2D("href2", "", 100, 0 ,20, 10, 1e-9, 1e-3);
	href2->GetYaxis()->SetTitle("Corrected Yield");
	href2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	TLegend* legendCorrectedYieldEta   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
 	legendCorrectedYieldEta->AddEntry(histoCorrectedYieldNormalEff_Eta_INT7_Data,"INT7","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendCorrectedYieldEta->AddEntry(histoCorrectedYieldNormalEff_Eta_EG1_Data,"EG1","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendCorrectedYieldEta->AddEntry(histoCorrectedYieldNormalEff_Eta_EG2_Data,"EG2","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	
	DrawGammaSetMarker(histoCorrectedYieldNormalEff_Eta_INT7_Data, 24, 1., 4,4 );
	DrawGammaSetMarker(histoCorrectedYieldNormalEff_Eta_EG1_Data, 25, 1.,2 ,2 );
	DrawGammaSetMarker(histoCorrectedYieldNormalEff_Eta_EG2_Data, 26, 1., 1, 1);

	href2->Draw();
	histoCorrectedYieldNormalEff_Eta_INT7_Data->Draw("samep");
	histoCorrectedYieldNormalEff_Eta_EG1_Data->Draw("samep");
	histoCorrectedYieldNormalEff_Eta_EG2_Data->Draw("samep");
	legendCorrectedYieldEta->Draw("samep");
	
	////////// Compare narrow/norm + wide/norm corrected yield ////////////////
	TH1D* historatioNarrowDefault_Pi0EtaBinning_Data;
	TH1D* historatioNarrowDefault_Pi0EtaBinning_MC;
	TH1D* historatioNarrowDefault_Eta_Data;
	TH1D* historatioNarrowDefault_Eta_MC;
	TH1D* historatioWideDefault_Pi0EtaBinning_Data;
	TH1D* historatioWideDefault_Pi0EtaBinning_MC;
	TH1D* historatioWideDefault_Eta_Data;
	TH1D* historatioWideDefault_Eta_MC;

	historatioNarrowDefault_Pi0EtaBinning_Data = (TH1D*) histoCorrectedYieldNormalEff_Narrow_Pi0EtaBinning_Data->Clone();
	historatioNarrowDefault_Pi0EtaBinning_Data->Divide(historatioNarrowDefault_Pi0EtaBinning_Data, histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_Data, 1., 1., "");	
	historatioNarrowDefault_Pi0EtaBinning_MC = (TH1D*) histoCorrectedYieldNormalEff_Narrow_Pi0EtaBinning_MC->Clone();
	historatioNarrowDefault_Pi0EtaBinning_MC->Divide(historatioNarrowDefault_Pi0EtaBinning_MC, histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_MC, 1., 1., "");	
	historatioNarrowDefault_Eta_Data = (TH1D*) histoCorrectedYieldNormalEff_Narrow_Eta_Data->Clone();
	historatioNarrowDefault_Eta_Data->Divide(historatioNarrowDefault_Eta_Data, histoCorrectedYieldNormalEff_Eta_INT7_Data, 1., 1., "");	
	historatioNarrowDefault_Eta_MC = (TH1D*) histoCorrectedYieldNormalEff_Narrow_Eta_MC->Clone();
	historatioNarrowDefault_Eta_MC->Divide(historatioNarrowDefault_Eta_MC, histoCorrectedYieldNormalEff_Eta_INT7_MC, 1., 1., "");		
	historatioWideDefault_Pi0EtaBinning_Data = (TH1D*) histoCorrectedYieldNormalEff_Wide_Pi0EtaBinning_Data->Clone();
	historatioWideDefault_Pi0EtaBinning_Data->Divide(historatioWideDefault_Pi0EtaBinning_Data, histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_Data, 1., 1., "");	
	historatioWideDefault_Pi0EtaBinning_MC = (TH1D*) histoCorrectedYieldNormalEff_Wide_Pi0EtaBinning_MC->Clone();
	historatioWideDefault_Pi0EtaBinning_MC->Divide(historatioWideDefault_Pi0EtaBinning_MC, histoCorrectedYieldNormalEff_Pi0EtaBinning_INT7_MC, 1., 1., "");	
	historatioWideDefault_Eta_Data = (TH1D*) histoCorrectedYieldNormalEff_Wide_Eta_Data->Clone();
	historatioWideDefault_Eta_Data->Divide(historatioWideDefault_Eta_Data, histoCorrectedYieldNormalEff_Eta_INT7_Data, 1., 1., "");	
	historatioWideDefault_Eta_MC = (TH1D*) histoCorrectedYieldNormalEff_Wide_Eta_MC->Clone();
	historatioWideDefault_Eta_MC->Divide(historatioWideDefault_Eta_MC, histoCorrectedYieldNormalEff_Eta_INT7_MC, 1., 1., "");
//
	TLegend* leg1   = GetAndSetLegend2(0.12+0.4, 0.14, 0.45+0.5, 0.30+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
//	
	leg1->AddEntry(historatioNarrowDefault_Pi0EtaBinning_Data,"#frac{Narrow}{Standard} Data","lp");
	leg1->AddEntry(historatioNarrowDefault_Pi0EtaBinning_MC,"#frac{Narrow}{Standard} MC","lp");
	leg1->AddEntry(historatioWideDefault_Pi0EtaBinning_Data,"#frac{Wide}{Standard} Data","lp");
	leg1->AddEntry(historatioWideDefault_Pi0EtaBinning_MC,"#frac{Wide}{Standard} MC","lp");
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
//	
	historatioNarrowDefault_Pi0EtaBinning_Data->SetAxisRange(0.9, 1.1, "Y");
	historatioNarrowDefault_Pi0EtaBinning_MC->SetAxisRange(0.9, 1.1, "Y");
	historatioWideDefault_Pi0EtaBinning_Data->SetAxisRange(0.9, 1.1, "Y");
	historatioWideDefault_Pi0EtaBinning_MC->SetAxisRange(0.9, 1.1, "Y");
	historatioNarrowDefault_Pi0EtaBinning_Data->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	historatioNarrowDefault_Pi0EtaBinning_MC->SetXTitle("p_{T}");
	historatioWideDefault_Pi0EtaBinning_Data->SetXTitle("p_{T}");
	historatioWideDefault_Pi0EtaBinning_MC->SetXTitle("p_{T}");
	historatioNarrowDefault_Pi0EtaBinning_Data->SetYTitle("#frac{yield_{diff}}{yield_{standard}}");
	
	historatioNarrowDefault_Pi0EtaBinning_Data->SetMarkerStyle(20);
	historatioNarrowDefault_Pi0EtaBinning_Data->SetMarkerSize(2);
	historatioNarrowDefault_Pi0EtaBinning_Data->SetMarkerColor(kBlack);
	historatioNarrowDefault_Pi0EtaBinning_Data->SetLineColor(kBlack);
	historatioNarrowDefault_Pi0EtaBinning_MC->SetMarkerStyle(25);
	historatioNarrowDefault_Pi0EtaBinning_MC->SetMarkerSize(2);
	historatioNarrowDefault_Pi0EtaBinning_MC->SetMarkerColor(kGreen+2);
	historatioNarrowDefault_Pi0EtaBinning_MC->SetLineColor(kGreen+2);
	historatioWideDefault_Pi0EtaBinning_Data->SetMarkerStyle(20);
	historatioWideDefault_Pi0EtaBinning_Data->SetMarkerSize(2);
	historatioWideDefault_Pi0EtaBinning_Data->SetMarkerColor(kRed);
	historatioWideDefault_Pi0EtaBinning_Data->SetLineColor(kRed);
	historatioWideDefault_Pi0EtaBinning_MC->SetMarkerStyle(25);
	historatioWideDefault_Pi0EtaBinning_MC->SetMarkerSize(2);
	historatioWideDefault_Pi0EtaBinning_MC->SetMarkerColor(kBlue);
	historatioWideDefault_Pi0EtaBinning_MC->SetLineColor(kBlue);


	canvas_ratio_DiffDefault_Pi0->cd();
	
	historatioNarrowDefault_Pi0EtaBinning_Data->Draw();
	historatioNarrowDefault_Pi0EtaBinning_MC->Draw("same");
	historatioWideDefault_Pi0EtaBinning_Data->Draw("same");
	historatioWideDefault_Pi0EtaBinning_MC->Draw("same");
	leg1->Draw("same");
	
	canvas_ratio_DiffDefault_Pi0->Write();
//
	////////// Ratio normal,narrow,wide/ right,left	/////////////////
	
	TH1D* historatioRightNormal_Pi0EtaBinning; 
	TH1D* historatioRightNarrow_Pi0EtaBinning; 
	TH1D* historatioRightWide_Pi0EtaBinning; 
	TH1D* historatioLeftNormal_Pi0EtaBinning; 
	TH1D* historatioLeftNarrow_Pi0EtaBinning; 
	TH1D* historatioLeftWide_Pi0EtaBinning; 

	historatioRightNormal_Pi0EtaBinning = (TH1D*) histoCorrectedYieldNormalEff_Right_Normal_Pi0EtaBinning -> Clone();
	historatioRightNormal_Pi0EtaBinning->Divide(historatioRightNormal_Pi0EtaBinning, histoCorrectedYieldNormalEff_Right_Normal_Pi0EtaBinning, 1., 1., "B");

	historatioRightNarrow_Pi0EtaBinning = (TH1D*) histoCorrectedYieldNormalEff_Right_Narrow_Pi0EtaBinning -> Clone();
	historatioRightNarrow_Pi0EtaBinning->Divide(historatioRightNarrow_Pi0EtaBinning, histoCorrectedYieldNormalEff_Right_Normal_Pi0EtaBinning, 1., 1., "B");
	historatioRightWide_Pi0EtaBinning = (TH1D*) histoCorrectedYieldNormalEff_Right_Wide_Pi0EtaBinning -> Clone();
	historatioRightWide_Pi0EtaBinning->Divide(historatioRightWide_Pi0EtaBinning, histoCorrectedYieldNormalEff_Right_Normal_Pi0EtaBinning, 1., 1., "B");
	historatioLeftNormal_Pi0EtaBinning = (TH1D*) histoCorrectedYieldNormalEff_Left_Normal_Pi0EtaBinning -> Clone();
	historatioLeftNormal_Pi0EtaBinning->Divide(historatioLeftNormal_Pi0EtaBinning, histoCorrectedYieldNormalEff_Left_Normal_Pi0EtaBinning, 1., 1., "B");
	historatioLeftNarrow_Pi0EtaBinning = (TH1D*) histoCorrectedYieldNormalEff_Left_Narrow_Pi0EtaBinning -> Clone();
	historatioLeftNarrow_Pi0EtaBinning->Divide(historatioLeftNarrow_Pi0EtaBinning, histoCorrectedYieldNormalEff_Left_Normal_Pi0EtaBinning, 1., 1., "B");
	historatioLeftWide_Pi0EtaBinning = (TH1D*) histoCorrectedYieldNormalEff_Left_Wide_Pi0EtaBinning -> Clone();
	historatioLeftWide_Pi0EtaBinning->Divide(historatioLeftWide_Pi0EtaBinning, histoCorrectedYieldNormalEff_Left_Normal_Pi0EtaBinning, 1., 1., "B");

	TLegend* leg2   = GetAndSetLegend2(0.12+0.4, 0.14, 0.45+0.5, 0.30+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
	leg2->AddEntry(historatioRightNormal_Pi0EtaBinning,"#frac{Normal}{Standard} right","lp");
	leg2->AddEntry(historatioRightNarrow_Pi0EtaBinning,"#frac{Narrow}{Standard} right","lp");
	leg2->AddEntry(historatioRightWide_Pi0EtaBinning,"#frac{Wide}{Standard} right","lp");
	leg2->AddEntry(historatioLeftNormal_Pi0EtaBinning,"#frac{Normal}{Standard} left","lp");
	leg2->AddEntry(historatioLeftNarrow_Pi0EtaBinning,"#frac{Narrow}{Standard} left","lp");
	leg2->AddEntry(historatioLeftWide_Pi0EtaBinning,"#frac{Wide}{Standard} left","lp");
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);

	canvas_ratio_diffright->cd();

	historatioRightNormal_Pi0EtaBinning->SetAxisRange(0.,2.,"Y");	
	
	historatioRightNormal_Pi0EtaBinning->Draw(); 
	historatioRightNarrow_Pi0EtaBinning->Draw("same"); 
	historatioRightWide_Pi0EtaBinning->Draw("same"); 
	historatioLeftNormal_Pi0EtaBinning->Draw("same"); 
	historatioLeftNarrow_Pi0EtaBinning->Draw("same"); 
	historatioLeftWide_Pi0EtaBinning->Draw("same"); 
	leg2->Draw("same");	
	canvas_ratio_diffright->Write();


	////////// Acceptance in jets/ Inclusive Jet //////////////////

	TH1D* histoAcceptance_ratio_Pi0EtaBinning;
	TH1D* histoAcceptance_ratio_Eta;

	histoAcceptance_ratio_Pi0EtaBinning = (TH1D*) histoAcceptance_Pi0EtaBinning -> Clone();
	histoAcceptance_ratio_Pi0EtaBinning->Divide(histoAcceptance_ratio_Pi0EtaBinning, histoAcceptance_Pi0EtaBinning_Inclusive, 1., 1., "");

	histoAcceptance_ratio_Eta = (TH1D*) histoAcceptance_Eta -> Clone();
	histoAcceptance_ratio_Eta->Divide(histoAcceptance_ratio_Eta, histoAcceptance_Eta_Inclusive, 1., 1., "");
	
	histoAcceptance_ratio_Pi0EtaBinning->SetTitle("Pi0 ratio Acceptance in jets / Inclusive Jet");
	histoAcceptance_ratio_Eta->SetTitle("Eta ratio Acceptance in jet / Inclusive");
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

	canvas_trigger_Eta->cd();	
	canvas_trigger_Eta->SetLogy();
	TH1D* histoCompareEG1INT7_Eta;
	TH1D* histoCompareEG2INT7_Eta;
	TH1D* histoCompareEG1EG2_Eta;
	TH2D *href6 = new TH2D("href6", "", 100, 0, 20, 100, 0, 1e+2);	
	href6->SetXTitle("#it{p}_{T} (GeV/#it{c})"); 
	href6->SetYTitle("Ratio of raw yield"); 
	
 	histoCompareEG1INT7_Eta = (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
	histoCompareEG2INT7_Eta = (TH1D*) histoYieldMeson_Eta_Data_EG2->Clone();
	histoCompareEG1EG2_Eta = (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
	histoCompareEG1INT7_Eta->Divide(histoCompareEG1INT7_Eta,histoYieldMeson_Eta_Data_INT7, 1, 1, "");
	histoCompareEG2INT7_Eta->Divide(histoCompareEG2INT7_Eta,histoYieldMeson_Eta_Data_INT7, 1, 1, "");
	histoCompareEG1EG2_Eta->Divide(histoCompareEG1EG2_Eta,histoYieldMeson_Eta_Data_EG2, 1, 1, "");
	DrawGammaSetMarker(histoCompareEG1INT7_Eta, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7_Eta, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2_Eta, 26, 1., 1, 1);
	TLegend* legendComparetrigger_Eta   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_Eta->AddEntry(histoCompareEG1INT7_Eta, "#frac{EG1}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Eta->AddEntry(histoCompareEG2INT7_Eta, "#frac{EG2}{INT7}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Eta->AddEntry(histoCompareEG1EG2_Eta, "#frac{EG1}{EG2}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	href6->Draw();
	histoCompareEG1INT7_Eta->Draw("samep");
	histoCompareEG2INT7_Eta->Draw("samep");
	histoCompareEG1EG2_Eta->Draw("samep");
	legendComparetrigger_Eta->Draw("same");
	histoCompareEG1INT7_Eta->Write();	
	
	///////// To compare each trigger (INT7, EG1, EG@) for raw yield ///////
	canvas_trigger_Pi0Yield->Divide(2,2);
	canvas_trigger_Pi0Yield->SetLogy();
	TH2D *href3 = new TH2D("href3","", 100, 0, 20, 100, 1e-8, 1e-4);
	href3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	href3->GetYaxis()->SetTitle("RAW Yield/N_{Evt}");
	TH1D* histoComparePi0YieldData_INT7;
	TH1D* histoComparePi0YieldData_EG1;
	TH1D* histoComparePi0YieldData_EG2;
	histoComparePi0YieldData_INT7 = (TH1D*) histoYieldMeson_Pi0_Data_INT7->Clone();
	histoComparePi0YieldData_EG1 = (TH1D*) histoYieldMeson_Pi0_Data_EG1->Clone();
	histoComparePi0YieldData_EG2 = (TH1D*) histoYieldMeson_Pi0_Data_EG2->Clone();

	TLegend* legendComparetrigger_Pi0YieldData   = GetAndSetLegend2(0.12+0.1, 0.24, 0.45+0.1, 0.24+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_Pi0YieldData->AddEntry(histoComparePi0YieldData_INT7, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData->AddEntry(histoComparePi0YieldData_EG1, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData->AddEntry(histoComparePi0YieldData_EG2, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...

	DrawGammaSetMarker(histoComparePi0YieldData_INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(histoComparePi0YieldData_EG1, 25, 1., 2, 2);
	DrawGammaSetMarker(histoComparePi0YieldData_EG2, 26, 1., 1, 1);
	canvas_trigger_Pi0Yield->cd(1);
	//canvas_trigger_Pi0Yield->SetTicks(0,1);
	gPad->SetLogy(1);
	href3->Draw();
 	histoComparePi0YieldData_INT7->Draw("samep");
	histoComparePi0YieldData_EG1->Draw("samep");
	histoComparePi0YieldData_EG2->Draw("samep");
	legendComparetrigger_Pi0YieldData->Draw("same");
	
	canvas_trigger_Pi0Yield->cd(3);
	gPad->SetLogy(1);
	TH2D *hhref3 = new TH2D("hhref3","", 100, 0, 20, 100, 7e-9, 1e-2);
	hhref3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hhref3->GetYaxis()->SetTitle("RAW Yield/N_{Evt}/RF_{Meson}");
	TH1D *hhistoComparePi0YieldData_INT7;
	TH1D *hhistoComparePi0YieldData_EG1;
	TH1D *hhistoComparePi0YieldData_EG2;
	hhistoComparePi0YieldData_INT7= (TH1D*) histoYieldMeson_Pi0_Data_INT7->Clone();
	hhistoComparePi0YieldData_EG1= (TH1D*) histoYieldMeson_Pi0_Data_EG1->Clone();
	hhistoComparePi0YieldData_EG2= (TH1D*) histoYieldMeson_Pi0_Data_EG2->Clone();
	hhistoComparePi0YieldData_EG1->Scale(1/RFEG1Pi0YieldData);
	hhistoComparePi0YieldData_EG2->Scale(1/RFEG2Pi0YieldData);
	DrawGammaSetMarker(hhistoComparePi0YieldData_INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(hhistoComparePi0YieldData_EG1, 25, 1., 2, 2);
	DrawGammaSetMarker(hhistoComparePi0YieldData_EG2, 26, 1., 1, 1);
	hhref3->Draw();
 	hhistoComparePi0YieldData_INT7->Draw("samep");
	hhistoComparePi0YieldData_EG1->Draw("samep");
	hhistoComparePi0YieldData_EG2->Draw("samep");
	//legendComparetrigger_Pi0YieldData->Draw("same");

	canvas_trigger_Pi0Yield->cd(2);	
	gPad->SetLogy(1);
	TH2D *hhhref3 = new TH2D("hhref3","", 100, 0, 20, 100, 1e-9, 1e-1);
	hhhref3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hhhref3->GetYaxis()->SetTitle("RAW Yield/N_{evt} for MC");
	TH1D *histoComparePi0YieldMC_INT7;
	TH1D *histoComparePi0YieldMC_EG1;
	TH1D *histoComparePi0YieldMC_EG2;
	histoComparePi0YieldMC_INT7= (TH1D*) histoYieldMeson_Pi0_MC_INT7->Clone();
	histoComparePi0YieldMC_EG1= (TH1D*) histoYieldMeson_Pi0_MC_EG1->Clone();
	histoComparePi0YieldMC_EG2= (TH1D*) histoYieldMeson_Pi0_MC_EG2->Clone();
	DrawGammaSetMarker(histoComparePi0YieldMC_INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(histoComparePi0YieldMC_EG1, 25, 1., 2, 2);
	DrawGammaSetMarker(histoComparePi0YieldMC_EG2, 26, 1., 1, 1);
	hhhref3->Draw();
 	histoComparePi0YieldMC_INT7->Draw("samep");
	histoComparePi0YieldMC_EG1->Draw("samep");
	histoComparePi0YieldMC_EG2->Draw("samep");
	//legendComparetrigger_Pi0Yield->Draw("same");
	canvas_trigger_Pi0Yield->cd(4);
	gPad->SetLogy(1);
	TH2D *hhhhref3 = new TH2D("hhhref3","", 100, 0, 20, 100, 1e-8, 1e-4);
	hhhhref3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hhhhref3->GetYaxis()->SetTitle("RAW Yield/N_{Evt}/RF_{Meson} for MC");
	TH1D *hhistoComparePi0YieldMC_INT7;
	TH1D *hhistoComparePi0YieldMC_EG1;
	TH1D *hhistoComparePi0YieldMC_EG2;
	hhistoComparePi0YieldMC_INT7= (TH1D*) histoYieldMeson_Pi0_MC_INT7->Clone();
	hhistoComparePi0YieldMC_EG1= (TH1D*) histoYieldMeson_Pi0_MC_EG1->Clone();
	hhistoComparePi0YieldMC_EG2= (TH1D*) histoYieldMeson_Pi0_MC_EG2->Clone();
//	hhistoComparePi0YieldMC_EG1->Scale(NEvtMCEG1);
//	hhistoComparePi0YieldMC_EG1->Scale(1/NEvtMCINT7);
//	hhistoComparePi0YieldMC_EG2->Scale(NEvtMCEG2);
//	hhistoComparePi0YieldMC_EG2->Scale(1/NEvtMCINT7);
//	hhistoComparePi0YieldMC_EG1->Scale(1/RFEG1Pi0ClusterMC);
//	hhistoComparePi0YieldMC_EG2->Scale(1/RFEG2Pi0ClusterMC);
	hhistoComparePi0YieldMC_EG1->Scale(1/RFEG1Pi0YieldMC); // N_{Event}/N_{EG1}
	hhistoComparePi0YieldMC_EG2->Scale(1/RFEG2Pi0YieldMC); // N_{Event}/N_{EG1}
	DrawGammaSetMarker(hhistoComparePi0YieldMC_INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(hhistoComparePi0YieldMC_EG1, 25, 1., 2, 2);
	DrawGammaSetMarker(hhistoComparePi0YieldMC_EG2, 26, 1., 1, 1);
	
	hhhhref3->Draw();
 	hhistoComparePi0YieldMC_INT7->Draw("samep");
	hhistoComparePi0YieldMC_EG1->Draw("samep");
	hhistoComparePi0YieldMC_EG2->Draw("samep");
	
	canvas_trigger_Pi0Yield->Write();
	canvas_trigger_Pi0Yield_Inclusive->cd();
	gPad->SetLogy(1);
	TH2D *hhref3_Inclusive = new TH2D("hhref3","", 100, 0, 20, 100, 7e-9, 1e-2);
	hhref3_Inclusive->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hhref3_Inclusive->GetYaxis()->SetTitle("RAW Yield/N_{Evt}/RF_{cluster}");
	TH1D *hhistoComparePi0YieldData_INT7_Inclusive;
	TH1D *hhistoComparePi0YieldData_EG1_Inclusive;
	TH1D *hhistoComparePi0YieldData_EG2_Inclusive;
	hhistoComparePi0YieldData_INT7_Inclusive= (TH1D*) histoYieldMeson_Pi0_Data_INT7_Inclusive->Clone();
	hhistoComparePi0YieldData_EG1_Inclusive= (TH1D*) histoYieldMeson_Pi0_Data_EG1_Inclusive->Clone();
	hhistoComparePi0YieldData_EG2_Inclusive= (TH1D*) histoYieldMeson_Pi0_Data_EG2_Inclusive->Clone();
	hhistoComparePi0YieldData_EG1_Inclusive->Scale(1/RFEG1_Data_Inclusive);
	hhistoComparePi0YieldData_EG2_Inclusive->Scale(1/RFEG2_Data_Inclusive);
	
	TLegend* legendComparetrigger_Pi0YieldData_Inclusive   = GetAndSetLegend2(0.12+0.1, 0.24, 0.45+0.1, 0.24+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_Pi0YieldData_Inclusive->AddEntry(hhistoComparePi0YieldData_INT7_Inclusive, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData_Inclusive->AddEntry(hhistoComparePi0YieldData_EG1_Inclusive, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData_Inclusive->AddEntry(hhistoComparePi0YieldData_EG2_Inclusive, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...

	DrawGammaSetMarker(hhistoComparePi0YieldData_INT7_Inclusive, 24, 1., 4, 4);
	DrawGammaSetMarker(hhistoComparePi0YieldData_EG1_Inclusive, 25, 1., 2, 2);
	DrawGammaSetMarker(hhistoComparePi0YieldData_EG2_Inclusive, 26, 1., 1, 1);
	hhref3_Inclusive->Draw();
 	hhistoComparePi0YieldData_INT7_Inclusive->Draw("samep");
	hhistoComparePi0YieldData_EG1_Inclusive->Draw("samep");
	hhistoComparePi0YieldData_EG2_Inclusive->Draw("samep");
	legendComparetrigger_Pi0YieldData_Inclusive->Draw("same");
	//legendComparetrigger_Pi0YieldData->Draw("same");

	// define histogram for CompareYieldOfTriggerForEtaWithRatio
	TH1D* histoCompareEG1INT7EtaYieldData;
	TH1D* histoCompareEG2INT7EtaYieldData;
	TH1D* histoCompareEG1EG2EtaYieldData;
	histoCompareEG1INT7EtaYieldData = (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
	histoCompareEG2INT7EtaYieldData = (TH1D*) histoYieldMeson_Eta_Data_EG2->Clone();
	histoCompareEG1EG2EtaYieldData = (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
 	histoCompareEG1INT7EtaYieldData->Divide(histoCompareEG1INT7EtaYieldData,histoYieldMeson_Eta_Data_INT7, 1, 1, "");
	histoCompareEG2INT7EtaYieldData->Divide(histoCompareEG2INT7EtaYieldData,histoYieldMeson_Eta_Data_INT7, 1, 1, "");
	histoCompareEG1EG2EtaYieldData->Divide(histoCompareEG1EG2EtaYieldData,histoYieldMeson_Eta_Data_EG2, 1, 1, "");
	TF1 *fYieldEG1INT7EtaData = new TF1("fYieldEG1INT7Eta","[0]", 0, 100);
	TF1 *fYieldEG2INT7EtaData = new TF1("fYieldEG2INT7Eta","[0]", 0, 100);
	TF1 *fYieldEG1EG2EtaData  = new TF1("fYieldEG1EG2Eta","[0]", 0, 100);
	histoCompareEG1INT7EtaYieldData->Fit("fYieldEG1INT7Eta","","",12,20);	
	histoCompareEG2INT7EtaYieldData->Fit("fYieldEG2INT7Eta","","",6,20);	
	histoCompareEG1EG2EtaYieldData ->Fit("fYieldEG1EG2Eta","","",12,20);	
	DrawGammaSetMarker(histoCompareEG1INT7EtaYieldData, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEG2INT7EtaYieldData, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEG1EG2EtaYieldData, 26, 1., 1, 1);
	Double_t RFEG1EtaYieldData = fYieldEG1INT7EtaData->GetParameter(0);
	Double_t RFEG2EtaYieldData = fYieldEG2INT7EtaData->GetParameter(0);
	Double_t RFEGX2EtaYieldData = fYieldEG1EG2EtaData->GetParameter(0);
	Double_t RFEG1EtaYieldDataError = fYieldEG1INT7EtaData->GetParError(0);
	Double_t RFEG2EtaYieldDataError = fYieldEG2INT7EtaData->GetParError(0);
	fYieldEG1INT7EtaData->Write("EG1INT7RejectionFactorOfYieldForEta");		
	fYieldEG2INT7EtaData->Write("EG2INT7RejectionFactorOfYieldForEta");		
	fYieldEG1EG2EtaData->Write("EG1EG2RejectionFactorOfYieldForEta");		


   
	TH1D *histoCompareEtaYieldData_INT7;
	TH1D *histoCompareEtaYieldData_EG1;
	TH1D *histoCompareEtaYieldData_EG2;
	histoCompareEtaYieldData_INT7= (TH1D*) histoYieldMeson_Eta_Data_INT7->Clone();
	histoCompareEtaYieldData_EG1= (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
	histoCompareEtaYieldData_EG2= (TH1D*) histoYieldMeson_Eta_Data_EG2->Clone();
	histoCompareEtaYieldData_EG1->Scale(1/RFEG1EtaYieldData);
	histoCompareEtaYieldData_EG2->Scale(1/RFEG2EtaYieldData);
	DrawGammaSetMarker(histoCompareEtaYieldData_INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEtaYieldData_EG1, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEtaYieldData_EG2, 26, 1., 1, 1);
	
	TLegend* legendComparetrigger_EtaYieldData_WithRatio   = GetAndSetLegend2(0.12+0.5, 0.24+0.5, 0.45+0.5, 0.24+0.5+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_EtaYieldData_WithRatio->AddEntry(histoCompareEtaYieldData_INT7, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_EtaYieldData_WithRatio->AddEntry(histoCompareEtaYieldData_EG1, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_EtaYieldData_WithRatio->AddEntry(histoCompareEtaYieldData_EG2, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	// define eta meson of inclusoive measurement 
//	TH1D* histoCompareEG1INT7EtaYieldData_Inclusive;
//	TH1D* histoCompareEG2INT7EtaYieldData_Inclusive;
//	TH1D* histoCompareEG1EG2EtaYieldData_Inclusive;
//	histoCompareEG1INT7EtaYieldData_Inclusive = (TH1D*) histoYieldMeson_Eta_Data_EG1_Inclusive->Clone();
//	histoCompareEG2INT7EtaYieldData_Inclusive = (TH1D*) histoYieldMeson_Eta_Data_EG2_Inclusive->Clone();
//	histoCompareEG1EG2EtaYieldData_Inclusive = (TH1D*) histoYieldMeson_Eta_Data_EG1_Inclusive->Clone();
// 	histoCompareEG1INT7EtaYieldData->Divide(histoCompareEG1INT7EtaYieldData,histoYieldMeson_Eta_Data_INT7, 1, 1, "");
//	histoCompareEG2INT7EtaYieldData->Divide(histoCompareEG2INT7EtaYieldData,histoYieldMeson_Eta_Data_INT7, 1, 1, "");
//	histoCompareEG1EG2EtaYieldData->Divide(histoCompareEG1EG2EtaYieldData,histoYieldMeson_Eta_Data_EG2, 1, 1, "");
//	TF1 *fYieldEG1INT7EtaData = new TF1("fYieldEG1INT7Eta","[0]", 0, 100);
//	TF1 *fYieldEG2INT7EtaData = new TF1("fYieldEG2INT7Eta","[0]", 0, 100);
//	TF1 *fYieldEG1EG2EtaData  = new TF1("fYieldEG1EG2Eta","[0]", 0, 100);
//	histoCompareEG1INT7EtaYieldData->Fit("fYieldEG1INT7Eta","","",12,20);	
//	histoCompareEG2INT7EtaYieldData->Fit("fYieldEG2INT7Eta","","",6,20);	
//	histoCompareEG1EG2EtaYieldData ->Fit("fYieldEG1EG2Eta","","",12,20);	
//	DrawGammaSetMarker(histoCompareEG1INT7EtaYieldData, 24, 1., 4, 4);
//	DrawGammaSetMarker(histoCompareEG2INT7EtaYieldData, 25, 1., 2, 2);
//	DrawGammaSetMarker(histoCompareEG1EG2EtaYieldData, 26, 1., 1, 1);
//	Double_t RFEG1EtaYieldData = fYieldEG1INT7EtaData->GetParameter(0);
//	Double_t RFEG2EtaYieldData = fYieldEG2INT7EtaData->GetParameter(0);
//	Double_t RFEGX2EtaYieldData = fYieldEG1EG2EtaData->GetParameter(0);
//	Double_t RFEG1EtaYieldDataError = fYieldEG1INT7EtaData->GetParError(0);
//	Double_t RFEG2EtaYieldDataError = fYieldEG2INT7EtaData->GetParError(0);
//   
//	TH1D *histoCompareEtaYieldData_INT7;
//	TH1D *histoCompareEtaYieldData_EG1;
//	TH1D *histoCompareEtaYieldData_EG2;
//	histoCompareEtaYieldData_INT7= (TH1D*) histoYieldMeson_Eta_Data_INT7->Clone();
//	histoCompareEtaYieldData_EG1= (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
//	histoCompareEtaYieldData_EG2= (TH1D*) histoYieldMeson_Eta_Data_EG2->Clone();
//	histoCompareEtaYieldData_EG1->Scale(1/RFEG1EtaYieldData);
//	histoCompareEtaYieldData_EG2->Scale(1/RFEG2EtaYieldData);
//	DrawGammaSetMarker(histoCompareEtaYieldData_INT7, 24, 1., 4, 4);
//	DrawGammaSetMarker(histoCompareEtaYieldData_EG1, 25, 1., 2, 2);
//	DrawGammaSetMarker(histoCompareEtaYieldData_EG2, 26, 1., 1, 1);
//	
//	TLegend* legendComparetrigger_EtaYieldData_WithRatio   = GetAndSetLegend2(0.12+0.5, 0.24+0.5, 0.45+0.5, 0.24+0.5+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
// 	legendComparetrigger_EtaYieldData_WithRatio->AddEntry(histoCompareEtaYieldData_INT7, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
// 	legendComparetrigger_EtaYieldData_WithRatio->AddEntry(histoCompareEtaYieldData_EG1, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
// 	legendComparetrigger_EtaYieldData_WithRatio->AddEntry(histoCompareEtaYieldData_EG2, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...


	// Canvas CompareYieldOfTriggerForPi0WithRatio
	TCanvas* canvasRawYield_Pi0 = new TCanvas("canvasRawYield","",1350,1500);// gives the page size
    DrawGammaCanvasSettings( canvasRawYield_Pi0, 0.13, 0.02, 0.02, 0.09);
    canvasRawYield_Pi0->SetLogy();

    TPad* padRawYieldHistosPi0 = new TPad("padRawYieldHistosPi0", "", 0., 0.4, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldHistosPi0, 0.12, 0.02, 0.02, 0.);
    padRawYieldHistosPi0->Draw();
	
    TPad* padRawYieldRatiosPi0 = new TPad("padRawYieldRatiosPi0", "", 0., 0.2, 1., 0.4,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldRatiosPi0, 0.12, 0.02, 0., 0.18);
    padRawYieldRatiosPi0->Draw();
	
	TPad* padClusterRatiosPi0 = new TPad("padClusterRatiosPi0", "", 0., 0., 1., 0.2,-1, -1, -2);
    DrawGammaPadSettings( padClusterRatiosPi0, 0.12, 0.02, 0., 0.18);
    padClusterRatiosPi0->Draw();


   	TLegend* legendComparetrigger_Pi0YieldData_WithRatio   = GetAndSetLegend2(0.12+0.5, 0.24+0.5, 0.45+0.5, 0.24+0.5+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_Pi0YieldData_WithRatio->AddEntry(hhistoComparePi0YieldData_INT7, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData_WithRatio->AddEntry(hhistoComparePi0YieldData_EG1, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData_WithRatio->AddEntry(hhistoComparePi0YieldData_EG2, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
    padRawYieldHistosPi0->cd();
    padRawYieldHistosPi0->SetLogy();
   	TH2D *histoRawYieldPi0Dummy = new TH2D("histoRawYieldPi0Dummy","", 100, 0, 20, 100, 7e-9, 1e-2);
	SetStyleHistoTH2ForGraphs(histoRawYieldPi0Dummy, "#it{p}_{T} (GeV/#it{c})", "RAW Yield/N_{Evt}/RF_{Meson}", 0.033,0.04, 0.033,0.04, 1,1.35);
	histoRawYieldPi0Dummy->Draw();
	hhistoComparePi0YieldData_INT7->Draw("same");
	hhistoComparePi0YieldData_EG1->Draw("same");
	hhistoComparePi0YieldData_EG2->Draw("same");
	legendComparetrigger_Pi0YieldData_WithRatio->Draw("same");
	latexComparetriggerPi0YieldData->DrawLatex(2, 5e-8, Form("EG1/INT7: %3.1f #pm %3.2f", RFEG1Pi0YieldData, RFEG1Pi0YieldDataError));
	latexComparetriggerPi0YieldData->DrawLatex(2, 3e-8, Form("EG2/INT7: %3.1f #pm %3.2f", RFEG2Pi0YieldData, RFEG2Pi0YieldDataError));

    padRawYieldRatiosPi0->cd();
    padRawYieldRatiosPi0->SetTickx();
    padRawYieldRatiosPi0->SetTicky();
    padRawYieldRatiosPi0->SetLogy(1);
	
	TH2D *histoRejectionFactorOfPi0Dummy = new TH2D("histoRejectionFactorOfPi0Dummy", "", 100, 0, 20, 100, 0, 3e+2);	
	SetStyleHistoTH2ForGraphs(histoRejectionFactorOfPi0Dummy, "#it{p}_{T} (GeV/#it{c})", "RF_{meson}",0.07,0.1,0.07,0.1,0.8,0.55,510,505 );
	histoRejectionFactorOfPi0Dummy->Draw();
	histoCompareEG2INT7Pi0YieldData->Draw("same");
	histoCompareEG1EG2Pi0YieldData->Draw("same");
	histoCompareEG1INT7Pi0YieldData->Draw("same");

    padClusterRatiosPi0->cd();
    padClusterRatiosPi0->SetTickx();
    padClusterRatiosPi0->SetTicky();
    padClusterRatiosPi0->SetLogy();
	TH2D *histoRejectionFactorOfClusterPi0Dummy = new TH2D("histoRejectionFactorOfClusterPi0Dummy", "", 100, 0, 100, 100, 1e-1, 3e+6);
	SetStyleHistoTH2ForGraphs(histoRejectionFactorOfClusterPi0Dummy, "#it{E} (GeV)", "RF_{Cluster}",0.07,0.1,0.07,0.1,0.8,0.55,510,505 );
	histoRejectionFactorOfClusterPi0Dummy->Draw();
	histoCompareEG1INT7->Draw("same");
	histoCompareEG2INT7->Draw("same");
	histoCompareEG1EG2->Draw("same");

    canvasRawYield_Pi0->Update();
	
	TCanvas* canvasRawYield_Pi0_Inclusive = new TCanvas("canvasRawYield_Pi0_Inclusive","",1350,1500);// gives the page size
    DrawGammaCanvasSettings( canvasRawYield_Pi0_Inclusive, 0.13, 0.02, 0.02, 0.09);
    canvasRawYield_Pi0_Inclusive->SetLogy();

    TPad* padRawYieldHistosPi0_Inclusive = new TPad("padRawYieldHistosPi0_Inclusive", "", 0., 0.4, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldHistosPi0_Inclusive, 0.12, 0.02, 0.02, 0.);
    padRawYieldHistosPi0_Inclusive->Draw();
	
    TPad* padRawYieldRatiosPi0_Inclusive = new TPad("padRawYieldRatiosPi0_Inclusive", "", 0., 0.2, 1., 0.4,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldRatiosPi0_Inclusive, 0.12, 0.02, 0., 0.18);
    padRawYieldRatiosPi0_Inclusive->Draw();
	
	TPad* padClusterRatiosPi0_Inclusive = new TPad("padClusterRatiosPi0_Inclusive", "", 0., 0., 1., 0.2,-1, -1, -2);
    DrawGammaPadSettings( padClusterRatiosPi0_Inclusive, 0.12, 0.02, 0., 0.18);
    padClusterRatiosPi0_Inclusive->Draw();


   	TLegend* legendComparetrigger_Pi0YieldData_WithRatio_Inclusive   = GetAndSetLegend2(0.12+0.5, 0.24+0.5, 0.45+0.5, 0.24+0.5+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_Pi0YieldData_WithRatio_Inclusive->AddEntry(hhistoComparePi0YieldData_INT7_Inclusive, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData_WithRatio_Inclusive->AddEntry(hhistoComparePi0YieldData_EG1_Inclusive, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_Pi0YieldData_WithRatio_Inclusive->AddEntry(hhistoComparePi0YieldData_EG2_Inclusive, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...

    padRawYieldHistosPi0_Inclusive->cd();
    padRawYieldHistosPi0_Inclusive->SetLogy();
	TH2D *histoRawYieldPi0Dummy_Inclusive = new TH2D("histoRawYieldPi0Dummy_Inclusive","", 100, 0, 20, 100, 7e-9, 1e-2);
	SetStyleHistoTH2ForGraphs(histoRawYieldPi0Dummy_Inclusive, "#it{p}_{T} (GeV/#it{c})", "RAW Yield/N_{Evt}/RF_{Meson}", 0.033,0.04, 0.033,0.04, 1,1.35);
	histoRawYieldPi0Dummy_Inclusive->Draw();
	hhistoComparePi0YieldData_INT7_Inclusive->Draw("same");
	hhistoComparePi0YieldData_EG1_Inclusive->Draw("same");
	hhistoComparePi0YieldData_EG2_Inclusive->Draw("same");
	legendComparetrigger_Pi0YieldData_WithRatio_Inclusive->Draw("same");
	latexComparetriggerPi0YieldData_Inclusive->DrawLatex(2, 5e-8, Form("EG1/INT7: %3.1f #pm %3.1f", RFEG1_Data_Inclusive, RFEG1_DataError_Inclusive));
	latexComparetriggerPi0YieldData_Inclusive->DrawLatex(2, 3e-8, Form("EG2/INT7: %3.1f #pm %3.1f", RFEG2_Data_Inclusive, RFEG2_DataError_Inclusive));

    padRawYieldRatiosPi0_Inclusive->cd();
    padRawYieldRatiosPi0_Inclusive->SetTickx();
    padRawYieldRatiosPi0_Inclusive->SetTicky();
    padRawYieldRatiosPi0_Inclusive->SetLogy(1);
	TH2D *histoRejectionFactorOfMesonPi0Dummy_Inclusive = new TH2D("histoRejectionFactorOfMesonPi0Dummy_Inclusive", "", 100, 0, 20, 100, 1e-1, 3e+6);	
	SetStyleHistoTH2ForGraphs(histoRejectionFactorOfMesonPi0Dummy_Inclusive, "#it{p}_{T} (GeV/#it{c})", "RF_{meson}",0.07,0.1,0.07,0.1,0.8,0.55,510,505 );
	histoRejectionFactorOfMesonPi0Dummy_Inclusive->Draw();
	histoCompareEG2INT7Pi0YieldData_Inclusive->Draw("same");
	histoCompareEG1EG2Pi0YieldData_Inclusive->Draw("same");
	histoCompareEG1INT7Pi0YieldData_Inclusive->Draw("same");

    padClusterRatiosPi0_Inclusive->cd();
    padClusterRatiosPi0_Inclusive->SetTickx();
    padClusterRatiosPi0_Inclusive->SetTicky();
    padClusterRatiosPi0_Inclusive->SetLogy();
	TH2D *histoRejectionFactorOfClusterPi0Dummy_Inclusive = new TH2D("histoRejectionFactorOfClusterPi0Dummy_Inclusive", "", 100, 0, 100, 100, 1e-1, 3e+6);
	SetStyleHistoTH2ForGraphs(histoRejectionFactorOfClusterPi0Dummy_Inclusive, "#it{E} (GeV)", "RF_{Cluster}",0.07,0.1,0.07,0.1,0.8,0.55,510,505 );
	histoRejectionFactorOfClusterPi0Dummy_Inclusive->Draw();
	histoCompareEG1INT7_Inclusive->Draw("same");
	histoCompareEG2INT7_Inclusive->Draw("same");
	histoCompareEG1EG2_Inclusive->Draw("same");

    canvasRawYield_Pi0_Inclusive->Update();

	// compare cluster with raw yield for eta meson's rejection factor 
   	
	TCanvas* canvasRawYield_Eta = new TCanvas("canvasRawYield_Eta","",1350,1500);// gives the page size
    DrawGammaCanvasSettings( canvasRawYield_Eta, 0.13, 0.02, 0.02, 0.09);
    canvasRawYield_Eta->SetLogy();

    TPad* padRawYieldHistosEta = new TPad("padRawYieldHistosEta", "", 0., 0.4, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldHistosEta, 0.12, 0.02, 0.02, 0.);
    padRawYieldHistosEta->Draw();
	
    TPad* padRawYieldRatiosEta = new TPad("padRawYieldRatiosEta", "", 0., 0.2, 1., 0.4,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldRatiosEta, 0.12, 0.02, 0., 0.18);
    padRawYieldRatiosEta->Draw();
	
	TPad* padClusterRatiosEta = new TPad("padClusterRatiosEta", "", 0., 0., 1., 0.2,-1, -1, -2);
    DrawGammaPadSettings( padClusterRatiosEta, 0.12, 0.02, 0., 0.18);
    padClusterRatiosEta->Draw();
 
	padRawYieldHistosEta->cd();
    padRawYieldHistosEta->SetLogy();
   	TH2D *histoRawYieldEtaDummy = new TH2D("histoRawYieldEtaDummy","", 100, 0, 20, 100, 7e-9, 1e-5);
	SetStyleHistoTH2ForGraphs(histoRawYieldEtaDummy, "#it{p}_{T} (GeV/#it{c})", "RAW Yield/N_{Evt}/RF_{Meson}", 0.033,0.04, 0.033,0.04, 1,1.35);
	histoRawYieldEtaDummy->Draw();
	histoCompareEtaYieldData_INT7->Draw("same");
	histoCompareEtaYieldData_EG1->Draw("same");
	histoCompareEtaYieldData_EG2->Draw("same");
	legendComparetrigger_EtaYieldData_WithRatio->Draw("same");
	latexComparetriggerEtaYieldData->DrawLatex(2, 5e-8, Form("EG1/INT7: %3.1f #pm %3.2f", RFEG1EtaYieldData, RFEG1EtaYieldDataError));
	latexComparetriggerEtaYieldData->DrawLatex(2, 3e-8, Form("EG2/INT7: %3.1f #pm %3.2f", RFEG2EtaYieldData, RFEG2EtaYieldDataError));

    padRawYieldRatiosEta->cd();
    padRawYieldRatiosEta->SetTickx();
    padRawYieldRatiosEta->SetTicky();
    padRawYieldRatiosEta->SetLogy(1);
	
	TH2D *histoRejectionFactorOfEtaDummy = new TH2D("histoRejectionFactorOfEtaDummy", "", 100, 0, 20, 100, 0, 7e+1);	
	SetStyleHistoTH2ForGraphs(histoRejectionFactorOfEtaDummy, "#it{p}_{T} (GeV/#it{c})", "RF_{meson}",0.07,0.1,0.07,0.1,0.8,0.55,510,505 );
	histoRejectionFactorOfEtaDummy->Draw();
	histoCompareEG2INT7EtaYieldData->Draw("same");
	//histoCompareEG1EG2EtaYieldData->Draw("same");
	histoCompareEG1INT7EtaYieldData->Draw("same");

    padClusterRatiosEta->cd();
    padRawYieldRatiosEta->SetTickx();
    padRawYieldRatiosEta->SetTicky();
    padClusterRatiosEta->SetLogy();
	TH2D *histoRejectionFactorOfClusterEtaDummy = new TH2D("histoRejectionFactorOfClusterEtaDummy", "", 100, 0, 100, 100, 1e-1, 3e+6);
	SetStyleHistoTH2ForGraphs(histoRejectionFactorOfClusterEtaDummy, "#it{E} (GeV)", "RF_{Cluster}",0.07,0.1,0.07,0.1,0.8,0.55,510,505 );
	histoRejectionFactorOfClusterEtaDummy->Draw();
	histoCompareEG1INT7->Draw("same");
	histoCompareEG2INT7->Draw("same");
	histoCompareEG1EG2->Draw("same");


    canvasRawYield_Eta->Update();

	canvas_trigger_EtaYield->cd();
	canvas_trigger_EtaYield->SetLogy();
	TH2D *href4 = new TH2D("href4", "", 100, 0, 40, 100, 1e-10, 1e-4);
	href4->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	href4->GetYaxis()->SetTitle("RAW Yield/N_{Evt}");
	TH1D* histoCompareEtaYield_INT7;
	TH1D* histoCompareEtaYield_EG1;
	TH1D* histoCompareEtaYield_EG2;
	histoCompareEtaYield_INT7 = (TH1D*) histoYieldMeson_Eta_Data_INT7->Clone();
	histoCompareEtaYield_EG1 = (TH1D*) histoYieldMeson_Eta_Data_EG1->Clone();
	histoCompareEtaYield_EG2 = (TH1D*) histoYieldMeson_Eta_Data_EG2->Clone();
	TLegend* legendComparetrigger_EtaYield   = GetAndSetLegend2(0.12, 0.74, 0.45, 0.74+(0.15*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendComparetrigger_EtaYield->AddEntry(histoCompareEtaYield_INT7, "INT7", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_EtaYield->AddEntry(histoCompareEtaYield_EG1, "EG1", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendComparetrigger_EtaYield->AddEntry(histoCompareEtaYield_EG2, "EG2", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	href4->Draw();
 	histoCompareEtaYield_INT7->Draw("samep");
	histoCompareEtaYield_EG1->Draw("samep");
	histoCompareEtaYield_EG2->Draw("samep");
	legendComparetrigger_EtaYield->Draw("same");
	DrawGammaSetMarker(histoCompareEtaYield_INT7, 24, 1., 4, 4);
	DrawGammaSetMarker(histoCompareEtaYield_EG1, 25, 1., 2, 2);
	DrawGammaSetMarker(histoCompareEtaYield_EG2, 26, 1., 1, 1);
	
	canvas_trigger_EtaYield->Write();

	/////////// Efficiency in jets / Inclusive Jet ///////////////
	

	TH1D* histoEfficiency_ratio_Pi0EtaBinning;
	TH1D* histoEfficiency_ratio_Pi0EtaBinning_True;
	TH1D* histoEfficiency_ratio_Eta;
	TH1D* histoEfficiency_ratio_Eta_True;
	
	cout << "Debug: " << __LINE__ <<endl;
	canvas_ratio_Efficiency_Pi0EtaBinning->cd(); // Draw the Efficiency ratio_Pi0 in jet / Inclusive jet in canvas

	histoEfficiency_ratio_Pi0EtaBinning = (TH1D*) histoEfficiency_Pi0EtaBinning -> Clone();
	histoEfficiency_ratio_Pi0EtaBinning->Divide(histoEfficiency_ratio_Pi0EtaBinning, histoEfficiency_Pi0EtaBinning_Inclusive, 1., 1., "");

	histoEfficiency_ratio_Pi0EtaBinning_True = (TH1D*) histoEfficiency_Pi0EtaBinning_True -> Clone();
	histoEfficiency_ratio_Pi0EtaBinning_True->Divide(histoEfficiency_ratio_Pi0EtaBinning_True, histoEfficiency_Pi0EtaBinning_True_Inclusive, 1., 1., "");

	DrawGammaSetMarker(histoEfficiency_ratio_Pi0EtaBinning, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoEfficiency_ratio_Pi0EtaBinning_True, 20, 1., kGreen, kGreen);
	histoEfficiency_ratio_Pi0EtaBinning->SetAxisRange(0.,2.,"Y");
	histoEfficiency_ratio_Pi0EtaBinning->SetYTitle("Ratio #varepsilon_{#pi^{0}}");
	histoEfficiency_ratio_Pi0EtaBinning->SetXTitle("p_{T}");
	
	TLegend* legendratioEfficiency_Pi0EtaBinning   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendratioEfficiency_Pi0EtaBinning->AddEntry(histoEfficiency_ratio_Pi0EtaBinning, "#frac{Rec #varepsilon_{jet}}{Rec #varepsilon_{Inclusive}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioEfficiency_Pi0EtaBinning->AddEntry(histoEfficiency_ratio_Pi0EtaBinning_True, "#frac{True #varepsilon_{Jet}}{True #varepsilon_{Inclusive}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	
	histoEfficiency_ratio_Pi0EtaBinning->Draw();
	histoEfficiency_ratio_Pi0EtaBinning_True->Draw("same");	
	legendratioEfficiency_Pi0EtaBinning->Draw("same");
	f1->Draw("same");
	canvas_ratio_Efficiency_Pi0EtaBinning->Write();

	canvas_ratio_Efficiency_Eta->cd(); // Draw the Efficiency ratio_Eta in jet / Inclusive Jet in canvas
	
	histoEfficiency_ratio_Eta = (TH1D*) histoEfficiency_Eta -> Clone();
	histoEfficiency_ratio_Eta->Divide(histoEfficiency_ratio_Eta, histoEfficiency_Eta_Inclusive, 1., 1., "");

	histoEfficiency_ratio_Eta_True = (TH1D*) histoEfficiency_Eta_True -> Clone();
	histoEfficiency_ratio_Eta_True->Divide(histoEfficiency_ratio_Eta_True, histoEfficiency_Eta_True_Inclusive, 1., 1., "");
	
	DrawGammaSetMarker(histoEfficiency_ratio_Eta, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoEfficiency_ratio_Eta_True, 20, 1., kGreen, kGreen);
	histoEfficiency_ratio_Eta->SetAxisRange(0.,2.,"Y");
	histoEfficiency_ratio_Eta->SetYTitle("Ratio #varepsilon_{#eta}");
	histoEfficiency_ratio_Eta->SetXTitle("p_{T}");

	TLegend* legendratioEfficiency_Eta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendratioEfficiency_Eta->AddEntry(histoEfficiency_ratio_Eta, "#frac{Rec #varepsilon_{jet}}{Rec #varepsilon_{Inclusive}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendratioEfficiency_Eta->AddEntry(histoEfficiency_ratio_Eta_True, "#frac{Rec #varepsilon_{jet}}{Rec #varepsilon_{Inclusive}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	histoEfficiency_ratio_Eta->Draw();
	histoEfficiency_ratio_Eta_True->Draw("same");
	legendratioEfficiency_Eta->Draw("same");
	f1->Draw("same");
	canvas_ratio_Efficiency_Eta->Write();

	histoEfficiency_ratio_Pi0EtaBinning->Write();
	histoEfficiency_ratio_Pi0EtaBinning_True->Write();
	histoEfficiency_ratio_Eta->Write();
	histoEfficiency_ratio_Eta_True->Write();

	canvas_Acceptance->cd();
	TLegend* legendAcceptance   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.14*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendAcceptance->AddEntry(histoAcceptance_Pi0EtaBinning, "A_{jet, #pi^{0}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendAcceptance->AddEntry(histoAcceptance_Eta, "A_{jet, #eta}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendAcceptance->AddEntry(histoAcceptance_Pi0EtaBinning_Inclusive, "A_{inclusive, #pi^{0}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendAcceptance->AddEntry(histoAcceptance_Eta_Inclusive, "A_{inclusive, #eta}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	DrawGammaSetMarker(histoAcceptance_Pi0EtaBinning, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoAcceptance_Eta, 20, 1., kRed, kRed);
	DrawGammaSetMarker(histoAcceptance_Pi0EtaBinning_Inclusive, 20, 1., kGreen, kGreen);
	DrawGammaSetMarker(histoAcceptance_Eta_Inclusive, 20, 1., kBlue, kBlue);
	histoAcceptance_Pi0EtaBinning->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	histoAcceptance_Pi0EtaBinning->GetYaxis()->SetTitle("Acceptance");
	histoAcceptance_Pi0EtaBinning->Draw();
	histoAcceptance_Eta->Draw("same");
	histoAcceptance_Pi0EtaBinning_Inclusive->Draw("same");
	histoAcceptance_Eta_Inclusive->Draw("same");
	legendAcceptance->Draw("same");
	canvas_Acceptance->Write();
	
	canvas_Efficiency->cd();
	TLegend* legendEfficiency   = GetAndSetLegend2(0.12, 0.8, 0.45, 0.8+(0.14*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendEfficiency->AddEntry(histoEfficiency_Pi0, "#varepsilon_{jet, #pi^{0}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendEfficiency->AddEntry(histoEfficiency_Eta, "#varepsilon_{jet, #eta}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendEfficiency->AddEntry(histoEfficiency_Pi0_Inclusive, "#varepsilon_{inclusive, #pi^{0}}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendEfficiency->AddEntry(histoEfficiency_Eta_Inclusive, "#varepsilon_{inclusive, #eta}", "p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	DrawGammaSetMarker(histoEfficiency_Pi0, 20, 1., kBlack, kBlack);
	DrawGammaSetMarker(histoEfficiency_Eta, 20, 1., kRed, kRed);
	DrawGammaSetMarker(histoEfficiency_Pi0_Inclusive, 20, 1., kGreen, kGreen);
	DrawGammaSetMarker(histoEfficiency_Eta_Inclusive, 20, 1., kBlue, kBlue);
	histoEfficiency_Pi0->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	histoEfficiency_Pi0->GetYaxis()->SetTitle("#it{#varepsilon}");
	histoEfficiency_Pi0->SetAxisRange(0.,0.6,"Y");
	histoEfficiency_Pi0->Draw();
	histoEfficiency_Eta->Draw("same");
	histoEfficiency_Pi0_Inclusive->Draw("same");
	histoEfficiency_Eta_Inclusive->Draw("same");
	legendEfficiency->Draw("same");
	canvas_Efficiency->Write();	
	
	//////////////////// test /////////////////////////////////
	ctest->cd();
	//htest->Draw();
	ctest->Write();
	
    //	Filipad *filip = new Filipad("c1",0.21, 0.13, 0.05, 0.1, 1./3.5);
//	TCanvas *c1 = new TCanvas("c1", "c1", 5);
//	c1->SetLogy(1);
//	filip->SetCanvas();
//	filip->Draw();
//	TPad* padTop=filip->GetPad(1);
//	//filip->GetPad(1)->Write();
//	TH2D *hCorr = new TH2D("hCorr","", 100, 0, 20, 100, 1e-8, 1e-4);
//	hCorr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
//	hCorr->GetYaxis()->SetTitle("RAW Yield/N_{Evt}");
//	hCorr->Draw();
//	//histoYieldMeson_Pi0_Data_INT7->Draw("same");
//	histoComparePi0YieldData_INT7->Draw("same");
//	histoComparePi0YieldData_EG1->Draw("same");
//	histoComparePi0YieldData_EG2->Draw("same");
//	//filip->GetPad(1)->Write();
//	//histoComparePi0YieldData_INT7->Draw("same");
//	TPad* padBottom=filip->GetPad(2);
//	href5->Draw();
//	histoCompareEG1INT7Pi0YieldData->Draw("same");
//	histoCompareEG2INT7Pi0YieldData->Draw("same");
//	histoCompareEG1EG2Pi0YieldData->Draw("same");
//	//filip->GetPad(2)->Write();
//	c1->Update();
//	c1->SaveAs("test2.eps");

	/////// Save the Canvas file /////////////////////////////////
	canvas_ratio_Pi0Eta->SaveAs("RatioPi0Eta.eps");
	canvas_ratio_Pi0Eta_trigger->SaveAs("RatioPi0EtaTrigger.eps");
	canvas_ratio_Acceptance_Pi0->SaveAs("RatioAcceptancePi0.eps");
	canvas_ratio_Acceptance_Eta->SaveAs("RatioAcceptanceEta.eps");
	canvas_ratio_Efficiency_Pi0EtaBinning->SaveAs("RatioEfficiencyPi0.eps");
	canvas_ratio_Efficiency_Eta->SaveAs("RatioEfficiencyEta.eps");
	canvas_ratio_DiffDefault_Pi0->SaveAs("RatioDiffDefault.eps");
	canvas_Acceptance->SaveAs("Acceptance.eps");
	canvas_Efficiency->SaveAs("Efficiency.eps");
	canvas_trigger_Pi0Yield->SaveAs("CompareYieldOfTriggerForPi0.eps");
	canvas_trigger_Pi0Yield_Inclusive->SaveAs("CompareYieldOfTriggerForPi0_Inclusive.eps");
	canvasRawYield_Pi0->SaveAs("CompareYieldOfTriggerForPi0WithRatio.eps");
	canvasRawYield_Pi0_Inclusive->SaveAs("CompareYieldOfTriggerForPi0WithRatioInclusive.eps");
	canvasRawYield_Eta->SaveAs("CompareYieldOfTriggerForEtaWithRatio.eps");
	canvas_trigger_EtaYield->SaveAs("CompareYieldOfTriggerForEta.eps");
	canvas_triggerClusterData->SaveAs("CompareTriggerClusterData.eps");
	canvas_triggerClusterData_Inclusive->SaveAs("CompareTriggerClusterDataInclusive.eps");
	canvas_triggerClusterMC->SaveAs("CompareTriggerClusterMC.eps");
	canvas_triggerPi0YieldDataRatio->SaveAs("CompareTrigger_Pi0YieldDataRatio.eps");
	canvas_triggerPi0YieldDataRatio_Inclusive->SaveAs("CompareTrigger_Pi0YieldDataRatioInclusive.eps");
	canvas_triggerPi0YieldMCRatio->SaveAs("CompareTrigger_Pi0YieldMCRatio.eps");
	canvas_trigger_Eta->SaveAs("CompareTrigger_Eta.eps");
	ctest->SaveAs("test.eps");
	canvas_CorrectedYield_Pi0->SaveAs("CorrectedYield_Pi0.eps");
	canvas_CorrectedYield_Eta->SaveAs("CorrectedYield_Eta.eps");
	MyAnalysis->Close();	

	TFile *MyRejectionFactor = new TFile("MyRejectionFactor.root","recreate");
	fYieldEG1INT7Pi0Data->Write("EG1INT7RejectionFactorOfYield");		
	fYieldEG2INT7Pi0Data->Write("EG2INT7RejectionFactorOfYield");		
	fYieldEG1EG2Pi0Data->Write("EG1EG2RejectionFactorOfYield");		
	fYieldEG1INT7EtaData->Write("EG1INT7RejectionFactorOfYieldForEta");		
	fYieldEG2INT7EtaData->Write("EG2INT7RejectionFactorOfYieldForEta");		
	fYieldEG1EG2EtaData->Write("EG1EG2RejectionFactorOfYieldForEta");		

	MyRejectionFactor->Close();
}
