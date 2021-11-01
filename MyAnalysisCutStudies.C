//////////////////////////////////////////////////////////////
///////// My Analysis ///////////////////////////////////////
/////////////////////////////////////////////////////////////

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/ExtractSignalBinning.h"
//#include "MyCommonHeaders/Filipad.h"

void MyAnalysisCutStudies(TString nameMeson = "", TString nameTrig = ""){
    const Int_t MaxNumberOfFiles    = 13;

	// input var //
	Int_t gMeson = 0;
	Int_t gTrig = 0;	
	if (nameMeson.CompareTo("Pi0") == 0) gMeson = 0;
	if (nameMeson.CompareTo("Eta") == 0) gMeson = 1;
	if (nameMeson.CompareTo("Pi0EtaBinning") == 0) gMeson = 2;
	if (nameTrig.CompareTo("INT7") == 0) gTrig = 0;
	if (nameTrig.CompareTo("EG1") == 0) gTrig = 1;
	if (nameTrig.CompareTo("EG2") == 0) gTrig = 2;

    Int_t NBins = 200;
    Double_t* BinEdges = new Double_t[300];
    Int_t MaxBin = GetBinning( BinEdges, NBins, nameMeson.Data(), "13TeV", 4, 0, kFALSE, "", kTRUE ); // speicaltirgg set is needed
    Int_t StartBin = GetStartBin(nameMeson.Data(), "13TeV", 4, 0, "", "", kTRUE); // specialtrigg set is needed
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
    TString Energy = "13TeV";			

	
	const Int_t NMESON = 3;
	const Int_t NTRIG = 3;
	const Int_t NCUTS = 8;
	const Int_t NCUTSCALO = 6;
	const Int_t NCUTSMESON = 3;
	const Int_t NVAR = 4;
	//const char *Meson[NMESON] = {"Pi0", "Eta", "Pi0EtaBinning"};
	const char *Trig[NTRIG] = {"INT7", "EG1", "EG2"};
	const char *Cuts[NCUTS] = {"Standard", "ClusterMinEnergy", "ClusterNCells", "ClusterTiming", "ClusterM02", "ClusterTrackMatching", "OpeningAngle", "Alpha"};
	const char *CutsCalo[NCUTSCALO] = {"Standard", "ClusterMinEnergy", "ClusterNCells", "ClusterTiming", "ClusterM02", "ClusterTrackMatching"};
	const char *CutsMeson[NCUTSMESON] = {"Standard", "OpeningAngle", "Alpha"};
	const char *Var[NVAR] = {"Standard", "Cutvariation1", "Cutvariation2", "CutVariation3"};
	Int_t ColorVAR[NVAR];
	for (Int_t icolor = 1; icolor < NVAR+1; icolor++){
	ColorVAR[icolor] = icolor;
	}
	
	TString Meson[NMESON];
	Meson[0] = "Pi0";
	Meson[1] = "Eta";
	Meson[2] = "Pi0EtaBinning";

	TString EventCuts[NTRIG];
	EventCuts[0] = "00010113";
	EventCuts[1] = "0008d113";
	EventCuts[2] = "0008e113";
	
	TString CaloCuts[NCUTSCALO][NVAR];
	// Standard 
	CaloCuts[0][0] = "411790607l032230000";

	// ClusterMinEnergy
	CaloCuts[1][0] = "411790607l032230000";	
	CaloCuts[1][1] = "411790607l022230000";
	CaloCuts[1][2] = "411790607l042230000";
	CaloCuts[1][3] = "411790607l052230000";

	// NCell
	CaloCuts[2][0] = "411790607l032230000";	
	CaloCuts[2][1] = "411790607l031230000"; 
	CaloCuts[2][2] = "411790607l033230000";
	CaloCuts[2][3] = "411790607l032230000"; 

	// ClusterTiming
	CaloCuts[3][0] = "411790607l032230000";	
	CaloCuts[3][1] = "411790605l032230000";
	CaloCuts[3][2] = "411790608l032230000";
	CaloCuts[3][3] = "411790609l032230000";
		
	// ClusterM02
	CaloCuts[4][0] = "411790607l032230000";	
	CaloCuts[4][1] = "411790607l032220000";
	CaloCuts[4][2] = "411790607l032250000"; 
	CaloCuts[4][3] = "411790607l032230000";

	// ClusterTrackMatching
	CaloCuts[5][0] = "411790607l032230000";	
	CaloCuts[5][1] = "411790607e032230000";
	CaloCuts[5][2] = "411790607g032230000";
	CaloCuts[5][3] = "411790607h032230000";
	//CaluCuts[5][4] = "4117906077032230000";

	TString MesonCuts[NCUTSMESON][NVAR];
	// Standard
	MesonCuts[0][0] = "2l631031000000d0";

	// OpeningAngle	
	MesonCuts[1][0] = "2l631031000000d0";
	MesonCuts[1][1] = "2l631031000000b0";
	MesonCuts[1][2] = "2l631031000000g0";
	MesonCuts[1][3] = "2l631031000000a0";

	// Alpha
	MesonCuts[2][0] = "2l631031000000d0";
	MesonCuts[2][1] = "2l631041000000d0"; 
	MesonCuts[2][2] = "2l631051000000d0"; 
	MesonCuts[2][3] = "2l631061000000d0"; 

	// summary cut variation
	TString CaloMesonCuts[NCUTS][NVAR];
		for (Int_t icutscalo = 0; icutscalo < NCUTSCALO; icutscalo++){
			for (Int_t ivar =0; ivar < NVAR; ivar++){
				CaloMesonCuts[icutscalo][ivar] = Form("%s_%s", CaloCuts[icutscalo][ivar].Data(), MesonCuts[0][0].Data()); 
			}
		}
		for (Int_t icutsmeson = 0; icutsmeson < NCUTSMESON-1; icutsmeson++){
			for (Int_t ivar = 0; ivar < NVAR; ivar++)
				CaloMesonCuts[icutsmeson+NCUTSCALO][ivar] = Form("%s_%s", CaloCuts[0][0].Data(), MesonCuts[icutsmeson+1][ivar].Data()); 
		}
	TString EventCaloMesonCuts[NTRIG][NCUTS][NVAR];
		for (Int_t itrig = 0; itrig < NTRIG; itrig++){
			for (Int_t icuts = 0; icuts < NCUTS; icuts++){
				for (Int_t ivar = 0; ivar < NVAR; ivar++){
				EventCaloMesonCuts[itrig][icuts][ivar] = Form("%s_%s", EventCuts[itrig].Data(), CaloMesonCuts[icuts][ivar].Data());
				}
			}
		}

	cout << "Cut: " << EventCaloMesonCuts[0][0][0].Data() << endl;
		
	// read File and then Corredted yield
	TFile *myfile[NMESON][NTRIG][NCUTS][NVAR];
	TH1D* histoCorrectedYieldNormalEff[NMESON][NTRIG][NCUTS][NVAR];
	for (Int_t imeson = 0; imeson < NMESON; imeson++){
		for (Int_t itrig = 0; itrig < NTRIG; itrig++){
		//for (Int_t itrig =0; itrig < 2; itrig++){
			myfile[imeson][itrig][0][0] = new TFile(Form("../%s/%s/%s_data_GammaConvV1Correction_%s.root", EventCaloMesonCuts[itrig][0][0].Data(), Energy.Data(), Meson[imeson].Data(), EventCaloMesonCuts[itrig][0][0].Data()), "read");
			histoCorrectedYieldNormalEff[imeson][itrig][0][0] = (TH1D*) myfile[imeson][itrig][0][0]->Get("CorrectedYieldNormEff");
			for (Int_t icuts = 1; icuts < NCUTS; icuts++){
				for (Int_t ivar = 1; ivar < NVAR; ivar++){
				myfile[imeson][itrig][icuts][ivar] = new TFile(Form("../%s/%s/%s_data_GammaConvV1Correction_%s.root", EventCaloMesonCuts[itrig][icuts][ivar].Data(), Energy.Data(), Meson[imeson].Data(), EventCaloMesonCuts[itrig][icuts][ivar].Data()), "read");
				histoCorrectedYieldNormalEff[imeson][itrig][icuts][ivar] = (TH1D*) myfile[imeson][itrig][icuts][ivar]->Get("CorrectedYieldNormEff");
				}
			}
		}
	}
	TH1D* histoCorrectedYieldNormalEffRatio[NMESON][NTRIG][NCUTS][NVAR];
	for (Int_t imeson = 0; imeson < NMESON; imeson++){
		for (Int_t itrig = 0; itrig < NTRIG; itrig++){
			for (Int_t icuts = 1; icuts < NCUTS; icuts++){
				for (Int_t ivar = 1; ivar < NVAR; ivar++){
				histoCorrectedYieldNormalEffRatio[imeson][itrig][icuts][ivar] = (TH1D*) histoCorrectedYieldNormalEff[imeson][itrig][icuts][ivar]->Clone();	
				histoCorrectedYieldNormalEffRatio[imeson][itrig][icuts][ivar]->Divide(histoCorrectedYieldNormalEffRatio[imeson][itrig][icuts][ivar], histoCorrectedYieldNormalEff[imeson][itrig][0][0],1.,1.,"");	
				}
			}
		}
	}
	TH2F* histo2DDummyPt[NMESON][NTRIG][NCUTS];
	TH2F* histo2DDummyPtRatio[NMESON][NMESON][NCUTS];
	TH2F* histo2DDummyPtRatioOnly[NMESON][NMESON][NCUTS];

//	for(Int_t imeson =0; imeson<NMESON; imeson++){
//		for(Int_t gTrig =0; gTrig<NTRIG; gTrig++){
			histo2DDummyPt[gMeson][gTrig][0]               = new TH2F("histo2DDummyPt","histo2DDummyPt",1000,0, histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->GetNbinsX()),
					10000, 0.01*FindSmallestBin1DHist(histoCorrectedYieldNormalEff[gMeson][gTrig][0][0],1e6 ), 3*histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->GetMaximum());
			SetStyleHistoTH2ForGraphs(histo2DDummyPt[gMeson][gTrig][0], "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.033,0.04, 0.033,0.04, 1,1.35);

			histo2DDummyPtRatio[gMeson][gTrig][0] = new TH2F("histo2DDummyPtRatio", "", 1000, 0, histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->GetNbinsX()), 10000, 0.5, 1.5);	
			SetStyleHistoTH2ForGraphs(histo2DDummyPtRatio[gMeson][gTrig][0], "#it{p}_{T} (GeV/#it{c})", "#frac{var}{standard}",0.07,0.04,0.07,0.1,0.8,0.55,510,505 );

			for (Int_t icuts= 1; icuts<NCUTS; icuts++){
				histo2DDummyPt[gMeson][gTrig][icuts]               = new TH2F("histo2DDummyPt","histo2DDummyPt",1000,0, histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetNbinsX()),
						10000, 0.01*FindSmallestBin1DHist(histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1],1e6 ), 3*histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetMaximum());
				SetStyleHistoTH2ForGraphs(histo2DDummyPt[gMeson][gTrig][icuts], "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.033,0.04, 0.033,0.04, 1,1.35);
				histo2DDummyPtRatio[gMeson][gTrig][icuts] = new TH2F("histo2DDummyPtRatio", "", 1000, 0, histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetNbinsX()), 10000, 0.5, 1.5);	
				SetStyleHistoTH2ForGraphs(histo2DDummyPtRatio[gMeson][gTrig][icuts], "#it{p}_{T} (GeV/#it{c})", "#frac{var}{standard}",0.07,0.04,0.07,0.1,0.8,0.55,510,505 );
				histo2DDummyPtRatioOnly[gMeson][gTrig][icuts] = new TH2F("histo2DDummyPtRatioOnly", "", 1000, 0, histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][1]->GetNbinsX()), 10000, 0.5, 1.5);	
				SetStyleHistoTH2ForGraphs(histo2DDummyPtRatioOnly[gMeson][gTrig][icuts], "#it{p}_{T} (GeV/#it{c})", "#frac{var}{standard}",0.035,0.04,0.035,0.04,0.9,0.8,510,505 );

			}
//		}
//	}
	TCanvas *ctest = new TCanvas("canvas","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(ctest, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	TPad *pad[NCUTS];	
	TPad *padRatio[NCUTS];
	pad[0] = new TPad("pad0", "", 0.,0.7,0.25,1.,-1,-1,-2);
	pad[1] = new TPad("pad1", "", 0.25,0.7,0.5,1.,-1,-1,-2);
	pad[2] = new TPad("pad2", "", 0.5,0.7,0.75,1.,-1,-1,-2);
	pad[3] = new TPad("pad3", "", 0.75,0.7,1.,1.,-1,-1,-2);
	pad[4] = new TPad("pad4", "", 0.,0.2,0.25,0.5,-1,-1,-2);
	pad[5] = new TPad("pad5", "", 0.25,0.2,0.5,0.5,-1,-1,-2);
	pad[6] = new TPad("pad6", "", 0.5,0.2,0.75,0.5,-1,-1,-2);
	pad[7] = new TPad("pad7", "", 0.75,0.2,1.,0.5,-1,-1,-2);

	padRatio[0] = new TPad("padRatio0", "", 0.,0.5,0.25,0.7,-1,-1,-2);
	padRatio[1] = new TPad("padRatio1", "", 0.25,0.5,0.5,0.7,-1,-1,-2);
	padRatio[2] = new TPad("padRatio2", "", 0.5,0.5,0.75,0.7,-1,-1,-2);
	padRatio[3] = new TPad("padRatio3", "", 0.75,0.5,1.,0.7,-1,-1,-2);
	padRatio[4] = new TPad("padRatio4", "", 0.,0.,0.25,0.2,-1,-1,-2);
	padRatio[5] = new TPad("padRatio5", "", 0.25,0.,0.5,0.2,-1,-1,-2);
	padRatio[6] = new TPad("padRatio6", "", 0.5,0.,0.75,0.2,-1,-1,-2);
	padRatio[7] = new TPad("padRatio7", "", 0.75,0.,1.,0.2,-1,-1,-2);

	for (Int_t icuts = 0; icuts<NCUTS; icuts++){
		DrawGammaPadSettings( pad[icuts], 0.14, 0.02, 0.10, 0.);
		DrawGammaPadSettings( padRatio[icuts], 0.14, 0.02, 0., 0.10);
		pad[icuts]->Draw();
		padRatio[icuts]->Draw();
	}

	Int_t textSizeLabelsPixel = 900*0.04;
	Int_t expectedLinesInLegend = 1;
	TLegend* legcuts = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	ctest->cd();
	pad[0]->cd();
	pad[0]->SetLogy();
	histo2DDummyPt[gMeson][gTrig][0]->SetTitle(Cuts[0]);
	histo2DDummyPt[gMeson][gTrig][0]->Draw();
	histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->SetMarkerColor(1);
	histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->SetMarkerSize(1);
	histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->SetMarkerStyle(22);
	histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->SetLineColor(1);
	histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->SetLineWidth(1);
	histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->Draw("e1,same");

	padRatio[0]->cd();
	histo2DDummyPtRatio[gMeson][gTrig][0]->Draw();

	for (Int_t icuts = 1; icuts<NCUTS; icuts++){
		if(icuts != 2 && icuts != 4){

			pad[icuts]->cd();
			pad[icuts]->SetLogy();
			histo2DDummyPt[gMeson][gTrig][icuts]->SetTitle(Cuts[icuts]);
			histo2DDummyPt[gMeson][gTrig][icuts]->Draw();
			histoCorrectedYieldNormalEff[gMeson][gTrig][0][0]->Draw("e1,same");
			for(Int_t ivar= 1; ivar< NVAR; ivar++){
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetMarkerColor(ivar+1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetMarkerSize(1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetMarkerStyle(22+ivar);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetLineColor(ivar+1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetLineWidth(1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->Draw("e1,same");
			}
			padRatio[icuts]->cd();
			histo2DDummyPtRatio[gMeson][gTrig][icuts]->Draw();
			for(Int_t ivar = 1; ivar < NVAR; ivar++){
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetMarkerColor(ivar+1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetMarkerSize(1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetMarkerStyle(22+ivar);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetLineColor(ivar+1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetLineWidth(1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->Draw("e1,same");
			}
		} else {
			pad[icuts]->cd();
			pad[icuts]->SetLogy();
			histo2DDummyPt[gMeson][gTrig][icuts]->SetTitle(Cuts[icuts]);
			histo2DDummyPt[gMeson][gTrig][icuts]->Draw();
			for(Int_t ivar= 1; ivar< NVAR-1; ivar++){
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetMarkerColor(ivar+1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetMarkerSize(1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetMarkerStyle(22+ivar);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetLineColor(ivar+1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->SetLineWidth(1);
				histoCorrectedYieldNormalEff[gMeson][gTrig][icuts][ivar]->Draw("e1,same");
			}
			padRatio[icuts]->cd();
			histo2DDummyPtRatio[gMeson][gTrig][icuts]->Draw();
			for(Int_t ivar = 1; ivar < NVAR-1; ivar++){
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetMarkerColor(ivar+1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetMarkerSize(1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetMarkerStyle(22+ivar);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetLineColor(ivar+1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->SetLineWidth(1);
				histoCorrectedYieldNormalEffRatio[gMeson][gTrig][icuts][ivar]->Draw("e1,same");
			}

		}
	}	

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

	// To measure combined bin

	// take te value to upper a for only positive value
	Int_t index=0;
	Int_t indexForToCheck=0;
	Int_t NumberOfBins[NVAR]={0};
	Int_t TotalBin=0;
	for(Int_t ivar =1; ivar<NVAR; ivar++){
	NumberOfBins[ivar] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->GetNbinsX();
	TotalBin +=NumberOfBins[ivar];
	}
	
	// Combined bin
	
	Double_t x_Combined[200] = {0};
	Double_t xerror_Combined[200] = {0};
	Double_t y_Combined[200] = {0};
	Double_t yerror_Combined[200] = {0};
	Double_t ToCheck[200] = {0};

		// ClusterMinEnergy
	cout << "############################" << endl;
	cout << "ClusterMinEnergy" << endl;
	for(Int_t ivar=1; ivar<NVAR; ivar++){
		for(Int_t a =0; a<200; a++){
			x_Combined[a] = 0;
			xerror_Combined[a] = 0;
			y_Combined[a] = 0;
			xerror_Combined[a] = 0;
		}
		for(Int_t ibin=0; ibin<NumberOfBins[ivar]; ibin++){
			x_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->GetBinCenter(ibin+1);
			xerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->GetBinWidth(ibin+1);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->GetBinContent(ibin+1);
			if(y_Combined[index]<1.0 && y_Combined[index]>0.01) histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->SetBinContent(ibin+1, 2-y_Combined[index]);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->GetBinContent(ibin+1);
			yerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->GetBinError(ibin+1);
			index++;
		}	
	}
	


		// graph
	TGraphErrors* graphCombinedRatioCME = new TGraphErrors(index, x_Combined, y_Combined, xerror_Combined, yerror_Combined); 
//	TF1 *fCombinedRatioCMEpol1 = new TF1 ("fCombinedRatioCMEpol1", "[0]+[1]*x+[2]*x*x",2,20);
//	fCombinedRatioCMEpol1->SetParameters(1.0,0.01, 0.01);
	TF1 *fCombinedRatioCMEpol1 = new TF1 ("fCombinedRatioCMEpol1", "[0]",2,20);
	fCombinedRatioCMEpol1->SetParameter(0,1.02);
	graphCombinedRatioCME->Fit("fCombinedRatioCMEpol1", "R");
	TCanvas *cCombinedRatioCME = new TCanvas("cCombinedRatioCME","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCombinedRatioCME, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	cCombinedRatioCME->cd();
	//graphCombinedRatioCME->Draw();
	histo2DDummyPtRatioOnly[gMeson][gTrig][1]->Draw();
	for(Int_t ivar=1; ivar<NVAR; ivar++){
	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->Draw("e1,same");
	}	
	fCombinedRatioCMEpol1->Draw("same");
	
		// NCells
	cout << "############################" << endl;
	cout << "NCells" << endl;
	index = 0;
	indexForToCheck = 0;
	for(Int_t ivar=1; ivar<NVAR-1; ivar++){
		for(Int_t a = 0; a<200; a++){
			ToCheck[a]=0;
			x_Combined[a] = 0;
			xerror_Combined[a] = 0;
			y_Combined[a] = 0;
			xerror_Combined[a] = 0;
		}
		for(Int_t ibin=0; ibin<NumberOfBins[ivar]; ibin++){
			x_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->GetBinCenter(ibin+1);
			xerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->GetBinWidth(ibin+1);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->GetBinContent(ibin+1);
			if(y_Combined[index]<1.0 && y_Combined[index]>0.01) histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->SetBinContent(ibin+1, 2-y_Combined[index]);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->GetBinContent(ibin+1);
			yerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->GetBinError(ibin+1);
			index++;
		}	
	}
		// graph
			// For changing range to fit function
	for(Int_t ibin=0; ibin<NumberOfBins[1]; ibin++){
	//Int_t i = NumberOfBins[1]; // only varition 1
	Int_t i = 0; // only varition 2
	x_Combined[i] =0;
	xerror_Combined[i] = 0;
	y_Combined[i] = 0;
	xerror_Combined[i] = 0;
	i++;
	}

	TGraphErrors* graphCombinedRatioNCells = new TGraphErrors(index, x_Combined, y_Combined, xerror_Combined, yerror_Combined); 
	//TF1 *fCombinedRatioNCellspol1 = new TF1 ("fCombinedRatioNCellspol1", "[0]+[1]*x+[2]*x*x",2,20);
	TF1 *fCombinedRatioNCellspol1 = new TF1 ("fCombinedRatioNCellspol1", "[0]",2,20);
	//TF1 *fCombinedRatioNCellspol1 = new TF1 ("fCombinedRatioNCellspol1", "[0]+[1]",2,8);
	fCombinedRatioNCellspol1->SetParameters(1.0,0.01,0.01);
	//fCombinedRatioNCellspol2->SetParameters(1.0,0.01,0.01);
	graphCombinedRatioNCells->Fit("fCombinedRatioNCellspol1", "R");
	//graphCombinedRatioNCells->Fit("fCombinedRatioNCellspol2", "R");
	TCanvas *cCombinedRatioNCells = new TCanvas("cCombinedRatioNCells","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCombinedRatioNCells, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	cCombinedRatioNCells->cd();
	//graphCombinedRatioNCells->Draw();
	histo2DDummyPtRatioOnly[gMeson][gTrig][2]->Draw();
	for(Int_t ivar=1; ivar<NVAR-1; ivar++){
	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][2][ivar]->Draw("e1,same");
	}
	fCombinedRatioNCellspol1->Draw("same");
	//fCombinedRatioNCellspol2->Draw("same");
//	cout << "############################" << endl;
//	cout << " Fitting function of NCells for sys error" << endl;
//	cout << fCombinedRatioNCellspol1->GetParameter(0) << "+" << fCombinedRatioNCellspol1->GetParameter(1) << "x+"<< fCombinedRatioNCellspol1->GetParameter(2) << "x^2" <<endl; 	
//	cout << " Constant function using above fitting function over 8GeV/c" << endl;
//	cout << fCombinedRatioNCellspol1->GetParameter(0)+fCombinedRatioNCellspol1->GetParameter(1)*8+fCombinedRatioNCellspol1->GetParameter(2)*8*8 <<endl;
//	cout << "############################" << endl;
		// Timing

	cout << "############################" << endl;
	cout << "Timing" << endl;
	index = 0;
	for(Int_t ivar=1; ivar<NVAR; ivar++){
		for(Int_t a =0; a<200; a++){
			x_Combined[a] = 0;
			xerror_Combined[a] = 0;
			y_Combined[a] = 0;
			xerror_Combined[a] = 0;
		}
		for(Int_t ibin=0; ibin<NumberOfBins[ivar]; ibin++){
			x_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->GetBinCenter(ibin+1);
			xerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->GetBinWidth(ibin+1);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->GetBinContent(ibin+1);
			if(y_Combined[index]<1.0 && y_Combined[index]>0.01) histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->SetBinContent(ibin+1, 2-y_Combined[index]);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->GetBinContent(ibin+1);
			yerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->GetBinError(ibin+1);
			index++;
		}	
	}
		// graph
	TGraphErrors* graphCombinedRatioTiming = new TGraphErrors(index, x_Combined, y_Combined, xerror_Combined, yerror_Combined); 
	TF1 *fCombinedRatioTimingpol1 = new TF1 ("fCombinedRatioTimingpol1", "[0]",2,20);
	fCombinedRatioTimingpol1->SetParameters(0,1.01);
	graphCombinedRatioTiming->Fit("fCombinedRatioTimingpol1", "R");
	TCanvas *cCombinedRatioTiming = new TCanvas("cCombinedRatioTiming","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCombinedRatioTiming, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	cCombinedRatioTiming->cd();
	//graphCombinedRatioTiming->Draw();
	histo2DDummyPtRatioOnly[gMeson][gTrig][3]->Draw();
	for(Int_t ivar=1; ivar<NVAR; ivar++){
	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][3][ivar]->Draw("e1,same");
	}	
	fCombinedRatioTimingpol1->Draw("same");
	

		// Shape
	
	cout << "############################" << endl;
	cout << "Shape" <<endl;
	index = 0;
	for(Int_t ivar=1; ivar<NVAR-1; ivar++){
		for(Int_t a =0; a<200; a++){
			x_Combined[a] = 0;
			xerror_Combined[a] = 0;
			y_Combined[a] = 0;
			xerror_Combined[a] = 0;
		}
		for(Int_t ibin=0; ibin<NumberOfBins[ivar]; ibin++){
			x_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->GetBinCenter(ibin+1);
			xerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->GetBinWidth(ibin+1);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->GetBinContent(ibin+1);
			if(y_Combined[index]<1.0 && y_Combined[index]>0.01) histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->SetBinContent(ibin+1, 2-y_Combined[index]);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->GetBinContent(ibin+1);
			yerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->GetBinError(ibin+1);
			index++;
		}	
	}
	// For changing range to fit function
//	for(Int_t ibin=0; ibin<NumberOfBins[1]; ibin++){
//	Int_t i = NumberOfBins[1]; // only varition 1
//	//Int_t i = 0; // only varition 2
//	x_Combined[i] =0;
//	xerror_Combined[i] = 0;
//	y_Combined[i] = 0;
//	xerror_Combined[i] = 0;
//	}	

		// graph
	TGraphErrors* graphCombinedRatioShape = new TGraphErrors(index, x_Combined, y_Combined, xerror_Combined, yerror_Combined); 
	TF1 *fCombinedRatioShapepol1 = new TF1 ("fCombinedRatioShapepol1", "[0]",2,20);
	fCombinedRatioShapepol1->SetParameter(0, 1.02);
	graphCombinedRatioShape->Fit("fCombinedRatioShapepol1", "R");
	TCanvas *cCombinedRatioShape = new TCanvas("cCombinedRatioShape","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCombinedRatioShape, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	cCombinedRatioShape->cd();
	//graphCombinedRatioShape->Draw();
	histo2DDummyPtRatioOnly[gMeson][gTrig][4]->Draw();
	for(Int_t ivar=1; ivar<NVAR-1; ivar++){
	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][4][ivar]->Draw("e1,same");
	}	
	fCombinedRatioShapepol1->Draw("same");
	
		// TM

	cout << "############################" << endl;
	cout << "TM" << endl;
	index = 0;
	for(Int_t ivar=1; ivar<NVAR; ivar++){
		for(Int_t a = 0; a<200; a++){
			x_Combined[a] = 0;
			xerror_Combined[a] = 0;
			y_Combined[a] = 0;
			xerror_Combined[a] = 0;
		}
		for(Int_t ibin=0; ibin<NumberOfBins[ivar]; ibin++){
			x_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->GetBinCenter(ibin+1);
			xerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->GetBinWidth(ibin+1);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->GetBinContent(ibin+1);
			if(y_Combined[index]<1.0 && y_Combined[index]>0.01) histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->SetBinContent(ibin+1, 2-y_Combined[index]);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->GetBinContent(ibin+1);
			yerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->GetBinError(ibin+1);
			index++;
		}	
	}
		// graph
	TGraphErrors* graphCombinedRatioTM = new TGraphErrors(index, x_Combined, y_Combined, xerror_Combined, yerror_Combined); 
	TF1 *fCombinedRatioTMpol1 = new TF1 ("fCombinedRatioTMpol1", "[0]",2,20);
	fCombinedRatioTMpol1->SetParameters(0,1.01);
	graphCombinedRatioTM->Fit("fCombinedRatioTMpol1", "R");
	TCanvas *cCombinedRatioTM = new TCanvas("cCombinedRatioTM","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCombinedRatioTM, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	cCombinedRatioTM->cd();
	//graphCombinedRatioTM->Draw();
	histo2DDummyPtRatioOnly[gMeson][gTrig][5]->Draw();
	for(Int_t ivar=1; ivar<NVAR; ivar++){
	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][5][ivar]->Draw("e1,same");
	}	
	fCombinedRatioTMpol1->Draw("same");
	
		// OpeningAngle 
	
	cout << "############################" << endl;
	cout << "OpeningAngle" << endl;
	index = 0;
	for(Int_t ivar=1; ivar<NVAR; ivar++){
		for(Int_t a =0; a<200; a++){
			x_Combined[a] = 0;
			xerror_Combined[a] = 0;
			y_Combined[a] = 0;
			xerror_Combined[a] = 0;
		}
		for(Int_t ibin=0; ibin<NumberOfBins[ivar]; ibin++){
			x_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->GetBinCenter(ibin+1);
			xerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->GetBinWidth(ibin+1);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->GetBinContent(ibin+1);
			if(y_Combined[index]<1.0 && y_Combined[index]>0.01) histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->SetBinContent(ibin+1, 2-y_Combined[index]);
			y_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->GetBinContent(ibin+1);
			yerror_Combined[index] = histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->GetBinError(ibin+1);
			index++;
		}	
	}
		// graph
	TGraphErrors* graphCombinedRatioOpeningAngle = new TGraphErrors(index, x_Combined, y_Combined, xerror_Combined, yerror_Combined); 
	TF1 *fCombinedRatioOpeningAnglepol1 = new TF1 ("fCombinedRatioOpeningAnglepol1", "[0]",2,20);
	fCombinedRatioOpeningAnglepol1->SetParameters(0,1.01);
	graphCombinedRatioOpeningAngle->Fit("fCombinedRatioOpeningAnglepol1", "R");
	TCanvas *cCombinedRatioOpeningAngle = new TCanvas("cCombinedRatioOpeningAngle","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCombinedRatioOpeningAngle, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	cCombinedRatioOpeningAngle->cd();
	//graphCombinedRatioOpeningAngle->Draw();
	histo2DDummyPtRatioOnly[gMeson][gTrig][6]->Draw();
	for(Int_t ivar=1; ivar<NVAR; ivar++){
	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][6][ivar]->Draw("e1,same");
	}	
	fCombinedRatioOpeningAnglepol1->Draw("same");
	
//	ctest2->cd();
//	histo2DDummyPtRatioOnly[gMeson][gTrig][1]->Draw();
//	for(Int_t ivar=1; ivar<NVAR; ivar++){
//	histoCorrectedYieldNormalEffRatio[gMeson][gTrig][1][ivar]->Draw("e1,same");
//	ftest->Draw("same");
//	}

	TFile *MyAnalysisCutStudies = new TFile("MyAnalysisCutStudies.root", "Update");
	ctest->Write(Form("%s_%s_CutStudiesWithRatio",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioCME->Write(Form("%s_%s_CombinedRatioCME",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioNCells->Write(Form("%s_%s_CombinedRatioNCells",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioTiming->Write(Form("%s_%s_CombinedRatioTiming",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioShape->Write(Form("%s_%s_CombinedRatioShape",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioTM->Write(Form("%s_%s_CombinedRatioTM",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioOpeningAngle->Write(Form("%s_%s_CombinedRatioOpeningAngle",Meson[gMeson].Data(),Trig[gTrig]));
	MyAnalysisCutStudies->Write();
	MyAnalysisCutStudies->Close();

	ctest->SaveAs(Form("%s_%s_CutStudiesWithRatio.eps",Meson[gMeson].Data(), Trig[gTrig]));		
	cCombinedRatioCME->SaveAs(Form("%s_%s_CombinedRatioCME.eps",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioNCells->SaveAs(Form("%s_%s_CombinedRatioNCells.eps",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioTiming->SaveAs(Form("%s_%s_CombinedRatioTiming.eps",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioShape->SaveAs(Form("%s_%s_CombinedRatioShape.eps",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioTM->SaveAs(Form("%s_%s_CombinedRatioTM.eps",Meson[gMeson].Data(),Trig[gTrig]));
	cCombinedRatioOpeningAngle->SaveAs(Form("%s_%s_CombinedRatioOpeningAngle.eps",Meson[gMeson].Data(),Trig[gTrig]));
}

