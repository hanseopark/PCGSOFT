
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"



void CompareInvMassBinTrigger(){
	Int_t textSizeLabelsPixel = 900*0.04;
	Int_t expectedLinesInLegend = 1;

	TCanvas* ctest = new TCanvas("ctest","",200,10,1800,1200);  // gives the page size
	DrawGammaCanvasSettings(ctest, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	ctest->Divide(3,2);
	TString EventCuts = "00010113";
	TString EventCuts_EG1 = "0008d113";
	TString EventCuts_EG2 = "0008e113";
	TString CaloCuts  = "411790607l032230000";
	TString MesonCuts = "2l631031000000d0";
	TString CaloCuts_MB = "411790607l032230000";
	//TString CaloCuts_Eta_MB = "411791106f032230000";
	TString MesonCuts_MB = "0l631031000000d0";
	TString Energy = "13TeV";

	const Int_t NTRIG = 3;
	const Int_t color[NTRIG] = {4, 2, 1};
	const char *Trig[NTRIG] = {"INT7", "EG1", "EG2"};

	TString CutSelection[NTRIG];
	CutSelection[0] = Form("%s_%s_%s",EventCuts.Data(),CaloCuts.Data(),MesonCuts.Data());
	CutSelection[1] = Form("%s_%s_%s",EventCuts_EG1.Data(),CaloCuts.Data(),MesonCuts.Data());
	CutSelection[2] = Form("%s_%s_%s",EventCuts_EG2.Data(),CaloCuts.Data(),MesonCuts.Data());
	
	TFile *Pi0_Data[NTRIG];
	for(Int_t itrig=0; itrig< NTRIG; itrig++){
	Pi0_Data[itrig] 	= new TFile(Form("../%s/%s/Pi0_data_GammaConvV1WithoutCorrection_%s.root",CutSelection[itrig].Data(),Energy.Data(),CutSelection[itrig].Data()),"read");
	}

	Int_t Bin = 3;
	TH1D *histoPi0InvMassBin[NTRIG][Bin];
   	TH2D* href = new TH2D("href", "",100,0,0.3, 100000, 0, 1000);
 	href->SetXTitle("M_{#gamma#gamma} (GeV/c)");
	href->SetStats(0);
   	TH2D* href2 = new TH2D("href2", "",100,0,0.3, 100, 0, 2);
 	href2->SetXTitle("M_{#gamma#gamma} (GeV/c)");
 	href2->SetYTitle("Ratio");
	href2->SetStats(0);
	for(Int_t itrig=0; itrig < NTRIG; itrig++){	
	histoPi0InvMassBin[itrig][0] = (TH1D*) Pi0_Data[itrig] -> Get("fHistoMappingSignalInvMass_in_Pt_Bin12");
	histoPi0InvMassBin[itrig][1] = (TH1D*) Pi0_Data[itrig] -> Get("fHistoMappingSignalInvMass_in_Pt_Bin13");
	histoPi0InvMassBin[itrig][2] = (TH1D*) Pi0_Data[itrig] -> Get("fHistoMappingSignalInvMass_in_Pt_Bin14");
	}	
	TLegend* leg[Bin];
	TLegend* leg2[Bin];
	TLegend* legtrig   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);

	for(Int_t bin =0; bin<Bin; bin++){
	leg[bin]   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	leg2[bin]   = GetAndSetLegend2(0.12,0.74, 0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
		for(Int_t itrig =0; itrig<NTRIG; itrig++){
		histoPi0InvMassBin[itrig][bin]->SetLineWidth(1);
		histoPi0InvMassBin[itrig][bin]->SetLineColor(color[itrig]);
		histoPi0InvMassBin[itrig][bin]->SetMarkerColor(color[itrig]);
		histoPi0InvMassBin[itrig][bin]->SetMarkerStyle(24);
		}
	}
	TRatioPlot* rp;
	Double_t xmin[3]; 
	Double_t xmax[3];
	for(Int_t bin= 0; bin <Bin; bin++){
	ctest->cd(bin+1);
	href->Clear();
	xmin[bin] = 5.6 + 0.8*(bin); 
	xmax[bin] = 6.4 + 0.8*(bin); 
	//href->SetTitle(Form("%f GeV/#it{c}< #it{p}_T < %f GeV/#it{c}",xmin[bin],xmax[bin]));
 	href->Draw();	
		for(Int_t itrig =0; itrig <NTRIG; itrig++){
		histoPi0InvMassBin[itrig][bin]->Draw("same p");
		leg[bin]->AddEntry(histoPi0InvMassBin[itrig][bin], Trig[itrig]);
		leg[bin]->Draw("same p");
		}
//	rp = new TRatioPlot(histoPi0InvMassBin[2][bin],histoPi0InvMassBin[1][bin]);
//	rp->Draw("same");
//	rp->Clear();
	}	
	TH1D *hratio[Bin];
	
	for(Int_t bin= 0; bin<Bin; bin++){
	ctest->cd(4+bin);
	href2->Draw("p");
		hratio[bin]= (TH1D*) histoPi0InvMassBin[2][bin]->Clone();
		hratio[bin]->Divide(hratio[bin] ,histoPi0InvMassBin[0][bin],1. ,1. ,"");
		leg2[bin]->AddEntry(hratio[bin], "EG2/INT7");
		hratio[bin]->Draw("SAME P");
		leg2[bin]->Draw("SAME P");
	}

	ctest->SaveAs("test.eps");
}



