///////////////////////////////////////////////////////
/////// Unfolding macro //////////////////////
/////////////////////////////////////////////////////////

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"



void MyUnfold(){

	TLatex *t = new TLatex();
	t->SetTextAlign(12);
	t->SetTextSize(0.05);

	Int_t textSizeLabelsPixel = 900*0.04;
	Int_t expectedLinesInLegend = 1;

	TCanvas* canvas_compare_Unfolding_Pi0 = new TCanvas("canvas_compare_Unfolding_Pi0","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_compare_Unfolding_Pi0, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	TCanvas* canvas_compare_Unfolding_Pi0_DataPoints = new TCanvas("canvas_compare_Unfolding_Pi0_DataPoints","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_compare_Unfolding_Pi0_DataPoints, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	TCanvas* canvas_compare_Unfolding_Eta = new TCanvas("canvas_compare_Unfolding_Eta","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_compare_Unfolding_Eta, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)
	TCanvas* canvas_compare_Unfolding_Eta_DataPoints = new TCanvas("canvas_compare_Unfolding_Eta_DataPoints","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings(canvas_compare_Unfolding_Eta_DataPoints, 0.08, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	const Int_t NTRIG = 3;
	TFile *Unfolding[NTRIG];
	Unfolding[0] = new TFile("../RooUnfold/Jet_Unfolding_Corrections_13TeV_4_INT7.root","read");
	Unfolding[1] = new TFile("../RooUnfold/Jet_Unfolding_Corrections_13TeV_4_EG1.root","read");
	Unfolding[2] = new TFile("../RooUnfold/Jet_Unfolding_Corrections_13TeV_4_EG2.root","read");
	TH1D* histoUnfolding_Pi0[NTRIG];
	TH1D* histoUnfolding_Eta[NTRIG];
	TF1* functionUnfolding_Pi0[NTRIG];
	TF1* functionUnfolding_Eta[NTRIG];
	Int_t numberOfBinsPi0[NTRIG];	
	Int_t numberOfBinsEta[NTRIG];

	Double_t x_Pi0_Unfolding[3][100] = {0};
	Double_t y_Pi0_Unfolding[3][100] = {0};
	Double_t xerror_Pi0_Unfolding[3][100] = {0};
	Double_t yerror_Pi0_Unfolding[3][100] = {0};

	Double_t x_Eta_Unfolding[3][100];
	Double_t y_Eta_Unfolding[3][100];
	Double_t xerror_Eta_Unfolding[3][100];
	Double_t yerror_Eta_Unfolding[3][100];

	Double_t x_Pi0_Combined[200] = {0};
	Double_t xerror_Pi0_Combined[200] = {0};
	Double_t y_Pi0_Combined[200] = {0};
	Double_t yerror_Pi0_Combined[200] = {0};

	Double_t x_Eta_Combined[200] = {0};
	Double_t xerror_Eta_Combined[200] = {0};
	Double_t y_Eta_Combined[200] = {0};
	Double_t yerror_Eta_Combined[200] = {0};

	for(Int_t itrig = 0; itrig < NTRIG; itrig++){
		for(Int_t itrig = 0; itrig < NTRIG; itrig++){
		histoUnfolding_Pi0[itrig] = (TH1D*) Unfolding[itrig]->Get("Ratio_Total_Pi0");
		histoUnfolding_Eta[itrig] = (TH1D*) Unfolding[itrig]->Get("Ratio_Total_Eta");
		functionUnfolding_Pi0[itrig] = (TF1*) Unfolding[itrig]->Get("FinalFit_Pi0");
		functionUnfolding_Eta[itrig] = (TF1*) Unfolding[itrig]->Get("FinalFit_Eta");
		numberOfBinsPi0[itrig] = histoUnfolding_Pi0[itrig]->GetNbinsX();
		numberOfBinsEta[itrig] = histoUnfolding_Eta[itrig]->GetNbinsX();
		}
	}	
	Int_t ind_Pi0 =0; //index to keep track
	Int_t ind_Eta =0; //index to keep track
	for(Int_t itrig = 0; itrig < NTRIG; itrig++){
		for(Int_t i =0; i < numberOfBinsPi0[itrig]; i++){
			x_Pi0_Combined[ind_Pi0]        = histoUnfolding_Pi0[itrig]->GetBinCenter(i+1);
			y_Pi0_Combined[ind_Pi0]        = histoUnfolding_Pi0[itrig]->GetBinContent(i+1);
			xerror_Pi0_Combined[ind_Pi0]   = histoUnfolding_Pi0[itrig]->GetBinWidth(i+1)/2;
			yerror_Pi0_Combined[ind_Pi0]   = histoUnfolding_Pi0[itrig]->GetBinError(i+1);
			ind_Pi0++;
		} 
		for(Int_t i =0; i < numberOfBinsEta[itrig]; i++){
			x_Eta_Combined[ind_Eta]        = histoUnfolding_Eta[itrig]->GetBinCenter(i+1);
			y_Eta_Combined[ind_Eta]        = histoUnfolding_Eta[itrig]->GetBinContent(i+1);
			xerror_Eta_Combined[ind_Eta]   = histoUnfolding_Eta[itrig]->GetBinWidth(i+1)/2;
			yerror_Eta_Combined[ind_Eta]   = histoUnfolding_Eta[itrig]->GetBinError(i+1);
			ind_Eta++;
		} 
		for(Int_t i =0; i < numberOfBinsPi0[itrig]; i++){
			x_Pi0_Unfolding[itrig][i]      = histoUnfolding_Pi0[itrig]->GetBinCenter(i+1);
			y_Pi0_Unfolding[itrig][i]      = histoUnfolding_Pi0[itrig]->GetBinContent(i+1);
			xerror_Pi0_Unfolding[itrig][i] = histoUnfolding_Pi0[itrig]->GetBinWidth(i+1)/2;
			yerror_Pi0_Unfolding[itrig][i] = histoUnfolding_Pi0[itrig]->GetBinError(i+1);
		} 
		for(Int_t i =0; i < numberOfBinsEta[itrig]; i++){
			x_Eta_Unfolding[itrig][i]      = histoUnfolding_Eta[itrig]->GetBinCenter(i+1);
			y_Eta_Unfolding[itrig][i]      = histoUnfolding_Eta[itrig]->GetBinContent(i+1);
			xerror_Eta_Unfolding[itrig][i] = histoUnfolding_Eta[itrig]->GetBinWidth(i+1)/2;
			yerror_Eta_Unfolding[itrig][i] = histoUnfolding_Eta[itrig]->GetBinError(i+1);
		}
	} 	
// Pi0 meson
	TGraphErrors* graphUnfoldingPi0[NTRIG];
	for (Int_t itrig =0; itrig < NTRIG; itrig++){
	graphUnfoldingPi0[itrig] = new TGraphErrors(numberOfBinsPi0[itrig],x_Pi0_Unfolding[itrig],y_Pi0_Unfolding[itrig], xerror_Pi0_Unfolding[itrig],yerror_Pi0_Unfolding[itrig]);
	}
	TGraphErrors* graphUnfoldingPi0_Combined = new TGraphErrors(ind_Pi0, x_Pi0_Combined, y_Pi0_Combined, xerror_Pi0_Combined, yerror_Pi0_Combined);
	for(Int_t itrig =0; itrig< NTRIG; itrig++){
	graphUnfoldingPi0[itrig]->SetMarkerSize(3);
	graphUnfoldingPi0[itrig]->SetMarkerStyle(20+itrig);
	graphUnfoldingPi0[itrig]->SetMarkerColor(kGray+itrig);
	graphUnfoldingPi0[itrig]->SetLineColor(kGray+itrig);
	graphUnfoldingPi0[itrig]->SetLineWidth(1);
	graphUnfoldingPi0[itrig]->SetFillStyle(0);
	}
		
	TMultiGraph *multigraphUnfoldingPi0 = new TMultiGraph();
	multigraphUnfoldingPi0->SetTitle(";#it{p}_{T} (GeV/#it{c});");
	multigraphUnfoldingPi0->Add(graphUnfoldingPi0[0]);
	multigraphUnfoldingPi0->Add(graphUnfoldingPi0[1]);
	multigraphUnfoldingPi0->Add(graphUnfoldingPi0[2]);
	multigraphUnfoldingPi0->GetHistogram()->GetXaxis()->SetRangeUser(1., 20.);
	multigraphUnfoldingPi0->GetHistogram()->GetYaxis()->SetRangeUser(1.,1.6);
	TLegend* legendMultiUnfolding_Pi0   = GetAndSetLegend2(0.12+0.6,0.74, 0.45+0.6, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
 	legendMultiUnfolding_Pi0 ->AddEntry(graphUnfoldingPi0[0],"INT7","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendMultiUnfolding_Pi0 ->AddEntry(graphUnfoldingPi0[1],"EG1","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendMultiUnfolding_Pi0 ->AddEntry(graphUnfoldingPi0[2],"EG2","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	//legendMultiUnfolding_Pi0 ->AddEntry(graphUnfoldingPi0_Combined,"Combined","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	//graphUnfoldingPi0_Combined->Draw();
// Eta meson
	TGraphErrors* graphUnfoldingEta[NTRIG];
	for (Int_t itrig =0; itrig < NTRIG; itrig++){
	graphUnfoldingEta[itrig] = new TGraphErrors(numberOfBinsEta[itrig],x_Eta_Unfolding[itrig],y_Eta_Unfolding[itrig], xerror_Eta_Unfolding[itrig],yerror_Eta_Unfolding[itrig]);
	}
	TGraphErrors* graphUnfoldingEta_Combined = new TGraphErrors(ind_Eta, x_Eta_Combined, y_Eta_Combined, xerror_Eta_Combined, yerror_Eta_Combined);
	for(Int_t itrig =0; itrig< NTRIG; itrig++){
	graphUnfoldingEta[itrig]->SetMarkerSize(3);
	graphUnfoldingEta[itrig]->SetMarkerStyle(20+itrig);
	graphUnfoldingEta[itrig]->SetMarkerColor(kGray+itrig);
	graphUnfoldingEta[itrig]->SetLineColor(kGray+itrig);
	graphUnfoldingEta[itrig]->SetLineWidth(1);
	graphUnfoldingEta[itrig]->SetFillStyle(0);
	}
		
	TMultiGraph *multigraphUnfoldingEta = new TMultiGraph();
	multigraphUnfoldingEta->SetTitle(";#it{p}_{T} (GeV/#it{c});");
	multigraphUnfoldingEta->Add(graphUnfoldingEta[0]);
	multigraphUnfoldingEta->Add(graphUnfoldingEta[1]);
	multigraphUnfoldingEta->Add(graphUnfoldingEta[2]);
	multigraphUnfoldingEta->GetHistogram()->GetXaxis()->SetRangeUser(1.,25.);
	multigraphUnfoldingEta->GetHistogram()->GetYaxis()->SetRangeUser(1.,1.6);
	TLegend* legendMultiUnfolding_Eta   = GetAndSetLegend2(0.12+0.6,0.74, 0.45+0.6, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
 	legendMultiUnfolding_Eta ->AddEntry(graphUnfoldingEta[0],"INT7","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendMultiUnfolding_Eta ->AddEntry(graphUnfoldingEta[1],"EG1","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendMultiUnfolding_Eta ->AddEntry(graphUnfoldingEta[2],"EG2","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	//legendMultiUnfolding_Eta ->AddEntry(graphUnfoldingEta_Combined,"Combined","p"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	//graphUnfoldingEta_Combined->Draw();

	canvas_compare_Unfolding_Pi0_DataPoints->cd();

	multigraphUnfoldingPi0->Draw("AP"); 
	legendMultiUnfolding_Pi0->Draw("same");
  	TF1 *functionUnfolding_Pi0_Combined = new TF1 ("functionUnfolding_Pi0_Combined", "[0]+TMath::Exp([1]/x)", 2, 25);
//  	functionUnfolding_Pi0_Combined -> SetParameters(0.5, 1.0);
 	functionUnfolding_Pi0_Combined -> SetParameters(0.5, 1.0);
	//graphUnfoldingPi0_Combined->Fit(functionUnfolding_Pi0_Combined, "RQ0", "", 2, 15);
	graphUnfoldingPi0_Combined->Fit("functionUnfolding_Pi0_Combined", "R");
	functionUnfolding_Pi0_Combined->Draw("same");
	t->DrawLatex(3, 1.55, Form("Pi0, Combined: ") );
	t->DrawLatex(3, 1.50, Form("%.3e + exp(%.2e/x)", functionUnfolding_Pi0_Combined->GetParameter(0),functionUnfolding_Pi0_Combined->GetParameter(1)) );
	
//	canvas_compare_Unfolding_Pi0_DataPoints->Write();

	canvas_compare_Unfolding_Eta_DataPoints->cd();

	multigraphUnfoldingEta->Draw("AP"); 
	legendMultiUnfolding_Eta->Draw("same");
	TF1 *functionUnfolding_Eta_Combined = new TF1 ("functionUnfolding_Eta_Combined", "[0]+TMath::Exp([1]/x)", 2, 25);
	functionUnfolding_Eta_Combined -> SetParameters(0.5, 1.0);
	//graphUnfoldingEta_Combined->Fit(functionUnfolding_Eta_Combined, "RQ0", "", 2, 15);
	graphUnfoldingEta_Combined->Fit("functionUnfolding_Eta_Combined", "R");
	functionUnfolding_Eta_Combined->Draw("same");
	t->DrawLatex(3, 1.90, Form("Eta, Combined: ") );
	t->DrawLatex(3, 1.80, Form("%.3e + exp(%.2e/x)", functionUnfolding_Eta_Combined->GetParameter(0),functionUnfolding_Eta_Combined->GetParameter(1)) );


	////////////// compare unfolding //////////////////////////
	canvas_compare_Unfolding_Pi0->cd();
	functionUnfolding_Pi0[0]->SetLineColor(2);
	functionUnfolding_Pi0[0]->GetYaxis()->SetRangeUser(1.,1.5);
	functionUnfolding_Pi0[0]->SetTitle();
	functionUnfolding_Pi0[1]->SetLineColor(3);
	functionUnfolding_Pi0[2]->SetLineColor(4);
	TLegend* legendUnfolding_Pi0   = GetAndSetLegend2(0.12, 0.8, 0.45, 0.8+(0.14*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendUnfolding_Pi0->AddEntry(functionUnfolding_Pi0[0], "INT7", "l"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendUnfolding_Pi0->AddEntry(functionUnfolding_Pi0[1], "EG1", "l"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendUnfolding_Pi0->AddEntry(functionUnfolding_Pi0[2], "EG2", "l"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	functionUnfolding_Pi0[0]->Draw();
	functionUnfolding_Pi0[1]->Draw("same");
	functionUnfolding_Pi0[2]->Draw("same");
	legendUnfolding_Pi0->Draw("same");	
	
	canvas_compare_Unfolding_Eta->cd();
	functionUnfolding_Eta[0]->SetLineColor(2);
	functionUnfolding_Eta[0]->GetYaxis()->SetRangeUser(1.,1.5);
	functionUnfolding_Eta[0]->SetTitle();
	functionUnfolding_Eta[1]->SetLineColor(3);
	functionUnfolding_Eta[2]->SetLineColor(4);
	TLegend* legendUnfolding_Eta   = GetAndSetLegend2(0.12, 0.8, 0.45, 0.8+(0.14*expectedLinesInLegend), textSizeLabelsPixel, 2, "", 43, 0);
 	legendUnfolding_Eta->AddEntry(functionUnfolding_Eta[0], "INT7", "l"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendUnfolding_Eta->AddEntry(functionUnfolding_Eta[1], "EG1", "l"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
 	legendUnfolding_Eta->AddEntry(functionUnfolding_Eta[2], "EG2", "l"); // last argument: p (point), l (line), f (fill), e (error), combinations possible, ie "pf", "lep", ...
	functionUnfolding_Eta[0]->Draw();
	functionUnfolding_Eta[1]->Draw("same");
	functionUnfolding_Eta[2]->Draw("same");
	legendUnfolding_Eta->Draw("same");	
	
	TFile *MyAnalysis = new TFile("../RooUnfold/Jet_Unfolding_Corrections_13TeV_4.root","recreate");
	functionUnfolding_Pi0_Combined->Write("FinalFit_Pi0");
	functionUnfolding_Eta_Combined->Write("FinalFit_Eta");
	canvas_compare_Unfolding_Pi0->Write("CompareUnfolding_Pi0");
	canvas_compare_Unfolding_Eta->Write("CompareUnfolding_Eta");
	canvas_compare_Unfolding_Pi0_DataPoints->Write("CompareUnfolding_Pi0_DataPoints");
	canvas_compare_Unfolding_Eta_DataPoints->Write("CompareUnfolding_Eta_DataPoints");

	MyAnalysis->Close();

	canvas_compare_Unfolding_Pi0->SaveAs("CompareUnfolding_Pi0.eps");
	canvas_compare_Unfolding_Eta->SaveAs("CompareUnfolding_Eta.eps");
	canvas_compare_Unfolding_Pi0_DataPoints->SaveAs("CompareUnfolding_Pi0_DataPoints.eps");
	canvas_compare_Unfolding_Eta_DataPoints->SaveAs("CompareUnfolding_Eta_DataPoints.eps");

}
