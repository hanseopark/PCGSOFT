#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"
#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include <string>
#include "TGaxis.h"
#include "TFractionFitter.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "THStack.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TAttAxis.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "TEllipse.h"
#include "TPaveText.h"
#include "TRandom3.h"

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"

#include "CommonHeaders/ExtractSignalBinning.h"
//#include "CommonHeaders/ExtractSignalPlotting.h"


#include <iostream>
using namespace std;

// 내일 할꺼 ! for 구문 안 x값이랑 y값이 일치안함. 또, maxbin을 조정해서 알맞게 해야함.   
// 오류 발견 !! 분모의 0이나와서 그렇다 !! 현재는 통계수가 적으니 어쩔수없음. 
// change fiiting method. For example, expenatial + linerarity and linearity1(0~10GeV) + linearity2(10GeV~20GeV)


void Jet_Unfolding_Macro_EG2(
    TString CutNumber = "",
    TString Energy = "",
    Int_t Mode = 0,
    Int_t NumberOfBinsPi0 = 0,
    Int_t NumberOfBinsEta = 0
)
{
//	NumberOfBinsPi0 = NumberOfBinsPi0-1; //-10;	//NuberOfBinsMeson - StartBins;
//	NumberOfBinsEta = NumberOfBinsEta-1; //-4;	//8;

   TH1D* Histo_AsData_Pi0;
   TH1D* Histo_Missed_Pi0;
   TH1D* Histo_Reject_Pi0;
   TH1D* Histo_AsData_Eta;
   TH1D* Histo_Missed_Eta;
   TH1D* Histo_Reject_Eta;
   TFile *filePi0Data = new TFile(Form("%s/%s_Unfolding_AsData/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",CutNumber.Data(),Energy.Data(),CutNumber.Data()));
   Histo_AsData_Pi0                 = (TH1D*)filePi0Data->Get("histoYieldMeson");
   TFile *filePi0Missed = new TFile(Form("%s/%s_Unfolding_Missed/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",CutNumber.Data(),Energy.Data(),CutNumber.Data()));
   Histo_Missed_Pi0                 = (TH1D*)filePi0Missed->Get("histoYieldMeson");
   TFile *filePi0Reject = new TFile(Form("%s/%s_Unfolding_Reject/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",CutNumber.Data(),Energy.Data(),CutNumber.Data()));
   Histo_Reject_Pi0                 = (TH1D*)filePi0Reject->Get("histoYieldMeson");

   TFile *fileEtaData = new TFile(Form("%s/%s_Unfolding_AsData/Eta_MC_GammaConvV1WithoutCorrection_%s.root",CutNumber.Data(),Energy.Data(),CutNumber.Data()));
   Histo_AsData_Eta                 = (TH1D*)fileEtaData->Get("histoYieldMeson");
   TFile *fileEtaMissed = new TFile(Form("%s/%s_Unfolding_Missed/Eta_MC_GammaConvV1WithoutCorrection_%s.root",CutNumber.Data(),Energy.Data(),CutNumber.Data()));
   Histo_Missed_Eta                 = (TH1D*)fileEtaMissed->Get("histoYieldMeson");
   TFile *fileEtaReject = new TFile(Form("%s/%s_Unfolding_Reject/Eta_MC_GammaConvV1WithoutCorrection_%s.root",CutNumber.Data(),Energy.Data(),CutNumber.Data()));
   Histo_Reject_Eta                 = (TH1D*)fileEtaReject->Get("histoYieldMeson");

   Int_t NBins = 200;
   Double_t* BinEdgesPi0 = new Double_t[300];
   Double_t* BinEdgesEta = new Double_t[300];
   Int_t MaxBinPi0 = GetBinning( BinEdgesPi0, NBins, "Pi0", Energy, Mode, 3, kFALSE, "", kTRUE ); // speicaltirgg set is needed
   Int_t MaxBinEta = GetBinning( BinEdgesEta, NBins, "Eta", Energy, Mode, 3, kFALSE, "", kTRUE ); // specialtrigg set is needed

        
//   Int_t StartBinPi0 = GetStartBin("Pi0", Energy, Mode, -1, "", "", kTRUE);
//   Int_t StartBinEta = GetStartBin("Eta", Energy, Mode, -1, "", "", kTRUE);
   Int_t StartBinPi0 = GetStartBin("Pi0", Energy, Mode, 3, "", "", kTRUE); // specialtrigg set is needed
   Int_t StartBinEta = GetStartBin("Eta", Energy, Mode, 3, "", "", kTRUE); // specialtrigg set is needed
    
   TH1D* Ratio_Added_Pi0  = new TH1D("Ratio_Added_Pi0","Ratio_Added_Pi0",NumberOfBinsPi0,BinEdgesPi0);
   TH1D* Ratio_Subtr_Pi0  = new TH1D("Ratio_Subtr_Pi0","Ratio_Subtr_Pi0",NumberOfBinsPi0,BinEdgesPi0);
   TH1D* Ratio_Total_Pi0  = new TH1D("Ratio_Total_Pi0","Ratio_Total_Pi0",NumberOfBinsPi0,BinEdgesPi0);
   
   TH1D* Ratio_Added_Eta  = new TH1D("Ratio_Added_Eta","Ratio_Added_Eta",NumberOfBinsEta,BinEdgesEta);
   TH1D* Ratio_Subtr_Eta  = new TH1D("Ratio_Subtr_Eta","Ratio_Subtr_Eta",NumberOfBinsEta,BinEdgesEta);
   TH1D* Ratio_Total_Eta  = new TH1D("Ratio_Total_Eta","Ratio_Total_Eta",NumberOfBinsEta,BinEdgesEta);
   
   Double_t Value_Data;
   Double_t Error_Data;
   Double_t Value_Miss;
   Double_t Error_Miss;
   Double_t Value_Rej;
   Double_t Error_Rej;

   for(Int_t i =StartBinPi0+1; i<=NumberOfBinsPi0 ; i++){
     Value_Data = Histo_AsData_Pi0->GetBinContent(i);
     Error_Data = Histo_AsData_Pi0->GetBinError(i);
     Value_Miss = Histo_Missed_Pi0->GetBinContent(i);
     Error_Miss = Histo_Missed_Pi0->GetBinError(i);
     Value_Rej = Histo_Reject_Pi0->GetBinContent(i);
     Error_Rej = Histo_Reject_Pi0->GetBinError(i);
        
     Double_t ErrorRatioAdded = (Value_Miss/Value_Data)*TMath::Sqrt(pow(Error_Miss/Value_Miss,2)+pow(Error_Data/Value_Data,2));
     Double_t ErrorRatioSubtr = (Value_Rej/Value_Data)*TMath::Sqrt(pow(Error_Rej/Value_Rej,2)+pow(Error_Data/Value_Data,2));

     Ratio_Added_Pi0->SetBinContent(i,Value_Miss/Value_Data + 1);
     Ratio_Added_Pi0->SetBinError(i,ErrorRatioAdded);
        
     Ratio_Subtr_Pi0->SetBinContent(i,1 - Value_Rej/Value_Data);
     Ratio_Subtr_Pi0->SetBinError(i,ErrorRatioSubtr);
        
     Ratio_Total_Pi0->SetBinContent(i,Value_Miss/Value_Data - Value_Rej/Value_Data + 1);
     Ratio_Total_Pi0->SetBinError(i,ErrorRatioAdded + ErrorRatioSubtr);
	
//	cout << "Value_Data[" << i <<"] : "<< Value_Data<< endl;   
//	cout << "Value_Miss[" << i <<"] : "<< Value_Miss<< endl;   
//	cout << "Value_Rej[" << i <<"] : "<< Value_Rej<< endl;   
 if(Value_Data == 0 || (Value_Miss == 0 || Value_Rej == 0)) cout<< "bin: "<< i << "Pi0: Value Data or Miss or Rej =0 Maybe make the null. low statistics" <<endl;

}
   for(Int_t i =StartBinEta+1; i<=NumberOfBinsEta ; i++){
     Value_Data = Histo_AsData_Eta->GetBinContent(i);
     Error_Data = Histo_AsData_Eta->GetBinError(i);
     Value_Miss = Histo_Missed_Eta->GetBinContent(i);
     Error_Miss = Histo_Missed_Eta->GetBinError(i);
     Value_Rej = Histo_Reject_Eta->GetBinContent(i);
     Error_Rej = Histo_Reject_Eta->GetBinError(i);
        
     Double_t ErrorRatioAdded = (Value_Miss/Value_Data)*TMath::Sqrt(pow(Error_Miss/Value_Miss,2)+pow(Error_Data/Value_Data,2));
     Double_t ErrorRatioSubtr = (Value_Rej/Value_Data)*TMath::Sqrt(pow(Error_Rej/Value_Rej,2)+pow(Error_Data/Value_Data,2));

     Ratio_Added_Eta->SetBinContent(i,Value_Miss/Value_Data + 1);
     Ratio_Added_Eta->SetBinError(i,ErrorRatioAdded);
        
     Ratio_Subtr_Eta->SetBinContent(i,1 - Value_Rej/Value_Data);
     Ratio_Subtr_Eta->SetBinError(i,1,ErrorRatioSubtr);
        
     Ratio_Total_Eta->SetBinContent(i,Value_Miss/Value_Data - Value_Rej/Value_Data + 1);
     Ratio_Total_Eta->SetBinError(i,ErrorRatioAdded + ErrorRatioSubtr);
   	
//	cout << "Value_Data[" << i <<"] : "<< Value_Data<< endl;   
//	cout << "Value_Miss[" << i <<"] : "<< Value_Miss<< endl;   
//	cout << "Value_Rej[" << i <<"] : "<< Value_Rej<< endl;   
 if(Value_Data == 0 || (Value_Miss == 0 || Value_Rej == 0)) cout<< "bin: "<< i <<"Eta: Value Data or Miss or Rej =0 Maybe make the null. low statistics" <<endl;
}
    
   TCanvas* c1 = new TCanvas("c1","",0,0,850,650);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetTitleFontSize(0.08);
   gPad->SetRightMargin(0.015);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.1);
   gPad->SetTickx();
   gPad->SetTicky();   
    
   TH1D* hEmpty1 = new TH1D("hEmpty1","",100,0,21);
// change the fitting. For example, pol1 -> exp()    
//   TF1* fit = new TF1("fit","pol1",0,20);
//   Ratio_Total_Pi0->Fit("fit","Q","R");
   TF1* fit = new TF1("fit", "[0]+exp(-[1]/x)",1.4,20);
   Ratio_Total_Pi0->Fit("fit","Q","R");

   cout<<"Parameter 0: "<<fit->GetParameter(0)<<endl;
   cout<<"Error Par 0: "<<fit->GetParError(0)<<endl;
   cout<<"Parameter 1: "<<fit->GetParameter(1)<<endl;
   cout<<"Error Par 1: "<<fit->GetParError(1)<<endl;
    
   hEmpty1->GetYaxis()->SetTitle("Correction");
   hEmpty1->GetYaxis()->SetLabelFont(65);
   hEmpty1->GetYaxis()->SetLabelSize(14);
   hEmpty1->GetYaxis()->SetTitleFont(63);
   hEmpty1->GetYaxis()->SetTitleSize(20);
   hEmpty1->GetYaxis()->SetTitleOffset(1);
   hEmpty1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   hEmpty1->GetXaxis()->SetLabelFont(65);
   hEmpty1->GetXaxis()->SetLabelSize(14);
   hEmpty1->GetXaxis()->SetTitleFont(63);
   hEmpty1->GetXaxis()->SetTitleSize(20);
   hEmpty1->GetXaxis()->SetTitleOffset(1);
   hEmpty1->GetXaxis()->SetRangeUser(0,21);
   hEmpty1->GetYaxis()->SetRangeUser(0.4,2.2);
    
   Ratio_Added_Pi0->SetMarkerStyle(20);
   Ratio_Added_Pi0->SetMarkerSize(1);
   Ratio_Added_Pi0->SetMarkerColor(kCyan+2);
   Ratio_Added_Pi0->SetLineColor(kCyan+2);
   Ratio_Subtr_Pi0->SetMarkerStyle(20);
   Ratio_Subtr_Pi0->SetMarkerSize(1);
   Ratio_Subtr_Pi0->SetMarkerColor(kGreen+2);
   Ratio_Subtr_Pi0->SetLineColor(kGreen+2);
   Ratio_Total_Pi0->SetMarkerStyle(20);
   Ratio_Total_Pi0->SetMarkerSize(1);
   Ratio_Total_Pi0->SetMarkerColor(kBlack);
   Ratio_Total_Pi0->SetLineColor(kBlack);
    
   hEmpty1->Draw();    
   Ratio_Added_Pi0->Draw("SAME P");
   Ratio_Subtr_Pi0->Draw("SAME P");
   Ratio_Total_Pi0->Draw("SAME P");
   TLine *line1 = new TLine(0,1,21,1);
   line1->SetLineWidth(2);
   line1->SetLineStyle(1);
   line1->SetLineColor(kBlack);
   line1->Draw();
   
   TString ModeString;
   if(Mode == 0) ModeString = "PCM";
   if(Mode == 2) ModeString = "PCM-EMCal";
   if(Mode == 3) ModeString = "PCM-PHOS";
   if(Mode == 4) ModeString = "EMCal";
   if(Mode == 5) ModeString = "PHOS";
    
   TLatex T1;
   T1.SetTextSize(0.025);
   T1.SetTextAlign(12);
   T1.SetNDC();
   T1.DrawLatex(0.3, 0.90, Form("#pi^{0} #rightarrow #gamma#gamma with %s in charged jets",ModeString.Data()));
   T1.DrawLatex(0.3, 0.85, "Jet #it{p}_{T} > 10 GeV/c");
   T1.DrawLatex(0.3, 0.80, "pp, #sqrt{S_{NN}} = 13 TeV");
    
   TLegend* leg1 = new TLegend(0.75,0.73,0.95,0.95);     
   leg1->AddEntry(Ratio_Added_Pi0,"#frac{Added}{data}","lp");
   leg1->AddEntry(Ratio_Subtr_Pi0,"#frac{Subtracted}{data}","lp");
   leg1->AddEntry(Ratio_Total_Pi0,"#frac{Total correction}{data}","lp");
   leg1->SetBorderSize(0);
   leg1->SetFillStyle(0);
   leg1->Draw();
    
   c1->Update();   
   
   TCanvas* c2 = new TCanvas("c2","",0,0,850,650);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetTitleFontSize(0.08);
   gPad->SetRightMargin(0.015);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.02);
   gPad->SetLeftMargin(0.1);
   gPad->SetTickx();
   gPad->SetTicky();   
    
   TH1D* hEmpty2 = new TH1D("hEmpty2","",100,0,21);
    
//   TF1* fit2 = new TF1("fit2","pol1",0,20);
//   Ratio_Total_Eta->Fit("fit2","Q","R");
   TF1* fit2 = new TF1("fit2","[0]+exp(-[1]/x)",2,25);
   Ratio_Total_Eta->Fit("fit2","Q","R");


   cout<<"Parameter 0: "<<fit2->GetParameter(0)<<endl;
   cout<<"Error Par 0: "<<fit2->GetParError(0)<<endl;
   cout<<"Parameter 1: "<<fit2->GetParameter(1)<<endl;
   cout<<"Error Par 1: "<<fit2->GetParError(1)<<endl;
    
   hEmpty2->GetYaxis()->SetTitle("Correction");
   hEmpty2->GetYaxis()->SetLabelFont(65);
   hEmpty2->GetYaxis()->SetLabelSize(14);
   hEmpty2->GetYaxis()->SetTitleFont(63);
   hEmpty2->GetYaxis()->SetTitleSize(20);
   hEmpty2->GetYaxis()->SetTitleOffset(1);
   hEmpty2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   hEmpty2->GetXaxis()->SetLabelFont(65);
   hEmpty2->GetXaxis()->SetLabelSize(14);
   hEmpty2->GetXaxis()->SetTitleFont(63);
   hEmpty2->GetXaxis()->SetTitleSize(20);
   hEmpty2->GetXaxis()->SetTitleOffset(1);
   hEmpty2->GetXaxis()->SetRangeUser(0,21);
   hEmpty2->GetYaxis()->SetRangeUser(0.0,2.2);
    
   Ratio_Added_Eta->SetMarkerStyle(20);
   Ratio_Added_Eta->SetMarkerSize(1);
   Ratio_Added_Eta->SetMarkerColor(kCyan+2);
   Ratio_Added_Eta->SetLineColor(kCyan+2);
   Ratio_Subtr_Eta->SetMarkerStyle(20);
   Ratio_Subtr_Eta->SetMarkerSize(1);
   Ratio_Subtr_Eta->SetMarkerColor(kGreen+2);
   Ratio_Subtr_Eta->SetLineColor(kGreen+2);
   Ratio_Total_Eta->SetMarkerStyle(20);
   Ratio_Total_Eta->SetMarkerSize(1);
   Ratio_Total_Eta->SetMarkerColor(kBlack);
   Ratio_Total_Eta->SetLineColor(kBlack);
    
   hEmpty2->Draw();    
   Ratio_Added_Eta->Draw("SAME P");
   Ratio_Subtr_Eta->Draw("SAME P");
   Ratio_Total_Eta->Draw("SAME P");
   TLine *line2 = new TLine(0,1,21,1);
   line2->SetLineWidth(2);
   line2->SetLineStyle(1);
   line2->SetLineColor(kBlack);
   line2->Draw();
    
   TLatex T2;
   T2.SetTextSize(0.025);
   T2.SetTextAlign(12);
   T2.SetNDC();
   T2.DrawLatex(0.3, 0.90, Form("#eta #rightarrow #gamma#gamma with %s in charged jets",ModeString.Data()));
   T2.DrawLatex(0.3, 0.85, "Jet #it{p}_{T} > 10 GeV/c");
   T2.DrawLatex(0.3, 0.80, "pp, #sqrt{S_{NN}} = 13 TeV");
    
   TLegend* leg2 = new TLegend(0.75,0.73,0.95,0.95);     
   leg2->AddEntry(Ratio_Added_Eta,"#frac{Added}{data}","lp");
   leg2->AddEntry(Ratio_Subtr_Eta,"#frac{Subtracted}{data}","lp");
   leg2->AddEntry(Ratio_Total_Eta,"#frac{Total correction}{data}","lp");
   leg2->SetBorderSize(0);
   leg2->SetFillStyle(0);
   leg2->Draw();
    
   c2->Update();
   
   TFile *outputFile = new TFile(Form("RooUnfold/Jet_Unfolding_Corrections_%s_%i_EG2.root",Energy.Data(),Mode),"RECREATE");
   Ratio_Added_Pi0->Write("Ratio_Added_Pi0");
   Ratio_Subtr_Pi0->Write("Ratio_Subtr_Pi0");
   Ratio_Total_Pi0->Write("Ratio_Total_Pi0");
   fit->Write("FinalFit_Pi0");
   Ratio_Added_Eta->Write("Ratio_Added_Eta");
   Ratio_Subtr_Eta->Write("Ratio_Subtr_Eta");
   Ratio_Total_Eta->Write("Ratio_Total_Eta");
   fit2->Write("FinalFit_Eta");
	 
   c1->Write();
   c2->Write();

   outputFile->Close();

 
}
