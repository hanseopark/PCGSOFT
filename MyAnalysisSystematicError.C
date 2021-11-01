//////////////////////////////////////////////////////////////
///////// My Analysis ///////////////////////////////////////
/////////////////////////////////////////////////////////////

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"
//#include "MyCommonHeaders/Filipad.h"

void MyAnalysisSystematicError(TString nameTrig = ""){

	Int_t textSizeLabelsPixel = 900*0.04;
	Int_t expectedLinesInLegend = 1;

	//Read root data
	Bool_t gTrigINT7 = kFALSE;
	Bool_t gTrigEG1  = kFALSE;
	Bool_t gTrigEG2  = kFALSE;


	const Int_t NISMC = 2;
	const Int_t NMESON = 3;
	//const Int_t NTRIG = 3;
	const Int_t NTRIG = 1;
	const Int_t NCUTS = 8;
	const Int_t NCUTSCALO = 6;
	const Int_t NCUTSMESON = 3;
	const Int_t NVAR = 4;
	Int_t NBIN = 99;
	//const char *Trig[NTRIG] = {"INT7", "EG1", "EG2"};
	const char *Cuts[NCUTS] = {"Standard", "ClusterMinEnergy", "ClusterNCells", "ClusterTiming", "ClusterM02", "ClusterTrackMatching", "OpeningAngle", "Alpha"};
	const char *CutsCalo[NCUTSCALO] = {"Standard", "ClusterMinEnergy", "ClusterNCells", "ClusterTiming", "ClusterM02", "ClusterTrackMatching"};
	const char *CutsMeson[NCUTSMESON] = {"Standard", "OpeningAngle", "Alpha"};
	const char *Var[NVAR] = {"Standard", "Cutvariation1", "Cutvariation2", "CutVariation3"};
	Int_t ColorVAR[NVAR];
	for (Int_t icolor = 1; icolor < NVAR+1; icolor++){
	ColorVAR[icolor] = icolor;
	}

	TString Ismc[NISMC];
	Ismc[0] = "MC";
	Ismc[1] = "Data";	

	TString Meson[NMESON];
	Meson[0] = "Pi0";
	Meson[1] = "Eta";
	Meson[2] = "Pi0EtaBinning";

	TString EventCuts;
	Int_t NSubWagon = 95;
	if (nameTrig.CompareTo("INT7") == 0) {EventCuts = "00010113"; NSubWagon = 95;}
	if (nameTrig.CompareTo("EG1") == 0) {EventCuts = "0008d113"; NSubWagon = 96;}
	if (nameTrig.CompareTo("EG2") == 0) {EventCuts = "0008e113"; NSubWagon = 97;}
	
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
	CaloCuts[2][3] = "411790607l022230000"; 

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

	TString EventCaloMesonCuts[NCUTS][NVAR];
	for (Int_t icuts = 0; icuts < NCUTS; icuts++){
		for (Int_t ivar = 0; ivar < NVAR; ivar++){
			EventCaloMesonCuts[icuts][ivar] = Form("%s_%s", EventCuts.Data(), CaloMesonCuts[icuts][ivar].Data());
		}
	}

	// read root file 
	TFile *GammaCalo[NISMC];
	TList *TopDir;
	TList *SubDir;
	TList *SubSubDir;
	TH1D* histoClusGamma_Pt[NISMC][NCUTS][NVAR];
	Int_t gBinContent[NISMC][NCUTS][NVAR][NBIN];
	Double_t gBin[NISMC][NCUTS][NVAR][NBIN];
	for (Int_t iismc=0; iismc<NISMC; iismc++){
		GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_950.root", Ismc[iismc].Data())); // it is included all trigger.
			TopDir = (TList*) GammaCalo[iismc]->Get("GammaCalo_950");
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[0][0].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[0][0].Data()));
					histoClusGamma_Pt[iismc][0][0] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][0][0][ibin] = histoClusGamma_Pt[iismc][0][0]->GetBinContent(ibin);	
				gBin[iismc][0][0][ibin] = histoClusGamma_Pt[iismc][0][0]->GetXaxis()->GetBinCenter(ibin);
				}
		delete TopDir;
//		delete SubDir;
//		delete SubSubDir;
		delete GammaCalo[iismc];
	}
	cout << "//////////////////" <<endl;
	cout << "STandard " << endl;
	cout << "//////////////////" <<endl;
	for(Int_t ibin=0; ibin<NBIN; ibin++){
	cout << gBin[0][0][0][ibin] << "    " << gBinContent[0][0][0][ibin] << endl;
	}
	
	// ClusterMinEnergy
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d1.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d1",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[1][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[1][ivar].Data()));
					histoClusGamma_Pt[iismc][1][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][1][ivar][ibin] = histoClusGamma_Pt[iismc][1][ivar]->GetBinContent(ibin);	
				gBin[iismc][1][ivar][ibin] = histoClusGamma_Pt[iismc][1][ivar]->GetXaxis()->GetBinCenter(ibin);
				}
			}
		delete TopDir;
		delete GammaCalo[iismc];
	}
	cout << "//////////////////" <<endl;
	cout << "ClusterMinEnergy " << endl;
	cout << "//////////////////" <<endl;
	for (Int_t ivar=1; ivar<NVAR; ivar++){
		cout << "//////////////////" <<endl;
		cout << "Cut Variation: " << ivar << endl;
		cout << "//////////////////" <<endl;
		for(Int_t ibin=0; ibin<NBIN; ibin++){
		cout << gBin[0][1][ivar][ibin] << "    " << gBinContent[0][1][ivar][ibin] << endl;
		}		
	}
	cout << gBin[0][1][1][10] << "    " << gBinContent[0][1][1][10] << endl;
	cout << "DEBUG: " << __LINE__ << endl;
	// ClusterNCells
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d1.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d1",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[2][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[2][ivar].Data()));
					histoClusGamma_Pt[iismc][2][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][2][ivar][ibin] = histoClusGamma_Pt[iismc][2][ivar]->GetBinContent(ibin+1);	
				gBin[iismc][2][ivar][ibin] = histoClusGamma_Pt[iismc][2][ivar]->GetXaxis()->GetBinCenter(ibin+1);
				}
			}
		delete TopDir;
		delete GammaCalo[iismc];
	}
	cout << gBin[0][2][1][10] << "    " << gBinContent[0][2][1][10] << endl;
	cout << "DEBUG: " << __LINE__ << endl;
	// ClusterTiming
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d2.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d2",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[3][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[3][ivar].Data()));
					histoClusGamma_Pt[iismc][3][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][3][ivar][ibin] = histoClusGamma_Pt[iismc][3][ivar]->GetBinContent(ibin+1);	
				gBin[iismc][3][ivar][ibin] = histoClusGamma_Pt[iismc][3][ivar]->GetXaxis()->GetBinCenter(ibin+1);
				}
			}
		delete TopDir;
		delete GammaCalo[iismc];
	}
	cout << gBin[0][3][1][10] << "    " << gBinContent[0][3][1][10] << endl;
	cout << "DEBUG: " << __LINE__ << endl;
	// ClusterM02
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d3.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d3",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[4][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[4][ivar].Data()));
					histoClusGamma_Pt[iismc][4][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][4][ivar][ibin] = histoClusGamma_Pt[iismc][4][ivar]->GetBinContent(ibin+1);	
				gBin[iismc][4][ivar][ibin] = histoClusGamma_Pt[iismc][4][ivar]->GetXaxis()->GetBinCenter(ibin+1);
				}
			}
		delete TopDir;
		delete GammaCalo[iismc];
	}
	cout << gBin[0][4][1][10] << "    " << gBinContent[0][4][1][10] << endl;
	cout << "DEBUG: " << __LINE__ << endl;
	// ClusterTrackMatching
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d3.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d3",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[5][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[5][ivar].Data()));
					histoClusGamma_Pt[iismc][5][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][5][ivar][ibin] = histoClusGamma_Pt[iismc][5][ivar]->GetBinContent(ibin+1);	
				gBin[iismc][5][ivar][ibin] = histoClusGamma_Pt[iismc][5][ivar]->GetXaxis()->GetBinCenter(ibin+1);
				}
			}
			delete TopDir;
			delete GammaCalo[iismc];
	}
	cout << gBin[0][5][1][10] << "    " << gBinContent[0][5][1][10] << endl;
cout << "DEBUG: " << __LINE__ << endl;
	// OpeningAngle
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d4.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d4",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[6][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[6][ivar].Data()));
					histoClusGamma_Pt[iismc][6][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][6][ivar][ibin] = histoClusGamma_Pt[iismc][6][ivar]->GetBinContent(ibin+1);	
				gBin[iismc][6][ivar][ibin] = histoClusGamma_Pt[iismc][6][ivar]->GetXaxis()->GetBinCenter(ibin+1);
				}
			}
			delete TopDir;
			delete GammaCalo[iismc];
	}	
	cout << gBin[0][6][1][10] << "    " << gBinContent[0][6][1][10] << endl;
	cout << "DEBUG: " << __LINE__ << endl;
	// Alpha	
	for (Int_t iismc=0; iismc<NISMC; iismc++){
			GammaCalo[iismc] = new TFile(Form("/home/alidock/alice/work/Data/pp_13TeV/Injet/%sMerged/GammaCalo_%d4.root", Ismc[iismc].Data(), NSubWagon), "read");
			TopDir = (TList*) GammaCalo[iismc]->Get(Form("GammaCalo_%d4",NSubWagon));
			for (Int_t ivar=1; ivar<NVAR; ivar++){
				SubDir = (TList*)TopDir->FindObject(Form("Cut Number %s",EventCaloMesonCuts[7][ivar].Data()));
					SubSubDir = (TList*) SubDir->FindObject(Form("%s ESD histograms",EventCaloMesonCuts[7][ivar].Data()));
					histoClusGamma_Pt[iismc][7][ivar] = (TH1D*) SubSubDir->FindObject("ClusGamma_Pt");
				for (Int_t ibin=0; ibin<NBIN; ibin++){
				gBinContent[iismc][7][ivar][ibin] = histoClusGamma_Pt[iismc][7][ivar]->GetBinContent(ibin+1);	
				gBin[iismc][7][ivar][ibin] = histoClusGamma_Pt[iismc][7][ivar]->GetXaxis()->GetBinCenter(ibin+1);
				}
			}
			delete TopDir;
			delete GammaCalo[iismc];
	}	
	cout << gBin[0][7][1][10] << "    " << gBinContent[0][7][1][10] << endl;
	cout << "DEBUG: " << __LINE__ << endl;
	cout << "DEBUG: "<< __LINE__ << endl;	

	// read histogram in root file
	cout << "DEBUG: "<< __LINE__ << endl;	


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

	// The First plot, Number of cout of cluster (%) vs variation

		// ClusterMinEnergy
 	// E(pT) 0.5-1.0, 1.0-1.5, 1.5-2.0
	const Int_t NEPtBin = 5;
	Int_t Count[NEPtBin][NVAR] = {0};
	Int_t TotalCount[NEPtBin] = {0};
	Double_t ReCount[NEPtBin][NVAR];
	Int_t EPtStart = 4;
	//Double_t EPtWidth = 2.5;
	Double_t EPtWidth =5.0; 
	Int_t EPtEnd   = EPtStart+EPtWidth;
	
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){	
		for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
			Count[iEPtBin][0] +=gBinContent[1][0][0][ibin];
			cout << "bin: " << gBinContent[1][0][0][ibin] << endl; 
			cout << "count: " << Count[iEPtBin][0] << endl;
		}
		for(Int_t ivar=1; ivar<NVAR; ivar++){
			cout << "ivar" << ivar <<endl;
			for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
				Count[iEPtBin][ivar] +=gBinContent[1][1][ivar][ibin];
				cout << "bin: " << gBinContent[1][1][ivar][ibin] <<endl;
				cout << "count: " << Count[iEPtBin][ivar] << endl;
			}
		}
		//TotalCount[iEPtBin] = Count[iEPtBin][0]+Count[iEPtBin][1]+Count[iEPtBin][2]+Count[iEPtBin][3]; // It is possibility to change the value of total count to all add like bin 1 to 200 (0~20GeV/c)
		for (Int_t ivar=0; ivar<NVAR; ivar++){
		TotalCount[iEPtBin] += Count[iEPtBin][ivar];
		}
	}
	cout << "Divide E(Pt)" << endl;
	cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
	cout << Count[1][0] << "    "<< Count[1][1] << "    " << Count[1][2] << "    " << Count[1][3] << endl;
	cout << Count[2][0] << "    "<< Count[2][1] << "    " << Count[2][2] << "    " << Count[2][3] << endl;
	// redefintion
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar=0; ivar<NVAR; ivar++){
			ReCount[iEPtBin][ivar] = static_cast<double>(Count[iEPtBin][ivar])/static_cast<double>(TotalCount[iEPtBin]);
		}	
	}
	//cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
	cout << "Redif of Divide E(Pt)" << endl;
	cout << ReCount[0][0] << "    "<< ReCount[0][1] << "    " << ReCount[0][2] << "    " << ReCount[0][3] << endl;
	cout << ReCount[1][0] << "    "<< ReCount[1][1] << "    " << ReCount[1][2] << "    " << ReCount[1][3] << endl;
	cout << ReCount[2][0] << "    "<< ReCount[2][1] << "    " << ReCount[2][2] << "    " << ReCount[2][3] << endl;
	
	TH2D* histoCountCutVariationCMEDummy = new TH2D("","",4,0,4,10000,0,1);
	histoCountCutVariationCMEDummy->SetStats(0);
	histoCountCutVariationCMEDummy->GetXaxis()->SetBinLabel(1, "E_{min} #geq 0.6 GeV/c");
	histoCountCutVariationCMEDummy->GetXaxis()->SetBinLabel(2, "E_{min} #geq 0.7 GeV/c (Standard)");
	histoCountCutVariationCMEDummy->GetXaxis()->SetBinLabel(3, "E_{min} #geq 0.8 GeV/c");
	histoCountCutVariationCMEDummy->GetXaxis()->SetBinLabel(4, "E_{min} #geq 0.9 GeV/c");
	//TH1D* histoCountCutVariationCMEStdPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationCMELowPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationCMELMPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationCMEMidPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationCMEMHPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationCMEHighPt = new TH1D("","",4,0,4);
	cout << "DEBUG: "<< __LINE__ << endl;	

	histoCountCutVariationCMELowPt->SetBinContent(1, ReCount[0][1]);
	histoCountCutVariationCMELMPt->SetBinContent(1, ReCount[1][1]);
	histoCountCutVariationCMEMidPt->SetBinContent(1, ReCount[2][1]);
	histoCountCutVariationCMEMHPt->SetBinContent(1, ReCount[3][1]);
	histoCountCutVariationCMEHighPt->SetBinContent(1, ReCount[4][1]);
	histoCountCutVariationCMELowPt->SetBinContent(2, ReCount[0][0]);
	histoCountCutVariationCMELMPt->SetBinContent(2, ReCount[1][0]);
	histoCountCutVariationCMEMidPt->SetBinContent(2, ReCount[2][0]);
	histoCountCutVariationCMEMHPt->SetBinContent(2, ReCount[3][0]);
	histoCountCutVariationCMEHighPt->SetBinContent(2, ReCount[4][0]);
	for(Int_t ivar=2; ivar<NVAR; ivar++){
	histoCountCutVariationCMELowPt->SetBinContent(ivar+1, ReCount[0][ivar]);
	histoCountCutVariationCMELMPt->SetBinContent(ivar+1, ReCount[1][ivar]);
	histoCountCutVariationCMEMidPt->SetBinContent(ivar+1, ReCount[2][ivar]);
	histoCountCutVariationCMEMHPt->SetBinContent(ivar+1, ReCount[3][ivar]);
	histoCountCutVariationCMEHighPt->SetBinContent(ivar+1, ReCount[4][ivar]);
	}

	//histoCountCutVariationCMELowPt->SetBinLabel(ivar+1);

	histoCountCutVariationCMELowPt->SetLineColor(kRed+1);
	histoCountCutVariationCMELMPt->SetLineColor(kGreen+1);
	histoCountCutVariationCMEMidPt->SetLineColor(kBlue+1);
	histoCountCutVariationCMEMHPt->SetLineColor(kYellow+1);
	histoCountCutVariationCMEHighPt->SetLineColor(kCyan+1);

	TLegend* legendCountCutVariationCME   = GetAndSetLegend2(0.12+0.45,0.74, 0.45+0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	legendCountCutVariationCME->AddEntry(histoCountCutVariationCMELowPt, "0.5 GeV/c < E <1.0 GeV/c", "l");
	legendCountCutVariationCME->AddEntry(histoCountCutVariationCMELMPt,  "1.0 GeV/c < E <1.5 GeV/c", "l");
	legendCountCutVariationCME->AddEntry(histoCountCutVariationCMEMidPt, "1.5 GeV/c < E <2.0 GeV/c", "l");
	legendCountCutVariationCME->AddEntry(histoCountCutVariationCMEMHPt,  "2.0 GeV/c < E <2.5 GeV/c", "l");
	legendCountCutVariationCME->AddEntry(histoCountCutVariationCMEHighPt,"2.5 GeV/c < E <3.0 GeV/c", "l");
	
	TCanvas *cCountCutVariationCME = new TCanvas("cCountCutVariationCME","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCountCutVariationCME, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	cCountCutVariationCME->cd();
	//gPad->SetLogy();
	histoCountCutVariationCMEDummy->Draw();	
	histoCountCutVariationCMELowPt->Draw("same");
	histoCountCutVariationCMELMPt->Draw("same");
	histoCountCutVariationCMEMidPt->Draw("same");
	histoCountCutVariationCMEMHPt->Draw("same");
	histoCountCutVariationCMEHighPt->Draw("same");
	legendCountCutVariationCME->Draw("same");
	cout << "DEBUG: "<< __LINE__ << endl;	


	// ClusterNCells

	// initialize
	for(Int_t iEPtBin =0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar =0; ivar<NVAR; ivar++){
			Count[iEPtBin][ivar] = 0;
			ReCount[iEPtBin][ivar] = 0;
		}
		TotalCount[iEPtBin] =0;
	}

	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){	
		for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
			Count[iEPtBin][0] +=gBinContent[1][0][0][ibin];
			cout << "bin: " << gBinContent[1][0][0][ibin] << endl; 
			cout << "count: " << Count[iEPtBin][0] << endl;
		}
		for(Int_t ivar=1; ivar<NVAR-1; ivar++){ // variattion4 in not defined
			cout << "ivar" << ivar <<endl;
			for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
				Count[iEPtBin][ivar] +=gBinContent[1][2][ivar][ibin];
				cout << "bin: " << gBinContent[1][2][ivar][ibin] <<endl;
				cout << "count: " << Count[iEPtBin][ivar] << endl;
			}
		}
		//TotalCount[iEPtBin] = Count[iEPtBin][0]+Count[iEPtBin][1]+Count[iEPtBin][2]+Count[iEPtBin][3]; // It is possibility to change the value of total count to all add like bin 1 to 200 (0~20GeV/c)
		for (Int_t ivar=0; ivar<NVAR-1; ivar++){
		TotalCount[iEPtBin] += Count[iEPtBin][ivar];
		}
	}
//	cout << "Divide E(Pt)" << endl;
//	cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
//	cout << Count[1][0] << "    "<< Count[1][1] << "    " << Count[1][2] << "    " << Count[1][3] << endl;
//	cout << Count[2][0] << "    "<< Count[2][1] << "    " << Count[2][2] << "    " << Count[2][3] << endl;
	// redefintion
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar=0; ivar<NVAR-1; ivar++){
			ReCount[iEPtBin][ivar] = static_cast<double>(Count[iEPtBin][ivar])/static_cast<double>(TotalCount[iEPtBin]);
		}	
	}
	//cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
//	cout << "Redif of Divide E(Pt)" << endl;
//	cout << ReCount[0][0] << "    "<< ReCount[0][1] << "    " << ReCount[0][2] << "    " << ReCount[0][3] << endl;
//	cout << ReCount[1][0] << "    "<< ReCount[1][1] << "    " << ReCount[1][2] << "    " << ReCount[1][3] << endl;
//	cout << ReCount[2][0] << "    "<< ReCount[2][1] << "    " << ReCount[2][2] << "    " << ReCount[2][3] << endl;
	
	TH2D* histoCountCutVariationNCellsDummy = new TH2D("","",3,0,3,10000,0,1);
	histoCountCutVariationNCellsDummy->SetStats(0);
	histoCountCutVariationNCellsDummy->GetXaxis()->SetBinLabel(1, "N_{Cell} #geq 1");
	histoCountCutVariationNCellsDummy->GetXaxis()->SetBinLabel(2, "N_{Cell} #geq 2 (Standard)");
	histoCountCutVariationNCellsDummy->GetXaxis()->SetBinLabel(3, "N_{Cell} #geq 3");
	TH1D* histoCountCutVariationNCellsStdPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationNCellsLowPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationNCellsLMPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationNCellsMidPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationNCellsMHPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationNCellsHighPt = new TH1D("","",4,0,4);
	cout << "DEBUG: "<< __LINE__ << endl;	

	histoCountCutVariationNCellsLowPt->SetBinContent(1, ReCount[0][1]);
	histoCountCutVariationNCellsLMPt->SetBinContent(1, ReCount[1][1]);
	histoCountCutVariationNCellsMidPt->SetBinContent(1, ReCount[2][1]);
	histoCountCutVariationNCellsMHPt->SetBinContent(1, ReCount[3][1]);
	histoCountCutVariationNCellsHighPt->SetBinContent(1, ReCount[4][1]);
	histoCountCutVariationNCellsLowPt->SetBinContent(2, ReCount[0][0]);
	histoCountCutVariationNCellsLMPt->SetBinContent(2, ReCount[1][0]);
	histoCountCutVariationNCellsMidPt->SetBinContent(2, ReCount[2][0]);
	histoCountCutVariationNCellsMHPt->SetBinContent(2, ReCount[3][0]);
	histoCountCutVariationNCellsHighPt->SetBinContent(2, ReCount[4][0]);
	for(Int_t ivar=2; ivar<NVAR-1; ivar++){
	histoCountCutVariationNCellsLowPt->SetBinContent(ivar+1, ReCount[0][ivar]);
	histoCountCutVariationNCellsLMPt->SetBinContent(ivar+1, ReCount[1][ivar]);
	histoCountCutVariationNCellsMidPt->SetBinContent(ivar+1, ReCount[2][ivar]);
	histoCountCutVariationNCellsMHPt->SetBinContent(ivar+1, ReCount[3][ivar]);
	histoCountCutVariationNCellsHighPt->SetBinContent(ivar+1, ReCount[4][ivar]);
	}
	//histoCountCutVariationNCellsLowPt->SetBinLabel(ivar+1);

	histoCountCutVariationNCellsLowPt->SetLineColor(kRed+1);
	histoCountCutVariationNCellsLMPt->SetLineColor(kGreen+1);
	histoCountCutVariationNCellsMidPt->SetLineColor(kBlue+1);
	histoCountCutVariationNCellsMHPt->SetLineColor(kYellow+1);
	histoCountCutVariationNCellsHighPt->SetLineColor(kCyan+1);

	TLegend* legendCountCutVariationNCells   = GetAndSetLegend2(0.12+0.45,0.74, 0.45+0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	legendCountCutVariationNCells->AddEntry(histoCountCutVariationNCellsLowPt, "0.5 GeV/c < E <1.0 GeV/c", "l");
	legendCountCutVariationNCells->AddEntry(histoCountCutVariationNCellsLMPt,  "1.0 GeV/c < E <1.5 GeV/c", "l");
	legendCountCutVariationNCells->AddEntry(histoCountCutVariationNCellsMidPt, "1.5 GeV/c < E <2.0 GeV/c", "l");
	legendCountCutVariationNCells->AddEntry(histoCountCutVariationNCellsMHPt,  "2.0 GeV/c < E <2.5 GeV/c", "l");
	legendCountCutVariationNCells->AddEntry(histoCountCutVariationNCellsHighPt,"2.5 GeV/c < E <3.0 GeV/c", "l");
	
	TCanvas *cCountCutVariationNCells = new TCanvas("cCountCutVariationNCells","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCountCutVariationNCells, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	cCountCutVariationNCells->cd();
	//gPad->SetLogy();
	histoCountCutVariationNCellsDummy->Draw();	
	histoCountCutVariationNCellsLowPt->Draw("same");
	histoCountCutVariationNCellsLMPt->Draw("same");
	histoCountCutVariationNCellsMidPt->Draw("same");
	histoCountCutVariationNCellsMHPt->Draw("same");
	histoCountCutVariationNCellsHighPt->Draw("same");
	legendCountCutVariationNCells->Draw("same");

	// ClusterTiming 
	
		// initialize
	for(Int_t iEPtBin =0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar =0; ivar<NVAR; ivar++){
			Count[iEPtBin][ivar] = 0;
			ReCount[iEPtBin][ivar] = 0;
		}
		TotalCount[iEPtBin] =0;
	}


	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){	
		for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
			Count[iEPtBin][0] +=gBinContent[1][0][0][ibin];
			cout << "bin: " << gBinContent[1][0][0][ibin] << endl; 
			cout << "count: " << Count[iEPtBin][0] << endl;
		}
		for(Int_t ivar=1; ivar<NVAR; ivar++){
			cout << "ivar" << ivar <<endl;
			for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
				Count[iEPtBin][ivar] +=gBinContent[1][3][ivar][ibin];
				cout << "bin: " << gBinContent[1][3][ivar][ibin] <<endl;
				cout << "count: " << Count[iEPtBin][ivar] << endl;
			}
		}
		//TotalCount[iEPtBin] = Count[iEPtBin][0]+Count[iEPtBin][1]+Count[iEPtBin][2]+Count[iEPtBin][3]; // It is possibility to change the value of total count to all add like bin 1 to 200 (0~20GeV/c)
		for (Int_t ivar=0; ivar<NVAR; ivar++){
		TotalCount[iEPtBin] += Count[iEPtBin][ivar];
		}
	}

	// redefintion
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar=0; ivar<NVAR; ivar++){
			ReCount[iEPtBin][ivar] = static_cast<double>(Count[iEPtBin][ivar])/static_cast<double>(TotalCount[iEPtBin]);
		}	
	}
	
	TH2D* histoCountCutVariationTimingDummy = new TH2D("","",4,0,4,10000,0,1);
	histoCountCutVariationTimingDummy->SetStats(0);
	histoCountCutVariationTimingDummy->GetXaxis()->SetBinLabel(1, "|t| #leq 50 ns");
	histoCountCutVariationTimingDummy->GetXaxis()->SetBinLabel(2, "|t| #leq 30 (Standard)");
	histoCountCutVariationTimingDummy->GetXaxis()->SetBinLabel(3, "-20 ns #leq t #leq 30 ns");
	histoCountCutVariationTimingDummy->GetXaxis()->SetBinLabel(4, "-20 ns #leq t #leq 25 ns");
	//TH1D* histoCountCutVariationTimingStdPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTimingLowPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTimingLMPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTimingMidPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTimingMHPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTimingHighPt = new TH1D("","",4,0,4);
	cout << "DEBUG: "<< __LINE__ << endl;	

	histoCountCutVariationTimingLowPt->SetBinContent(1, ReCount[0][1]);
	histoCountCutVariationTimingLMPt->SetBinContent(1, ReCount[1][1]);
	histoCountCutVariationTimingMidPt->SetBinContent(1, ReCount[2][1]);
	histoCountCutVariationTimingMHPt->SetBinContent(1, ReCount[3][1]);
	histoCountCutVariationTimingHighPt->SetBinContent(1, ReCount[4][1]);
	histoCountCutVariationTimingLowPt->SetBinContent(2, ReCount[0][0]);
	histoCountCutVariationTimingLMPt->SetBinContent(2, ReCount[1][0]);
	histoCountCutVariationTimingMidPt->SetBinContent(2, ReCount[2][0]);
	histoCountCutVariationTimingMHPt->SetBinContent(2, ReCount[3][0]);
	histoCountCutVariationTimingHighPt->SetBinContent(2, ReCount[4][0]);
	for(Int_t ivar=2; ivar<NVAR; ivar++){
	histoCountCutVariationTimingLowPt->SetBinContent(ivar+1, ReCount[0][ivar]);
	histoCountCutVariationTimingLMPt->SetBinContent(ivar+1, ReCount[1][ivar]);
	histoCountCutVariationTimingMidPt->SetBinContent(ivar+1, ReCount[2][ivar]);
	histoCountCutVariationTimingMHPt->SetBinContent(ivar+1, ReCount[3][ivar]);
	histoCountCutVariationTimingHighPt->SetBinContent(ivar+1, ReCount[4][ivar]);
	}

	//histoCountCutVariationTimingLowPt->SetBinLabel(ivar+1);

	histoCountCutVariationTimingLowPt->SetLineColor(kRed+1);
	histoCountCutVariationTimingLMPt->SetLineColor(kGreen+1);
	histoCountCutVariationTimingMidPt->SetLineColor(kBlue+1);
	histoCountCutVariationTimingMHPt->SetLineColor(kYellow+1);
	histoCountCutVariationTimingHighPt->SetLineColor(kCyan+1);

	TLegend* legendCountCutVariationTiming   = GetAndSetLegend2(0.12+0.45,0.74, 0.45+0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	legendCountCutVariationTiming->AddEntry(histoCountCutVariationTimingLowPt, "0.5 GeV/c < E <1.0 GeV/c", "l");
	legendCountCutVariationTiming->AddEntry(histoCountCutVariationTimingLMPt,  "1.0 GeV/c < E <1.5 GeV/c", "l");
	legendCountCutVariationTiming->AddEntry(histoCountCutVariationTimingMidPt, "1.5 GeV/c < E <2.0 GeV/c", "l");
	legendCountCutVariationTiming->AddEntry(histoCountCutVariationTimingMHPt,  "2.0 GeV/c < E <2.5 GeV/c", "l");
	legendCountCutVariationTiming->AddEntry(histoCountCutVariationTimingHighPt,"2.5 GeV/c < E <3.0 GeV/c", "l");
	
	TCanvas *cCountCutVariationTiming = new TCanvas("cCountCutVariationTiming","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCountCutVariationTiming, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	cCountCutVariationTiming->cd();
	//gPad->SetLogy();
	histoCountCutVariationTimingDummy->Draw();	
	histoCountCutVariationTimingLowPt->Draw("same");
	histoCountCutVariationTimingLMPt->Draw("same");
	histoCountCutVariationTimingMidPt->Draw("same");
	histoCountCutVariationTimingMHPt->Draw("same");
	histoCountCutVariationTimingHighPt->Draw("same");
	legendCountCutVariationTiming->Draw("same");
	

	// Cluster shape

		// initialize
	for(Int_t iEPtBin =0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar =0; ivar<NVAR; ivar++){
			Count[iEPtBin][ivar] = 0;
			ReCount[iEPtBin][ivar] = 0;
		}
		TotalCount[iEPtBin] =0;
	}

	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){	
		for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
			Count[iEPtBin][0] +=gBinContent[1][0][0][ibin];
			cout << "bin: " << gBinContent[1][0][0][ibin] << endl; 
			cout << "count: " << Count[iEPtBin][0] << endl;
		}
		for(Int_t ivar=1; ivar<NVAR-1; ivar++){ // variattion4 in not defined
			cout << "ivar" << ivar <<endl;
			for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
				Count[iEPtBin][ivar] +=gBinContent[1][4][ivar][ibin];
				cout << "bin: " << gBinContent[1][4][ivar][ibin] <<endl;
				cout << "count: " << Count[iEPtBin][ivar] << endl;
			}
		}
		//TotalCount[iEPtBin] = Count[iEPtBin][0]+Count[iEPtBin][1]+Count[iEPtBin][2]+Count[iEPtBin][3]; // It is possibility to change the value of total count to all add like bin 1 to 200 (0~20GeV/c)
		for (Int_t ivar=0; ivar<NVAR-1; ivar++){
		TotalCount[iEPtBin] += Count[iEPtBin][ivar];
		}
	}
//	cout << "Divide E(Pt)" << endl;
//	cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
//	cout << Count[1][0] << "    "<< Count[1][1] << "    " << Count[1][2] << "    " << Count[1][3] << endl;
//	cout << Count[2][0] << "    "<< Count[2][1] << "    " << Count[2][2] << "    " << Count[2][3] << endl;
	// redefintion
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar=0; ivar<NVAR-1; ivar++){
			ReCount[iEPtBin][ivar] = static_cast<double>(Count[iEPtBin][ivar])/static_cast<double>(TotalCount[iEPtBin]);
		}	
	}
	//cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
//	cout << "Redif of Divide E(Pt)" << endl;
//	cout << ReCount[0][0] << "    "<< ReCount[0][1] << "    " << ReCount[0][2] << "    " << ReCount[0][3] << endl;
//	cout << ReCount[1][0] << "    "<< ReCount[1][1] << "    " << ReCount[1][2] << "    " << ReCount[1][3] << endl;
//	cout << ReCount[2][0] << "    "<< ReCount[2][1] << "    " << ReCount[2][2] << "    " << ReCount[2][3] << endl;
	
	TH2D* histoCountCutVariationShapeDummy = new TH2D("","",3,0,3,10000,0,1);
	histoCountCutVariationShapeDummy->SetStats(0);
	histoCountCutVariationShapeDummy->GetXaxis()->SetBinLabel(1, "0.1 #leq M #leq 0.7");
	histoCountCutVariationShapeDummy->GetXaxis()->SetBinLabel(2, "0.1 #leq M #leq 0.5 (Standard)");
	histoCountCutVariationShapeDummy->GetXaxis()->SetBinLabel(3, "0.1 #leq M #leq 0.3");
	TH1D* histoCountCutVariationShapeStdPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationShapeLowPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationShapeLMPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationShapeMidPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationShapeMHPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationShapeHighPt = new TH1D("","",4,0,4);
	cout << "DEBUG: "<< __LINE__ << endl;	

	histoCountCutVariationShapeLowPt->SetBinContent(1, ReCount[0][1]);
	histoCountCutVariationShapeLMPt->SetBinContent(1, ReCount[1][1]);
	histoCountCutVariationShapeMidPt->SetBinContent(1, ReCount[2][1]);
	histoCountCutVariationShapeMHPt->SetBinContent(1, ReCount[3][1]);
	histoCountCutVariationShapeHighPt->SetBinContent(1, ReCount[4][1]);
	histoCountCutVariationShapeLowPt->SetBinContent(2, ReCount[0][0]);
	histoCountCutVariationShapeLMPt->SetBinContent(2, ReCount[1][0]);
	histoCountCutVariationShapeMidPt->SetBinContent(2, ReCount[2][0]);
	histoCountCutVariationShapeMHPt->SetBinContent(2, ReCount[3][0]);
	histoCountCutVariationShapeHighPt->SetBinContent(2, ReCount[4][0]);
	for(Int_t ivar=2; ivar<NVAR-1; ivar++){
	histoCountCutVariationShapeLowPt->SetBinContent(ivar+1, ReCount[0][ivar]);
	histoCountCutVariationShapeLMPt->SetBinContent(ivar+1, ReCount[1][ivar]);
	histoCountCutVariationShapeMidPt->SetBinContent(ivar+1, ReCount[2][ivar]);
	histoCountCutVariationShapeMHPt->SetBinContent(ivar+1, ReCount[3][ivar]);
	histoCountCutVariationShapeHighPt->SetBinContent(ivar+1, ReCount[4][ivar]);
	}
	//histoCountCutVariationShapeLowPt->SetBinLabel(ivar+1);

	histoCountCutVariationShapeLowPt->SetLineColor(kRed+1);
	histoCountCutVariationShapeLMPt->SetLineColor(kGreen+1);
	histoCountCutVariationShapeMidPt->SetLineColor(kBlue+1);
	histoCountCutVariationShapeMHPt->SetLineColor(kYellow+1);
	histoCountCutVariationShapeHighPt->SetLineColor(kCyan+1);

	TLegend* legendCountCutVariationShape   = GetAndSetLegend2(0.12+0.45,0.74, 0.45+0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	legendCountCutVariationShape->AddEntry(histoCountCutVariationShapeLowPt, "0.5 GeV/c < E <1.0 GeV/c", "l");
	legendCountCutVariationShape->AddEntry(histoCountCutVariationShapeLMPt,  "1.0 GeV/c < E <1.5 GeV/c", "l");
	legendCountCutVariationShape->AddEntry(histoCountCutVariationShapeMidPt, "1.5 GeV/c < E <2.0 GeV/c", "l");
	legendCountCutVariationShape->AddEntry(histoCountCutVariationShapeMHPt,  "2.0 GeV/c < E <2.5 GeV/c", "l");
	legendCountCutVariationShape->AddEntry(histoCountCutVariationShapeHighPt,"2.5 GeV/c < E <3.0 GeV/c", "l");
	
	TCanvas *cCountCutVariationShape = new TCanvas("cCountCutVariationShape","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCountCutVariationShape, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	cCountCutVariationShape->cd();
	//gPad->SetLogy();
	histoCountCutVariationShapeDummy->Draw();	
	histoCountCutVariationShapeLowPt->Draw("same");
	histoCountCutVariationShapeLMPt->Draw("same");
	histoCountCutVariationShapeMidPt->Draw("same");
	histoCountCutVariationShapeMHPt->Draw("same");
	histoCountCutVariationShapeHighPt->Draw("same");
	legendCountCutVariationShape->Draw("same");


	// clusterTrackMatching

		// initialize
	for(Int_t iEPtBin =0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar =0; ivar<NVAR; ivar++){
			Count[iEPtBin][ivar] = 0;
			ReCount[iEPtBin][ivar] = 0;
		}
		TotalCount[iEPtBin] =0;
	}


	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){	
		for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
			Count[iEPtBin][0] +=gBinContent[1][0][0][ibin];
			cout << "bin: " << gBinContent[1][0][0][ibin] << endl; 
			cout << "count: " << Count[iEPtBin][0] << endl;
		}
		for(Int_t ivar=1; ivar<NVAR; ivar++){
			cout << "ivar" << ivar <<endl;
			for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
				Count[iEPtBin][ivar] +=gBinContent[1][5][ivar][ibin];
				cout << "bin: " << gBinContent[1][5][ivar][ibin] <<endl;
				cout << "count: " << Count[iEPtBin][ivar] << endl;
			}
		}
		//TotalCount[iEPtBin] = Count[iEPtBin][0]+Count[iEPtBin][1]+Count[iEPtBin][2]+Count[iEPtBin][3]; // It is possibility to change the value of total count to all add like bin 1 to 200 (0~20GeV/c)
		for (Int_t ivar=0; ivar<NVAR; ivar++){
		TotalCount[iEPtBin] += Count[iEPtBin][ivar];
		}
	}

	// redefintion
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar=0; ivar<NVAR; ivar++){
			ReCount[iEPtBin][ivar] = static_cast<double>(Count[iEPtBin][ivar])/static_cast<double>(TotalCount[iEPtBin]);
		}	
	}
	
	TH2D* histoCountCutVariationTMDummy = new TH2D("","",4,0,4,10000,0,1);
	histoCountCutVariationTMDummy->SetStats(0);
	histoCountCutVariationTMDummy->GetXaxis()->SetBinLabel(1, "E/p #leq 2");
	histoCountCutVariationTMDummy->GetXaxis()->SetBinLabel(2, "Std+Secondary TM (Standard)");
	histoCountCutVariationTMDummy->GetXaxis()->SetBinLabel(3, "E/p #leq 1.75");
	histoCountCutVariationTMDummy->GetXaxis()->SetBinLabel(4, "E/p #leq 1.5");
	//TH1D* histoCountCutVariationTMStdPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTMLowPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTMLMPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTMMidPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTMMHPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationTMHighPt = new TH1D("","",4,0,4);
	cout << "DEBUG: "<< __LINE__ << endl;	

	histoCountCutVariationTMLowPt->SetBinContent(1, ReCount[0][1]);
	histoCountCutVariationTMLMPt->SetBinContent(1, ReCount[1][1]);
	histoCountCutVariationTMMidPt->SetBinContent(1, ReCount[2][1]);
	histoCountCutVariationTMMHPt->SetBinContent(1, ReCount[3][1]);
	histoCountCutVariationTMHighPt->SetBinContent(1, ReCount[4][1]);
	histoCountCutVariationTMLowPt->SetBinContent(2, ReCount[0][0]);
	histoCountCutVariationTMLMPt->SetBinContent(2, ReCount[1][0]);
	histoCountCutVariationTMMidPt->SetBinContent(2, ReCount[2][0]);
	histoCountCutVariationTMMHPt->SetBinContent(2, ReCount[3][0]);
	histoCountCutVariationTMHighPt->SetBinContent(2, ReCount[4][0]);
	for(Int_t ivar=2; ivar<NVAR; ivar++){
	histoCountCutVariationTMLowPt->SetBinContent(ivar+1, ReCount[0][ivar]);
	histoCountCutVariationTMLMPt->SetBinContent(ivar+1, ReCount[1][ivar]);
	histoCountCutVariationTMMidPt->SetBinContent(ivar+1, ReCount[2][ivar]);
	histoCountCutVariationTMMHPt->SetBinContent(ivar+1, ReCount[3][ivar]);
	histoCountCutVariationTMHighPt->SetBinContent(ivar+1, ReCount[4][ivar]);
	}

	//histoCountCutVariationTMLowPt->SetBinLabel(ivar+1);

	histoCountCutVariationTMLowPt->SetLineColor(kRed+1);
	histoCountCutVariationTMLMPt->SetLineColor(kGreen+1);
	histoCountCutVariationTMMidPt->SetLineColor(kBlue+1);
	histoCountCutVariationTMMHPt->SetLineColor(kYellow+1);
	histoCountCutVariationTMHighPt->SetLineColor(kCyan+1);

	TLegend* legendCountCutVariationTM   = GetAndSetLegend2(0.12+0.45,0.74, 0.45+0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	legendCountCutVariationTM->AddEntry(histoCountCutVariationTMLowPt, "0.5 GeV/c < E <1.0 GeV/c", "l");
	legendCountCutVariationTM->AddEntry(histoCountCutVariationTMLMPt,  "1.0 GeV/c < E <1.5 GeV/c", "l");
	legendCountCutVariationTM->AddEntry(histoCountCutVariationTMMidPt, "1.5 GeV/c < E <2.0 GeV/c", "l");
	legendCountCutVariationTM->AddEntry(histoCountCutVariationTMMHPt,  "2.0 GeV/c < E <2.5 GeV/c", "l");
	legendCountCutVariationTM->AddEntry(histoCountCutVariationTMHighPt,"2.5 GeV/c < E <3.0 GeV/c", "l");
	
	TCanvas *cCountCutVariationTM = new TCanvas("cCountCutVariationTM","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCountCutVariationTM, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	cCountCutVariationTM->cd();
	//gPad->SetLogy();
	histoCountCutVariationTMDummy->Draw();	
	histoCountCutVariationTMLowPt->Draw("same");
	histoCountCutVariationTMLMPt->Draw("same");
	histoCountCutVariationTMMidPt->Draw("same");
	histoCountCutVariationTMMHPt->Draw("same");
	histoCountCutVariationTMHighPt->Draw("same");
	legendCountCutVariationTM->Draw("same");
	


	// OpeningAngle

		// initialize
	for(Int_t iEPtBin =0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar =0; ivar<NVAR; ivar++){
			Count[iEPtBin][ivar] = 0;
			ReCount[iEPtBin][ivar] = 0;
		}
		TotalCount[iEPtBin] =0;
	}

	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){	
		for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
			Count[iEPtBin][0] +=gBinContent[1][0][0][ibin];
			cout << "bin: " << gBinContent[1][0][0][ibin] << endl; 
			cout << "count: " << Count[iEPtBin][0] << endl;
		}
		for(Int_t ivar=1; ivar<NVAR-1; ivar++){ // variattion4 in not defined
			cout << "ivar" << ivar <<endl;
			for(Int_t ibin=EPtStart+EPtWidth*iEPtBin; ibin<EPtEnd+EPtWidth*iEPtBin; ibin++){
				Count[iEPtBin][ivar] +=gBinContent[1][6][ivar][ibin];
				cout << "bin: " << gBinContent[1][6][ivar][ibin] <<endl;
				cout << "count: " << Count[iEPtBin][ivar] << endl;
			}
		}
		//TotalCount[iEPtBin] = Count[iEPtBin][0]+Count[iEPtBin][1]+Count[iEPtBin][2]+Count[iEPtBin][3]; // It is possibility to change the value of total count to all add like bin 1 to 200 (0~20GeV/c)
		for (Int_t ivar=0; ivar<NVAR-1; ivar++){
		TotalCount[iEPtBin] += Count[iEPtBin][ivar];
		}
	}
//	cout << "Divide E(Pt)" << endl;
//	cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
//	cout << Count[1][0] << "    "<< Count[1][1] << "    " << Count[1][2] << "    " << Count[1][3] << endl;
//	cout << Count[2][0] << "    "<< Count[2][1] << "    " << Count[2][2] << "    " << Count[2][3] << endl;
	// redefintion
	for(Int_t iEPtBin=0; iEPtBin<NEPtBin; iEPtBin++){
		for(Int_t ivar=0; ivar<NVAR-1; ivar++){
			ReCount[iEPtBin][ivar] = static_cast<double>(Count[iEPtBin][ivar])/static_cast<double>(TotalCount[iEPtBin]);
		}	
	}
	//cout << Count[0][0] << "    "<< Count[0][1] << "    " << Count[0][2] << "    " << Count[0][3] << endl;
//	cout << "Redif of Divide E(Pt)" << endl;
//	cout << ReCount[0][0] << "    "<< ReCount[0][1] << "    " << ReCount[0][2] << "    " << ReCount[0][3] << endl;
//	cout << ReCount[1][0] << "    "<< ReCount[1][1] << "    " << ReCount[1][2] << "    " << ReCount[1][3] << endl;
//	cout << ReCount[2][0] << "    "<< ReCount[2][1] << "    " << ReCount[2][2] << "    " << ReCount[2][3] << endl;
	
	TH2D* histoCountCutVariationOpeningAngleDummy = new TH2D("","",3,0,3,10000,0,1);
	histoCountCutVariationOpeningAngleDummy->SetStats(0);
	histoCountCutVariationOpeningAngleDummy->GetXaxis()->SetBinLabel(1, "#theta > 0.0152");
	histoCountCutVariationOpeningAngleDummy->GetXaxis()->SetBinLabel(2, "#theta > 0.017 (Standard)");
	histoCountCutVariationOpeningAngleDummy->GetXaxis()->SetBinLabel(3, "#theta > 0.0202");
	TH1D* histoCountCutVariationOpeningAngleStdPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationOpeningAngleLowPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationOpeningAngleLMPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationOpeningAngleMidPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationOpeningAngleMHPt = new TH1D("","",4,0,4);
	TH1D* histoCountCutVariationOpeningAngleHighPt = new TH1D("","",4,0,4);
	cout << "DEBUG: "<< __LINE__ << endl;	

	histoCountCutVariationOpeningAngleLowPt->SetBinContent(1, ReCount[0][1]);
	histoCountCutVariationOpeningAngleLMPt->SetBinContent(1, ReCount[1][1]);
	histoCountCutVariationOpeningAngleMidPt->SetBinContent(1, ReCount[2][1]);
	histoCountCutVariationOpeningAngleMHPt->SetBinContent(1, ReCount[3][1]);
	histoCountCutVariationOpeningAngleHighPt->SetBinContent(1, ReCount[4][1]);
	histoCountCutVariationOpeningAngleLowPt->SetBinContent(2, ReCount[0][0]);
	histoCountCutVariationOpeningAngleLMPt->SetBinContent(2, ReCount[1][0]);
	histoCountCutVariationOpeningAngleMidPt->SetBinContent(2, ReCount[2][0]);
	histoCountCutVariationOpeningAngleMHPt->SetBinContent(2, ReCount[3][0]);
	histoCountCutVariationOpeningAngleHighPt->SetBinContent(2, ReCount[4][0]);
	for(Int_t ivar=2; ivar<NVAR-1; ivar++){
	histoCountCutVariationOpeningAngleLowPt->SetBinContent(ivar+1, ReCount[0][ivar]);
	histoCountCutVariationOpeningAngleLMPt->SetBinContent(ivar+1, ReCount[1][ivar]);
	histoCountCutVariationOpeningAngleMidPt->SetBinContent(ivar+1, ReCount[2][ivar]);
	histoCountCutVariationOpeningAngleMHPt->SetBinContent(ivar+1, ReCount[3][ivar]);
	histoCountCutVariationOpeningAngleHighPt->SetBinContent(ivar+1, ReCount[4][ivar]);
	}
	//histoCountCutVariationOpeningAngleLowPt->SetBinLabel(ivar+1);

	histoCountCutVariationOpeningAngleLowPt->SetLineColor(kRed+1);
	histoCountCutVariationOpeningAngleLMPt->SetLineColor(kGreen+1);
	histoCountCutVariationOpeningAngleMidPt->SetLineColor(kBlue+1);
	histoCountCutVariationOpeningAngleMHPt->SetLineColor(kYellow+1);
	histoCountCutVariationOpeningAngleHighPt->SetLineColor(kCyan+1);

	TLegend* legendCountCutVariationOpeningAngle   = GetAndSetLegend2(0.12+0.45,0.74, 0.45+0.45, 0.74+(0.17*expectedLinesInLegend), textSizeLabelsPixel, 1, "", 43, 0);
	legendCountCutVariationOpeningAngle->AddEntry(histoCountCutVariationOpeningAngleLowPt, "0.5 GeV/c < E <1.0 GeV/c", "l");
	legendCountCutVariationOpeningAngle->AddEntry(histoCountCutVariationOpeningAngleLMPt,  "1.0 GeV/c < E <1.5 GeV/c", "l");
	legendCountCutVariationOpeningAngle->AddEntry(histoCountCutVariationOpeningAngleMidPt, "1.5 GeV/c < E <2.0 GeV/c", "l");
	legendCountCutVariationOpeningAngle->AddEntry(histoCountCutVariationOpeningAngleMHPt,  "2.0 GeV/c < E <2.5 GeV/c", "l");
	legendCountCutVariationOpeningAngle->AddEntry(histoCountCutVariationOpeningAngleHighPt,"2.5 GeV/c < E <3.0 GeV/c", "l");
	
	TCanvas *cCountCutVariationOpeningAngle = new TCanvas("cCountCutVariationOpeningAngle","",1600,800);  // gives the page size
	DrawGammaCanvasSettings(cCountCutVariationOpeningAngle, 0.1, 0.01, 0.01, 0.09); // sets margins (left, right, top, bottom)

	cCountCutVariationOpeningAngle->cd();
	//gPad->SetLogy();
	histoCountCutVariationOpeningAngleDummy->Draw();	
	histoCountCutVariationOpeningAngleLowPt->Draw("same");
	histoCountCutVariationOpeningAngleLMPt->Draw("same");
	histoCountCutVariationOpeningAngleMidPt->Draw("same");
	histoCountCutVariationOpeningAngleMHPt->Draw("same");
	histoCountCutVariationOpeningAngleHighPt->Draw("same");
	legendCountCutVariationOpeningAngle->Draw("same");


	// Alpha
	

	// The Second Plot, Ratio of variaion/standard of sytematic error vs pT
	// Adding of Fitting function like pol1 and pol2 and so on.



	// The Third Plot, all cuts of systemaitic uncertainty vs pT

	// Save the canvas
	TFile *MyAnalysisSystematicError = new TFile("MyAnalysisSystematicError.root", "Update");
	cCountCutVariationCME->Write(Form("Cluster_ClusterMinEnergy_%s", nameTrig.Data()));
	cCountCutVariationNCells->Write(Form("Cluster_NCells_%s", nameTrig.Data()));
	cCountCutVariationTiming->Write(Form("Cluster_Timing_%s", nameTrig.Data()));
	cCountCutVariationShape->Write(Form("Cluster_Shape_%s", nameTrig.Data()));
	cCountCutVariationTM->Write(Form("Cluster_TM_%s", nameTrig.Data()));
	cCountCutVariationOpeningAngle->Write(Form("Cluster_TM_%s", nameTrig.Data()));
	MyAnalysisSystematicError->Write();
	MyAnalysisSystematicError->Close();
	
	cCountCutVariationCME->SaveAs(Form("Cluster_ClusterMinEnergy_%s.eps", nameTrig.Data()));
	cCountCutVariationNCells->SaveAs(Form("Cluster_NCells_%s.eps", nameTrig.Data()));
	cCountCutVariationTiming->SaveAs(Form("Cluster_Timing_%s.eps", nameTrig.Data()));
	cCountCutVariationShape->SaveAs(Form("Cluster_Shape_%s.eps", nameTrig.Data()));
	cCountCutVariationTM->SaveAs(Form("Cluster_TM_%s.eps", nameTrig.Data()));
	cCountCutVariationOpeningAngle->SaveAs(Form("Cluster_OpeningAngle_%s.eps", nameTrig.Data()));
	
}

