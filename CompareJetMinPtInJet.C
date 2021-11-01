
void CompareJetMinPtInJet(
    TString pwd5GeV = "211019_0_pp_13TeV_EMCal_5GeV",
    TString pwd10GeV = "211020_0_pp_13TeV_EMCal_10GeV",
    TString pwd10GeVCorr = "211022_0_pp_13TeV_EMCal_10GeVCorr",
    TString pwdMB = "211025_0_pp_13TeV_EMCal_InclusiveJet",
	TString Energy = "13TeV"
    ){


    const Int_t NMeson = 3;
    const char *Meson[NMeson] = {"Pi0", "Eta", "Pi0EtaBinning"};
    
    TString MesonSelection[NMeson];

    const Int_t NDataSample = 2;
	const Int_t colorData[NDataSample] = {1,8};
	const Int_t markerStyleData[NDataSample] = {20,24};
    const char *DataSample[NDataSample] = {"data", "MC"};
    TString DataSampleSelection[NMeson];

    //const Int_t NEnergy = 4;
	//const Int_t colorEnergy[NEnergy] = {4,3,2,1};
    //const char *energy[NEnergy] = {"5 GeV", "10 GeV", "20 GeV", "10 GeV Cor"};
    const Int_t NEnergy = 4;
	const Int_t colorEnergy[NEnergy] = {4,2,8,1};
    const char *energy[NEnergy] = {"5 GeV", "10 GeV", "10 GeV Cor", "Inclusive jet"};

    TString EnergySelection[NEnergy];
    EnergySelection[0] = pwd5GeV;
    EnergySelection[1] = pwd10GeV;
    //EnergySelection[2] = pwd20GeV;
    EnergySelection[2] = pwd10GeVCorr;
    EnergySelection[3] = pwdMB;

	TString EventCuts = "00010113";
	TString EventCuts_EG1 = "0008d113";
	TString EventCuts_EG2 = "0008e113";
	TString CaloCuts  = "411790607l032230000";
	TString MesonCuts = "2l631031000000d0";
	TString MesonCuts_MB = "0l631031000000d0";

	const Int_t NTRIG = 3;
	const Int_t colorTrig[NTRIG] = {4, 2, 1};
	const char *Trig[NTRIG] = {"INT7", "EG1", "EG2"};

	TString TrigSelection[NTRIG];
	TrigSelection[0] = Form("%s_%s_%s",EventCuts.Data(),CaloCuts.Data(),MesonCuts.Data());
	TrigSelection[1] = Form("%s_%s_%s",EventCuts_EG1.Data(),CaloCuts.Data(),MesonCuts.Data());
	TrigSelection[2] = Form("%s_%s_%s",EventCuts_EG2.Data(),CaloCuts.Data(),MesonCuts.Data());

	TString MBSelection = Form("%s_%s_%s",EventCuts.Data(),CaloCuts.Data(),MesonCuts_MB.Data());

    // Reading corrected root file
    TFile *CorrData[NMeson][NDataSample][NEnergy][NTRIG];
    TH1D *histoYieldMesonPerEvent[NMeson][NDataSample][NEnergy][NTRIG];
    TH1D *histoCorrectedYield[NMeson][NDataSample][NEnergy][NTRIG];
    for(Int_t iMeson = 0; iMeson < NMeson; iMeson++){
        MesonSelection[iMeson] = Meson[iMeson];
        for(Int_t iDataSample = 0; iDataSample < NDataSample; iDataSample++){
            DataSampleSelection[iDataSample] = DataSample[iDataSample];
            for(Int_t iEnergy = 0; iEnergy < NEnergy; iEnergy++){
				if (iEnergy == 3){
					CorrData[iMeson][iDataSample][iEnergy][0] = new TFile(Form("%s/%s/%s/%s_%s_GammaConvV1Correction_%s.root",EnergySelection[iEnergy].Data(), MBSelection.Data(), Energy.Data(), MesonSelection[iMeson].Data(), DataSampleSelection[iDataSample].Data(), MBSelection.Data()), "read");
					histoYieldMesonPerEvent[iMeson][iDataSample][iEnergy][0] = (TH1D*) CorrData[iMeson][iDataSample][iEnergy][0]->Get("histoYieldMesonPerEvent");
					histoCorrectedYield[iMeson][iDataSample][iEnergy][0] = (TH1D*) CorrData[iMeson][iDataSample][iEnergy][0]->Get("CorrectedYieldNormEff");
				} else {
					for (Int_t itrig = 0; itrig < NTRIG; itrig++){
						CorrData[iMeson][iDataSample][iEnergy][itrig] = new TFile(Form("%s/%s/%s/%s_%s_GammaConvV1Correction_%s.root",EnergySelection[iEnergy].Data(), TrigSelection[itrig].Data(), Energy.Data(), MesonSelection[iMeson].Data(), DataSampleSelection[iDataSample].Data(), TrigSelection[itrig].Data()), "read");
						histoYieldMesonPerEvent[iMeson][iDataSample][iEnergy][itrig] = (TH1D*) CorrData[iMeson][iDataSample][iEnergy][itrig]->Get("histoYieldMesonPerEvent");
						histoCorrectedYield[iMeson][iDataSample][iEnergy][itrig] = (TH1D*) CorrData[iMeson][iDataSample][iEnergy][itrig]->Get("CorrectedYieldNormEff");
					}
				}
                        
            }
        }
    }

	TH1D *histoRatioPi0Eta[NDataSample][NEnergy][NTRIG];
	for(Int_t iDataSample = 0; iDataSample < NDataSample; iDataSample++){
		for(Int_t iEnergy = 0; iEnergy < NEnergy; iEnergy++){
			if (iEnergy == 3){
				histoRatioPi0Eta[iDataSample][iEnergy][0] = (TH1D*) histoCorrectedYield[1][iDataSample][iEnergy][0]->Clone(); 
				histoRatioPi0Eta[iDataSample][iEnergy][0]->Divide(histoRatioPi0Eta[iDataSample][iEnergy][0], histoCorrectedYield[2][iDataSample][iEnergy][0],1.,1.,"");
			} else {
				for (Int_t itrig = 0; itrig < NTRIG; itrig++){
					histoRatioPi0Eta[iDataSample][iEnergy][itrig] = (TH1D*) histoCorrectedYield[1][iDataSample][iEnergy][itrig]->Clone(); 
					histoRatioPi0Eta[iDataSample][iEnergy][itrig]->Divide(histoRatioPi0Eta[iDataSample][iEnergy][itrig], histoCorrectedYield[2][iDataSample][iEnergy][itrig],1.,1.,"");
				}
			}

		}
	}

	TH1D *histoRatioInclusiveToInJet[NMeson][NDataSample][NEnergy];
	TH1D *histoHardFactor[NMeson][NDataSample];
	for(Int_t iMeson = 0; iMeson<NMeson-1; iMeson++){
		for(Int_t iDataSample = 0; iDataSample < NDataSample; iDataSample++){
			for(Int_t iEnergy = 0; iEnergy<NEnergy-1; iEnergy++){
				histoRatioInclusiveToInJet[iMeson][iDataSample][iEnergy] = (TH1D*) histoCorrectedYield[iMeson][iDataSample][iEnergy][0]->Clone();
				histoRatioInclusiveToInJet[iMeson][iDataSample][iEnergy]->Divide(histoRatioInclusiveToInJet[iMeson][iDataSample][iEnergy], histoCorrectedYield[iMeson][iDataSample][3][0],1.,1.,"");
			}
			histoHardFactor[iMeson][iDataSample] = (TH1D*) histoRatioInclusiveToInJet[iMeson][iDataSample][0]->Clone();
			histoHardFactor[iMeson][iDataSample]->Divide(histoHardFactor[iMeson][iDataSample],histoRatioInclusiveToInJet[iMeson][iDataSample][1],1.,1.,"");
		}
	}

	auto legEnergy = new TLegend(0.7,0.7,0.9,0.9);

   	TH2D* hrefPi0 = new TH2D("hrefPi0", "",100,0,20, 100, 10e-10, 0.1*10e-2);
	hrefPi0->SetStats(0);
	hrefPi0->SetXTitle("#bf{#it{p}_{T}}");
	
   	TH2D* hrefEta = new TH2D("hrefEta", "",100,0,20, 100, 10e-10, 0.1*10e-3);
	hrefEta->SetStats(0);
	hrefEta->SetXTitle("#bf{#it{p}_{T}}");

   	TH2D* hrefPi0Eta = new TH2D("hrefPi0Eta", "",100,0,20, 100, 0, 1);
	hrefPi0Eta->SetStats(0);
	hrefPi0Eta->SetXTitle("#bf{#it{p}_{T}}");

   	TH2D* hrefRatioPi0Injet = new TH2D("hrefRatioPi0Injet", "",100,0,20, 100, 0, 1);
	hrefRatioPi0Injet->SetStats(0);
	hrefRatioPi0Injet->SetXTitle("#bf{#it{p}_{T}}");
	hrefRatioPi0Injet->SetYTitle("#frac{#pi_{0, in jet}}{#pi_{0, Inclusive}}");

	TH2D* hrefRatioEtaInjet = new TH2D("hrefRatioEtaInjet", "",100,0,20, 100, 0, 1);
	hrefRatioEtaInjet->SetStats(0);
	hrefRatioEtaInjet->SetXTitle("#bf{#it{p}_{T}}");
	hrefRatioEtaInjet->SetYTitle("#frac{#eta_{0, in jet}}{#eta_{0, Inclusive}}");

   	TH2D* hrefHardFactor = new TH2D("hrefHardFactor", "",100,0,20, 100, 0, 10);
	hrefHardFactor->SetStats(0);
	hrefHardFactor->SetXTitle("#bf{#it{p}_{T}}");
	hrefHardFactor->SetYTitle("#frac{InJet (5 GeV)}{InJet (10 GeV)}");

	TCanvas* ctest = new TCanvas("ctest","",200,10,1800,1200);  // gives the page size
	ctest->Divide(3,2);
	Int_t i = 0;
	Int_t j = 0;
	Int_t k = 0;
	while(kTRUE) {
	legEnergy->Clear();
	i++;
	ctest->cd(i); // Combined Pi0, Eta, EtaToPi0
	ctest->SetLogy();

	if (j==0) hrefPi0->Draw();
	else if (j==1) hrefEta->Draw();
	else hrefPi0Eta->Draw();

	if (j==2) {
		for(Int_t iEnergy = 0; iEnergy < NEnergy; iEnergy++){
			for (Int_t itrig = 0; itrig < NTRIG-2; itrig++){
				histoRatioPi0Eta[k][iEnergy][itrig]->SetLineWidth(1);
				histoRatioPi0Eta[k][iEnergy][itrig]->SetLineColor(colorEnergy[iEnergy]);
				histoRatioPi0Eta[k][iEnergy][itrig]->SetMarkerColor(colorEnergy[iEnergy]);
				histoRatioPi0Eta[k][iEnergy][itrig]->SetMarkerStyle(24);
				histoRatioPi0Eta[k][iEnergy][itrig]->SetMarkerSize(1);

				legEnergy ->AddEntry(histoRatioPi0Eta[k][iEnergy][itrig],energy[iEnergy]);
				histoRatioPi0Eta[k][iEnergy][itrig]->Draw("SAME P");
				legEnergy->Draw("SAME P");
			}        
		}
	} else {
		for(Int_t iEnergy = 0; iEnergy < NEnergy; iEnergy++){
			for (Int_t itrig = 0; itrig < NTRIG-2; itrig++){
				histoCorrectedYield[j][k][iEnergy][itrig]->SetLineWidth(1);
				histoCorrectedYield[j][k][iEnergy][itrig]->SetLineColor(colorEnergy[iEnergy]);
				histoCorrectedYield[j][k][iEnergy][itrig]->SetMarkerColor(colorEnergy[iEnergy]);
				histoCorrectedYield[j][k][iEnergy][itrig]->SetMarkerStyle(24);
				histoCorrectedYield[j][k][iEnergy][itrig]->SetMarkerSize(1);

				legEnergy ->AddEntry(histoCorrectedYield[j][k][iEnergy][itrig],energy[iEnergy]);
				histoCorrectedYield[j][k][iEnergy][itrig]->Draw("SAME P");
				legEnergy->Draw("SAME P");
			}        
		}
	}
	j++;
	if (j>2){
	j=0;
	k++;
	}
	if (k>1) break;
	}
	ctest->SaveAs("Result/CombinedMesonCorrectedYieldJetMinEnergy.eps");

	TCanvas* ctest2 = new TCanvas("ctest2","",200,10,1800,1200);  // gives the page size
	auto legData = new TLegend(0.7,0.7,0.9,0.9);
	ctest2->Divide(2,1); // Comparing in jet analysis 5 GeV with 10 GeV
	
	for (Int_t iMeson=0; iMeson<NMeson-1; iMeson++){
		legData->Clear();
		ctest2->cd(iMeson+1);
		hrefHardFactor->Draw();
		for (Int_t iDataSample=0; iDataSample<NDataSample; iDataSample++){
			legData->AddEntry(histoHardFactor[iMeson][iDataSample],DataSample[iDataSample]);
			histoHardFactor[iMeson][iDataSample]->Draw("SAME P");
			legData->Draw("SAME P");
		}
	}
	ctest2->SaveAs("Result/CombinedMesonRatioJetMinEnergy.eps");

	
	TCanvas* ctest3 = new TCanvas("ctest3","",200,10,1800,1200);  // gives the page size
	auto legPi0Injet = new TLegend(0.8,0.8,0.9,0.9);
	ctest3->cd(); // about Pi0 meson Comparing inclusive jet with in jet analysis
	hrefRatioPi0Injet->Draw();
	for (Int_t iDataSample =0; iDataSample<NDataSample; iDataSample++){
		histoRatioInclusiveToInJet[0][iDataSample][0]->SetLineWidth(1);
		histoRatioInclusiveToInJet[0][iDataSample][0]->SetLineColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[0][iDataSample][0]->SetMarkerColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[0][iDataSample][0]->SetMarkerStyle(markerStyleData[iDataSample]);
		histoRatioInclusiveToInJet[0][iDataSample][0]->SetMarkerSize(1);

		histoRatioInclusiveToInJet[0][iDataSample][1]->SetLineWidth(1);
		histoRatioInclusiveToInJet[0][iDataSample][1]->SetLineColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[0][iDataSample][1]->SetMarkerColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[0][iDataSample][1]->SetMarkerStyle(markerStyleData[iDataSample]+1);
		histoRatioInclusiveToInJet[0][iDataSample][1]->SetMarkerSize(1);

	
		histoRatioInclusiveToInJet[0][iDataSample][0]->Draw("SAME P");
		histoRatioInclusiveToInJet[0][iDataSample][1]->Draw("SAME P");

	}
	legPi0Injet->AddEntry(histoRatioInclusiveToInJet[0][0][0],energy[0]);
	legPi0Injet->AddEntry(histoRatioInclusiveToInJet[0][0][1],energy[1]);
	legPi0Injet->Draw("SAME P");
	ctest3->SaveAs("Result/CombinedPi0RatioInclusiveToInJet.eps");

	TCanvas* ctest4 = new TCanvas("ctest4","",200,10,1800,1200);  // gives the page size
	auto legEtaInjet = new TLegend(0.8,0.8,0.9,0.9);
	ctest4->cd(); // about Eta meson Comparing inclusive jet with in jet analysis
	hrefRatioEtaInjet->Draw();
	for (Int_t iDataSample =0; iDataSample<NDataSample; iDataSample++){
		histoRatioInclusiveToInJet[1][iDataSample][0]->SetLineWidth(1);
		histoRatioInclusiveToInJet[1][iDataSample][0]->SetLineColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[1][iDataSample][0]->SetMarkerColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[1][iDataSample][0]->SetMarkerStyle(markerStyleData[iDataSample]);
		histoRatioInclusiveToInJet[1][iDataSample][0]->SetMarkerSize(1);

		histoRatioInclusiveToInJet[1][iDataSample][1]->SetLineWidth(1);
		histoRatioInclusiveToInJet[1][iDataSample][1]->SetLineColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[1][iDataSample][1]->SetMarkerColor(colorData[iDataSample]);
		histoRatioInclusiveToInJet[1][iDataSample][1]->SetMarkerStyle(markerStyleData[iDataSample]+1);
		histoRatioInclusiveToInJet[1][iDataSample][1]->SetMarkerSize(1);

	
		histoRatioInclusiveToInJet[1][iDataSample][0]->Draw("SAME P");
		histoRatioInclusiveToInJet[1][iDataSample][1]->Draw("SAME P");

	}
	legEtaInjet->AddEntry(histoRatioInclusiveToInJet[1][0][0],energy[0]);
	legEtaInjet->AddEntry(histoRatioInclusiveToInJet[1][0][1],energy[1]);
	legEtaInjet->Draw("SAME P");
	ctest3->SaveAs("Result/CombinedEtaRatioInclusiveToInJet.eps");

	TCanvas* cref[100];

	// Pi0
	cref[0] = new TCanvas("","",200,10,1800,1200);
	cref[0]->cd(); // Multi Energy, Pi0 Data INT7
	legEnergy->Clear();
	cref[0]->SetLogy();
	hrefPi0->Draw();
	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
		histoCorrectedYield[0][0][iEnergy][0]->Draw("SAME P");
		legEnergy->AddEntry(histoCorrectedYield[0][0][iEnergy][0], energy[iEnergy]);
		legEnergy->Draw("SAME P");
	}
	cref[0]->SaveAs("Result/Pi0_Data_CorrectedYield_MultiEnergy.eps");

	cref[1] = new TCanvas("","",200,10,1800,1200);
	cref[1]->cd(); // Multi Energy, Eta Data INT7
	legEnergy->Clear();
	cref[1]->SetLogy();
	hrefEta->Draw();
	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
		histoCorrectedYield[1][0][iEnergy][0]->Draw("SAME P");
		legEnergy->AddEntry(histoCorrectedYield[1][0][iEnergy][0], energy[iEnergy]);
		legEnergy->Draw("SAME P");
	}
	cref[1]->SaveAs("Result/Eta_Data_CorrectedYield_MultiEnergy.eps");

	cref[2] = new TCanvas("","",200,10,1800,1200);
	cref[2]->cd(); // Multi Energy, Pi0 MC INT7
	legEnergy->Clear();
	cref[2]->SetLogy();
	hrefPi0->Draw();
	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
		histoCorrectedYield[0][1][iEnergy][0]->Draw("SAME P");
		legEnergy->AddEntry(histoCorrectedYield[0][1][iEnergy][0], energy[iEnergy]);
		legEnergy->Draw("SAME P");
	}
	cref[2]->SaveAs("Result/Pi0_MC_CorrectedYield_MultiEnergy.eps");

	cref[3] = new TCanvas("","",200,10,1800,1200);
	cref[3]->cd(); // Multi Energy, Eta MC INT7
	legEnergy->Clear();
	cref[3]->SetLogy();
	hrefEta->Draw();
	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
		histoCorrectedYield[1][1][iEnergy][0]->Draw("SAME P");
		legEnergy->AddEntry(histoCorrectedYield[1][1][iEnergy][0], energy[iEnergy]);
		legEnergy->Draw("SAME P");
	}
	cref[3]->SaveAs("Result/Eta_MC_CorrectedYield_MultiEnergy.eps");

	cref[4] = new TCanvas("","",200,10,1800,1200);
	cref[4]->cd(); // Multi Energy, EtaToPi0 Data INT7
	legEnergy->Clear();
	hrefPi0Eta->Draw();
	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
		histoRatioPi0Eta[0][iEnergy][0]->Draw("SAME P");
		legEnergy->AddEntry(histoRatioPi0Eta[0][iEnergy][0], energy[iEnergy]);
		legEnergy->Draw("SAME P");
	}
	cref[4]->SaveAs("Result/EtaToPi0_Data_MultiEnergy.eps");

	cref[5] = new TCanvas("","",200,10,1800,1200);
	cref[5]->cd(); // Multi Energy, EtaToPi0 Data INT7
	legEnergy->Clear();
	hrefPi0Eta->Draw();
	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
		histoRatioPi0Eta[1][iEnergy][0]->Draw("SAME P");
		legEnergy->AddEntry(histoRatioPi0Eta[1][iEnergy][0], energy[iEnergy]);
		legEnergy->Draw("SAME P");
	}
	cref[5]->SaveAs("Result/EtaToPi0_MC_MultiEnergy.eps");


//	cref[4] = new TCanvas("","",200,10,1800,1200);
//	cref[4]->cd(); // Multi Energy, Pi0 MC INT7
//	legEnergy->Clear();
//	cref[4]->SetLogy();
//	//hrefPi0->Draw();
//	for(Int_t iEnergy =0; iEnergy<NEnergy; iEnergy++){
//		histoCorrectedYield[0][0][iEnergy][1]->Draw("SAME P");
//		legEnergy->AddEntry(histoCorrectedYield[0][0][iEnergy][1], energy[iEnergy]);
//		legEnergy->Draw("SAME P");
//	}
//	cref[4]->SaveAs("Result/Pi0_Data_EG1_CorrectedYield_MultiEnergy.eps");











}


















