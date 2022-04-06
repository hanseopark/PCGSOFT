///////////////////////////////////////////////////////
/////// Unfolding macro //////////////////////
/////////////////////////////////////////////////////////

//#include "CommonHeaders/PlottingGammaConversionHistos.h"
//#include "CommonHeaders/PlottingGammaConversionAdditional.h"
//#include "CommonHeaders/FittingGammaConversion.h"
//#include "CommonHeaders/ConversionFunctions.h"


void MyUnfoldDefault(){

	TF1 *Pi0 = new TF1("pi0", "[0]+[1]*x", 2, 25);
	Pi0->SetParameters(1,0);
	TF1 *Eta = new TF1("Eta", "[0]+[1]*x", 2, 25);
	Eta->SetParameters(1,0);

	TFile *MyAnalysis = new TFile("Jet_Unfolding_Corrections_13TeV_4_Default.root","recreate");
	Pi0->Write("FinalFit_Pi0");
	Eta->Write("FinalFit_Eta");
	MyAnalysis->Close();
}
