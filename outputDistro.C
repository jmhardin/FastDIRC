void outputDistro(TString ifile_pref)
{
	double xmin = -1500;
	double xmax = -xmin;
	double ymin = -350;
	double ymax = -ymin;	

	TFile *f1 = new TFile(ifile_pref + ".root");
	TH2F *pi_dist = (TH2F*) f1->Get("pion_dist_xy");
	TH2F *k_dist = (TH2F*) f1->Get("kaon_dist_xy");

	pi_dist->GetXaxis()->SetRangeUser(xmin,xmax);
	pi_dist->GetYaxis()->SetRangeUser(ymin,ymax);
	pi_dist->SetStats(false);

   	const UInt_t Number = 3;
   	//Double_t Red[Number]    = { 0.016, 0.016, 0.336};
  	//Double_t Green[Number]  = { 0.157, 0.157, 0.625};
   	//Double_t Blue[Number]   = { 0.218, 0.218, 0.824};
   	Double_t Red[Number]    = { 0.743, 0.336, 0.016};
  	Double_t Green[Number]  = { 0.854, 0.625, 0.157};
   	Double_t Blue[Number]   = { 0.930, 0.824, 0.218};
   	Double_t Length[Number] = { 0.00, 0.30, 1.00 };
   	Int_t nb=50;
   	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   	pi_dist->Draw("colz");

	c1->Print(ifile_pref + "_pion_dist.gif");	

}
