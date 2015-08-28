void outputDistro(TString ifile_pref)
{
	double xmin = -1500;
	double xmax = 1500;
	double ymin = -150;
	double ymax = 350;	

	double tmin = 0;
	double tmax = 250;

	TFile *f1 = new TFile(ifile_pref + ".root");
	TH2F *pi_dist = (TH2F*) f1->Get("pion_dist_xy");
	TH2F *pi_dist_yt = (TH2F*) f1->Get("pion_dist_yt");
	TH2F *pi_dist_xt = (TH2F*) f1->Get("pion_dist_xt");
	TH1F *pi_dist_t = (TH1F*) f1->Get("pion_dist_t");
	TH1F *pion_lambda = (TH1F*) f1->Get("pion_lambda");
	TH2F *pi_geant_dist = (TH2F*) f1->Get("pion_dist_geant_xy");
	TH2F *k_dist = (TH2F*) f1->Get("kaon_dist_xy");

	pi_dist->GetXaxis()->SetRangeUser(xmin,xmax);
	pi_dist->GetYaxis()->SetRangeUser(ymin,ymax);
	pi_geant_dist->GetXaxis()->SetRangeUser(xmin,xmax);
	pi_geant_dist->GetYaxis()->SetRangeUser(ymin,ymax);
	pi_dist_yt->GetXaxis()->SetRangeUser(ymin,ymax);
	pi_dist_yt->GetXaxis()->SetTitle("Y pos (mm)");
	pi_dist_yt->GetYaxis()->SetRangeUser(tmin,tmax);
	pi_dist_yt->GetYaxis()->SetTitle("Time (ns)");
	pi_dist_xt->GetXaxis()->SetRangeUser(xmin,xmax);
	pi_dist_xt->GetXaxis()->SetTitle("X pos (mm)");
	pi_dist_xt->GetYaxis()->SetRangeUser(tmin,tmax);
	pi_dist_xt->GetYaxis()->SetTitle("Time (ns)");
	pi_dist_t->GetXaxis()->SetRangeUser(tmin,tmax);
	pi_dist_t->GetXaxis()->SetTitle("Time (ns)");
	pi_dist->SetStats(false);
	pi_geant_dist->SetStats(false);
	pi_dist_yt->SetStats(false);
	pi_dist_xt->SetStats(false);
	pi_dist_t->SetStats(false);
	pion_lambda->SetStats(false);

	
	pi_dist->SetTitle("");
	pi_dist_yt->SetTitle("");
	pi_dist_xt->SetTitle("");
	pi_dist_t->SetTitle("");
	pi_geant_dist->SetTitle("");
	pion_lambda->SetTitle("");

	pion_lambda->GetXaxis()->SetTitle("Wavelength (nm)");

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

	c1->SetWindowSize(3600,2400);

	c1->Print(ifile_pref + "_pion_dist.gif");	


   	pi_geant_dist->Draw("colz");
	c1->Print(ifile_pref + "_pion_geant_dist.gif");	
   	
	pi_dist_yt->Draw("colz");
	c1->Print(ifile_pref + "_pion_dist_yt.gif");	
	
	pi_dist_xt->Draw("colz");
	c1->Print(ifile_pref + "_pion_dist_xt.gif");	

	pi_dist_t->Draw();
	c1->Print(ifile_pref + "_pion_dist_t.gif");	

	pion_lambda->SetTitle("");
	pion_lambda->GetXaxis()->SetTitle("\\lambda (nm)");
	pion_lambda->SetStats(false);
	
	pion_lambda->Draw();
	c1->Print(ifile_pref + "_pion_lambda.pdf");	


}
