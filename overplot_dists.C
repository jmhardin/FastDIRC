void overplot_dists(TString ifile)
{
	TFile *f1 = new TFile(ifile);
	TH2F *xypion = (TH2F*) f1->Get("pion_dist_xy");
	TH2F *xykaon = (TH2F*) f1->Get("kaon_dist_xy");

	xypion->SetTitle(ifile);

	int redc[3];
	int bluec[3];
	int redblue[2];
	redblue[0] = 2;
	redblue[1] = 9;

	redc[0] = 45;
	redc[1] = 46;
	redc[2] = 2;

	bluec[0] = 33;
	bluec[1] = 38;
	bluec[2] = 9;


	gStyle->SetPalette(3,bluec);
	
	xypion->Draw("COL");

	c1->SetWindowSize(1000,800);
	c1->Print("pion_col_xy.gif");

	
	gStyle->SetPalette(3,bluec);
	xykaon->Draw("col");
	c1->SetWindowSize(1000,800);
	c1->Print("kaon_col_xy.gif");
}
