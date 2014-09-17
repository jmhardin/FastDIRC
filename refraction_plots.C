void refraction_plots(TString ifile)
{
	double peak_zoom_before_a = 15;
	double peak_zoom_before_b = 20;
	double peak_zoom_after_a = 15;
	double peak_zoom_after_b = 22;
	TFile *f = new TFile(ifile);

	TH1F *pion_before = (TH1F*) f->Get("ref_pion_before");
	TH1F *kaon_before = (TH1F*) f->Get("ref_kaon_before");
	TH1F *pion_after = (TH1F*) f->Get("ref_pion_after");
	TH1F *kaon_after = (TH1F*) f->Get("ref_kaon_after");
	
	pion_before->SetStats(false);

	pion_before->SetAxisRange(10,60);
	kaon_before->SetAxisRange(10,60);

	pion_before->SetTitle("Pion (blue) versus Kaon (red) distributions before interface");

	pion_before->SetLineColor(kBlue);
	kaon_before->SetLineColor(kRed);

	pion_before->Draw();
	kaon_before->Draw("SAME");

	c1->SetWindowSize(1000,800);

	c1->Print("Refraction_before.gif");

	pion_before->SetAxisRange(peak_zoom_before_a,peak_zoom_before_b);
	kaon_before->SetAxisRange(peak_zoom_before_a,peak_zoom_before_b);

	pion_before->Draw();
	kaon_before->Draw("SAME");

	c1->SetWindowSize(1000,800);

	c1->Print("Refraction_before_zoom.gif");

	pion_after->SetStats(false);

	pion_after->SetAxisRange(10,60);
	kaon_after->SetAxisRange(10,60);

	pion_after->SetTitle("Pion (blue) versus Kaon (red) distributions after interface");

	pion_after->SetLineColor(kBlue);
	kaon_after->SetLineColor(kRed);

	pion_after->Draw();
	kaon_after->Draw("SAME");

	c1->SetWindowSize(1000,800);

	c1->Print("Refraction_after.gif");

	pion_after->SetAxisRange(peak_zoom_after_a,peak_zoom_after_b);
	kaon_after->SetAxisRange(peak_zoom_after_a,peak_zoom_after_b);

	pion_after->Draw();
	kaon_after->Draw("SAME");

	c1->SetWindowSize(1000,800);

	c1->Print("Refraction_after_zoom.gif");
}

