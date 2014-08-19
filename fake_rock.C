{
double spread = 3.8;
double seperation = 5.23;
double mean_close = 0;
double mean_far = mean_close + seperation;

//denominator pf spread
double spreadsq2 = 2*spread*spread;

TRandom3* randgen = new TRandom3();


double hmin = -40;
double hmax = 40;

TFile *f1 = new TFile("./fitdirc45.root");
TH1F *fhclose = (TH1F*) f1->Get("ll_diff_close");
TH1F *fhfar = (TH1F*) f1->Get("ll_diff_far");

fhclose->Reset();
fhfar->Reset();

double close_obs, far_obs;
double close_ll_diff, far_ll_diff;

for (int ii = 0; ii < 100000; ii++)
{
	close_obs = randgen->Gaus(mean_close,spread);
	far_obs = randgen->Gaus(mean_far,spread);
	
	close_ll_diff = -1*(close_obs - mean_close)*(close_obs - mean_close);
	close_ll_diff += (close_obs - mean_far)*(close_obs - mean_far);
	close_ll_diff /= spreadsq2;

	far_ll_diff = - (far_obs - mean_close)*(far_obs - mean_close);
	far_ll_diff += (far_obs - mean_far)*(far_obs - mean_far);
	far_ll_diff /= spreadsq2;
	
	fhclose->Fill(close_ll_diff);
	fhfar->Fill(far_ll_diff);
}


fhclose->SetAxisRange(hmin,hmax);
fhfar->SetAxisRange(hmin,hmax);

fhclose->SetLineColor(kRed);
//hclose->SetFillColorAlpha(kRed,.5);

fhfar->SetLineColor(kBlue);
//hfar->SetFillColorAlpha(kBlue,.5);

fhfar->SetTitle("Fake Seperation distance = 5.23 mrad.  spread = 2.5mrad");



TH1F *fpion_veto_eff = new TH1F(*hclose);
TH1F *fkaon_missid = new TH1F(*hfar);

fpion_veto_eff->SetName("pion_veto_eff");
fpion_veto_eff->SetTitle("Pion Veto Efficiency");

fkaon_missid->SetName("kaon_missid");
fkaon_missid->SetTitle("Kaon Miss ID");

for (int i = 0; i < fpion_veto_eff->GetNbinsX(); i++)
{
        fpion_veto_eff->SetBinContent(i,hclose->Integral(i,pion_veto_eff->GetNbinsX()));
        fkaon_missid->SetBinContent(i,hfar->Integral(0,i));
}

fpion_veto_eff->SetAxisRange(0,10000,"Y");

double scale_int = 1/fhclose->Integral(0,fpion_veto_eff->GetNbinsX());
fpion_veto_eff->Scale(scale_int);
scale_int = 1/fhfar->Integral(0,fkaon_missid->GetNbinsX());
fkaon_missid->Scale(fscale_int);

fhfar->Draw();
fhclose->Draw("SAME H");
c1->SetWindowSize(1000,800);

c1->Print("overlap_fake.gif");


fpion_veto_eff->Draw("");
fkaon_missid->Draw("SAME H");

c1->Print("overlap_integral_fake.gif");


TGraph* frock_graph;
int frock_n = fpion_veto_eff->GetNbinsX();
double fxr[rock_n];
double fyr[rock_n];

double fival = 0;
double flast_x = fpion_veto_eff->GetBinContent(0);
for (int i = 0; i < fpion_veto_eff->GetNbinsX(); i++)
{
        fxr[i] = fpion_veto_eff->GetBinContent(i);
        fyr[i] = fkaon_missid->GetBinContent(i);

	fival -= fyr[i]*(fxr[i] - flast_x);
        flast_x = fxr[i];
}

printf("Fake ROC integral: %12.04f\n",fival);

frock_graph = new TGraph(frock_n,fxr,fyr);
frock_graph->SetLineColor(2);
frock_graph->SetLineWidth(4);
frock_graph->SetTitle("Fake ROC graph");
frock_graph->GetXaxis()->SetTitle("\"Pion Veto Efficiency\"");
frock_graph->GetYaxis()->SetTitle("\"Kaon Missid rate\"");
frock_graph->GetXaxis()->SetLimits(0,1.01);
frock_graph->SetMinimum(0);
frock_graph->SetMaximum(1.01);

rock_graph->Draw("ACP");
c1->Print("roc_curve_fake.gif");
}
