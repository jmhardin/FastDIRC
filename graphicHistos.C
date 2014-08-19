{
double hmin = -20;
double hmax = 20;


double spread = 1.9;
double seperation = 5.23;
double mean_pion = 0;
double mean_kaon = mean_pion + seperation;

//denominator pf spread
double spreadsq2 = 2*spread*spread;

TRandom3* randgen = new TRandom3();

TFile *f1 = new TFile("./tmpfitdirc.root");
TH1F *hpion = (TH1F*) f1->Get("ll_diff_pion");
TH1F *hkaon = (TH1F*) f1->Get("ll_diff_kaon");

hpion->SetAxisRange(hmin,hmax);
hkaon->SetAxisRange(hmin,hmax);

hpion->SetLineColor(kRed);
//hpion->SetFillColorAlpha(kRed,.5);

hkaon->SetLineColor(kBlue);
//hkaon->SetFillColorAlpha(kBlue,.5);

hkaon->SetTitle("log(P(Pi)/P(K)) for actual Pi (red) and K (blue) 4.5 GeV");



TH1F *pion_veto_eff = new TH1F(*hpion);
TH1F *kaon_missid = new TH1F(*hkaon);

pion_veto_eff->SetName("pion_veto_eff");
pion_veto_eff->SetTitle("Pion Veto Efficiency");

kaon_missid->SetName("kaon_missid");
kaon_missid->SetTitle("Kaon Miss ID");

for (int i = 0; i < pion_veto_eff->GetNbinsX(); i++)
{
	pion_veto_eff->SetBinContent(i,hpion->Integral(i,pion_veto_eff->GetNbinsX()));
	kaon_missid->SetBinContent(i,hkaon->Integral(0,i));
}

pion_veto_eff->SetAxisRange(0,10000,"Y");

double scale_int = 1/hpion->Integral(0,pion_veto_eff->GetNbinsX());
pion_veto_eff->Scale(scale_int);
scale_int = 1/hkaon->Integral(0,kaon_missid->GetNbinsX());
kaon_missid->Scale(scale_int);

hkaon->Draw();
hpion->Draw("SAME H");
c1->SetWindowSize(1000,800);

c1->Print("overlap.gif");


pion_veto_eff->Draw("");
kaon_missid->Draw("SAME H");

c1->SetWindowSize(1000,800);
c1->Print("overlap_integral.gif");


TGraph* rock_graph;
int rock_n = pion_veto_eff->GetNbinsX();
double xr[rock_n], yr[rock_n];
double last_x = pion_veto_eff->GetBinContent(0);
double ival = 0;

for (int i = 0; i < pion_veto_eff->GetNbinsX(); i++)
{
	xr[i] = pion_veto_eff->GetBinContent(i);
	yr[i] = kaon_missid->GetBinContent(i);

	ival -= yr[i]*(xr[i] - last_x);
	last_x = xr[i];
//	printf("%8.04f %8.04f\n",xr[i],yr[i]);	
}
printf("ROC integral: %12.04f\n",ival);
rock_graph = new TGraph(rock_n,xr,yr);
rock_graph->SetLineColor(2);
rock_graph->SetLineWidth(4);
//rock_graph->SetMarkerColor(4);
//rock_graph->SetMarkerStyle(21);
rock_graph->SetTitle("ROC Curve");
rock_graph->GetXaxis()->SetTitle("Pion Veto Efficiency");
rock_graph->GetYaxis()->SetTitle("Kaon Missid rate");
rock_graph->GetXaxis()->SetLimits(0,1.01);
rock_graph->SetMinimum(0);
rock_graph->SetMaximum(1.01);

rock_graph->Draw("ACP");
//c1->Print("roc_curve.gif");




//FAKE version stuff below
/*---------------------------------------------------------------------------------------------------------------------------*/




TH1F *fhpion = (TH1F*) f1->Get("ll_diff_pion");
TH1F *fhkaon = (TH1F*) f1->Get("ll_diff_kaon");

fhpion->Reset();
fhkaon->Reset();

double pion_obs, kaon_obs;
double pion_ll_diff, kaon_ll_diff;

for (int ii = 0; ii < 100000; ii++)
{
        pion_obs = randgen->Gaus(mean_pion,spread);
        kaon_obs = randgen->Gaus(mean_kaon,spread);

        pion_ll_diff = -1*(pion_obs - mean_pion)*(pion_obs - mean_pion);
        pion_ll_diff += (pion_obs - mean_kaon)*(pion_obs - mean_kaon);
        pion_ll_diff /= spreadsq2;

        kaon_ll_diff = - (kaon_obs - mean_pion)*(kaon_obs - mean_pion);
        kaon_ll_diff += (kaon_obs - mean_kaon)*(kaon_obs - mean_kaon);
        kaon_ll_diff /= spreadsq2;

        fhpion->Fill(pion_ll_diff);
        fhkaon->Fill(kaon_ll_diff);
}


fhpion->SetAxisRange(hmin,hmax);
fhkaon->SetAxisRange(hmin,hmax);

fhpion->SetLineColor(kRed);
//hpion->SetFillColorAlpha(kRed,.5);

fhkaon->SetLineColor(kBlue);
//hkaon->SetFillColorAlpha(kBlue,.5);

fhkaon->SetTitle("Fake Seperation distance = 5.23 mrad.  spread = 2.5mrad");

TH1F *fpion_veto_eff = new TH1F(*fhpion);
TH1F *fkaon_missid = new TH1F(*fhkaon);

fpion_veto_eff->SetName("pion_veto_eff");
fpion_veto_eff->SetTitle("Pion Veto Efficiency");

fkaon_missid->SetName("kaon_missid");
fkaon_missid->SetTitle("Kaon Miss ID");

for (int i = 0; i < fpion_veto_eff->GetNbinsX(); i++)
{
        fpion_veto_eff->SetBinContent(i,fhpion->Integral(i,fpion_veto_eff->GetNbinsX()));
        fkaon_missid->SetBinContent(i,fhkaon->Integral(0,i));
}

fpion_veto_eff->SetAxisRange(0,10000,"Y");

double fscale_int = 1/fhpion->Integral(0,fpion_veto_eff->GetNbinsX());
fpion_veto_eff->Scale(fscale_int);
fscale_int = 1/fhkaon->Integral(0,fkaon_missid->GetNbinsX());
fkaon_missid->Scale(fscale_int);
/*
Move down so the ROC Curves can be overlayed
fhkaon->Draw();
fhpion->Draw("SAME H");
c1->SetWindowSize(1000,800);

c1->Print("overlap_fake.gif");


fpion_veto_eff->Draw("");
fkaon_missid->Draw("SAME H");

c1->Print("overlap_integral_fake.gif");
*/

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
frock_graph->SetLineColor(4);
frock_graph->SetLineWidth(4);
frock_graph->SetTitle("Fake ROC graph");
frock_graph->GetXaxis()->SetTitle("\"Pion Veto Efficiency\"");
frock_graph->GetYaxis()->SetTitle("\"Kaon Missid rate\"");
frock_graph->GetXaxis()->SetLimits(0,1.01);
frock_graph->SetMinimum(0);
frock_graph->SetMaximum(1.01);

frock_graph->Draw("CP");
c1->Print("roc_curve_fake.gif");

}

