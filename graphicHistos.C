double gaussian_ROC_integral(double sep, double sig)
{
	int npoints = 10000;
	double min = sep - 10*sig;
	double max = sep + 10*sig;
	double step = (max - min)/npoints;
	double ival = 0;
	double x = 0;
	double y = 0;	

	for (int i = 0; i < npoints; i++)
	{
		x = min + i*step;
		y = TMath::Gaus(x,sep,sig,true)*(.5+TMath::Erf(x/(sqrt(2)*sig))/2);
		ival += y*step;
	}
	return ival;
}
double find_sig_val(double seperation, double roc_integral, double sig_start, double stop_y = .00001, int max_iter = 1000)
{
	//newton's method - quick and dirty
	
	double dsig = .01;

	double cur_sig = sig_start;
	double cur_y = 0;

	double dydsig = 0;

	double cur_fval = 0;

	for (int i = 0; i < max_iter; i++)
	{
		cur_y = gaussian_ROC_integral(seperation,cur_sig);
		cur_fval = cur_y - roc_integral;

		if (fabs(cur_fval) < stop_y) break;
		
		dydsig = (gaussian_ROC_integral(seperation,cur_sig-dsig/2) - gaussian_ROC_integral(seperation,cur_sig+dsig/2))/dsig;
		
		cur_sig = cur_sig + cur_fval/dydsig;
	}

	return cur_sig;
}
double runGraphicHistos(TString ifile = "tmpfitdirc.root", bool verbose_out=true, double ienergy = 5, int iupdown = 0)
{
double hmin = -100;
double hmax = 100;

hmin = -25;
hmax = 25;

double energy = ienergy;

double pi_mass = .13957;
double k_mass = .49367;

double pi_beta = sqrt(1-pi_mass*pi_mass/(energy*energy));
double k_beta = sqrt(1-k_mass*k_mass/(energy*energy));

double quartz_index = 1.47;

double pi_mrad = 1000*acos(1/(pi_beta*quartz_index));
double k_mrad = 1000*acos(1/(k_beta*quartz_index));

double spread = 2;
//double seperation = 5.2;
double seperation = pi_mrad - k_mrad;
double mean_pion = 0;
double mean_kaon = mean_pion + seperation;
if (verbose_out == true)
{
	printf("Energy: %4.02f\n",energy);
	printf("Mrad Seperation: %8.03f\n",seperation);
}
//denominator pf spread
double spreadsq2 = 2*spread*spread;

TRandom3* randgen = new TRandom3();

TFile *f1 = new TFile(ifile);
TH1F *hpion;
TH1F *hkaon;
TH1F *phots_pion;
if (iupdown == 0)
{
	hpion = (TH1F*) f1->Get("ll_diff_pion");
	hkaon = (TH1F*) f1->Get("ll_diff_kaon");
	phots_pion = (TH1F*) f1->Get("phot_found_pion");
}
else if (iupdown == 1)
{
	hpion = (TH1F*) f1->Get("ll_diff_pion_up");
	hkaon = (TH1F*) f1->Get("ll_diff_kaon_up");
	phots_pion = (TH1F*) f1->Get("phot_found_pion_up");
}
else if (iupdown == -1)
{
	hpion = (TH1F*) f1->Get("ll_diff_pion_down");
	hkaon = (TH1F*) f1->Get("ll_diff_kaon_down");
	phots_pion = (TH1F*) f1->Get("phot_found_pion_down");
}
else
{
	printf("Unrecognize updown arguement: %d \nFailing....\n",iupdown);
	return -1;
}

printf("pion_ll mean, spread: %12.04f, %12.04f\n",hpion->GetMean(),hpion->GetRMS());
printf("kaon_ll mean, spread: %12.04f, %12.04f\n",hkaon->GetMean(),hkaon->GetRMS());

//Swap kaon and pion numbers

for (int i = 1; i < hpion->GetNbinsX()/2; i++)
{
	double t_swap = hpion->GetBinContent(i);
	hpion->SetBinContent(i,hpion->GetBinContent(hpion->GetNbinsX() - i + 1));
	hpion->SetBinContent(hpion->GetNbinsX() - i + 1, t_swap);
}
for (int i = 1; i < hkaon->GetNbinsX()/2; i++)
{
	double t_swap = hkaon->GetBinContent(i);
	hkaon->SetBinContent(i,hkaon->GetBinContent(hkaon->GetNbinsX() - i + 1));
	hkaon->SetBinContent(hkaon->GetNbinsX() - i + 1, t_swap);
}

int rebin = 20;
hpion->Rebin(rebin);
hkaon->Rebin(rebin);

hpion->SetAxisRange(hmin,hmax);
hkaon->SetAxisRange(hmin,hmax);

hkaon->GetXaxis()->SetTitle("Loglikelihood difference");

hpion->SetStats(false);
hkaon->SetStats(false);

hpion->SetLineColor(kCyan);
//hpion->SetFillColorAlpha(kRed,.5);

hkaon->SetLineColor(kBlue);
//hkaon->SetFillColorAlpha(kBlue,.5);

TLegend *leg_ll = new TLegend(.6,.6,.8,.8);
leg_ll->AddEntry(hpion,"Pion");
leg_ll->AddEntry(hkaon,"Kaon");

hkaon->SetTitle("log(P(Pi)/P(K)) for actual Pi (red) and K (blue) at 5 GeV");



TH1F *pion_veto_eff = new TH1F(*hpion);
TH1F *kaon_missid = new TH1F(*hkaon);

pion_veto_eff->SetName("pion_veto_eff");
pion_veto_eff->SetTitle("");

kaon_missid->SetName("kaon_missid");
kaon_missid->SetTitle("");



for (int i = 0; i < pion_veto_eff->GetNbinsX(); i++)
{
	pion_veto_eff->SetBinContent(i,hpion->Integral(0,i));
	kaon_missid->SetBinContent(i,hkaon->Integral(i,kaon_missid->GetNbinsX()));
//	printf("%12.04f %12.04f %d\n",1,1,i);
}

pion_veto_eff->SetAxisRange(0,10000,"Y");

double scale_int = 1/hpion->Integral(0,pion_veto_eff->GetNbinsX());
pion_veto_eff->Scale(scale_int);
scale_int = 1/hkaon->Integral(0,kaon_missid->GetNbinsX());
kaon_missid->Scale(scale_int);

hkaon->SetTitle("");
hpion->SetTitle("");


if (verbose_out == true)
{
    hkaon->Draw();
    hpion->Draw("SAME H");
    leg_ll->Draw("SAME");
    c1->SetWindowSize(1000,800);

    c1->Print("overlap.pdf");
}

if (verbose_out == true)
{
    pion_veto_eff->Draw("");
    kaon_missid->Draw("SAME H");

    c1->SetWindowSize(1000,800);
    c1->Print("overlap_integral.pdf");
}

TGraph* roc_graph;
int roc_n = pion_veto_eff->GetNbinsX();
TVectorF xr(roc_n);//gross
TVectorF yr(roc_n);
double ival = 0;

for (int i = 0; i < pion_veto_eff->GetNbinsX(); i++)
{
	xr[i] = pion_veto_eff->GetBinContent(i);
	yr[i] = kaon_missid->GetBinContent(i);

//	printf("%8.04f %8.04f\n",xr[i],yr[i]);	
}

double y1,y2,x1,x2;
x1 = pion_veto_eff->GetBinContent(0);
double last_x = pion_veto_eff->GetBinContent(0);
double last_y = kaon_missid->GetBinContent(0);

for (int i = 0; i < pion_veto_eff->GetNbinsX()-1; i++)
{
	ival += (yr[i]+last_y)*(xr[i] - last_x)/2;
	last_x = xr[i];
	last_y = yr[i];
}

if (verbose_out == true)
{
	printf("ROC integral: %12.04f\n",ival);
}

roc_graph = new TGraph(xr,yr);
roc_graph->SetLineColor(2);
roc_graph->SetLineWidth(4);
//roc_graph->SetMarkerColor(4);
//roc_graph->SetMarkerStyle(21);
roc_graph->SetTitle("");
roc_graph->GetXaxis()->SetTitle("Kaon Efficiency");
roc_graph->GetYaxis()->SetTitle("Pion Rejection");
roc_graph->GetXaxis()->SetLimits(0,1.01);
roc_graph->SetMinimum(0);
roc_graph->SetMaximum(1.01);

if (verbose_out == true)
{
    roc_graph->Draw("ACP");
    c1->Print("roc_curve.gif");
}
spread = find_sig_val(seperation,ival,spread); 


//FAKE version stuff below
/*---------------------------------------------------------------------------------------------------------------------------*/




TH1F *fhpion = (TH1F*) f1->Get("ll_diff_pion");
TH1F *fhkaon = (TH1F*) f1->Get("ll_diff_kaon");

fhpion->Reset();
fhkaon->Reset();

fhpion->SetBins(1000,hmin,hmax);
fhkaon->SetBins(1000,hmin,hmax);

double pion_obs, kaon_obs;
double pion_ll_diff, kaon_ll_diff;

for (int ii = 0; ii < 10000; ii++)
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

fhkaon->SetTitle("");

TH1F *fpion_veto_eff = new TH1F(*fhpion);
TH1F *fkaon_missid = new TH1F(*fhkaon);

fpion_veto_eff->SetName("pion_veto_eff");
fpion_veto_eff->SetTitle("");

fkaon_missid->SetName("kaon_missid");
fkaon_missid->SetTitle("");

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




TGraph* froc_graph;
int froc_n = fpion_veto_eff->GetNbinsX();
TVectorF fxr(froc_n);
TVectorF fyr(froc_n);

double fival = 0;
double flast_x = fpion_veto_eff->GetBinContent(0);
double flast_y = fkaon_missid->GetBinContent(0);
for (int i = 0; i < fpion_veto_eff->GetNbinsX(); i++)
{
        fxr[i] = fpion_veto_eff->GetBinContent(i);
        fyr[i] = fkaon_missid->GetBinContent(i);


        fival -= (fyr[i]+flast_y)*(fxr[i] - flast_x)/2;
        flast_x = fxr[i];
	flast_y = fyr[i];
}
ival = 0;
flast_x = fpion_veto_eff->GetBinContent(0);
flast_y = fkaon_missid->GetBinContent(0);
for (int i = 0; i < fpion_veto_eff->GetNbinsX(); i++)
{

//Why oh why is Erf not the standard definition
	double t = hmin + i*(hmax-hmin)/fpion_veto_eff->GetNbinsX();
	fxr[i] = .5 + TMath::Erf(t/(sqrt(2)*spread))/2;
	fyr[i] = .5 - TMath::Erf((t-seperation)/(sqrt(2)*spread))/2;

        fival -= (fyr[i]+flast_y)*(fxr[i] - flast_x)/2;
        flast_x = fxr[i];
	flast_y = fyr[i];
}	
//printf("Fake ROC integral: %12.04f\n",fival);

froc_graph = new TGraph(fxr,fyr);
froc_graph->SetLineColor(4);
froc_graph->SetLineWidth(4);
froc_graph->SetLineStyle(2);
froc_graph->SetTitle("");
froc_graph->GetXaxis()->SetTitle("\"Kaon Efficiency\"");
froc_graph->GetYaxis()->SetTitle("\"Pion Rejection\"");
froc_graph->GetXaxis()->SetLimits(0,1.01);
froc_graph->SetMinimum(0);
froc_graph->SetMaximum(1.01);


roc_graph->SetFillColorAlpha(kWhite,1);
froc_graph->SetFillColorAlpha(kWhite,1);
TLegend *leg_roc = new TLegend(.3,.5,.7,.7);
leg_roc->AddEntry(roc_graph,"ROC Curve");
leg_roc->AddEntry(froc_graph,"Matched Gaussian ROC Curve");


if (verbose_out == true)
{
	froc_graph->Draw("SAME");
	leg_roc->Draw("SAME");
	c1->Print("roc_curve_overlay.pdf");
}



if (verbose_out == true)
{
	printf("Matching resolution: %6.03f\n",spread);
	printf("Matching resolution per photon: %6.03f\n",spread*sqrt(phots_pion->GetMean()));
}
else
{
	printf("%6.04f\n",spread);
}
return spread;

}

void graphicHistos(TString ifile = "tmpfitdirc.root", bool verbose_out = true, double ienergy = 5, int updown=0)
{
	runGraphicHistos(ifile, verbose_out,ienergy,updown);
}
