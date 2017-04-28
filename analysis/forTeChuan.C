//================================================
void forTeChuan()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  theory();
  //pp500();
}

//================================================
void pp500()
{
  // systematics for MTD response efficiency
  TFile *f1 = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","read");
  TH1F *hSys1 = (TH1F*)f1->Get("hSysMtdRespEff_1");
  TH1F *hSys2 = (TH1F*)f1->Get("hSysMtdRespEff_2");


  // sysematics for smearing
  TFile *f2 = TFile::Open("Rootfiles/Run13.pp500.TrkResScan.root","read");
  TH1F *hdef = (TH1F*)f2->Get("JpsiSmearEff_def_cent00100_rebin");
  TH1F *hmax = (TH1F*)f2->Get("JpsiSmearEff_max_cent00100_rebin");
  TH1F *hmin = (TH1F*)f2->Get("JpsiSmearEff_min_cent00100_rebin");
  TH1F *hSys3 = (TH1F*)hdef->Clone("SysMcSmearing");
  hSys3->Reset();
  for(int bin=1; bin<=hSys3->GetNbinsX(); bin++)
    {
      double def = hdef->GetBinContent(bin);
      double max = hmax->GetBinContent(bin);
      double min = hmin->GetBinContent(bin);
      double sys = (fabs(max-def)+fabs(min-def))/2;
      hSys3->SetBinContent(bin,1);
      hSys3->SetBinError(bin,sys);
    }

  // my pp 500 results
  TFile *f3 = TFile::Open("Rootfiles/Pico.Run13.pp500.jpsi.VtxCut.pt1.5.pt1.0.xsec.root","read");
  TH1F *hyield = (TH1F*)f3->Get("Jpsi_InvYield_cent00100");

  TFile *fout = TFile::Open("Rootfiles/Run13.pp500.forTeChuan.root","recreate");
  hSys1->Write("SysMtdRespEff_stat");
  hSys2->Write("SysMtdRespEff_sys");
  hSys3->Write();
  hyield->Write("Jpsi_InvYield");
}

//================================================
void theory()
{
  TCanvas *c2 = new TCanvas("c2","c2", 700, 700);
  SetPadMargin(gPad,0.12, 0.14, 0.03,0.03);
  gPad->SetLogy();
  TH1F *h = new TH1F("h2",";p_{T} (GeV/c);B#times1/(2#pip_{T})#timesd^{2}#sigma/(dp_{T}dy)   (nb/GeV/c)^{2}",10,0,25);
  ScaleHistoTitle(h,0.045,1,0.035,0.045,1.4,0.035,62);
  h->GetYaxis()->SetRangeUser(1e-6,100);
  h->GetXaxis()->SetRangeUser(0,22.5);
  h->GetYaxis()->CenterTitle(1);
  h->Draw();

  // BNL + PKU
  const double br = 0.0594;
  const double phase = 4*pi;
  const int cgc_low_npoints = 59;
  double cgc_lowPt[cgc_low_npoints][3] = {{0.2, 312.223, 614.502}, {0.3, 467.386, 919.656}, {0.4, 625.062,1229.52}, {0.5, 773.646, 1521.35}, {0.6, 922.245, 1812.75}, {0.7,1059.23, 2081.11}, {0.8, 1177.25, 2312.14}, {0.9, 1272.2,2498.04}, {1., 1342.55, 2635.53}, {1.1, 1391.37, 2730.39}, {1.2,1414.01, 2773.97}, {1.3, 1410.09, 2765.83}, {1.4, 1382.57,2711.7}, {1.5, 1336.26, 2620.82}, {1.6, 1274.51, 2499.62}, {1.7,1203.05, 2359.43}, {1.8, 1125.74, 2207.75}, {1.9, 1045.77,2050.91}, {2., 965.624, 1893.87}, {2.1, 886.599, 1739.33}, {2.2,811.005, 1591.53}, {2.3, 739.973, 1452.56}, {2.4, 673.882,1323.15}, {2.5, 612.682, 1203.32}, {2.6, 556.144, 1092.67}, {2.7,504.512, 991.62}, {2.8, 457.631, 899.854}, {2.9, 415.23,816.838}, {3., 376.98, 741.938}, {3.1, 342.763, 674.932}, {3.2,311.949, 614.579}, {3.3, 284.155, 560.133}, {3.4, 259.053,510.95}, {3.5, 236.355, 466.466}, {3.6, 215.774, 426.126}, {3.7,197.142, 389.595}, {3.8, 180.264, 356.487}, {3.9, 164.962,326.456}, {4., 151.075, 299.19}, {4.1, 138.491, 274.461}, {4.2,127.037, 251.942}, {4.3, 116.6, 231.412}, {4.4, 107.077,212.674}, {4.5, 98.3808, 195.555}, {4.6, 90.373, 179.793}, {4.7,83.0709, 165.412}, {4.8, 76.4167, 152.298}, {4.9, 70.355,140.345}, {5., 64.834, 129.448}, {5.1, 59.8364, 119.57}, {5.2,55.2726, 110.543}, {5.3, 51.0982, 102.279}, {5.4, 47.2741,94.703}, {5.5, 43.7658, 87.7481}, {5.6, 40.5324, 81.3353}, {5.7,37.5627, 75.4416}, {5.8, 34.8342, 70.0233}, {5.9, 32.3263,65.0398}, {6., 30.0201, 60.4543}};
  double cgc_low_pt[cgc_low_npoints];
  double cgc_low_pt_err[cgc_low_npoints];
  double cgc_low_y[cgc_low_npoints];
  double cgc_low_yh[cgc_low_npoints];
  double cgc_low_yl[cgc_low_npoints];
  for(int i=0; i<cgc_low_npoints; i++)
    {
      double scale = br * 1./phase * 1./cgc_lowPt[i][0];
      cgc_low_pt[i] = cgc_lowPt[i][0];
      cgc_low_pt_err[i] = 0;
      cgc_low_y[i] = (cgc_lowPt[i][1]+cgc_lowPt[i][2])/2 * scale;
      cgc_low_yh[i] = (cgc_lowPt[i][1]* scale - cgc_low_y[i]) ;
      cgc_low_yl[i] = (cgc_low_y[i] - cgc_lowPt[i][2] * scale);
    }
  TGraphAsymmErrors *gCgcLowPt = new TGraphAsymmErrors(cgc_low_npoints,cgc_low_pt,cgc_low_y,cgc_low_pt_err,cgc_low_pt_err,cgc_low_yl,cgc_low_yh);
  gCgcLowPt->SetFillStyle(1001);
  gCgcLowPt->SetLineColor(kOrange+1);
  gCgcLowPt->SetFillColor(gCgcLowPt->GetLineColor());

  const int cgc_high_npoints = 33;
  double cgc_highPt[cgc_high_npoints][3] = {{4., 212.673, 454.546}, {4.5, 114.815, 245.635}, {5., 63.9102,136.874}, {5.5, 36.7717, 78.852}, {6., 21.7079, 46.6017}, {6.5,13.1481, 28.2551}, {7., 8.16417, 17.5623}, {7.5, 5.17138,11.1338}, {8., 3.36145, 7.24488}, {8.5, 2.23414, 4.82106}, {9.,1.51397, 3.27119}, {9.5, 1.04376, 2.25809}, {10., 0.730681,1.58275}, {10.5, 0.518107, 1.12357}, {11., 0.372124,0.80788}, {11.5, 0.270555, 0.588004}, {12., 0.199, 0.432953}, {12.5,0.148088, 0.32254}, {13., 0.111343, 0.242775}, {13.5, 0.0845518,0.184553}, {14., 0.0647922, 0.141571}, {14.5, 0.0500633,0.109506}, {15., 0.0389818, 0.0853615}, {15.5, 0.0305712,0.0670226}, {16., 0.0241333, 0.0529745}, {16.5, 0.01917,0.0421348}, {17., 0.0153213, 0.0337202}, {17.5, 0.0123179,0.027146}, {18., 0.00996177, 0.0219814}, {18.5, 0.00810017,0.0178955}, {19., 0.00661969, 0.0146425}, {19.5, 0.00543558,0.0120378}, {20., 0.00448337, 0.00994106}};
  double cgc_high_pt[cgc_high_npoints];
  double cgc_high_pt_err[cgc_high_npoints];
  double cgc_high_y[cgc_high_npoints];
  double cgc_high_yh[cgc_high_npoints];
  double cgc_high_yl[cgc_high_npoints];
  for(int i=0; i<cgc_high_npoints; i++)
    {
      double scale = br * 1./phase * 1./cgc_highPt[i][0];
      cgc_high_pt[i] = cgc_highPt[i][0];
      cgc_high_pt_err[i] = 0;
      cgc_high_y[i] = (cgc_highPt[i][1]+cgc_highPt[i][2])/2 * scale;
      cgc_high_yh[i] = (cgc_highPt[i][1] * scale - cgc_high_y[i]);
      cgc_high_yl[i] = (cgc_high_y[i] - cgc_highPt[i][2] * scale);
    }
  TGraphAsymmErrors *gCgcHighPt = new TGraphAsymmErrors(cgc_high_npoints,cgc_high_pt,cgc_high_y,cgc_high_pt_err,cgc_high_pt_err,cgc_high_yl,cgc_high_yh);
  gCgcHighPt->SetFillStyle(1001);
  gCgcHighPt->SetLineColor(kCyan+1);
  gCgcHighPt->SetFillColor(gCgcHighPt->GetLineColor());
  gCgcLowPt->Draw("sames E3");
  gCgcHighPt->Draw("sames E3");

  // data
  TFile *fdata = TFile::Open(Form("sptrum.root"),"read");	
  TGraphErrors *gData = (TGraphErrors*)fdata->Get("Jpsi_pp500");
  gData->Draw("sames PEZ");
  for(int i=0; i<19; i++)
    {
      TBox *box = (TBox*)fdata->Get(Form("sys_uncert_%d",i));
      box->Draw();
    }

  leg = new TLegend(0.5,0.7,0.7,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 500 GeV");
  leg->AddEntry(gData,"STAR  J/#psi#rightarrowe^{+}e^{-}, |y|<1","P");
  leg->AddEntry(gCgcLowPt,"CGC+NRQCD","F");
  leg->AddEntry(gCgcHighPt,"NLO NRQCD","F");
  leg->Draw();

}
