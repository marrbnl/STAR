const char *run_type = "Run14_AuAu200";
TString run_cfg_name = "2016sQM";
const double luminosity = 14.2;
const double ppInelastic = 42.; // mb
const double ppInelasticErr = 3.; // mb

//================================================
void plot_2016sQM()
{  
  //rawSignal();
  //xsec();
  //ppRef();
  nPart();
  //compareQM2015();
}

//================================================
void nPart(const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptStat(0);
  const double ncoll[6] = {64, 124, 224, 377, 609, 964}; 
  const double ncollErr[6] = {18, 25, 30, 33, 31,  28};
  const double npart[6] =    {48, 76, 116, 168, 236, 325};
  const double npartErr[6] = {9,  11, 11,  11,  9,   5};
  
  // Ncoll uncertainty
  double pub_npart[6]     = {325, 236, 168, 116, 76, 48};
  double pub_npart_err[6] = {5,   9,   11,  11,  11, 9};
  double pub_ncoll[6]     = {964, 609, 377, 224, 124, 64};
  double pub_ncoll_err[6] = {28,  31,  33,  30,  25,  18};
  double pub_npart_xerr[6];
  double pub_npart_xsys[6];
  double pub_npart_relerr[6];
  double pub_npart_y[6];
  for(int i=0; i<6; i++)
    {
      pub_npart_xerr[i] = 0;
      pub_npart_xsys[i] = 6;
      pub_npart_y[i] = 1;
      pub_npart_relerr[i] = pub_ncoll_err[i]/pub_ncoll[i];
    }
  TGraphErrors *gPubNpartErr = new TGraphErrors(6, pub_npart, pub_npart_y, pub_npart_err, pub_npart_relerr);
  gPubNpartErr->SetFillColor(kGreen-10);
  gPubNpartErr->SetLineColor(kGreen-10);

  // template
  TH1F *hRaaVsNpart = new TH1F("hRaaVsNpart",";N_{part};J/#psi R_{AA}",100,-100,420);
  ScaleHistoTitle(hRaaVsNpart,28,1,24,28,1,24,63);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,2);
  hRaaVsNpart->GetYaxis()->CenterTitle();

  const int npart_color[3] = {1, 4, 2};

  // STAR published
  double pub_raa_vs_npart[5] = {0.64, 0.68, 0.72, 0.95, 1.08};
  double pub_raa_vs_npart_err[5] = {0.14, 0.14, 0.14, 0.16, 0.14};
  double pub_raa_vs_npart_sys[5] = {0.03, 0.03, 0.06, 0.06, 0.08};
  TGraphErrors *pubRaaVsNpart = new TGraphErrors(5, pub_npart, pub_raa_vs_npart, pub_npart_xerr, pub_raa_vs_npart_err);
  pubRaaVsNpart->SetMarkerStyle(20);
  pubRaaVsNpart->SetMarkerColor(npart_color[0]);
  pubRaaVsNpart->SetMarkerSize(1.8);
  pubRaaVsNpart->SetLineColor(npart_color[0]);
  TGraphErrors *pubRaaVsNpartSys = new TGraphErrors(5, pub_npart, pub_raa_vs_npart, pub_npart_xsys, pub_raa_vs_npart_sys);
  pubRaaVsNpartSys->SetMarkerStyle(24);
  pubRaaVsNpartSys->SetMarkerColor(npart_color[0]);
  pubRaaVsNpartSys->SetMarkerSize(0);
  pubRaaVsNpartSys->SetLineColor(npart_color[0]);
  pubRaaVsNpartSys->SetFillStyle(0);
  TBox *star_pub = new TBox(360,1-0.14,365,1+0.14);
  star_pub->SetLineColor(kGray+3);
  star_pub->SetFillColor(kGray+3);
  star_pub->SetLineWidth(2.);
  star_pub->SetFillStyle(1001);

 // MTD results
  TFile *fref = TFile::Open(Form("Rootfiles/%s/Publication.Jpsi.200GeV.root",run_cfg_name.Data()),"read");
  TH1F *hpp = (TH1F*)fref->Get("pp200_Jpsi_Integrated");

  const char *name[2] = {"0-15","5-15"};
  TFile *fdata = TFile::Open("Rootfiles/Run14_AuAu200.Npart.pt1.5.pt1.2.root","read");
  TFile *fsys = TFile::Open("Rootfiles/Run14_AuAu200.Npart.pt1.5.pt1.2.sys.root","read");
  TH1F *hdata[2], *hsys[2];
  TH1F *hRaa[2];
  TGraphErrors *raaVsNpart[2], *raaVsNpartSys[2];
  TBox *globalSys[2];
  for(int i=0; i<2; i++)
    {
      hdata[i] = (TH1F*)fdata->Get(Form("Jpsi_InvYield_pt%s",name[i]));
      hsys[i] = (TH1F*)fsys->Get(Form("Sys_all_pt%s",name[0]));
      int npoints = hdata[i]->GetNbinsX();
      raaVsNpart[i] = new TGraphErrors(npoints);
      raaVsNpart[i]->SetName(Form("Run14_raaVsNpart_pt%s",name[i]));
      raaVsNpartSys[i] = new TGraphErrors(npoints);
      raaVsNpartSys[i]->SetName(Form("Run14_raaVsNpartSys_pt%s",name[i]));

      hRaa[i] = (TH1F*)hdata[i]->Clone(Form("hRaaVsNpart_pt%s",name[i]));
      hRaa[i]->Reset();

      double pp_yield = hpp->GetBinContent(i+1);
      double pp_err = hpp->GetBinError(i+1);
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  double x = npart[ipoint];
	  double auau_yield = hdata[i]->GetBinContent(ipoint+1) * 2 * pi;
	  double auau_err   = hdata[i]->GetBinError(ipoint+1) * 2 * pi;
	  double raa = ppInelastic/ncoll[ipoint] * 1e6 * auau_yield/pp_yield;
	  double raa_err =  auau_err / auau_yield * raa;
	  double raa_sys = hsys[i]->GetBinContent(ipoint+1) * raa;
	  raaVsNpart[i]->SetPoint(ipoint,x,raa);
	  raaVsNpart[i]->SetPointError(ipoint,0,raa_err);
	  raaVsNpartSys[i]->SetPoint(ipoint,x,raa);
	  raaVsNpartSys[i]->SetPointError(ipoint,5,raa_sys);
	  hRaa[i]->SetBinContent(ipoint+1, raa);
	  hRaa[i]->SetBinError(ipoint+1, raa_err);
	  cout << auau_yield << "  " << pp_yield << "  " << raa << "  " << raa_err << "  " << raa_sys << endl;
	}

      raaVsNpart[i]->SetMarkerStyle(29);
      raaVsNpart[i]->SetMarkerColor(npart_color[i+1]);
      raaVsNpart[i]->SetMarkerSize(2.5);
      raaVsNpart[i]->SetLineColor(npart_color[i+1]);
      raaVsNpartSys[i]->SetMarkerStyle(20);
      raaVsNpartSys[i]->SetMarkerColor(npart_color[i+1]);
      raaVsNpartSys[i]->SetMarkerSize(0);
      raaVsNpartSys[i]->SetLineColor(npart_color[i+1]);
      raaVsNpartSys[i]->SetLineWidth(1);
      raaVsNpartSys[i]->SetFillStyle(0);

      double gSys = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err*pp_err/pp_yield/pp_yield);
      globalSys[i] = new TBox(350+i*5,1-gSys,355+i*5,1+gSys);
      globalSys[i]->SetLineColor(npart_color[i+1]);
      globalSys[i]->SetFillColor(npart_color[i+1]);
      globalSys[i]->SetLineWidth(2.);
      globalSys[i]->SetFillStyle(1001);
      cout << gSys << endl;
    }

  //==============================================
  // Raa vs Npart dimuon
  //==============================================
  TCanvas *c0 = new TCanvas("Jpsi_raa_vs_npart_dimuon","Jpsi_raa_vs_npart_dimuon",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->GetXaxis()->SetRangeUser(0,420);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,1.5);
  hRaaVsNpart->DrawCopy();
  //gPubNpartErr->Draw("e3sames");
  TLine *line = GetLine(0,1,370,1,1);
  //line->Draw();
  for(int i=1; i<2; i++)
    {
      raaVsNpart[i]->Draw("samesPEZ");
      raaVsNpartSys[i]->Draw("samesE5");
      //globalSys[i]->Draw("fsame");
    }
  TLegend *leg0 = new TLegend(0.16,0.75,0.4,0.95);
  leg0->SetBorderSize(0);
  leg0->SetFillColor(0);
  leg0->SetTextFont(62);
  leg0->SetTextSize(0.04);
  leg0->SetHeader("Au+Au @ 200 GeV");
  leg0->AddEntry(raaVsNpart[0],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 0 GeV/c","P");
  leg0->AddEntry(raaVsNpart[1],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 5 GeV/c","P");
  //leg0->Draw();
  TLegend *leg_ncoll = new TLegend(0.15,0.17,0.35,0.22);
  leg_ncoll->SetBorderSize(0);
  leg_ncoll->SetFillColor(0);
  leg_ncoll->SetTextFont(62);
  leg_ncoll->SetTextSize(0.035);
  leg_ncoll->AddEntry(gPubNpartErr,"N_{coll} uncertainty","f");
  //leg_ncoll->Draw();
  TPaveText *star = GetPaveText(0.75,0.85,0.92,0.97,26,23);
  //star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c0->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_dimuon.pdf",run_type,run_cfg_name.Data()));
      c0->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_dimuon.png",run_type,run_cfg_name.Data()));
      c0->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_dimuon.jpg",run_type,run_cfg_name.Data()));
      c0->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_dimuon.eps",run_type,run_cfg_name.Data()));
    }
  return;

  // Add RpA data points
  TGraphErrors *rpaVsNpart[2];
  TGraphErrors *rpaVsNpartSys[2];
  double rpa[2] = {0.689, 0.814};
  double rpa_err[2] = {0.0285, 0.0687};
  double rpa_sys[2] = {0.0745, 0.0880};
  double rpa_gsys[2] = {0.1269, 0.1431};
  double rpa_npart[2] = {5.7, 5.7};
  double rpa_npart_err[2] = {0.3, 0.3};
  for(int i=0; i<2; i++)
    {
      rpaVsNpart[i] = new TGraphErrors(1);
      rpaVsNpart[i]->SetPoint(0, rpa_npart[i], rpa[i]);
      rpaVsNpart[i]->SetPointError(0, rpa_npart_err[i], rpa_err[i]);
      rpaVsNpartSys[i] = new TGraphErrors(1);
      rpaVsNpartSys[i]->SetPoint(0, rpa_npart[i], rpa[i]);
      rpaVsNpartSys[i]->SetPointError(0, 5, rpa_sys[i]);

      rpaVsNpart[i]->SetMarkerStyle(24);
      rpaVsNpart[i]->SetMarkerColor(npart_color[i+1]);
      rpaVsNpart[i]->SetMarkerSize(1.5);
      rpaVsNpart[i]->SetLineColor(npart_color[i+1]);
      rpaVsNpartSys[i]->SetMarkerColor(npart_color[i+1]);
      rpaVsNpartSys[i]->SetMarkerSize(0);
      rpaVsNpartSys[i]->SetLineColor(npart_color[i+1]);
      rpaVsNpartSys[i]->SetLineWidth(1);
      rpaVsNpartSys[i]->SetFillStyle(0);
    }


  TCanvas *c00 = new TCanvas("Jpsi_rpaa_vs_npart_dimuon","Jpsi_rpaa_vs_npart_dimuon",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->GetXaxis()->SetRangeUser(-10,375);
  hRaaVsNpart->SetYTitle("J/#psi R_{pA,AA}");
  hRaaVsNpart->DrawCopy();
  hRaaVsNpart->SetYTitle("J/#psi R_{AA}");
  gPubNpartErr->Draw("e3sames");
  TBox *npart_global_sys = new TBox(rpa_npart[0]-5,1-0.054,rpa_npart[0]+5,1+0.054);
  npart_global_sys->SetLineColor(kGreen-10);
  npart_global_sys->SetFillColor(kGreen-10);
  npart_global_sys->SetFillStyle(1001);
  npart_global_sys->Draw("samef");
  line = GetLine(-10,1,375,1,1);
  line->Draw();
  for(int i=0; i<2; i++)
    {
      raaVsNpart[i]->Draw("samesPEZ");
      raaVsNpartSys[i]->Draw("samesE5");
      globalSys[i]->Draw("fsame");

      rpaVsNpart[i]->Draw("samesPEZ");
      rpaVsNpartSys[i]->Draw("samesE5");

      TBox *rpa_global_sys = new TBox(360+i*6.8,1-rpa_gsys[i],365+i*6.8,1+rpa_gsys[i]);
      rpa_global_sys->SetLineColor(npart_color[i+1]);
      rpa_global_sys->SetFillColor(npart_color[i+1]);
      rpa_global_sys->SetLineWidth(2.);
      rpa_global_sys->SetFillStyle(0);
      rpa_global_sys->Draw("samef");
    }

  TLegend *leg00 = new TLegend(0.16,0.7,0.4,0.95);
  leg00->SetBorderSize(0);
  leg00->SetFillColor(0);
  leg00->SetTextFont(62);
  leg00->SetTextSize(0.04);
  leg00->SetHeader("STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg00->AddEntry(raaVsNpart[0],"Au+Au @ 200 GeV, p_{T} > 0 GeV/c","P");
  leg00->AddEntry(raaVsNpart[1],"Au+Au @ 200 GeV, p_{T} > 5 GeV/c","P");
  leg00->AddEntry(rpaVsNpart[0],"p+Au @ 200 GeV, p_{T} > 1 GeV/c","P");
  leg00->AddEntry(rpaVsNpart[1],"p+Au @ 200 GeV, p_{T} > 5 GeV/c","P");
  leg00->Draw();
  leg_ncoll->Draw();
  star->Draw();
  if(savePlot)
    {
      c00->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaRpaNpart_dimuon.pdf",run_type,run_cfg_name.Data()));
      c00->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaRpaNpart_dimuon.png",run_type,run_cfg_name.Data()));
      c00->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaRpaNpart_dimuon.jpg",run_type,run_cfg_name.Data()));
      c00->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaRpaNpart_dimuon.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // Raa vs Npart dimuon vs dielectron
  //==============================================
  TCanvas *c1 = new TCanvas("Jpsi_raa_vs_npart","Jpsi_raa_vs_npart",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->GetXaxis()->SetRangeUser(0,370);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  TLine *line = GetLine(0,1,370,1,1);
  line->Draw();
  pubRaaVsNpart->Draw("samesPEZ");
  pubRaaVsNpartSys->Draw("samesE5");
  star_pub->Draw("fsame");
  for(int i=0; i<2; i++)
    {
      raaVsNpart[i]->Draw("samesPEZ");
      raaVsNpartSys[i]->Draw("samesE5");
      globalSys[i]->Draw("fsame");
    }
  TLegend *leg1 = new TLegend(0.16,0.75,0.4,0.95);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.035);
  leg1->SetHeader("Au+Au @ 200 GeV");
  leg1->AddEntry(raaVsNpart[0],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 0 GeV/c","P");
  leg1->AddEntry(raaVsNpart[1],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 5 GeV/c","P");
  leg1->AddEntry(pubRaaVsNpart,"STAR J/#psi#rightarrowe^{+}e^{-}, |y| < 1, p_{T} > 5 GeV/c","P");
  leg1->Draw();
  leg_ncoll->Draw();
  star->Draw();
  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_pub.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_pub.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_pub.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_pub.eps",run_type,run_cfg_name.Data()));
    }


  //==============================================
  // Raa vs model
  //==============================================
  // tsinghua group
  char tmp[256];
  double dtmp;
  ifstream ifin;
  ifin.open("Rootfiles/Published/Jpsi_Raa_200/model/data_and_figures_from_Yunpeng/data_pt_5-10_GeV.txt");
  ifin.getline(tmp,256);
  ifin.getline(tmp,256);
  double raa_5_tsu[11];
  double npart_tsu[11];
  for(int i=0; i<11; i++)
    {
      ifin >> npart_tsu[i] >> raa_5_tsu[i];
      cout << npart_tsu[i] << "  " <<   raa_5_tsu[i] << endl;
    }
  ifin.close();
  TGraph *gTsuRaaVsNpart[2];
  gTsuRaaVsNpart[0] = new TGraph("Rootfiles/Published/Zebo/QWG2013/Jpsi_Cent/Jpsi_Zhuang_RHIC_mid.dat");
  gTsuRaaVsNpart[1] = new TGraph(11,npart_tsu,raa_5_tsu);
  for(int i=0; i<2; i++)
    {
      gTsuRaaVsNpart[i]->SetMarkerStyle(21);
      gTsuRaaVsNpart[i]->SetMarkerColor(npart_color[i+1]);
      gTsuRaaVsNpart[i]->SetLineColor(npart_color[i+1]);
      gTsuRaaVsNpart[i]->SetLineWidth(2);
      gTsuRaaVsNpart[i]->SetLineStyle(1);
    }
  
  // tamu group
  ifin.open("Rootfiles/Published/Jpsi_Raa_200/model/data_from_Xingbo/raa_centra_rhic_pt4.5to10_tozebo_110513.dat");
  double raa_5_tamu[15];
  double npart_tamu[15];
  for(int i=0; i<15; i++)
    {
      ifin >> npart_tamu[i] >> raa_5_tamu[i];
    }
  ifin.close();
  TGraph *gTamuRaaVsNpart[2];
  gTamuRaaVsNpart[1] = new TGraph(15,npart_tamu,raa_5_tamu);

  TFile *ftamu = TFile::Open("Rootfiles/2016sQM/raa200theory1.root","read");
  TGraph *gtmp = (TGraph*)ftamu->Get("gr1");
  gTamuRaaVsNpart[0] = new TGraph(gtmp->GetN(), gtmp->GetX(), gtmp->GetY());
  for(int i=0; i<2; i++)
    {
      gTamuRaaVsNpart[i]->SetMarkerStyle(20);
      gTamuRaaVsNpart[i]->SetMarkerColor(npart_color[i+1]);
      gTamuRaaVsNpart[i]->SetLineColor(npart_color[i+1]);
      gTamuRaaVsNpart[i]->SetLineWidth(2);
      gTamuRaaVsNpart[i]->SetLineStyle(7);
    }

  TCanvas *c2 = new TCanvas("Jpsi_Raa_vs_model","Jpsi_Raa_vs_model",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  line->Draw();
  gTsuRaaVsNpart[0]->Draw("lsames");
  gTsuRaaVsNpart[1]->Draw("lsames");
  gTamuRaaVsNpart[0]->Draw("lsames");
  gTamuRaaVsNpart[1]->Draw("lsames");
  pubRaaVsNpart->Draw("samesPEZ");
  pubRaaVsNpartSys->Draw("samesE5");
  star_pub->Draw("fsame");
  for(int i=0; i<2; i++)
    {
      raaVsNpart[i]->Draw("samesPEZ");
      raaVsNpartSys[i]->Draw("samesE5");
      globalSys[i]->Draw("fsame");
    }
  leg1->Draw();
  leg_ncoll->Draw();
  star->Draw();
  
  TPaveText *tmodel = GetPaveText(0.213,0.35,0.64,0.74);
  tmodel->SetTextFont(62);
  tmodel->SetTextAlign(11);
  tmodel->SetTextSize(0.035);
  tmodel->AddText("Tsinghua Model");
  tmodel->AddText("TAMU Model");
  tmodel->Draw();

  TLegend *leg21 = new TLegend(0.45,0.65,0.67,0.75);
  leg21->SetBorderSize(0);
  leg21->SetFillColor(0);
  leg21->SetTextFont(62);
  leg21->SetTextSize(0.035);
  leg21->AddEntry(gTsuRaaVsNpart[0],"p_{T} > 0 GeV/c","L");
  leg21->AddEntry(gTamuRaaVsNpart[0],"p_{T} > 0 GeV/c","L");
  leg21->Draw();
  TLegend *leg22 = new TLegend(0.67,0.65,0.85,0.75);
  leg22->SetBorderSize(0);
  leg22->SetFillColor(0);
  leg22->SetTextFont(62);
  leg22->SetTextSize(0.035);
  leg22->AddEntry(gTsuRaaVsNpart[1],"p_{T} > 5 GeV/c","L");
  leg22->AddEntry(gTamuRaaVsNpart[1],"p_{T} > 5 GeV/c","L");
  leg22->Draw();

  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_model.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_model.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_model.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpart_vs_model.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // STAR vs LHC low pT
  //==============================================
  
  // PHENIX low pT
  // http://www.phenix.bnl.gov/phenix/WWW/info/data/ppg068/auau_electrons.txt
  // PRL 98 (2007) 232301
  // global uncertainty = 12%
  const int nPhenix_lowpT = 8;
  double phenixLowPt_npart[nPhenix_lowpT] = {351.4, 299.0, 253.9, 215.3, 166.6, 114.2, 58.4, 14.5};
  double phenixLowPt_npart_stat[nPhenix_lowpT] = {0, 0, 0, 0, 0, 0, 0, 0};
  double phenixLowPt_npart_sys[nPhenix_lowpT] = {5, 5, 5, 5, 5, 5, 5, 5};
  double phenixLowPt_raa[nPhenix_lowpT] = {0.26, 0.34, 0.36, 0.45, 0.58, 0.58, 0.65, 0.74};
  double phenixLowPt_raa_stat[nPhenix_lowpT] = {0.05, 0.06, 0.06, 0.07, 0.07, 0.08, 0.07, 0.12};
  double phenixLowPt_raa_sys[nPhenix_lowpT] = {0.04, 0.05, 0.05, 0.07, 0.08, 0.08, 0.10, 0.21};
  TGraphErrors *grPhenixLowPt = new TGraphErrors(nPhenix_lowpT, phenixLowPt_npart, phenixLowPt_raa, phenixLowPt_npart_stat, phenixLowPt_raa_stat);
  TGraphErrors *grPhenixLowPtSys  = new TGraphErrors(nPhenix_lowpT, phenixLowPt_npart, phenixLowPt_raa, phenixLowPt_npart_sys, phenixLowPt_raa_sys);
  int tmp_color = kBlack;
  grPhenixLowPt->SetLineColor(tmp_color);
  grPhenixLowPt->SetMarkerStyle(24);
  grPhenixLowPt->SetMarkerSize(1.5);
  grPhenixLowPt->SetMarkerColor(tmp_color);
  grPhenixLowPtSys->SetMarkerColor(tmp_color);
  grPhenixLowPtSys->SetLineColor(tmp_color);
  grPhenixLowPtSys->SetFillStyle(0);
  TBox *bPhenixLowPt = new TBox(355,1-0.12,360,1+0.12);
  bPhenixLowPt->SetLineColor(kGray+3);
  bPhenixLowPt->SetFillColor(kGray+3);
  bPhenixLowPt->SetFillStyle(1001);

  // ALICE Low pT
  // http://hepdata.cedar.ac.uk/view/ins1263062
  // Phys. Lett. B 734 (2014) 314-327
  const int nAlice_lowpT = 3;
  double aliceLowPt_npart[nAlice_lowpT] = {356.0, 191.5, 37.9};
  double aliceLowPt_npart_stat[nAlice_lowpT] = {0, 0, 0};
  double aliceLowPt_npart_sys[nAlice_lowpT] = {5, 5, 5};
  double aliceLowPt_raa[nAlice_lowpT] = {0.73, 0.70, 0.79};
  double aliceLowPt_raa_stat[nAlice_lowpT] = {0.09, 0.08, 0.15};
  double aliceLowPt_raa_sys[nAlice_lowpT] = {0.06, 0.05, 0.09};
  TGraphErrors *grAliceLowPt = new TGraphErrors(nAlice_lowpT, aliceLowPt_npart, aliceLowPt_raa, aliceLowPt_npart_stat, aliceLowPt_raa_stat);
  TGraphErrors *grAliceLowPtSys  = new TGraphErrors(nAlice_lowpT, aliceLowPt_npart, aliceLowPt_raa, aliceLowPt_npart_sys, aliceLowPt_raa_sys);
  int tmp_color = kBlue;
  grAliceLowPt->SetLineColor(tmp_color);
  grAliceLowPt->SetMarkerStyle(21);
  grAliceLowPt->SetMarkerSize(1.5);
  grAliceLowPt->SetMarkerColor(tmp_color);
  grAliceLowPtSys->SetMarkerColor(tmp_color);
  grAliceLowPtSys->SetLineColor(tmp_color);
  grAliceLowPtSys->SetFillStyle(0);
  TBox *bAliceLowPt = new TBox(360,1-0.13,365,1+0.13);
  bAliceLowPt->SetLineColor(tmp_color);
  bAliceLowPt->SetFillColor(tmp_color);
  bAliceLowPt->SetFillStyle(1001);
  
  TCanvas *c3 = new TCanvas("Jpsi_raa_vs_npart_LowPt","Jpsi_raa_vs_npart_LowPt",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->SetMaximum(2.0);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  //TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  grAliceLowPtSys->Draw("samesE5");
  grAliceLowPt->Draw("samesPEZ");
  bAliceLowPt->Draw("fsame");

  grPhenixLowPtSys->Draw("samesE5");
  grPhenixLowPt->Draw("samesPEZ");
  bPhenixLowPt->Draw("fsame");

  TGraphErrors *grStarLowPt = (TGraphErrors*)raaVsNpart[0]->Clone("grStarLowPt");
  TGraphErrors *grStarLowPtSys = (TGraphErrors*)raaVsNpartSys[0]->Clone("grStarLowPtSys");
  TBox *bStarLowPt = (TBox*)globalSys[0]->Clone("bStarLowPt");
  int tmp_color = kRed;
  grStarLowPt->SetMarkerColor(tmp_color);
  grStarLowPt->SetLineColor(tmp_color);
  grStarLowPt->Draw("samesPEZ");
  grStarLowPtSys->SetMarkerColor(tmp_color);
  grStarLowPtSys->SetLineColor(tmp_color);
  grStarLowPtSys->Draw("samesE5");
  bStarLowPt->SetLineColor(tmp_color-7);
  bStarLowPt->SetFillColor(tmp_color-7);
  bStarLowPt->SetFillStyle(1001);
  bStarLowPt->Draw("fsame");

  TLegend *leg_ncoll_2 = new TLegend(0.15,0.17,0.35,0.22);
  leg_ncoll_2->SetBorderSize(0);
  leg_ncoll_2->SetFillColor(0);
  leg_ncoll_2->SetTextFont(62);
  leg_ncoll_2->SetTextSize(0.035);;
  leg_ncoll_2->AddEntry(gPubNpartErr,"STAR N_{coll} uncertainty","f");
  leg_ncoll_2->Draw();
  star->Draw();

  TLegend *leg3 = new TLegend(0.16,0.76,0.4,0.96);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->SetTextFont(62);
  leg3->SetTextSize(0.035);
  leg3->SetHeader("p_{T,J/#psi} > 0 GeV/c");
  leg3->AddEntry(grStarLowPt,"STAR: Au+Au, #sqrt{s_{NN}} = 200 GeV, J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg3->AddEntry(grPhenixLowPt,"PHENIX: Au+Au, #sqrt{s_{NN}} = 200 GeV, J/#psi#rightarrowe^{+}e^{-}, |y| < 0.35","P");
  leg3->AddEntry(grAliceLowPt,"ALICE: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, J/#psi#rightarrowe^{+}e^{-}, |y| < 0.8","P");
  leg3->Draw();

  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC.eps",run_type,run_cfg_name.Data()));
    }

  // Add in models
  TGraph *gTamuRaaVsNpartRHIC[2], *gTamuRaaVsNpartLHC[2];
  for(int i=0; i<2; i++) gTamuRaaVsNpartRHIC[i] = (TGraph*)gTamuRaaVsNpart[i]->Clone(Form("gTamuRaaVsNpart_%d",i));
  gTamuRaaVsNpartLHC[0] = new TGraph("Rootfiles/Published/Zebo/QWG2013/Jpsi_Cent/Jpsi_Rapp_mid.dat");
  gTamuRaaVsNpartLHC[1] = new TGraph("Rootfiles/Published/Zebo/QWG2013/Jpsi_Cent/Jpsi_Rapp_LHC_mid_highPt.dat");
  for(int i=0; i<2; i++)
    {
      gTamuRaaVsNpartRHIC[i]->SetLineColor(2);
      gTamuRaaVsNpartRHIC[i]->SetLineStyle(7);
      gTamuRaaVsNpartRHIC[i]->SetLineWidth(2);

      gTamuRaaVsNpartLHC[i]->SetFillStyle(1001);
      gTamuRaaVsNpartLHC[i]->SetFillColor(kCyan-9);
      gTamuRaaVsNpartLHC[i]->SetLineColor(kCyan-9);
      gTamuRaaVsNpartLHC[i]->SetLineWidth(1);
      gTamuRaaVsNpartLHC[i]->SetLineStyle(1);
    }

  TGraph *gTsuRaaVsNpartRHIC[2];
  TGraphErrors *gTsuRaaVsNpartLHC[2];
  for(int i=0; i<2; i++) gTsuRaaVsNpartRHIC[i] = (TGraph*)gTsuRaaVsNpart[i]->Clone(Form("gTsuRaaVsNpartRHIC_%d",i));
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/Mid_Y_RAA_Npart.dat");
  double tsu_lhc_lowPt_npart[26];
  double tsu_lhc_lowPt_npart_err[26];
  double tsu_lhc_lowPt_raa[26];
  double tsu_lhc_lowPt_raa_err[26];
  double tsu_lhc_lowPt_raa_low[26];
  double tsu_lhc_lowPt_raa_high[26];
  for(int i=0; i<11; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<26; i++)
    {
      ifin >> tsu_lhc_lowPt_npart[i] >> dtmp >> dtmp >> tsu_lhc_lowPt_raa_low[i] >> dtmp >> tsu_lhc_lowPt_raa_high[i];
      tsu_lhc_lowPt_npart_err[i] = 0;
      tsu_lhc_lowPt_raa[i] = (tsu_lhc_lowPt_raa_low[i] + tsu_lhc_lowPt_raa_high[i])/2;
      tsu_lhc_lowPt_raa_err[i] = (tsu_lhc_lowPt_raa_high[i] - tsu_lhc_lowPt_raa_low[i])/2;
      cout << tsu_lhc_lowPt_raa[i] << endl;
    }
  ifin.close();
  gTsuRaaVsNpartLHC[0] = new TGraphErrors(26, tsu_lhc_lowPt_npart, tsu_lhc_lowPt_raa, tsu_lhc_lowPt_npart_err, tsu_lhc_lowPt_raa_err);

  //ifin.open("Rootfiles/Published/Model/Jpsi_Raa/Mid_Y_RAA_Npart_highpt.dat");
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/cms_raa_npart");
  double tsu_lhc_highPt_npart[24];
  double tsu_lhc_highPt_npart_err[24];
  double tsu_lhc_highPt_raa[24];
  double tsu_lhc_highPt_raa_err[24];
  for(int i=0; i<4; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<24; i++)
    {
      ifin >> tsu_lhc_highPt_npart[i] >> tsu_lhc_highPt_raa[i];
      tsu_lhc_highPt_npart_err[i] = 0;
      tsu_lhc_highPt_raa_err[i] = 0;
      tsu_lhc_highPt_raa[i] = tsu_lhc_highPt_raa[i]*1.2;
    }
  ifin.close();
  
  gTsuRaaVsNpartLHC[1] = new TGraphErrors(24, tsu_lhc_highPt_npart, tsu_lhc_highPt_raa, tsu_lhc_highPt_npart_err, tsu_lhc_highPt_raa_err);

  for(int i=0; i<2; i++)
    {
      gTsuRaaVsNpartRHIC[i]->SetLineColor(2);
      gTsuRaaVsNpartRHIC[i]->SetLineStyle(1);
      gTsuRaaVsNpartRHIC[i]->SetLineWidth(2);

      gTsuRaaVsNpartLHC[i]->SetFillStyle(1001);
      gTsuRaaVsNpartLHC[i]->SetFillColor(kViolet-9);
      gTsuRaaVsNpartLHC[i]->SetLineColor(kViolet-9);
      gTsuRaaVsNpartLHC[i]->SetLineWidth(1);
      gTsuRaaVsNpartLHC[i]->SetLineStyle(1);
    }

  TCanvas *c31 = new TCanvas("Jpsi_raa_vs_npart_vs_model_LowPt","Jpsi_raa_vs_npart_vs_model_LowPt",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  line->Draw();

  gTamuRaaVsNpartLHC[0]->Draw("fsames");
  //gTamuRaaVsNpartLHC[0]->Draw("lsames");
  gTsuRaaVsNpartLHC[0]->Draw("samesE4");
  // int tmp_color = gTsuRaaVsNpartLHC[0]->GetLineColor();
  // TGraph *gTsuRaaVsNpartLHC_up = new TGraph(26,tsu_lhc_lowPt_npart,tsu_lhc_lowPt_raa_high);
  // gTsuRaaVsNpartLHC_up->SetLineColor(tmp_color);
  // gTsuRaaVsNpartLHC_up->SetLineWidth(1);
  // gTsuRaaVsNpartLHC_up->SetLineStyle(1);
  // //gTsuRaaVsNpartLHC_up->Draw("samesl");
  // TGraph *gTsuRaaVsNpartLHC_low = new TGraph(26,tsu_lhc_lowPt_npart,tsu_lhc_lowPt_raa_low);
  // gTsuRaaVsNpartLHC_low->SetLineColor(tmp_color);
  // gTsuRaaVsNpartLHC_low->SetLineWidth(1);
  // gTsuRaaVsNpartLHC_low->SetLineStyle(1);
  // //gTsuRaaVsNpartLHC_low->Draw("samesl");

  gTamuRaaVsNpartRHIC[0]->Draw("lsames");
  gTsuRaaVsNpartRHIC[0]->Draw("lsames");

  grAliceLowPtSys->Draw("samesE5");
  grAliceLowPt->Draw("samesPEZ");
  bAliceLowPt->Draw("fsame");

  grPhenixLowPtSys->Draw("samesE5");
  grPhenixLowPt->Draw("samesPEZ");
  bPhenixLowPt->Draw("fsame");
  grStarLowPt->Draw("samesPEZ");
  grStarLowPtSys->Draw("samesE5");
  bStarLowPt->Draw("fsame");
  leg_ncoll_2->Draw();
  star->Draw();
  leg3->Draw();

  tmodel->Draw();

  TLegend *leg31 = new TLegend(0.42,0.65,0.6,0.75);
  leg31->SetBorderSize(0);
  leg31->SetFillColor(0);
  leg31->SetTextFont(62);
  leg31->SetTextSize(0.035);
  leg31->AddEntry(gTsuRaaVsNpartRHIC[0],"RHIC","L");
  leg31->AddEntry(gTamuRaaVsNpartRHIC[0],"RHIC","L");
  leg31->Draw();
  TLegend *leg32 = new TLegend(0.55,0.65,0.7,0.75);
  leg32->SetBorderSize(0);
  leg32->SetFillColor(0);
  leg32->SetTextFont(62);
  leg32->SetTextSize(0.035);
  leg32->AddEntry(gTsuRaaVsNpartLHC[0],"LHC","f");
  leg32->AddEntry(gTamuRaaVsNpartLHC[0],"LHC","f");
  leg32->Draw();

  if(savePlot)
    {
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC_vs_model.pdf",run_type,run_cfg_name.Data()));
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC_vs_model.png",run_type,run_cfg_name.Data()));
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC_vs_model.jpg",run_type,run_cfg_name.Data()));
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartLowPt_vs_LHC_vs_model.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // STAR vs LHC high pT
  //==============================================
  double cms_npart[6] = {355.4, 261.4, 187.2, 130.0, 86.3, 22.1};
  double cms_npart_err[6] = {0.,0.,0.,0.,0.,0.};
  double cms_npart_sys[6] = {6,6,6,6,6,6};
  double cms_raa[6] = {0.24, 0.26, 0.31, 0.5, 0.7, 0.62};
  double cms_raa_err[6] = {0.03, 0.03, 0.04, 0.07, 0.11, 0.11};
  double cms_raa_sys[6] = {0.02, 0.02, 0.02, 0.05, 0.08, 0.10};
  TGraphErrors *cmsRaaVsNpart = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_err, cms_raa_err);
  TGraphErrors *cmsRaaVsNpartSys = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_sys, cms_raa_sys);
  int tmp_color = kBlue;
  cmsRaaVsNpart->SetMarkerStyle(21);
  cmsRaaVsNpart->SetMarkerColor(tmp_color);
  cmsRaaVsNpart->SetMarkerSize(1.5);
  cmsRaaVsNpart->SetLineColor(tmp_color);
  cmsRaaVsNpartSys->SetMarkerColor(tmp_color);
  cmsRaaVsNpartSys->SetMarkerSize(0);
  cmsRaaVsNpartSys->SetLineColor(tmp_color);
  cmsRaaVsNpartSys->SetFillStyle(0);
  TBox *bCmsHighPt = new TBox(360,1-0.082,365,1+0.082);
  bCmsHighPt->SetLineColor(tmp_color);
  bCmsHighPt->SetFillColor(tmp_color);
  bCmsHighPt->SetFillStyle(1001);

  TGraphErrors *grStarHighPt = (TGraphErrors*)raaVsNpart[1]->Clone("grStarHighPt");
  TGraphErrors *grStarHighPtSys = (TGraphErrors*)raaVsNpartSys[1]->Clone("grStarHighPtSys");
  TBox *bStarHighPt = (TBox*)globalSys[1]->Clone("bStarHighPt");
  int tmp_color = kRed;
  grStarHighPt->SetMarkerColor(tmp_color);
  grStarHighPt->SetLineColor(tmp_color);
  grStarHighPtSys->SetMarkerColor(tmp_color);
  grStarHighPtSys->SetLineColor(tmp_color);
  bStarHighPt->SetLineColor(tmp_color);
  bStarHighPt->SetFillColor(tmp_color);

  TCanvas *c4 = new TCanvas("Jpsi_raa_vs_npart_HighPt","Jpsi_raa_vs_npart_HighPt",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->SetMaximum(2);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  //TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  cmsRaaVsNpartSys->Draw("samesE5");
  cmsRaaVsNpart->Draw("samesPEZ");
  bCmsHighPt->Draw("fsame");

  grStarHighPtSys->Draw("samesE5");
  grStarHighPt->Draw("samesPEZ");
  bStarHighPt->Draw("fsame");

  leg_ncoll_2->Draw();

  TLegend *leg4 = new TLegend(0.16,0.8,0.4,0.96);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(0);
  leg4->SetTextFont(62);
  leg4->SetTextSize(0.035);
  leg4->SetHeader("J/#psi#rightarrow#mu^{+}#mu^{-}");
  leg4->AddEntry(grStarHighPt,"STAR: Au+Au, #sqrt{s_{NN}} = 200 GeV, |y| < 0.5, p_{T} > 5 GeV/c","P");
  leg4->AddEntry(cmsRaaVsNpart,"CMS: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 2.4, p_{T} > 6.5 GeV/c","P");
  leg4->Draw();
  // TPaveText *star_2 = GetPaveText(0.75,0.85,0.68,0.73,26,23);
  // star_2->AddText("STAR preliminary");
  // star_2->SetTextColor(2);
  // star_2->Draw();
  star->Draw();
  if(savePlot)
    {
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC.pdf",run_type,run_cfg_name.Data()));
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC.png",run_type,run_cfg_name.Data()));
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC.jpg",run_type,run_cfg_name.Data()));
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC.eps",run_type,run_cfg_name.Data()));
    }

  // add in models
  TCanvas *c41 = new TCanvas("Jpsi_raa_vs_npart_vs_model_HighPt","Jpsi_raa_vs_npart_vs_model_HighPt",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  line->Draw();
  leg_ncoll_2->Draw();


  gTamuRaaVsNpartLHC[1]->SetLineColor(4);
  gTamuRaaVsNpartLHC[1]->SetLineWidth(2);
  gTamuRaaVsNpartLHC[1]->SetLineStyle(7);
  gTsuRaaVsNpartLHC[1]->SetLineWidth(2);
  gTsuRaaVsNpartLHC[1]->SetLineColor(4);
  gTsuRaaVsNpartLHC[1]->SetLineStyle(1);
  // version 2
  const int ver2 = 0;
  if(ver2)
    {
      int color2 = kMagenta+2;
      gTamuRaaVsNpartLHC[1]->SetLineColor(color2);
      gTsuRaaVsNpartLHC[1]->SetLineColor(color2);
      cmsRaaVsNpartSys->SetLineColor(color2);
      cmsRaaVsNpart->SetLineColor(color2);
      cmsRaaVsNpart->SetMarkerColor(color2);
      bCmsHighPt->SetLineColor(color2);
      bCmsHighPt->SetFillColor(color2);
    }



  gTamuRaaVsNpartLHC[1]->Draw("lsames");
  gTsuRaaVsNpartLHC[1]->Draw("lsames");

  gTamuRaaVsNpartRHIC[1]->Draw("lsames");
  gTsuRaaVsNpartRHIC[1]->Draw("lsames");


  cmsRaaVsNpartSys->Draw("samesE5");
  cmsRaaVsNpart->Draw("samesPEZ");
  bCmsHighPt->Draw("fsame");

  grStarHighPtSys->Draw("samesE5");
  grStarHighPt->Draw("samesPEZ");
  bStarHighPt->Draw("fsame");

  leg4->Draw();
  star->Draw();

  TPaveText *tmodel2 = GetPaveText(0.213,0.35,0.68,0.78);
  tmodel2->SetTextFont(62);
  tmodel2->SetTextAlign(11);
  tmodel2->SetTextSize(0.035);
  tmodel2->AddText("Tsinghua Model");
  tmodel2->AddText("TAMU Model");
  tmodel2->Draw();

  TLegend *leg41 = new TLegend(0.42,0.69,0.6,0.79);
  leg41->SetBorderSize(0);
  leg41->SetFillColor(0);
  leg41->SetTextFont(62);
  leg41->SetTextSize(0.035);
  leg41->AddEntry(gTsuRaaVsNpartRHIC[1],"RHIC","L");
  leg41->AddEntry(gTamuRaaVsNpartRHIC[1],"RHIC","L");
  leg41->Draw();
  TLegend *leg42 = new TLegend(0.55,0.69,0.7,0.79);
  leg42->SetBorderSize(0);
  leg42->SetFillColor(0);
  leg42->SetTextFont(62);
  leg42->SetTextSize(0.035);
  leg42->AddEntry(gTsuRaaVsNpartLHC[1],"LHC","L");
  leg42->AddEntry(gTamuRaaVsNpartLHC[1],"LHC","L");
  leg42->Draw();

  if(savePlot)
    {
      if(ver2==0)
	{
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model.pdf",run_type,run_cfg_name.Data()));
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model.png",run_type,run_cfg_name.Data()));
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model.jpg",run_type,run_cfg_name.Data()));
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model.eps",run_type,run_cfg_name.Data()));
	}
      else
	{
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model_2.pdf",run_type,run_cfg_name.Data()));
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model_2.png",run_type,run_cfg_name.Data()));
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model_2.jpg",run_type,run_cfg_name.Data()));
	  c41->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaNpartHighPt_vs_LHC_vs_model_2.eps",run_type,run_cfg_name.Data()));
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("2016sQM.JpsiRaaVsNpart.root","recreate");
      hpp->Write("pp200_Jpsi_Integrated");
      for(int i=0; i<2; i++)
	{
	  raaVsNpart[i]->Write(Form("MTD_Run14AuAu_JpsiRaaVsNpart%dGeV",i*5));
	  raaVsNpartSys[i]->Write(Form("MTD_Run14AuAu_JpsiRaaVsNpart%dGeV_Sys",i*5));
	}
      pubRaaVsNpart->Write("STAR_JpsiRaaVsNpart5GeV_dielec");
      pubRaaVsNpartSys->Write("STAR_JpsiRaaVsNpart5GeV_dielec_Sys");
      grPhenixLowPt->Write("PHENIX_JpsiRaaVsNpart_LowPt");
      grPhenixLowPtSys->Write("PHENIX_JpsiRaaVsNpart_LowPt_Sys");
      grAliceLowPt->Write("ALICE_JpsiRaaVsNpart_LowPt");
      grAliceLowPtSys->Write("ALICE_JpsiRaaVsNpart_LowPt_Sys");
      cmsRaaVsNpart->Write("CMS_JpsiRaaVsNpart_HighPt");
      cmsRaaVsNpartSys->Write("CMS_JpsiRaaVsNpart_HighPt_Sys");

      gTsuRaaVsNpart[0]->Write("Tsinghua_RHIC_JpsiRaaVsNpart_LowPt");
      gTsuRaaVsNpart[1]->Write("Tsinghua_RHIC_JpsiRaaVsNpart_HighPt");
      gTsuRaaVsNpartLHC[0]->Write("Tsinghua_LHC_JpsiRaaVsNpart_LowPt");
      gTsuRaaVsNpartLHC[1]->Write("Tsinghua_LHC_JpsiRaaVsNpart_HighPt");
      gTamuRaaVsNpart[0]->Write("TAMU_RHIC_JpsiRaaVsNpart_LowPt");
      gTamuRaaVsNpart[1]->Write("TAMU_RHIC_JpsiRaaVsNpart_HighPt"); 
      gTamuRaaVsNpartLHC[0]->Write("TAMU_LHC_JpsiRaaVsNpart_LowPt");
      gTamuRaaVsNpartLHC[1]->Write("TAMU_LHC_JpsiRaaVsNpart_HighPt");
    }
}


//================================================
void xsec(const bool savePlot = 0, const bool saveHisto = 1)
{
  gStyle->SetOptStat(0);
  const int nCentBins = 4;
  const double ncoll[nCentBins] = {393., 785., 300., 95.}; 
  const double ncollErr[nCentBins] = {27., 29., 31., 21.};
  const double npart[nCentBins] = {161, 280, 142, 62};
  const int marker_style[nCentBins] = {kFullCircle, kFullStar, kFullSquare, kFullCross};
  const double marker_size[nCentBins] = {1.5,2,1.5,2};
  const int marker_color[nCentBins] = {1,2,4,6};
  const int pub_marker_style[nCentBins] = {kOpenCircle, kOpenStar, kOpenSquare, kOpenCross};
  const int pub_marker_color_low[nCentBins] = {kGray+2, kRed+2, kBlue+2, kMagenta+2};
  const int pub_marker_color_high[nCentBins] = {1,2,4,6};
  const double scale_factor[nCentBins] = {10,1,0.2,0.1};
  const int lowpt_color = 9;
  const int highpt_color = 12;

  // TBW fits
  TH1F *hTBW[nCentBins];
  TFile *fmodel = TFile::Open("Rootfiles/models.root","read");
  for(int k=0; k<nCentBins; k++)
    {
      hTBW[k] = (TH1F*)fmodel->Get(Form("TBW_JpsiInvYield_AuAu200_cent%s",cent_Title[k]));
    }

  // published data
  TFile *fpub = TFile::Open(Form("Rootfiles/%s/Publication.Jpsi.200GeV.root",run_cfg_name.Data()),"read");
  TGraphAsymmErrors *gAuAuLowPt[nCentBins];
  TGraphAsymmErrors *gAuAuLowPtSys[nCentBins];
  TGraphAsymmErrors *gAuAuHighPt[nCentBins];
  TGraphAsymmErrors *gAuAuHighPtSys[nCentBins];
  TGraphAsymmErrors *gRaaLowPt[nCentBins];
  TGraphAsymmErrors *gRaaLowPtSys[nCentBins];
  TGraphAsymmErrors *gRaaHighPt[nCentBins];
  TGraphAsymmErrors *gRaaHighPtSys[nCentBins];
  double x,y;
  for(int k=0; k<nCentBins; k++)
    {
      gAuAuLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[k]));
      gAuAuLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[k]));
      gAuAuHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[k]));
      gAuAuHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[k]));
      scaleGraph(gAuAuLowPt[k], scale_factor[k]);
      scaleGraph(gAuAuLowPtSys[k], scale_factor[k]);
      scaleGraph(gAuAuHighPt[k], scale_factor[k]);
      scaleGraph(gAuAuHighPtSys[k], scale_factor[k]);
      hTBW[k]->Scale(scale_factor[k]);

      gRaaLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",cent_Title[k]));
      gRaaLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",cent_Title[k]));
      gRaaHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",cent_Title[k]));
      gRaaHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",cent_Title[k]));
    }

  TH1F *hJpsipp = (TH1F*)fpub->Get("hPPJpsiFinal");
  // data
  TFile *fSys = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hSys = (TH1F*)fSys->Get(Form("Sys_all_%s",cent_Title[0]));
  TFile *fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiInvYield[nCentBins];
  TGraphErrors *hJpsiXsec[nCentBins];
  TGraphErrors *hJpsiXsecSys[nCentBins];
  TGraphErrors *hJpsiRaa[nCentBins];
  TGraphErrors *hJpsiRaaSys[nCentBins];
  TGraphErrors *hJpsiRaaSys2[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)fdata->Get(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
      int npoints = hJpsiInvYield[k]->GetNbinsX();
      hJpsiXsec[k] = new TGraphErrors(npoints);
      hJpsiXsec[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiXsecSys[k] = new TGraphErrors(npoints);
      hJpsiXsecSys[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));
      for(int i=0; i<npoints; i++)
	{
	  if(k>1 && i==npoints-1) continue;
	  x = hJpsiInvYield[k]->GetBinCenter(i+1);
	  y = hJpsiInvYield[k]->GetBinContent(i+1);
	  hJpsiXsec[k]->SetPoint(i,x,y*scale_factor[k]);
	  hJpsiXsec[k]->SetPointError(i,0,hJpsiInvYield[k]->GetBinError(i+1)*scale_factor[k]);

	  hJpsiXsecSys[k]->SetPoint(i,x,y*scale_factor[k]);
	  hJpsiXsecSys[k]->SetPointError(i,0.2,hSys->GetBinContent(i+1)*y*scale_factor[k]);
	}

      hJpsiRaa[k] = new TGraphErrors();
      hJpsiRaa[k]->SetName(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
      hJpsiRaaSys[k] = new TGraphErrors();
      hJpsiRaaSys[k]->SetName(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
      hJpsiRaaSys2[k] = new TGraphErrors();
      hJpsiRaaSys2[k]->SetName(Form("Graph_Jpsi_Raa_cent%s_sys_2",cent_Title[k]));
      for(int i=0; i<npoints; i++)
	{
	  if(k>1 && i==npoints-1) continue;
	  hJpsiXsec[k]->GetPoint(i,x,y);
	  double AuAu_val = y/scale_factor[k];
	  double AuAu_err = hJpsiXsec[k]->GetErrorY(i)/scale_factor[k];
	  double AuAu_sys = hSys->GetBinContent(i+1) * AuAu_val;
	  double pp_val = hJpsipp->GetBinContent(i+1);
	  double pp_err = hJpsipp->GetBinError(i+1);
	  double prefix = ppInelastic/ncoll[k] * 1e6;
	  double val = prefix * AuAu_val / pp_val;
	  double err = prefix * AuAu_err / pp_val;
	  double sys = prefix * AuAu_sys / pp_val;
	  double sys2 = val * sqrt(AuAu_sys*AuAu_sys/AuAu_val/AuAu_val + pp_err*pp_err/pp_val/pp_val);
	  hJpsiRaa[k]->SetPoint(i,x,val);
	  hJpsiRaa[k]->SetPointError(i,0,err);
	  //printf("%s: raa -> %2.2f +/- %2.2f%%, AuAu -> %2.2e +/- %2.2f%%, pp -> %2.2e +/-%2.2f%%\n",pt_Name[i+1],val,err/val*100,AuAu_val,AuAu_err/AuAu_val*100,pp_val,pp_err/pp_val*100);
	  //printf("%s: raa -> %2.2f +/- %2.2f%%, AuAu -> %2.2e +/- %2.2f%%, pp -> %2.2e +/-%2.2f%%\n",pt_Name[i+1],val,sys2/val*100,AuAu_val,AuAu_sys/AuAu_val*100,pp_val,pp_err/pp_val*100);

	  hJpsiRaaSys[k]->SetPoint(i,x,val);
	  hJpsiRaaSys[k]->SetPointError(i,0.2,sys2);

	  hJpsiRaaSys2[k]->SetPoint(i,x,val);
	  hJpsiRaaSys2[k]->SetPointError(i,0.2,sys2);
	}
    }


  TCanvas *c1 = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",800,700);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);B_{ll}d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{-2}]",15,0,15);
  hAuAu->GetYaxis()->SetRangeUser(1e-13,5e-3);
  ScaleHistoTitle(hAuAu,0.05,1,0.04,0.05,1.2,0.04,62);
  gPad->SetLogy();
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hAuAu->Draw();
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiXsec[k]->SetMarkerStyle(marker_style[k]);
      hJpsiXsec[k]->SetMarkerColor(marker_color[k]);
      hJpsiXsec[k]->SetLineColor(marker_color[k]);
      hJpsiXsec[k]->SetMarkerSize(marker_size[k]+0.5);
      hJpsiXsecSys[k]->SetFillStyle(0);
      hJpsiXsecSys[k]->SetLineColor(hJpsiXsec[k]->GetLineColor());
      hJpsiXsecSys[k]->Draw("samesE5");
      hJpsiXsec[k]->Draw("samesPEZ");
    }
  TLegend *leg = new TLegend(0.7,0.65,0.9,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hJpsiXsec[0],"0-60%#times10","P");
  leg->AddEntry(hJpsiXsec[1],"0-20%","P");
  leg->AddEntry(hJpsiXsec[2],"20-40%/5","P");
  leg->AddEntry(hJpsiXsec[3],"40-60%/10","P");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.15,0.3,0.85,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  t1->Draw();
  TPaveText *star = GetPaveText(0.25,0.5,0.25,0.3,0.04);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c2 = new TCanvas("AuAu200_Jpsi_vs_pub","AuAu200_Jpsi_vs_pub",800,700);
  gPad->SetLogy();
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hAuAu->Draw();
  TGraphErrors *gJpsi[nCentBins], *gJpsiSys[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      gAuAuLowPt[k]->SetMarkerStyle(pub_marker_style[k]);
      gAuAuLowPt[k]->SetMarkerSize(marker_size[k]);
      gAuAuLowPt[k]->SetMarkerColor(pub_marker_color_low[k]);
      gAuAuLowPt[k]->SetLineColor(pub_marker_color_low[k]);
      gAuAuLowPt[k]->Draw("sames PE");
      gAuAuLowPtSys[k]->SetMarkerColor(pub_marker_color_low[k]);
      gAuAuLowPtSys[k]->SetLineColor(pub_marker_color_low[k]);
      gAuAuLowPtSys[k]->Draw("sameE5");
      gAuAuHighPt[k]->SetMarkerStyle(pub_marker_style[k]);
      gAuAuHighPt[k]->SetMarkerSize(marker_size[k]);
      gAuAuHighPt[k]->SetMarkerColor(pub_marker_color_high[k]);
      gAuAuHighPt[k]->SetLineColor(pub_marker_color_high[k]);
      gAuAuHighPt[k]->Draw("sames PE");
      gAuAuHighPtSys[k]->SetMarkerColor(pub_marker_color_high[k]);
      gAuAuHighPtSys[k]->SetLineColor(pub_marker_color_high[k]);
      gAuAuHighPtSys[k]->Draw("sameE5");
      hTBW[k]->Draw("sames");
      hJpsiXsecSys[k]->Draw("samesE5");
      hJpsiXsec[k]->Draw("samesP");
    }
  TLegend *leg = new TLegend(0.7,0.65,0.9,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hJpsiXsec[0],"0-60%#times10","P");
  leg->AddEntry(hJpsiXsec[1],"0-20%","P");
  leg->AddEntry(hJpsiXsec[2],"20-40%/5","P");
  leg->AddEntry(hJpsiXsec[3],"40-60%/10","P");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.15,0.3,0.9,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV"));
  t1->Draw();

  TPaveText *t1 = GetPaveText(0.29,0.39,0.3,0.35,24,63);
  t1->AddText("Open: J/#psi#rightarrowe^{+}e^{-}, |y| < 1");
  t1->Draw();

  TPaveText *t1 = GetPaveText(0.31,0.41,0.25,0.3,24,63);
  t1->AddText("Filled: J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  t1->Draw();

  TLegend *leg = new TLegend(0.16,0.2,0.3,0.25);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(24);
  leg->SetTextFont(63);
  leg->AddEntry(hTBW[0],"TBW fit (#beta=0)","L");
  leg->Draw();
  

  TPaveText *star = GetPaveText(0.4,0.6,0.8,0.85,0.04);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // RAA plot: dimuon
  //==============================================
  TCanvas *c30 = new TCanvas("Raa_Jpsi_dimuon","Raa_Jpsi_dimuon",1100,700);
  TPad *pad = GetSinglePad("jet_trigger",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();

  TH1F *hRaa = new TH1F("Raa_Jpsi",";p_{T} (GeV/c);R_{AA}",10,0,15);
  hRaa->GetYaxis()->SetRangeUser(0.05,1.95);
  hRaa->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hRaa,24,1.9,20,24,1.9,20,63);
  TBox *box_ncoll[nCentBins];
  TPaveText *centLabel[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();

      hJpsiRaa[k]->SetMarkerStyle(29);
      hJpsiRaa[k]->SetMarkerColor(2);
      hJpsiRaa[k]->SetLineColor(2);
      hJpsiRaa[k]->SetMarkerSize(2.5);
      hJpsiRaaSys[k]->SetFillStyle(0);
      hJpsiRaaSys[k]->SetLineColor(hJpsiRaa[k]->GetLineColor());
      hJpsiRaaSys2[k]->SetLineColor(hJpsiRaa[k]->GetLineColor());

      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      if(k%2==0) centLabel[k] = GetPaveText(0.85,0.9,0.9,0.95);
      else       centLabel[k] = GetPaveText(0.82,0.87,0.9,0.95);
      centLabel[k]->SetTextFont(63);
      centLabel[k]->SetTextSize(24);
      centLabel[k]->AddText(Form("%s%%",cent_Name[k]));
      centLabel[k]->Draw();

      // Global systematics
      double gerr = sqrt(pow(ncollErr[k]/ncoll[k],2)+pow(ppInelasticErr/ppInelastic,2));
      box_ncoll[k] = new TBox(14.5,1-gerr,14.75,1+gerr);
      box_ncoll[k]->SetFillStyle(1001);
      box_ncoll[k]->SetFillColor(kRed-7);
      box_ncoll[k]->Draw();
    }
  pad->cd(1);
  TLegend *leg30 = new TLegend(0.2,0.78,0.45,0.97);
  leg30->SetBorderSize(0);
  leg30->SetFillColor(0);
  leg30->SetTextFont(63);
  leg30->SetTextSize(20);
  leg30->SetHeader("STAR Au+Au @ 200 GeV");
  leg30->AddEntry(hJpsiRaa[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg30->Draw();

  pad->cd(2);
  TPaveText *star = GetPaveText(0.2,0.4,0.85,0.9,24,23);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c30->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.pdf",run_type,run_cfg_name.Data()));
      c30->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.png",run_type,run_cfg_name.Data()));
      c30->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.jpg",run_type,run_cfg_name.Data()));
      c30->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // RAA plot: dimuon vs dielectron
  //==============================================
  TCanvas *c31 = new TCanvas("Raa_Jpsi_vs_pub","Raa_Jpsi_vs_pub",1100,700);
  TPad *pad = GetSinglePad("jet_trigger",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();

  TH1F *hRaa = new TH1F("Raa_Jpsi",";p_{T} (GeV/c);R_{AA}",10,0,15);
  hRaa->GetYaxis()->SetRangeUser(0.05,1.95);
  hRaa->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hRaa,24,1.9,20,24,1.9,20,63);
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();
      gRaaLowPt[k]->SetMarkerStyle(kOpenStar);
      gRaaLowPt[k]->SetMarkerSize(2);
      gRaaLowPt[k]->SetMarkerColor(lowpt_color);
      gRaaLowPt[k]->SetLineColor(lowpt_color);
      gRaaLowPtSys[k]->SetMarkerColor(lowpt_color);
      gRaaLowPtSys[k]->SetLineColor(lowpt_color);

      gRaaHighPt[k]->SetMarkerStyle(kOpenStar);
      gRaaHighPt[k]->SetMarkerSize(2);
      gRaaHighPt[k]->SetMarkerColor(highpt_color);
      gRaaHighPt[k]->SetLineColor(highpt_color);
      gRaaHighPtSys[k]->SetMarkerColor(highpt_color);
      gRaaHighPtSys[k]->SetLineColor(highpt_color);

      gRaaLowPtSys[k]->Draw("sameE5");
      gRaaHighPtSys[k]->Draw("sameE5");
      gRaaLowPt[k]->Draw("sames PEZ");
      gRaaHighPt[k]->Draw("sames PEZ");
      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");
      centLabel[k]->Draw();

      // Global systematics
      box_ncoll[k]->Draw();
    }
  pad->cd(1);
  TLegend *leg31 = new TLegend(0.2,0.6,0.45,0.78);
  leg31->SetBorderSize(0);
  leg31->SetFillColor(0);
  leg31->SetTextFont(63);
  leg31->SetTextSize(20);
  leg31->AddEntry(gRaaLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y| < 1 (MB)","P");
  leg31->AddEntry(gRaaHighPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y| < 1 (HT)","P");
  leg31->Draw();
  leg30->Draw();

  pad->cd(2);
  star->Draw();
  if(savePlot)
    {
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.pdf",run_type,run_cfg_name.Data()));
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.png",run_type,run_cfg_name.Data()));
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.jpg",run_type,run_cfg_name.Data()));
      c31->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  //RAA plot: dimuon vs dielectron vs model
  //==============================================
  ifstream ifin;
  char tmp[256];
  double dtmp;

  // model calculations
  // Tsinghua group
  ifin.open("Rootfiles/Published/Jpsi_Raa_200/model/data_and_figures_from_Yunpeng/data_Raa_Pt.dat");
  double pt_thu[11];
  double v2_init[4][11];
  double v2_gen[4][11];
  double v2_all[4][11];
  ifin.getline(tmp,256);
  cout << tmp << endl;
  ifin.getline(tmp,256);
  cout << tmp << endl;
  for(int i=0; i<11; i++)
    {
      ifin >> pt_thu[i] 
	   >> v2_init[1][i] >> v2_gen[1][i] >> v2_all[1][i]
	   >> v2_init[2][i] >> v2_gen[2][i] >> v2_all[2][i]
	   >> v2_init[3][i] >> v2_gen[3][i] >> v2_all[3][i]
	   >> v2_init[0][i] >> v2_gen[0][i] >> v2_all[0][i];
    }
  TGraph *grTHU[nCentBins]; 
  for(int k=0; k<nCentBins; k++)
    {
      grTHU[k] = new TGraph(11,pt_thu,v2_all[k]);
      grTHU[k]->SetName(Form("THU_v2_all_cent%s",cent_Title[k]));
      grTHU[k]->SetLineColor(1);
      grTHU[k]->SetLineWidth(2);
      grTHU[k]->SetLineStyle(1);
    }
  ifin.close();

  // TAMU group
  TGraph *grTAMU[nCentBins];
  const char *fname[nCentBins] = {"060","020","2040","4060"};
  double pt_tamu[51];
  double v2_tamu[51];
  for(int k=0; k<nCentBins; k++)
    {
      ifin.open(Form("Rootfiles/Published/Jpsi_Raa_200/model/data_from_Xingbo/raa_pt_rhic_%s_tozebo_110507.dat",fname[k]));
      for(int i=0; i<51; i++)
	{
	  ifin >> pt_tamu[i] >> v2_tamu[i];
	  //cout << pt_tamu[i] << endl;
	}
      grTAMU[k] = new TGraph(51,pt_tamu,v2_tamu);
      grTAMU[k]->SetName(Form("TAMU_v2_all_cent%s",cent_Title[k]));
      grTAMU[k]->SetLineColor(1);
      grTAMU[k]->SetLineWidth(2);
      grTAMU[k]->SetLineStyle(7);
      ifin.close();
    }

  TCanvas *c32 = new TCanvas("JpsiRaaVspt_dimuon_vs_dielectron_vs_model","JpsiRaaVspt_dimuon_vs_dielectron_vs_model",1100,700);
  TPad *pad = GetSinglePad("c32",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();

      grTAMU[k]->Draw("samesC");
      grTHU[k]->Draw("samesC");

      gRaaLowPtSys[k]->Draw("sameE5");
      gRaaHighPtSys[k]->Draw("sameE5");
      gRaaLowPt[k]->Draw("sames PEZ");
      gRaaHighPt[k]->Draw("sames PEZ");
      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");
      centLabel[k]->Draw();
      box_ncoll[k]->Draw();
    }
  pad->cd(1);
  leg30->Draw();
  leg31->Draw();
  pad->cd(3);
  TLegend *leg32 = new TLegend(0.2,0.7,0.45,0.85);
  leg32->SetBorderSize(0);
  leg32->SetFillColor(0);
  leg32->SetTextFont(63);
  leg32->SetTextSize(20);
  leg32->AddEntry(grTHU[0],"Transport Model I (Y.-P. Liu #it{et al.})","L");
  leg32->AddEntry(grTAMU[0],"Transport Model II (X. Zhao #it{et al.})","L");
  leg32->Draw();
  pad->cd(2);
  star->Draw();
  if(savePlot)
    {
      c32->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_dielectron_vs_Model.pdf",run_type,run_cfg_name.Data()));
      c32->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_dielectron_vs_Model.png",run_type,run_cfg_name.Data()));
      c32->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_dielectron_vs_Model.jpg",run_type,run_cfg_name.Data()));
      c32->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_dielectron_vs_Model.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // RAA plot: STAR vs PHENIX 
  //==============================================
  double phenix_pt[4][5];
  double phenix_pt_err[4][5];
  double phenix_pt_sys[4][5];
  double phenix_raa[4][5];
  double phenix_raa_stat[4][5];
  double phenix_raa_sys[4][5];
  const char *phenix_name[4] = {"060","020","2040","4060"};
  TGraphErrors *gPhenixRaa[4];
  TGraphErrors *gPhenixRaaSys[4];
  TBox *box_phenix[nCentBins];
  for(int k=0; k<4; k++)
    {
      ifin.open(Form("Rootfiles/Published/Jpsi_Raa_200/model/raa_phenix_%s.dat",phenix_name[k]));
      for(int j=0; j<5; j++)
	{
	  ifin >> dtmp >> dtmp >> dtmp >> dtmp >> phenix_pt[k][j] >> dtmp >> phenix_raa[k][j] >> phenix_raa_stat[k][j] >> phenix_raa_sys[k][j] >> dtmp;
	  phenix_pt_err[k][j] = 0;
	  phenix_pt_sys[k][j] = 0.2;
	}
      ifin.close();
      gPhenixRaa[k] = new TGraphErrors(5,phenix_pt[k], phenix_raa[k],phenix_pt_err[k],phenix_raa_stat[k]);
      gPhenixRaaSys[k] = new TGraphErrors(5,phenix_pt[k], phenix_raa[k],phenix_pt_sys[k],phenix_raa_sys[k]);
      gPhenixRaa[k]->SetMarkerStyle(20);
      gPhenixRaa[k]->SetMarkerSize(2);
      gPhenixRaa[k]->SetMarkerColor(kGreen+2);
      gPhenixRaa[k]->SetLineColor(kGreen+2);
      gPhenixRaaSys[k]->SetMarkerColor(kGreen+2);
      gPhenixRaaSys[k]->SetLineColor(kGreen+2);
      gPhenixRaaSys[k]->SetFillStyle(0);

      double gerr = 0.1;
      if(k==3) gerr = 0.13;
      box_phenix[k] = new TBox(14.25,1-gerr,14.5,1+gerr);
      box_phenix[k]->SetFillStyle(1001);
      box_phenix[k]->SetFillColor(gPhenixRaa[k]->GetMarkerColor());
    }
  
  TCanvas *c4 = new TCanvas("JpsiRaaVspt_STAR_vs_PHENIX","JpsiRaaVspt_STAR_vs_PHENIX",1100,700);
  TPad *pad = GetSinglePad("c4",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();

      gPhenixRaaSys[k]->Draw("sameE5");
      gPhenixRaa[k]->Draw("sames PEZ");

      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      centLabel[k]->Draw();
      box_ncoll[k]->Draw();
      if(k>0) box_phenix[k]->Draw();
    }
  pad->cd(1);
  TLegend *leg4 = new TLegend(0.2,0.72,0.45,0.97);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(0);
  leg4->SetTextFont(63);
  leg4->SetTextSize(20);
  leg4->SetHeader("Au+Au @ 200 GeV");
  leg4->AddEntry(hJpsiRaa[0],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg4->AddEntry(gPhenixRaa[0],"PHENIX: J/#psi#rightarrowe^{+}e^{-}, |y| < 0.35","P");
  leg4->Draw();

  pad->cd(2);
  star->Draw();

  if(savePlot)
    {
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX.pdf",run_type,run_cfg_name.Data()));
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX.png",run_type,run_cfg_name.Data()));
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX.jpg",run_type,run_cfg_name.Data()));
      c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // Raa vs pT: STAR vs PHENIX vs model
  //==============================================
  TCanvas *c5 = new TCanvas("JpsiRaaVspt_STAR_vs_PHENIX_vs_model","JpsiRaaVspt_STAR_vs_PHENIX_vs_model",1100,700);
  TPad *pad = GetSinglePad("c5",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();

      grTAMU[k]->Draw("samesC");
      grTHU[k]->Draw("samesC");

      gPhenixRaaSys[k]->Draw("sameE5");
      gPhenixRaa[k]->Draw("sames PEZ");

      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      centLabel[k]->Draw();
      box_ncoll[k]->Draw();
      if(k>0) box_phenix[k]->Draw();
    }
  pad->cd(1);
  leg4->Draw();
  TLegend *leg5 = new TLegend(0.2,0.55,0.45,0.72);
  leg5->SetBorderSize(0);
  leg5->SetFillColor(0);
  leg5->SetTextFont(63);
  leg5->SetTextSize(20);
  leg5->AddEntry(grTHU[0],"Transport Model I (Y.-P. Liu #it{et al.})","L");
  leg5->AddEntry(grTAMU[0],"Transport Model II (X. Zhao #it{et al.})","L");
  leg5->Draw();
  pad->cd(2);
  star->Draw();
  if(savePlot)
    {
      c5->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_vs_Model.pdf",run_type,run_cfg_name.Data()));
      c5->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_vs_Model.png",run_type,run_cfg_name.Data()));
      c5->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_vs_Model.jpg",run_type,run_cfg_name.Data()));
      c5->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_vs_Model.eps",run_type,run_cfg_name.Data()));
    }

  //==============================================
  // Raa vs pT: STAR mid-rapidity vs PHENIX forward
  //==============================================
  ifin.open("Rootfiles/2016sQM/Publications/phenix_raa_vs_pt_forward.txt");
  double phenix_forward_pt[5];
  double phenix_forward_pt_err[5];
  double phenix_forward_raa[5];
  double phenix_forward_raa_err[5];
  double phenix_forward_raa_sys1[5];
  double phenix_forward_raa_stat[5];
  double phenix_forward_raa_sys2[5];
  double phenix_forward_pt_err2[5];
  ifin.getline(tmp,256);
  cout << tmp << endl;
  ifin.getline(tmp,256);
  cout << tmp << endl;
  ifin.getline(tmp,256);
  cout << tmp << endl;

  ifin.getline(tmp,256);
  cout << tmp << endl;
  ifin.getline(tmp,256);
  cout << tmp << endl;
  for(int i=0; i<5; i++)
    {
      ifin >> phenix_forward_pt[i] 
	   >> phenix_forward_raa[i] 
	   >> phenix_forward_raa_err[i]
	   >> phenix_forward_raa_sys1[i]
	   >> phenix_forward_raa_sys2[i];

      phenix_forward_raa_stat[i] = sqrt(pow(phenix_forward_raa_err[i],2)+pow(phenix_forward_raa_sys1[i],2));
      phenix_forward_pt_err[i] = 0;
      phenix_forward_pt_err2[i] = 0.2;
    }
  TGraphErrors *grPhenixForward[nCentBins]; 
  TGraphErrors *grPhenixForwardSys[nCentBins]; 
  for(int k=0; k<1; k++)
    {
      int tmp_color = kGreen + 2;
      grPhenixForward[k] = new TGraphErrors(5,phenix_forward_pt,phenix_forward_raa,phenix_forward_pt_err,phenix_forward_raa_stat);
      grPhenixForward[k]->SetName(Form("PhenixForward_cent%s",cent_Title[k]));
      grPhenixForward[k]->SetLineColor(tmp_color);
      grPhenixForward[k]->SetMarkerStyle(20);
      grPhenixForward[k]->SetMarkerSize(2);
      grPhenixForward[k]->SetMarkerColor(tmp_color);

      grPhenixForwardSys[k] = new TGraphErrors(5,phenix_forward_pt,phenix_forward_raa,phenix_forward_pt_err2,phenix_forward_raa_sys2);
      grPhenixForwardSys[k]->SetName(Form("PhenixForwardSys_cent%s",cent_Title[k]));
      grPhenixForwardSys[k]->SetMarkerColor(tmp_color);
      grPhenixForwardSys[k]->SetLineColor(tmp_color);
      grPhenixForwardSys[k]->SetFillStyle(0);
    }
  ifin.close();
  
  TCanvas *c6 = new TCanvas("JpsiRaaVsPt_StarMid_PhenixForward","JpsiRaaVsPt_StarMid_PhenixForward",800,600);
  TH1F *hRaa2 = new TH1F("hRaa2_Jpsi",";p_{T} (GeV/c);R_{AA}",15,0,15);
  hRaa2->GetYaxis()->SetRangeUser(0.05,1.95);
  hRaa2->GetYaxis()->CenterTitle();
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  ScaleHistoTitle(hRaa2,28,1,24,28,1,24,63);
  hRaa2->SetMaximum(1.6);
  hRaa2->SetYTitle("J/#psi R_{AA}");
  hRaa2->DrawCopy();
  TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  hJpsiRaaSys[1]->Draw("sameE5");
  hJpsiRaa[1]->Draw("samesPEZ");
  grPhenixForwardSys[0]->Draw("sameE5");
  grPhenixForward[0]->Draw("samesPEZ");

  TLegend *leg5 = new TLegend(0.18,0.75,0.41,0.97);
  leg5->SetBorderSize(0);
  leg5->SetFillColor(0);
  leg5->SetTextFont(63);
  leg5->SetTextSize(22);
  leg5->SetHeader("Au+Au @ 200 GeV, 0-20%");
  leg5->AddEntry(hJpsiRaa[1],Form("STAR: J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5"),"P");
  leg5->AddEntry(grPhenixForward[0],"PHENIX: J/#psi#rightarrow#mu^{+}#mu^{-}, 1.2 < |y| < 2.2","P");
  leg5->Draw();

  TPaveText *star2 = GetPaveText(0.75,0.9,0.9,0.95,24,23);
  star2->AddText("STAR preliminary");
  star2->SetTextColor(2);
  star2->Draw();

  double gerr = sqrt(pow(ncollErr[1]/ncoll[1],2)+pow(ppInelasticErr/ppInelastic,2));
  TBox *box_star = new TBox(14.5,1-gerr,14.75,1+gerr);
  box_star->SetFillStyle(1001);
  box_star->SetFillColor(kRed-7);
  box_star->Draw();

  TBox *box_phenix_forward = new TBox(14.25,1-0.1,14.5,1+0.1);
  box_phenix_forward->SetFillStyle(1001);
  box_phenix_forward->SetFillColor(kGreen+2);
  box_phenix_forward->Draw();

  if(savePlot)
    {
      c6->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_forward.pdf",run_type,run_cfg_name.Data()));
      c6->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_forward.png",run_type,run_cfg_name.Data()));
      c6->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_forward.jpg",run_type,run_cfg_name.Data()));
      c6->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_PHENIX_forward.eps",run_type,run_cfg_name.Data()));
    }


  //==============================================
  // Raa vs pT: STAR mid-rapidity vs LHC mid-rapidity 0-40%
  //==============================================

  // ALICE JHEP 07 (2015) 051
  double alice_pt_040[2] = {1.56, 3.33};
  double alice_pt_040_stat[2] = {0, 0};
  double alice_pt_040_sys[2] = {0.3, 0.3};
  double alice_raa_040[2] = {0.82, 0.58};
  double alice_raa_040_stat[2] = {0.11, 0.06};
  double alice_raa_040_sys[2] = {0.10,0.08};
  TGraphErrors *grAlice0040 = new TGraphErrors(2, alice_pt_040, alice_raa_040, alice_pt_040_stat, alice_raa_040_stat);
  TGraphErrors *grAlice0040Sys  = new TGraphErrors(2, alice_pt_040, alice_raa_040, alice_pt_040_sys, alice_raa_040_sys);
  int tmp_color = kBlue;
  grAlice0040->SetLineColor(tmp_color);
  grAlice0040->SetMarkerStyle(21);
  grAlice0040->SetMarkerSize(2);
  grAlice0040->SetMarkerColor(tmp_color);
  grAlice0040Sys->SetMarkerColor(tmp_color);
  grAlice0040Sys->SetLineColor(tmp_color);
  grAlice0040Sys->SetFillStyle(0);

  // CMS JHEP 05 (2012) 063
  double cms_pt_040[1] = {9.35};
  double cms_pt_040_stat[1] = {0};
  double cms_pt_040_sys[1] = {0.3};
  double cms_raa_040[1] = {0.285};
  double cms_raa_040_stat[1] = {0.01};
  double cms_raa_040_sys[1] = {0.011};
  TGraphErrors *grCms0040 = new TGraphErrors(1, cms_pt_040, cms_raa_040, cms_pt_040_stat, cms_raa_040_stat);
  TGraphErrors *grCms0040Sys  = new TGraphErrors(1, cms_pt_040, cms_raa_040, cms_pt_040_sys, cms_raa_040_sys);
  int tmp_color = kBlack;
  grCms0040->SetLineColor(tmp_color);
  grCms0040->SetMarkerStyle(20);
  grCms0040->SetMarkerSize(2);
  grCms0040->SetMarkerColor(tmp_color);
  grCms0040Sys->SetMarkerColor(tmp_color);
  grCms0040Sys->SetLineColor(tmp_color);
  grCms0040Sys->SetFillStyle(0);

  // STAR 0-40%
  const int kNpoints = 9;
  double star_pt_040[kNpoints];
  double star_pt_040_stat[kNpoints];
  double star_pt_040_sys[kNpoints];
  double star_raa_040[kNpoints];
  double star_raa_040_stat[kNpoints];
  double star_raa_040_sys[kNpoints];
  double x1, y1;
  double ncoll_0040 = 543;
  double ncoll_0040_err = 30;
  for(int i=0; i<kNpoints; i++)
    {
      hJpsiXsec[1]->GetPoint(i,x,y);
      hJpsiXsec[2]->GetPoint(i,x1,y1);
      double AuAu_val = (y/scale_factor[1]+y1/scale_factor[2])/2;
      double AuAu_err = sqrt(pow(hJpsiXsec[1]->GetErrorY(i)/scale_factor[1],2)+pow(hJpsiXsec[2]->GetErrorY(i)/scale_factor[2],2))/2;
      double AuAu_sys = hSys->GetBinContent(i+1) * AuAu_val;
      double pp_val = hJpsipp->GetBinContent(i+1);
      double pp_err = hJpsipp->GetBinError(i+1);
      double prefix = ppInelastic/ncoll_0040 * 1e6;
      double val = prefix * AuAu_val / pp_val;
      double err = prefix * AuAu_err / pp_val;
      double sys = val * sqrt(AuAu_sys*AuAu_sys/AuAu_val/AuAu_val + pp_err*pp_err/pp_val/pp_val);
      star_pt_040[i] = x;
      star_pt_040_stat[i] = 0;
      star_pt_040_sys[i] = 0.2;
      star_raa_040[i] = val;
      star_raa_040_stat[i] = err;
      star_raa_040_sys[i] = sys;
      //cout << hJpsiXsec[1]->GetErrorY(i)/y << "  " << hJpsiXsec[2]->GetErrorY(i)/y1 << "  " << AuAu_err/AuAu_val << endl;
      cout << val << "  " << err << endl;
    }

  TGraphErrors *hStarRaa0040 = new TGraphErrors(kNpoints, star_pt_040, star_raa_040, star_pt_040_stat, star_raa_040_stat);
  TGraphErrors *hStarRaa0040Sys = new TGraphErrors(kNpoints, star_pt_040, star_raa_040, star_pt_040_sys, star_raa_040_sys);
  hStarRaa0040->SetLineColor(2);
  hStarRaa0040->SetMarkerStyle(29);
  hStarRaa0040->SetMarkerSize(2.5);
  hStarRaa0040->SetMarkerColor(2);
  hStarRaa0040Sys->SetMarkerColor(2);
  hStarRaa0040Sys->SetLineColor(2);
  hStarRaa0040Sys->SetFillStyle(0);

  TCanvas *c7 = new TCanvas("JpsiRaaVsPt_STAR_LHC","JpsiRaaVsPt_STAR_LHC",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  ScaleHistoTitle(hRaa2,28,1,24,28,1,24,63);
  hRaa2->GetXaxis()->SetRangeUser(0,10);
  hRaa2->SetMaximum(1.6);
  hRaa2->SetYTitle("J/#psi R_{AA}");
  hRaa2->DrawCopy();
  TLine *line = GetLine(hRaa2->GetXaxis()->GetXmin(),1,10,1,1);
  line->Draw();
  grAlice0040Sys->Draw("sameE5");
  grAlice0040->Draw("samesPEZ");
  grCms0040Sys->Draw("sameE5");
  grCms0040->Draw("samesPEZ");
  hStarRaa0040Sys->Draw("sameE5");
  hStarRaa0040->Draw("samesPEZ");

  TLegend *leg7 = new TLegend(0.18,0.7,0.41,0.95);
  leg7->SetBorderSize(0);
  leg7->SetFillColor(0);
  leg7->SetTextFont(63);
  leg7->SetTextSize(22);
  leg7->SetHeader("0-40% centrality");
  leg7->AddEntry(hStarRaa0040,"STAR:  Au+Au, #sqrt{s_{NN}} = 0.2 TeV,  |y| < 0.5","P");
  leg7->AddEntry(grAlice0040, "ALICE: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 0.8","P");
  leg7->AddEntry(grCms0040,   "CMS:   Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 2.4","P");
  leg7->Draw();

  star2->Draw();

  double gerr = sqrt(pow(ncoll_0040_err/ncoll_0040,2)+pow(ppInelasticErr/ppInelastic,2));
  TBox *box_star = new TBox(9.6,1-gerr,9.75,1+gerr);
  box_star->SetFillStyle(1001);
  box_star->SetFillColor(kRed-7);
  box_star->Draw();

  TBox *box_alice = new TBox(9.45,1-0.13,9.6,1+0.13);
  box_alice->SetFillStyle(1001);
  box_alice->SetFillColor(kBlue);
  box_alice->Draw();

  TBox *box_cms = new TBox(9.3,1-0.082,9.45,1+0.082);
  box_cms->SetFillStyle(1001);
  box_cms->SetFillColor(kGray+3);
  box_cms->Draw();

  if(savePlot)
    {
      c7->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040.pdf",run_type,run_cfg_name.Data()));
      c7->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040.png",run_type,run_cfg_name.Data()));
      c7->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040.jpg",run_type,run_cfg_name.Data()));
      c7->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040.eps",run_type,run_cfg_name.Data()));
    }


  //==============================================
  // Raa vs pT: STAR mid-rapidity vs LHC mid-rapidity 0-40% vs model
  //==============================================

  // Tsinghua model
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/MidYcharmoniumRAApT");
  double tsu_lhc_pt[21];
  double tsu_lhc_pt_err[21];
  double tsu_lhc_raa[21];
  double tsu_lhc_raa_err[21];
  double tsu_lhc_raa_low[21];
  double tsu_lhc_raa_high[21];
  for(int i=0; i<11; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<21; i++)
    {
      ifin >> tsu_lhc_pt[i] >> dtmp >> dtmp >> dtmp >> tsu_lhc_raa_low[i] >> tsu_lhc_raa_high[i];
      tsu_lhc_pt_err[i] = 0;
      tsu_lhc_raa[i] = (tsu_lhc_raa_low[i] + tsu_lhc_raa_high[i])/2;
      tsu_lhc_raa_err[i] = (tsu_lhc_raa_high[i] - tsu_lhc_raa_low[i])/2;
    }
  ifin.close();
  TGraphErrors *tsu_JpsiRaa_LHC_0040 = new TGraphErrors(21,tsu_lhc_pt,tsu_lhc_raa,tsu_lhc_pt_err,tsu_lhc_raa_err);
  int tmp_color = kViolet-9;
  tsu_JpsiRaa_LHC_0040->SetFillStyle(1001);
  tsu_JpsiRaa_LHC_0040->SetFillColor(tmp_color);
  tsu_JpsiRaa_LHC_0040->SetLineColor(tmp_color);
  tsu_JpsiRaa_LHC_0040->SetLineWidth(1);
  tsu_JpsiRaa_LHC_0040->SetLineStyle(1);
  TGraph *tsu_JpsiRaa_LHC_0040_up = new TGraph(21,tsu_lhc_pt,tsu_lhc_raa_high);
  tsu_JpsiRaa_LHC_0040_up->SetLineColor(tmp_color);
  tsu_JpsiRaa_LHC_0040_up->SetLineWidth(1);
  tsu_JpsiRaa_LHC_0040_up->SetLineStyle(1);
  TGraph *tsu_JpsiRaa_LHC_0040_low = new TGraph(21,tsu_lhc_pt,tsu_lhc_raa_low);
  tsu_JpsiRaa_LHC_0040_low->SetLineColor(tmp_color);
  tsu_JpsiRaa_LHC_0040_low->SetLineWidth(1);
  tsu_JpsiRaa_LHC_0040_low->SetLineStyle(1);
  double x1, y1, x2, y2;
  TGraph *tsu_JpsiRaa_RHIC_0040 = new TGraph(grTHU[0]->GetN());
  for(int i=0; i<grTHU[0]->GetN(); i++)
    {
      grTHU[1]->GetPoint(i,x1,y1);
      grTHU[2]->GetPoint(i,x2,y2);
      tsu_JpsiRaa_RHIC_0040->SetPoint(i,(x1+x2)/2,(y1*ncoll[1]+y2*ncoll[2])/(ncoll[1]+ncoll[2]));
    }
  tsu_JpsiRaa_RHIC_0040->SetLineWidth(2);
  tsu_JpsiRaa_RHIC_0040->SetLineStyle(1);
  tsu_JpsiRaa_RHIC_0040->SetLineColor(2);

  // tamu model
  double tamu_lhc_pt[76];
  double tamu_lhc_pt_err[76];
  double tamu_lhc_raa[76];
  double tamu_lhc_raa_err[76];
  double tamu_lhc_raa_low[76];
  double tamu_lhc_raa_high[76];
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/pt_lhc_u_b5.7_noshdw_ftBels_toJulianBook_140415.dat");
  for(int i=0; i<3; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<76; i++)
    {
      ifin >> tamu_lhc_pt[i] >> dtmp >> dtmp >> tamu_lhc_raa_low[i];
    }
  ifin.close();
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/pt_lhc_u_b5.7_shdw_ftBels_toJulianBook_140415.dat");
  for(int i=0; i<3; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<76; i++)
    {
      ifin >> dtmp >> dtmp >> dtmp >> tamu_lhc_raa_high[i];
      tamu_lhc_pt_err[i] = 0;
      tamu_lhc_raa[i] = (tamu_lhc_raa_low[i] + tamu_lhc_raa_high[i])/2;
      tamu_lhc_raa_err[i] = (tamu_lhc_raa_high[i] - tamu_lhc_raa_low[i])/2;
    }
  ifin.close();
  TGraphErrors *tamu_JpsiRaa_LHC_0040 = new TGraphErrors(67,tamu_lhc_pt,tamu_lhc_raa,tamu_lhc_pt_err,tamu_lhc_raa_err);
  int tmp_color = kCyan-9;
  tamu_JpsiRaa_LHC_0040->SetFillStyle(1001);
  tamu_JpsiRaa_LHC_0040->SetFillColor(tmp_color);
  tamu_JpsiRaa_LHC_0040->SetLineColor(tmp_color);
  tamu_JpsiRaa_LHC_0040->SetLineWidth(1);
  tamu_JpsiRaa_LHC_0040->SetLineStyle(1);
  TGraph *tamu_JpsiRaa_LHC_0040_up = new TGraph(67,tamu_lhc_pt,tamu_lhc_raa_high);
  tamu_JpsiRaa_LHC_0040_up->SetLineColor(tmp_color);
  tamu_JpsiRaa_LHC_0040_up->SetLineWidth(1);
  tamu_JpsiRaa_LHC_0040_up->SetLineStyle(1);
  TGraph *tamu_JpsiRaa_LHC_0040_low = new TGraph(67,tamu_lhc_pt,tamu_lhc_raa_low);
  tamu_JpsiRaa_LHC_0040_low->SetLineColor(tmp_color);
  tamu_JpsiRaa_LHC_0040_low->SetLineWidth(1);
  tamu_JpsiRaa_LHC_0040_low->SetLineStyle(1);
  TGraph *tamu_JpsiRaa_RHIC_0040 = new TGraph(grTAMU[0]->GetN());
  for(int i=0; i<grTAMU[0]->GetN(); i++)
    {
      grTAMU[1]->GetPoint(i,x1,y1);
      grTAMU[2]->GetPoint(i,x2,y2);
      tamu_JpsiRaa_RHIC_0040->SetPoint(i,(x1+x2)/2,(y1*ncoll[1]+y2*ncoll[2])/(ncoll[1]+ncoll[2]));
    }
  tamu_JpsiRaa_RHIC_0040->SetLineWidth(2);
  tamu_JpsiRaa_RHIC_0040->SetLineStyle(7);
  tamu_JpsiRaa_RHIC_0040->SetLineColor(2);


  TCanvas *c8 = new TCanvas("JpsiRaaVsPt_STAR_LHC_model","JpsiRaaVsPt_STAR_LHC_model",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  ScaleHistoTitle(hRaa2,28,1,24,28,1,24,63);
  hRaa2->SetMaximum(1.8);
  hRaa2->DrawCopy();
  TLine *line = GetLine(hRaa2->GetXaxis()->GetXmin(),1,10,1,1);
  line->Draw();
  grAlice0040Sys->Draw("sameE5");
  grAlice0040->Draw("samesPEZ");
  grCms0040Sys->Draw("sameE5");
  grCms0040->Draw("samesPEZ");
  hStarRaa0040Sys->Draw("sameE5");
  hStarRaa0040->Draw("samesPEZ");

  star2->Draw();
  box_star->Draw();
  box_alice->Draw();
  box_cms->Draw();

  TLegend *leg8 = new TLegend(0.18,0.75,0.41,0.95);
  leg8->SetBorderSize(0);
  leg8->SetFillColor(0);
  leg8->SetTextFont(63);
  leg8->SetTextSize(22);
  leg8->SetHeader("0-40% centrality");
  leg8->AddEntry(hStarRaa0040,"STAR:  Au+Au, #sqrt{s_{NN}} = 0.2 TeV,  |y| < 0.5","P");
  leg8->AddEntry(grAlice0040, "ALICE: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 0.8","P");
  leg8->AddEntry(grCms0040,   "CMS:   Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 2.4","P");
  leg8->Draw();

  if(savePlot)
    {
      c8->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040_2.pdf",run_type,run_cfg_name.Data()));
      c8->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040_2.png",run_type,run_cfg_name.Data()));
      c8->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040_2.jpg",run_type,run_cfg_name.Data()));
      c8->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_0040_2.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c81 = new TCanvas("JpsiRaaVsPt_STAR_LHC_model2","JpsiRaaVsPt_STAR_LHC_model2",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  //gPad->SetGrid(1,1);
  ScaleHistoTitle(hRaa2,28,1,24,28,1,24,63);
  hRaa2->DrawCopy();
  line->Draw();

  tamu_JpsiRaa_LHC_0040->Draw("samesE3");

  tsu_JpsiRaa_LHC_0040->Draw("samesE3");
  //tsu_JpsiRaa_LHC_0040_up->Draw("samesl");
  //tsu_JpsiRaa_LHC_0040_low->Draw("samesl");


  //tamu_JpsiRaa_LHC_0040_up->Draw("samesl");
  //tamu_JpsiRaa_LHC_0040_low->Draw("samesl");

  tsu_JpsiRaa_RHIC_0040->Draw("samesl");
  tamu_JpsiRaa_RHIC_0040->Draw("samesl");

  grAlice0040Sys->Draw("sameE5");
  grAlice0040->Draw("samesPEZ");
  grCms0040Sys->Draw("sameE5");
  grCms0040->Draw("samesPEZ");
  hStarRaa0040Sys->Draw("sameE5");
  hStarRaa0040->Draw("samesPEZ");

  star2->Draw();
  box_star->Draw();
  box_alice->Draw();
  box_cms->Draw();

  leg8->Draw();

  TPaveText *tmodel = GetPaveText(0.21,0.38,0.64,0.74);
  tmodel->SetTextFont(62);
  tmodel->SetTextAlign(11);
  tmodel->SetTextSize(0.035);
  tmodel->AddText("Tsinghua Model");
  tmodel->AddText("TAMU Model");
  tmodel->Draw();

  TLegend *leg21 = new TLegend(0.42,0.65,0.6,0.75);
  leg21->SetBorderSize(0);
  leg21->SetFillColor(0);
  leg21->SetTextFont(62);
  leg21->SetTextSize(0.035);
  leg21->AddEntry(tsu_JpsiRaa_RHIC_0040,"RHIC","L");
  leg21->AddEntry(tamu_JpsiRaa_RHIC_0040,"RHIC","L");
  leg21->Draw();
  TLegend *leg22 = new TLegend(0.55,0.65,0.7,0.75);
  leg22->SetBorderSize(0);
  leg22->SetFillColor(0);
  leg22->SetTextFont(62);
  leg22->SetTextSize(0.035);
  leg22->AddEntry(tsu_JpsiRaa_LHC_0040,"LHC","f");
  leg22->AddEntry(tamu_JpsiRaa_LHC_0040,"LHC","f");
  leg22->Draw();
  
  if(savePlot)
    {
      c81->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_vs_model_0040.pdf",run_type,run_cfg_name.Data()));
      c81->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_vs_model_0040.png",run_type,run_cfg_name.Data()));
      c81->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_vs_model_0040.jpg",run_type,run_cfg_name.Data()));
      c81->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_LHC_vs_model_0040.eps",run_type,run_cfg_name.Data()));
    }



 if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.root",run_cfg_name.Data()),"recreate");
      hJpsipp->Write();
      for(int k=0; k<nCentBins; k++)
	{
	  for(int i=0; i<hJpsiXsec[k]->GetN(); i++)
	    {
	      hJpsiXsec[k]->GetPoint(i,x,y);
	      hJpsiXsec[k]->SetPoint(i,x,y/scale_factor[k]);
	      hJpsiXsec[k]->SetPointError(i,0,hJpsiXsec[k]->GetErrorY(i)/scale_factor[k]);
	      
	      hJpsiXsecSys[k]->SetPoint(i,x,y/scale_factor[k]);
	      hJpsiXsecSys[k]->SetPointError(i,0.2,hSys->GetBinContent(i+1)*y/scale_factor[k]);
	    }
	  hJpsiXsec[k]->Write();
	  hJpsiXsecSys[k]->Write();
	  hJpsiRaa[k]->Write();
	  hJpsiRaaSys[k]->Write();
	}
      grAlice0040->Write("ALICE_JpsiRaaVsPt_0040");
      grAlice0040Sys->Write("ALICE_JpsiRaaVsPt_0040_sys");

      ofstream out;
      for(int k=0; k<nCentBins; k++)
	{
	  out.open(Form("Rootfiles/2016sQM/%s_JpsiXsec_%s.dat",run_type,cent_Title[k]),std::ofstream::out);
	  for(int i=0; i<hJpsiXsec[k]->GetN(); i++)
	    {
	      hJpsiXsec[k]->GetPoint(i,x,y);
	      out << setw(5) << x << setw(5)<< hJpsiInvYield[k]->GetBinWidth(i+1) << setw(15) << y << setw(15) << hJpsiXsec[k]->GetErrorY(i) << endl;
	    }
	  out.close();
	}
    }
}

//================================================
void rawSignal(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(kFALSE);
  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut));
  TH1F *hJpsiSignal[nPtBins];
  TH1F *hSeUL[nPtBins];
  TH1F *hSeLS[nPtBins];
  TH1F *hMixBkg[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiSignal[i] = (TH1F*)fin->Get(Form("Jpsi_Signal_pt%s_cent0060",pt_Name[i]));
      hSeUL[i] = (TH1F*)fin->Get(Form("DataUL_pt%s_cent0060",pt_Name[i]));
      hSeLS[i] = (TH1F*)fin->Get(Form("DataLS_pt%s_cent0060",pt_Name[i]));
      hMixBkg[i] = (TH1F*)fin->Get(Form("MEbkg_pt%s_cent0060",pt_Name[i]));
    }
  TH1F *hJpsiCount = (TH1F*)fin->Get("Jpsi_BinCountYield_cent0060");
  TH1F *hSigToBkg = (TH1F*)fin->Get("Jpsi_SigToBkg_cent0060");

  // raw signal
  TGaxis::SetMaxDigits(3); 
  int ptCuts[2] = {0,5};
  TH1F *hhSeUL, *hhSeLS, *hhMixBkg, *hhSignal;
  for(int i=0; i<2; i++)
    {
      TCanvas *c = new TCanvas(Form("Jpsi_All_pt%d",ptCuts[i]),Form("Jpsi_All_pt%d",ptCuts[i]),800,650);
      SetPadMargin(gPad,0.12,0.12,0.05,0.05);
      hhSeUL = (TH1F*)hSeUL[0]->Clone(Form("hhSeUL_pt%d",ptCuts[i]));
      hhSeLS = (TH1F*)hSeLS[0]->Clone(Form("hhSeLS_pt%d",ptCuts[i]));
      hhMixBkg = (TH1F*)hMixBkg[0]->Clone(Form("hhMixBkg_pt%d",ptCuts[i]));
      hhSignal = (TH1F*)hJpsiSignal[0]->Clone(Form("hhSignal_pt%d",ptCuts[i]));
      
      if(i==1)
	{
	  hhSeUL->Reset();
	  hhSeLS->Reset();
	  hhMixBkg->Reset();
	  hhSignal->Reset();
	  for(int j=6; j<nPtBins; j++)
	    {
	      hhSeUL->Add(hSeUL[j]);
	      hhSeLS->Add(hSeLS[j]);
	      hhMixBkg->Add(hMixBkg[j]);
	      hhSignal->Add(hJpsiSignal[j]);
	    }
	}
    
      ScaleHistoTitle(hhSeUL,0.05,1,0.04,0.05,1.1,0.04,62);
      hhSeUL->GetXaxis()->SetRangeUser(2.5,3.8);
      TH1F *hplot = (TH1F*)hhSeUL->Clone(Form("hplot_%d",i));
      if(i==0) 
	{
	  hplot->Scale(0.15);
	  hhSeLS->Scale(0.15);
	  hhMixBkg->Scale(0.15);
	  hplot->SetMaximum(4*hhSignal->GetMaximum());
	  hplot->SetMinimum(-5e3);
	}
      if(i==1) 
	{
	  hplot->SetMaximum(1.45*hplot->GetMaximum());
	  hplot->SetMinimum(-1e2);
	}
      hplot->GetYaxis()->SetNdivisions(505);
      hplot->SetYTitle("Counts");
      hplot->SetXTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
      hplot->Draw();
      hplot->SetMarkerStyle(25);
      hhSeLS->Draw("samesP");
      hhMixBkg->Draw("samesHIST");
      hhSignal->SetMarkerColor(2);
      hhSignal->SetLineColor(2);
      hhSignal->Draw("sames");

      TF1 *funcSignal = new TF1(Form("Fit_JpsiSignal_pt%d",ptCuts[i]),Form("gausn(0)+pol%d(3)",3),2.45,3.8);
      funcSignal->SetParameter(1,3.09);
      funcSignal->SetParameter(2,0.06);
      TFitResultPtr ptr = hhSignal->Fit(funcSignal,"IR0");
      funcSignal->SetLineColor(1);
      funcSignal->SetLineStyle(1);
      funcSignal->Draw("same");



      TPaveText *star = GetPaveText(0.18,0.4,0.85,0.92,0.045);
      star->AddText("STAR preliminary");
      star->SetTextFont(32);
      star->SetTextColor(2);
      star->Draw();

      TLegend *leg = new TLegend(0.15,0.3,0.3,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(63);
      leg->SetTextSize(18);
      if(i==0)
	{
	  leg->AddEntry(hplot,"Unlike-sign pairs (#times 0.15)","P");
	  leg->AddEntry(hhSeLS,"Like-sign pairs (#times 0.15)","PL");
	  leg->AddEntry(hhMixBkg,"Mixed-event (#times 0.15)","L");
	}
      else
	{
	  leg->AddEntry(hplot,"Unlike-sign pairs","P");
	  leg->AddEntry(hhSeLS,"Like-sign pairs","PL");
	  leg->AddEntry(hhMixBkg,"Mixed-event","L");
	}
      leg->Draw();

      double nAll = hhSeUL->Integral(hhSeUL->FindFixBin(3.025),hhSeUL->FindFixBin(3.175));
      double nJpsi = hJpsiCount->GetBinContent(0);
      double error = hJpsiCount->GetBinError(0);
      if(i==1) nJpsi = hJpsiCount->IntegralAndError(6,9,error);
      double nJpsiFit = funcSignal->GetParameter(0)/hhSignal->GetBinWidth(1);
      TPaveText *t1 = GetPaveText(0.6,0.75,0.62,0.92,0.038,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("Au+Au @ 200 GeV"));
      t1->AddText(Form("#it{L} ~ %1.1f nb^{-1}",luminosity));
      t1->AddText(Form("|y_{J/#psi}|<0.5, p_{T,J/#psi} > %d GeV/c",ptCuts[i]));
      if(i==0) t1->AddText(Form("N_{J/#psi} = %2.0f, S/B = 1:%2.1f",nJpsi, (nAll-nJpsi)/nJpsi));
      if(i==1) t1->AddText(Form("N_{J/#psi} = %2.0f, S/B = 1:%2.1f",nJpsiFit, (nAll-nJpsiFit)/nJpsiFit));
      t1->AddText(Form("Significance = %2.1f",nJpsi/error));
      t1->Draw();
      cout << funcSignal->GetParameter(0)/funcSignal->GetParError(0) << endl;

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt%d_cent0060.pdf",run_type,run_cfg_name.Data(),ptCuts[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt%d_cent0060.png",run_type,run_cfg_name.Data(),ptCuts[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt%d_cent0060.jpg",run_type,run_cfg_name.Data(),ptCuts[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt%d_cent0060.eps",run_type,run_cfg_name.Data(),ptCuts[i]));
	}
    }
  return;

  // all signal
  TCanvas *c1 = new TCanvas("Fit_Jpsi_All","Fit_Jpsi_All",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hJpsiSignal[0],0.05,1,0.04,0.05,1.1,0.04,62);
  hJpsiSignal[0]->SetMaximum(1.5*hJpsiSignal[0]->GetMaximum());
  hJpsiSignal[0]->SetYTitle("Counts");
  hJpsiSignal[0]->SetMarkerSize(1.5);
  hJpsiSignal[0]->Draw();
  TF1 *funcSignal = new TF1("Fit_JpsiSignal_AllPt","gausn(0)+pol3(3)",2.3,4);
  funcSignal->SetParameter(1,3.09);
  funcSignal->SetParameter(2,0.06);
  TFitResultPtr ptr = hJpsiSignal[0]->Fit(funcSignal,"IR0");
  funcSignal->SetLineColor(2);
  funcSignal->Draw("same");

  TPaveText *t1 = GetPaveText(0.55,0.75,0.5,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText("Au+Au @ 200 GeV 0-60%");
  t1->AddText(Form("#it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}");
  t1->AddText("|y|<0.5, p_{T} > 1 GeV/c");
  t1->AddText(Form("S/B = 1:%2.1f",1./hSigToBkg[0]->GetBinContent(0)));
  t1->AddText(Form("N_{J/#psi} = %2.0f",hJpsiCount[0]->GetBinContent(0)));
  t1->AddText(Form("Significance = %2.1f#sigma",hJpsiCount[0]->GetBinContent(0)/hJpsiCount[0]->GetBinError(0)));
  t1->Draw();

  TPaveText *star = GetPaveText(0.22,0.32,0.85,0.9,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.eps",run_type,run_cfg_name.Data()));
    }

  return;
  TCanvas *c2 = new TCanvas("Fit_jpsi_ptbins","Fit_jpsi_ptbins",1100,650);
  c2->Divide(3,2);
  for(int i=1; i<nPtBins; i++)
    {
      c2->cd(i);
      SetPadMargin(gPad,0.13,0.13,0.05,0.05);
      ScaleHistoTitle(hJpsiSignal[i],0.06,1,0.05,0.06,1.1,0.05,62);
      hJpsiSignal[i]->SetMaximum(1.8*hJpsiSignal[i]->GetMaximum());
      hJpsiSignal[i]->SetYTitle("Counts");
      hJpsiSignal[i]->SetMarkerSize(1.5);
      hJpsiSignal[i]->Draw();
      TF1 *funcSignal = new TF1(Form("Fit_JpsiSignal_pt%s",pt_Name[i]),"gausn(0)+pol1(3)",2.3,4);
      funcSignal->SetParameter(1,3.09);
      funcSignal->SetParameter(2,0.06);
      hJpsiSignal[i]->Fit(funcSignal,"IR0");
      funcSignal->SetLineColor(2);
      funcSignal->Draw("same");
      if(i==1)
	{
	  TPaveText *t1 = GetPaveText(0.17,0.37,0.65,0.92,0.06,62);
	  t1->SetTextAlign(11);
	  t1->AddText("Au+Au @ 200 GeV");
	  t1->AddText("0-20%");
	  t1->AddText(Form("#it{L} ~ %1.1f nb^{-1}",luminosity));
	  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}");
	  t1->Draw();
	}

      t1 = GetPaveText(0.53,0.75,0.5,0.83,0.06,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]));
      t1->AddText(Form("S/B = 1:%2.1f",1./hSigToBkg->GetBinContent(i)));
      t1->AddText(Form("N_{J/#psi} = %2.0f",hJpsiCount->GetBinContent(i)));
      t1->AddText(Form("Significance: %2.1f",hJpsiCount->GetBinContent(i)/hJpsiCount->GetBinError(i)));
      t1->Draw();
    }
  /*
  c2->cd(2);
  TPaveText *star = GetPaveText(0.25,0.5,0.85,0.93,0.06);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  */
  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0020.pdf",run_type,run_cfg_name.Data()));
      /*
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0060.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0060.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0060.eps",run_type,run_cfg_name.Data()));
      */
    }

}


//================================================
void compareQM2015(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *fout = TFile::Open(Form("Rootfiles/%s/Publication.Jpsi.200GeV.root",run_cfg_name.Data()),"read");

  // ===== pp reference

  // sQM2016
  TH1F *hpp_2016 = (TH1F*)fout->Get("hPPJpsiFinal");
  hpp_2016->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hpp_2016,"J/#Psi cross section in p+p collisions at #sqrt{s} = 200 GeV;p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",kTRUE);

  // QM2015
  TH1F *hpp_2015 = (TH1F*)fout->Get("Jpsi_xsec_pp200_combined");
  hpp_2015->Draw("sames");

  TLegend *leg = new TLegend(0.5,0.6,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p reference");
  leg->AddEntry(hpp_2016,"2016 sQM","P");
  leg->AddEntry(hpp_2015,"2015  QM","P");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Comp_pp_2015vs2016.pdf",run_type,run_cfg_name.Data()));
  
  TH1F *hppRatio = (TH1F*)hpp_2015->Clone("hppRatio");
  hppRatio->Reset();
  for(int bin=1; bin<=6; bin++)
    {
      double y15 = hpp_2015->GetBinContent(bin);
      double e15 = hpp_2015->GetBinError(bin);
      double y16 = 0;
      if(bin<=3) y16 = hpp_2016->GetBinContent(bin+1);
      else if(bin==4)
	{
	  y16 = hpp_2016->GetBinContent(5) * hpp_2016->GetBinCenter(5) + hpp_2016->GetBinContent(6) * hpp_2016->GetBinCenter(6);
	  y16 = y16/(hpp_2015->GetBinCenter(bin)*hpp_2015->GetBinWidth(bin));
	}
      else y16 = hpp_2016->GetBinContent(bin+2);
      hppRatio->SetBinContent(bin, y15/y16);
      hppRatio->SetBinError(bin, e15/y16);
    }
  hppRatio->GetYaxis()->SetRangeUser(0,2);
  TCanvas *c = draw1D(hppRatio,"Ratio of p+p reference at #sqrt{s} = 200 GeV;p_{T} (GeV/c);QM2015/sQM2016",kFALSE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Comp_pp_2015to2016.pdf",run_type,run_cfg_name.Data()));


  TFile *fin[2];
  fin[0] = TFile::Open("Rootfiles/2016sQM.root","read");
  fin[1] = TFile::Open("Rootfiles/2015QM.root","read");

  // ===== AuAu yield
  TGraphAsymmErrors *hJpsiAuAuYield[4][nCentBins];
  TGraphAsymmErrors *hJpsiAuAuYieldSys[4][nCentBins];
  TH1F *hTBW[4];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiAuAuYield[0][k] = (TGraphAsymmErrors*)fout->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[k]));
      hJpsiAuAuYieldSys[0][k] = (TGraphAsymmErrors*)fout->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[k]));

      hJpsiAuAuYield[1][k] = (TGraphAsymmErrors*)fout->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[k]));
      hJpsiAuAuYieldSys[1][k] = (TGraphAsymmErrors*)fout->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[k]));

      for(int i=0; i<2; i++)
	{
	  hJpsiAuAuYield[2+i][k] = (TGraphAsymmErrors*)fin[i]->Get(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
	  hJpsiAuAuYield[2+i][k]->SetName(Form("%s_%d",hJpsiAuAuYield[2+i][k]->GetName(),i));
	  hJpsiAuAuYieldSys[2+i][k] = (TGraphAsymmErrors*)fin[i]->Get(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));
	  hJpsiAuAuYieldSys[2+i][k]->SetName(Form("%s_%d",hJpsiAuAuYieldSys[2+i][k]->GetName(),i));
	}
      hTBW[k] = (TH1F*)fout->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[k]));
    }

  TGraphErrors *hJpsiAuAuYieldRatio[4][nCentBins];
  TGraphErrors *hJpsiAuAuYieldRatioSys[4][nCentBins];
  double x,y;
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<4; i++)
	{
	  hJpsiAuAuYieldRatio[i][k] = new TGraphErrors(hJpsiAuAuYield[i][k]->GetN());
	  hJpsiAuAuYieldRatio[i][k]->SetName(Form("Ratio_%d_%d",i,k));

	  hJpsiAuAuYieldRatioSys[i][k] = new TGraphErrors(hJpsiAuAuYieldSys[i][k]->GetN());
	  hJpsiAuAuYieldRatioSys[i][k]->SetName(Form("RatioSys_%d_%d",i,k));

	  for(int ipoint=0; ipoint<hJpsiAuAuYield[i][k]->GetN(); ipoint++)
	    {
	      hJpsiAuAuYield[i][k]->GetPoint(ipoint, x, y);
	      double dpT = 0;
	      if(i==0) dpT = 1;
	      if(i==1)
		{
		  if(ipoint<5) dpT = 1;
		  else         dpT = 2;
		}
	      if(i==2)
		{
		  if(ipoint<3) dpT = 1;
		  else         dpT = 2;
		}
	      if(i==3)
		{
		  if(ipoint<6) dpT = 1;
		  else if(ipoint<8) dpT = 2;
		  else dpT = 5;
		}
	      double pT = x;
	      if(i==0) pT = 0.5 + ipoint;
	      int start_bin = hTBW[k]->FindFixBin(pT-dpT/2);
	      int end_bin   = hTBW[k]->FindFixBin(pT+dpT/2);
	      double scale = 0;
	      for(int ibin=start_bin; ibin<=end_bin; ibin++)
		{
		  scale += hTBW[k]->GetBinContent(ibin) * hTBW[k]->GetBinCenter(ibin) * hTBW[k]->GetBinWidth(ibin);
		}
	      scale = scale / pT / dpT;
	      hJpsiAuAuYieldRatio[i][k]->SetPoint(ipoint,x,y/scale);
	      hJpsiAuAuYieldRatio[i][k]->SetPointError(ipoint, 0, hJpsiAuAuYield[i][k]->GetErrorYhigh(ipoint)/scale);
	      hJpsiAuAuYieldRatio[i][k]->SetMarkerStyle(24 - i/2*4);
	      hJpsiAuAuYieldRatio[i][k]->SetMarkerSize(1.5);
	      hJpsiAuAuYieldRatio[i][k]->SetMarkerColor(pow(2,i));
	      hJpsiAuAuYieldRatio[i][k]->SetLineColor(pow(2,i));
	      hJpsiAuAuYieldRatioSys[i][k]->SetPoint(ipoint,x,y/scale);
	      hJpsiAuAuYieldRatioSys[i][k]->SetPointError(ipoint, 0.2, hJpsiAuAuYieldSys[i][k]->GetErrorYhigh(ipoint)/scale);
	      hJpsiAuAuYieldRatioSys[i][k]->SetMarkerSize(0);
	      hJpsiAuAuYieldRatioSys[i][k]->SetLineColor(pow(2,i));
	      hJpsiAuAuYieldRatioSys[i][k]->SetFillStyle(0);
	    }
	}
    }

  TCanvas *c = new TCanvas("Yield_Jpsi_vs_pub","Yield_Jpsi_vs_pub",1100,700);
  TPad *pad = GetSinglePad("Jpsi",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();

  TH1F *hJpsi = new TH1F("Yield_Jpsi",";p_{T} (GeV/c);Data/TBW",10,-0.1,10);
  hJpsi->GetYaxis()->SetRangeUser(0,2);
  hJpsi->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hJpsi,24,1.9,20,24,1.9,20,63);
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hJpsi->DrawCopy();
      TLine *line = GetLine(hJpsi->GetXaxis()->GetXmin(),1,hJpsi->GetXaxis()->GetXmax(),1,1);
      line->Draw();
      for(int i=0; i<4; i++)
	{
	  hJpsiAuAuYieldRatioSys[i][k]->Draw("sameE5");
	  hJpsiAuAuYieldRatio[i][k]->Draw("samesPEZ");
	}
      if(k%2==0) t1 = GetPaveText(0.22,0.25,0.9,0.95);
      else       t1 = GetPaveText(0.1,0.12,0.9,0.95);
      t1->SetTextFont(63);
      t1->SetTextSize(24);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  pad->cd(1);
  TLegend *leg = new TLegend(0.65,0.05,0.85,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05);
  leg->AddEntry(hJpsiAuAuYieldRatio[0][0],"Pub: low p_{T}","P");
  leg->AddEntry(hJpsiAuAuYieldRatio[1][0],"Pub: high p_{T}","P");
  leg->AddEntry(hJpsiAuAuYieldRatio[2][0],"2016 sQM","P");
  leg->AddEntry(hJpsiAuAuYieldRatio[3][0],"2015  QM","P");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Comp_Yield_2015vs2016.pdf",run_type,run_cfg_name.Data()));

  // ===== RAA 
  TGraphErrors *hJpsiRaa[2][nCentBins];
  TGraphErrors *hJpsiRaaSys[2][nCentBins];
  for(int i=0; i<2; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiRaa[i][k] = (TGraphErrors*)fin[i]->Get(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
	  hJpsiRaa[i][k]->SetName(Form("%s_%d",hJpsiRaa[i][k]->GetName(),i));
	  hJpsiRaa[i][k]->SetMarkerStyle(20+i);
	  hJpsiRaa[i][k]->SetMarkerColor(i+1);
	  hJpsiRaa[i][k]->SetLineColor(i+1);
	  hJpsiRaa[i][k]->SetMarkerSize(2);

	  hJpsiRaaSys[i][k] = (TGraphErrors*)fin[i]->Get(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
	  hJpsiRaaSys[i][k]->SetName(Form("%s_%d",hJpsiRaaSys[i][k]->GetName(),i));
	  hJpsiRaaSys[i][k]->SetMarkerColor(i+1);
	  hJpsiRaaSys[i][k]->SetLineColor(i+1);
	}
    }

  TCanvas *c = new TCanvas("Raa_Jpsi_vs_pub","Raa_Jpsi_vs_pub",1100,700);
  TPad *pad = GetSinglePad("jet_trigger",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();

  TH1F *hRaa = new TH1F("Raa_Jpsi",";p_{T} (GeV/c);R_{AA}",10,-0.1,15);
  hRaa->GetYaxis()->SetRangeUser(0.05,1.5);
  hRaa->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hRaa,24,1.9,20,24,1.9,20,63);
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();
      for(int i=0; i<2; i++)
	{
	  hJpsiRaaSys[i][k]->Draw("sameE5");
	  hJpsiRaa[i][k]->Draw("samesPEZ");
	}
      if(k%2==0) t1 = GetPaveText(0.22,0.25,0.9,0.95);
      else       t1 = GetPaveText(0.1,0.12,0.9,0.95);
      t1->SetTextFont(63);
      t1->SetTextSize(24);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  pad->cd(1);
  TLegend *leg = new TLegend(0.75,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.06);
  leg->AddEntry(hJpsiRaa[0][0],"2016 sQM","P");
  leg->AddEntry(hJpsiRaa[1][0],"2015  QM","P");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Comp_Raa_2015vs2016.pdf",run_type,run_cfg_name.Data()));
}

//================================================
void upsilon(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("Rootfiles/UpsFinal.root","read");
  TH1F *hUL = (TH1F*)f->Get("SE_UL");
  TH1F *hLS = (TH1F*)f->Get("SE_LS");
  TH1F *hSub = (TH1F*)f->Get("SE_Ups");
  TF1 *func[4];
  func[0] = (TF1*)f->Get("fun");
  func[1] = (TF1*)f->Get("r1s");
  func[2] = (TF1*)f->Get("r2s");
  func[3] = (TF1*)f->Get("r3s");

  TCanvas *c1 = new TCanvas("Fit_Jpsi_All","Fit_Jpsi_All",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hUL,0.05,1,0.04,0.05,1.1,0.04,62);
  hUL->SetMaximum(120);
  hUL->SetTitle("");
  hUL->SetLineColor(2);
  hUL->Draw();
  hLS->Draw("sames HIST");

  int lowbin = hUL->FindBin(9);
  int highbin = hUL->FindBin(12);
  printf("[i] %d counts in unlike-sign\n",hUL->Integral(lowbin,highbin));
  printf("[i] %d counts in like-sign\n",hLS->Integral(lowbin,highbin));
  TPaveText *t1 = GetPaveText(0.15,0.25,0.85,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->Draw();
  TLegend *leg = new TLegend(0.5,0.55,0.75,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("#Upsilon#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hUL,"Unlike-sign pairs","P");
  leg->AddEntry(hLS,"Like-sign pairs","L");
  leg->Draw();
  TPaveText *star = GetPaveText(0.15,0.4,0.15,0.2,0.04,20);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c2 = new TCanvas("Fit_Upsilon","Fit_Upsilon",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  //gPad->SetGrid(1,1);
  ScaleHistoTitle(hSub,0.05,1,0.04,0.05,1.1,0.04,62);
  hSub->SetTitle("");
  hSub->SetMarkerColor(2);
  hSub->SetLineColor(1);
  hSub->GetYaxis()->SetRangeUser(-20,40);
  hSub->Draw();
  for(int i=0; i<4; i++)
    {
      func[i]->SetLineStyle(i+1);
      func[i]->Draw("sames");
    }
  TF1 *funcBkg = new TF1("funcBkg","[0]*expo(1)",8,12);
  for(int i=0; i<3; i++)
    {
      funcBkg->SetParameter(i,func[0]->GetParameter(9+i));
    }
  funcBkg->SetLineColor(kGreen+1);
  funcBkg->SetLineStyle(5);
  funcBkg->Draw("sames");
  printf("Fit: chi2/ndf = %1.1f/%1.0f\n",func[0]->GetChisquare(), func[0]->GetNDF());

  TPaveText *t1 = GetPaveText(0.15,0.25,0.85,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->Draw();
  TLegend *leg = new TLegend(0.65,0.63,0.85,0.92);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("#Upsilon#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hSub,"UL-LS","P");
  leg->AddEntry(func[0],"Combined fit","L");
  leg->AddEntry(func[1],"#Upsilon(1S)","L");
  leg->AddEntry(func[2],"#Upsilon(2S)","L");
  leg->AddEntry(func[3],"#Upsilon(3S)","L");
  leg->AddEntry(funcBkg,"Drell-Yan+b#bar{b}","L");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.65,0.85,0.55,0.60,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("#chi^{2}/NDF = %1.1f/%d\n",func[0]->GetChisquare(), func[0]->GetNDF()));
  t1->Draw();
  TPaveText *star = GetPaveText(0.7,0.85,0.15,0.2,0.04,20);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.eps",run_type,run_cfg_name.Data()));
    }
 
}

//================================================
void ppRef(const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TFile *fout = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","update");

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TFile *f = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.pp200GeV.root","read");

  // pp data
  TH1F *hpp = new TH1F("pp200_Jpsi",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hpp->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hpp,"",kTRUE);
  const int npp = 11;
  double xpp[npp] = {2.25, 2.75, 3.25, 3.75, 4.5, 5.5, 6.5, 7.5, 9, 11, 13};
  double exlpp[npp];
  double sxlpp[npp];
  for(int i=0; i<npp; i++) { exlpp[i] = 0; sxlpp[i] = 0.2; }
  double ypp[npp] = {0.68, 0.318, 0.187, 0.1032, 0.0334, 0.00905, 0.00154, 0.00084, 2.005e-4, 4.55e-5, 9.71e-6};
  double eylpp[npp] = {0.14, 0.050, 0.023, 0.0118, 0.0034, 0.00117, 0.00042, 0.00017, 2.42e-5, 7.2e-6, 2.41e-6};
  double sylpp[npp] = {0.06, 0.028, 0.018, 0.0089, 0.0028, 0.00077, 0.00013, 0.00025, 5.21e-5, 1.18e-5, 2.62e-6};
  double syhpp[npp] = {0.03, 0.007, 0.009, 0.0018, 0.0004, 0.00011, 0.00003, 0.00013, 1.6e-6, 4e-7, 1.6e-7};

  double xsec_5_pub = 0;
  for(int i=0; i<6; i++)
    {
      double pt = xpp[i+5];
      double dpt = 1;
      if(i>2) dpt = 2;
      xsec_5_pub += ypp[i+5] * 2 * pi * pt * dpt;
    }
  printf("[i] published pT > 5, xsec = %4.2e \n",xsec_5_pub);
  

  TGraphAsymmErrors *gHighPtPP = new TGraphAsymmErrors(npp, xpp, ypp, exlpp, exlpp, eylpp, eylpp);
  gHighPtPP->SetName("Jpsi_xsec_pp200_highPt");
  gHighPtPP->SetMarkerStyle(20);
  gHighPtPP->SetMarkerSize(1.5);
  gHighPtPP->Draw("sames PE");

  TGraphAsymmErrors *gHighPtPPSys = new TGraphAsymmErrors(npp, xpp, ypp, sxlpp, sxlpp, sylpp, syhpp);
  gHighPtPPSys->SetName("Jpsi_xsec_pp200_highPt_systematics");
  gHighPtPPSys->SetMarkerStyle(20);
  gHighPtPPSys->SetMarkerSize(0);
  gHighPtPPSys->SetFillStyle(0);
  gHighPtPPSys->Draw("sameE5");

  // Run 12 data
  TFile *fjpsi = TFile::Open(Form("Rootfiles/%s/jpsi_xsec_pp200_run12.root",run_cfg_name.Data()),"read");
  TGraphErrors *gRun12ppSys = (TGraphErrors*)fjpsi->Get("gJpsiXsecCombSys_run12");
  gRun12ppSys->SetMarkerStyle(20);
  gRun12ppSys->SetMarkerSize(0);
  gRun12ppSys->SetFillStyle(0);
  gRun12ppSys->Draw("sames E5");

  TGraphErrors *gRun12pp = (TGraphErrors*)fjpsi->Get("gJpsiXsecComb_run12");
  gRun12pp->SetMarkerStyle(29);
  gRun12pp->SetMarkerSize(2);
  gRun12pp->Draw("sames PE");

  // PHENIX measurements
  TGraphAsymmErrors *gPhenixSys = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix_Systematics");
  gPhenixSys->SetMarkerSize(0);
  gPhenixSys->SetFillStyle(0);
  gPhenixSys->SetLineColor(4);
  gPhenixSys->Draw("sameE5");

  TGraphAsymmErrors *gPhenix = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix");
  gPhenix->SetMarkerStyle(24);
  gPhenix->SetMarkerSize(1.5);
  gPhenix->SetMarkerColor(4);
  gPhenix->SetLineColor(4);
  gPhenix->Draw("sames PE");

  leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  leg->AddEntry(gHighPtPP,"STAR 2009 HT |y|<1","PL");
  leg->AddEntry(gRun12pp,"STAR 2012 |y|<1","PL");
  leg->AddEntry(gPhenix,"PHENIX |y|<0.35","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Compare_pp_ref.pdf",run_type,run_cfg_name.Data()));

  TGraphAsymmErrors *hJpsiRebin[3];
  TGraphAsymmErrors *hJpsiRebinSys[3];
  const char *name[3] = {"STAR_2009","STAR_2012","PHENIX"};

  TCanvas *c = new TCanvas("fit_xsec","fit_xsec",1100,600);
  c->Divide(2,1);
  TF1 *funcJpsiXsec[3];

  TGraphAsymmErrors *graph = 0x0, *gSys = 0x0;
  double pT, y1, xh1, xl1, yh1, yl1, syh1, syl1;
  for(int i=0; i<3; i++)
    {
      hJpsiRebin[i] = new TGraphAsymmErrors(nbins);
      hJpsiRebin[i]->SetName(Form("Jpsi_xsec_%s",name[i]));
      hJpsiRebinSys[i] = new TGraphAsymmErrors(nbins);
      hJpsiRebinSys[i]->SetName(Form("Jpsi_xsec_%s_Sys",name[i]));
      if(i==0)
	{
	  graph = gHighPtPP;
	  gSys = gHighPtPPSys;
	}
      if(i==1)
	{
	  graph = gRun12pp;
	  gSys = gRun12ppSys;
	}
      if(i==2)
	{
	  graph = gPhenix;
	  gSys = gPhenixSys;
	}
      cout << name[i] << endl;
      double fraction = 1;
      if(i<2)
	{
	  c->cd(i+1);
	  hpp->DrawCopy("");
	  gPad->SetLogy();
	  graph->GetXaxis()->SetRangeUser(0,15);
	  graph->GetYaxis()->SetRangeUser(1e-6,10);
	  funcJpsiXsec[i] = new TF1(Form("Func_Jpsi_xsec_%s",name[i]),"exp([0]+[1]*x+[2]*x*x)",4,15);
	  graph->Fit(funcJpsiXsec[i], "IR0Q");
	  graph->Draw("samesPE");
	  funcJpsiXsec[i]->Draw("sames");
	  fraction = funcJpsiXsec[i]->Integral(10,14)/funcJpsiXsec[i]->Integral(10,15);
	  TPaveText *t1 = GetTitleText(name[i],0.05);
	  t1->Draw();
	}
      for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
	{
	  graph->GetPoint(ipoint, pT, y1);
	  xh1 = graph->GetErrorXhigh(ipoint);
	  xl1 = graph->GetErrorXlow(ipoint);
	  yh1 = graph->GetErrorYhigh(ipoint);
	  yl1 = graph->GetErrorYlow(ipoint);
	  syh1 = gSys->GetErrorYhigh(ipoint);
	  syl1 = gSys->GetErrorYlow(ipoint);
	  if(y1<=0) continue;
	  //printf("[i] pT = %4.3f, stat = (-%4.3f%%, +%4.3f%%), sys = (-%4.3f%%, +%4.3f%%)\n",pT,yl1/y1*100,yh1/y1*100,syl1/y1*100,syh1/y1*100);
	}
      
      for(int bin=1; bin<=nbins; bin++)
	{
	  double x = 0, y = 0, xh = 0, xl = 0, yh = 0, yl = 0, syh = 0, syl = 0;
	  double dpT = 0;
	  x  = (xbins[bin]+xbins[bin-1])/2.;
	  xh = (xbins[bin]-xbins[bin-1])/2.;
	  xl = (xbins[bin]-xbins[bin-1])/2.;
	  for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
	    {
	      graph->GetPoint(ipoint, pT, y1);
	      if(pT<xbins[bin-1] || pT>xbins[bin]) continue;
	      xh1 = graph->GetErrorXhigh(ipoint);
	      xl1 = graph->GetErrorXlow(ipoint);
	      yh1 = graph->GetErrorYhigh(ipoint);
	      yl1 = graph->GetErrorYlow(ipoint);
	      syh1 = gSys->GetErrorYhigh(ipoint);
	      syl1 = gSys->GetErrorYlow(ipoint);
	      if(i<=1)
		{
		  if(pT<4) dpT = 0.5;
		  else if(pT<8) dpT = 1;
		  else dpT = 2;
		}
	      else
		{
		  if(pT<5) dpT = 0.25;
		  else     dpT = 1;
		}
	      double yield = y1 * dpT * pT;
	      y += yield;
	      yh += pow(yh1*dpT*pT, 2);
	      yl += pow(yl1*dpT*pT, 2);
	      syh += syh1/y1 * yield;
	      syl += syl1/y1 * yield;
	    }
	  syh /= y;
	  syl /= y;
	  y = y / (xbins[bin]-xbins[bin-1]) / x;
	  yh = sqrt(yh) / (xbins[bin]-xbins[bin-1]) / x;
	  yl = sqrt(yl) / (xbins[bin]-xbins[bin-1]) / x;
	  syh *= y;
	  syl *= y;
	  if(y==0) continue;
	  if(i==2 && x>8) continue;
	  if(xbins[bin]>14) 
	    {
	      y/ = fraction;
	      yl/ = fraction;
	      yh/ = fraction;
	      syl/ = fraction;
	      syh/ = fraction;
	      cout << fraction << endl;
	    }
	  hJpsiRebin[i]->SetPoint(bin-1, x, y);
	  hJpsiRebin[i]->SetPointError(bin-1, 0, 0, yl, yh);
	  hJpsiRebinSys[i]->SetPoint(bin-1, x, y);
	  hJpsiRebinSys[i]->SetPointError(bin-1, xl, xh, syl, syh);
	  //printf("[i] Rebin: pT = %4.3f, y = %4.3e, stat = (-%4.3f%%, +%4.3f%%), sys = (-%4.3f%%, +%4.3f%%)\n",x,y,yl/y*100,yh/y*100,syl/y*100,syh/y*100);
	}

      hJpsiRebinSys[i]->SetMarkerSize(0);
      hJpsiRebinSys[i]->SetFillStyle(0);
      hJpsiRebinSys[i]->SetLineColor(pow(2,i));

      hJpsiRebin[i]->SetMarkerColor(pow(2,i));
      hJpsiRebin[i]->SetMarkerSize(1.5);
      hJpsiRebin[i]->SetLineColor(pow(2,i));
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Fit_pp_ref.pdf",run_type,run_cfg_name.Data()));
  hJpsiRebin[0]->SetMarkerStyle(20);
  hJpsiRebin[1]->SetMarkerStyle(29);
  hJpsiRebin[1]->SetMarkerSize(2);
  hJpsiRebin[2]->SetMarkerStyle(24);

  TH1F *hppRebin = new TH1F("pp200_Jpsi_Rebin",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hppRebin->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hppRebin,"",kTRUE);
  for(int i=0; i<3; i++)
    {
      hJpsiRebinSys[i]->Draw("sameE5");
      hJpsiRebin[i]->Draw("samePE");
    }

  // take the average of three
  TH1F *hPPJpsiFinal = new TH1F("hPPJpsiFinal","",nbins,xbins);

  TGraphAsymmErrors *hPPJpsiFinalSys = new TGraphAsymmErrors(nbins);
  hPPJpsiFinalSys->SetName("hPPJpsiFinalSys");

  double xsec_0 = 0, xsec_0_err = 0;
  double xsec_5 = 0, xsec_5_err = 0;
  printf("+++ Final avergae +++\n");
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      x = 0, y = 0, xh = 0, xl = 0, yh = 0, yl = 0, syh = 0, syl = 0;
      x  = (xbins[ipoint+1]+xbins[ipoint])/2.;
      xh = (xbins[ipoint+1]-xbins[ipoint])/2.;
      xl = (xbins[ipoint+1]-xbins[ipoint])/2.;

      double total_weight = 0;
      for(int i=0; i<3; i++)
	{
	  hJpsiRebin[i]->GetPoint(ipoint, pT, y1);
	  if(y1<=0) continue;
	  yh1 = hJpsiRebin[i]->GetErrorYhigh(ipoint);
	  yl1 = hJpsiRebin[i]->GetErrorYlow(ipoint);
	  syh1 = hJpsiRebinSys[i]->GetErrorYhigh(ipoint);
	  syl1 = hJpsiRebinSys[i]->GetErrorYlow(ipoint);
	  double weight = 1./(pow((yh1+yl1)/2,2) + pow((syh1+syl1)/2,2));
	  total_weight += weight;
	  y += y1 * weight;
	  yh += weight;
	}
      y /= total_weight;
      yh = 1./sqrt(yh);
      hPPJpsiFinal->SetBinContent(ipoint+1, y);
      hPPJpsiFinal->SetBinError(ipoint+1, yh);
      double yield = y * 2 * pi * x * (xbins[ipoint+1]-xbins[ipoint]);
      double yield_err = yh/y * yield;
      if(x>0) 
	{
	  xsec_0 += yield;
	  xsec_0_err += pow(yield_err, 2);
	}
      if(x>5) 
	{
	  xsec_5 += yield;
	  xsec_5_err += pow(yield_err, 2);
	}
      printf("[i] Rebin: pT = %4.3f, y = %4.3e, stat = (-%4.3f%%, +%4.3f%%), sys = (-%4.3f%%, +%4.3f%%)\n",x,y,yh/y*100,yh/y*100,syl/y*100,syh/y*100);
    }
  hPPJpsiFinal->SetMarkerStyle(21);
  hPPJpsiFinal->SetMarkerSize(1.5);
  hPPJpsiFinal->SetMarkerColor(kGreen+2);
  hPPJpsiFinal->SetLineColor(kGreen+2);
  hPPJpsiFinal->Draw("samesPE");

  leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  leg->AddEntry(hJpsiRebin[0],"STAR 2009 HT |y|<1","PL");
  leg->AddEntry(hJpsiRebin[1],"STAR 2012 |y|<1","PL");
  leg->AddEntry(hJpsiRebin[2],"PHENIX |y|<0.35","PL");
  leg->AddEntry(hPPJpsiFinal,"Combined","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Combine_pp_ref.pdf",run_type,run_cfg_name.Data()));

  xsec_0_err = sqrt(xsec_0_err);
  xsec_5_err = sqrt(xsec_5_err);
  printf("[i] pT > 0, xsec = %4.2f #pm %4.2f\n",xsec_0,xsec_0_err);
  printf("[i] pT > 5, xsec = %4.2e #pm %4.2e\n",xsec_5,xsec_5_err);
  TH1F *hJpsiIntXsec = new TH1F("pp200_Jpsi_Integrated",";p_{T} (GeV/c);Bd#sigma/dy [nb]",2,0,2);
  hJpsiIntXsec->SetBinContent(1, xsec_0);
  hJpsiIntXsec->SetBinError(1, xsec_0_err);
  hJpsiIntXsec->SetBinContent(2, xsec_5);
  hJpsiIntXsec->SetBinError(2, xsec_5_err);
  
  if(saveHisto)
    {
      fout->cd();
      hJpsiIntXsec->Write("",TObject::kOverwrite);
      hPPJpsiFinal->Write("",TObject::kOverwrite);
      hPPJpsiFinalSys->Write("",TObject::kOverwrite);
    }
}


void scaleGraph(TGraphAsymmErrors *gr, double scale)
{
  int npoints = gr->GetN();
  Double_t x,y;
  for(int i=0; i<npoints; i++)
    {
      gr->GetPoint(i,x,y);
      double err_xlow = gr->GetErrorXlow(i);
      double err_xhigh = gr->GetErrorXhigh(i);
      double err_ylow = gr->GetErrorYlow(i);
      double err_yhigh = gr->GetErrorYhigh(i);
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,err_xlow,err_xhigh,err_ylow*scale,err_yhigh*scale);
    }
}

void drawBoxes(Int_t n, Double_t* x, Double_t dx, Double_t* y, Double_t* dy, Int_t lineWidth, Int_t lineColor){
  for(int i=0;i<n;i++){
    TBox *box = new TBox(x[i]-dx,y[i]-dy[i],x[i]+dx,y[i]+dy[i]);
    box->SetLineWidth(lineWidth);
    box->SetLineColor(lineColor);
    box->SetFillStyle(0);
    box->Draw("lsame");
  }
}

