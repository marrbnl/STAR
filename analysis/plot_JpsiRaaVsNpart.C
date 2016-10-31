const double ppInelastic = 42.; // mb
const double ppInelasticErr = 3.; // mb

//================================================
void plot_JpsiRaaVsNpart()
{  
  TFile *fin = TFile::Open("2016sQM.JpsiRaaVsNpart.root","read");

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
  TH1F *hRaaVsNpart = new TH1F("hRaaVsNpart",";N_{part};J/#psi R_{AA}",100,-1,370);
  ScaleHistoTitle(hRaaVsNpart,28,1,24,28,1,24,63);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,2);
  hRaaVsNpart->GetYaxis()->CenterTitle();

  const int npart_color[3] = {1, 4, 2};

  // STAR published
  TGraphErrors *pubRaaVsNpart = (TGraphErrors*)fin->Get("STAR_JpsiRaaVsNpart5GeV_dielec");
  pubRaaVsNpart->SetMarkerStyle(20);
  pubRaaVsNpart->SetMarkerColor(npart_color[0]);
  pubRaaVsNpart->SetMarkerSize(1.8);
  pubRaaVsNpart->SetLineColor(npart_color[0]);
  TGraphErrors *pubRaaVsNpartSys = (TGraphErrors*)fin->Get("STAR_JpsiRaaVsNpart5GeV_dielec_Sys");
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
  TH1F *hpp = (TH1F*)fin->Get("pp200_Jpsi_Integrated");
  TGraphErrors *raaVsNpart[2], *raaVsNpartSys[2];
  TBox *globalSys[2];
  for(int i=0; i<2; i++)
    {
      raaVsNpart[i] = (TGraphErrors*)fin->Get(Form("MTD_Run14AuAu_JpsiRaaVsNpart%dGeV",i*5));
      raaVsNpartSys[i] = (TGraphErrors*)fin->Get(Form("MTD_Run14AuAu_JpsiRaaVsNpart%dGeV_Sys",i*5));

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

      double pp_yield = hpp->GetBinContent(i+1);
      double pp_err = hpp->GetBinError(i+1);
      double gSys = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err*pp_err/pp_yield/pp_yield);
      globalSys[i] = new TBox(350+i*5,1-gSys,355+i*5,1+gSys);
      globalSys[i]->SetLineColor(npart_color[i+1]);
      globalSys[i]->SetFillColor(npart_color[i+1]);
      globalSys[i]->SetLineWidth(2.);
      globalSys[i]->SetFillStyle(1001);
    }

  //==============================================
  // Raa vs Npart dimuon
  //==============================================
  TCanvas *c0 = new TCanvas("Jpsi_raa_vs_npart_dimuon","Jpsi_raa_vs_npart_dimuon",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
  line->Draw();
  for(int i=0; i<2; i++)
    {
      raaVsNpart[i]->Draw("samesPEZ");
      raaVsNpartSys[i]->Draw("samesE5");
      globalSys[i]->Draw("fsame");
    }
  TLegend *leg0 = new TLegend(0.16,0.75,0.4,0.95);
  leg0->SetBorderSize(0);
  leg0->SetFillColor(0);
  leg0->SetTextFont(62);
  leg0->SetTextSize(0.04);
  leg0->SetHeader("Au+Au @ 200 GeV");
  leg0->AddEntry(raaVsNpart[0],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 0 GeV/c","P");
  leg0->AddEntry(raaVsNpart[1],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 5 GeV/c","P");
  leg0->Draw();
  TLegend *leg_ncoll = new TLegend(0.15,0.17,0.35,0.22);
  leg_ncoll->SetBorderSize(0);
  leg_ncoll->SetFillColor(0);
  leg_ncoll->SetTextFont(62);
  leg_ncoll->SetTextSize(0.035);
  leg_ncoll->AddEntry(gPubNpartErr,"N_{coll} uncertainty","f");
  leg_ncoll->Draw();
  TPaveText *star = GetPaveText(0.75,0.85,0.92,0.97,26,23);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();


  //==============================================
  // Raa vs Npart dimuon vs dielectron
  //==============================================
  TCanvas *c1 = new TCanvas("Jpsi_raa_vs_npart","Jpsi_raa_vs_npart",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
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


  //==============================================
  // Raa vs model
  //==============================================
  // tsinghua group
  TGraph *gTsuRaaVsNpart[2];
  gTsuRaaVsNpart[0] = (TGraph*)fin->Get("Tsinghua_RHIC_JpsiRaaVsNpart_LowPt");
  gTsuRaaVsNpart[1] = (TGraph*)fin->Get("Tsinghua_RHIC_JpsiRaaVsNpart_HighPt");
  for(int i=0; i<2; i++)
    {
      gTsuRaaVsNpart[i]->SetMarkerStyle(21);
      gTsuRaaVsNpart[i]->SetMarkerColor(npart_color[i+1]);
      gTsuRaaVsNpart[i]->SetLineColor(npart_color[i+1]);
      gTsuRaaVsNpart[i]->SetLineWidth(2);
      gTsuRaaVsNpart[i]->SetLineStyle(1);
    }
  
  // tamu group
  TGraph *gTamuRaaVsNpart[2];
  gTamuRaaVsNpart[0] = (TGraph*)fin->Get("TAMU_RHIC_JpsiRaaVsNpart_LowPt");
  gTamuRaaVsNpart[1] = (TGraph*)fin->Get("TAMU_RHIC_JpsiRaaVsNpart_HighPt");

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
  
  TPaveText *tmodel = GetPaveText(0.21,0.38,0.64,0.74);
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

  //==============================================
  // STAR vs LHC low pT
  //==============================================
  
  // PHENIX low pT
  // http://www.phenix.bnl.gov/phenix/WWW/info/data/ppg068/auau_electrons.txt
  // PRL 98 (2007) 232301
  // global uncertainty = 12%
  TGraphErrors *grPhenixLowPt = (TGraphErrors*)fin->Get("PHENIX_JpsiRaaVsNpart_LowPt");
  TGraphErrors *grPhenixLowPtSys  = (TGraphErrors*)fin->Get("PHENIX_JpsiRaaVsNpart_LowPt_Sys");
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
  TGraphErrors *grAliceLowPt = (TGraphErrors*)fin->Get("ALICE_JpsiRaaVsNpart_LowPt");
  TGraphErrors *grAliceLowPtSys  = (TGraphErrors*)fin->Get("ALICE_JpsiRaaVsNpart_LowPt_Sys");
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
  TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
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

  TLegend *leg3 = new TLegend(0.16,0.75,0.4,0.95);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->SetTextFont(62);
  leg3->SetTextSize(0.035);
  leg3->SetHeader("p_{T,J/#psi} > 0 GeV/c");
  leg3->AddEntry(grStarLowPt,"STAR: Au+Au, #sqrt{s_{NN}} = 200 GeV |y| < 0.5","P");
  leg3->AddEntry(grPhenixLowPt,"PHENIX: Au+Au, #sqrt{s_{NN}} = 200 GeV |y| < 0.35","P");
  leg3->AddEntry(grAliceLowPt,"ALICE: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV |y| < 0.8","P");
  leg3->Draw();

  // Add in models
  TGraph *gTamuRaaVsNpartRHIC[2], *gTamuRaaVsNpartLHC[2];
  for(int i=0; i<2; i++) gTamuRaaVsNpartRHIC[i] = (TGraph*)gTamuRaaVsNpart[i]->Clone(Form("gTamuRaaVsNpart_%d",i));

  //Zebo/QWG2013/Jpsi_Cent/Jpsi_Rapp_mid.dat
  gTamuRaaVsNpartLHC[0] = (TGraph*)fin->Get("TAMU_LHC_JpsiRaaVsNpart_LowPt");

  //Zebo/QWG2013/Jpsi_Cent/Jpsi_Rapp_LHC_mid_highPt.dat
  gTamuRaaVsNpartLHC[1] = (TGraph*)fin->Get("TAMU_LHC_JpsiRaaVsNpart_HighPt");

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
  gTsuRaaVsNpartLHC[0] = (TGraphErrors*)fin->Get("Tsinghua_LHC_JpsiRaaVsNpart_LowPt");
  gTsuRaaVsNpartLHC[1] =(TGraphErrors*)fin->Get("Tsinghua_LHC_JpsiRaaVsNpart_HighPt");

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
  gTsuRaaVsNpartLHC[0]->Draw("samesE4");

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


  //==============================================
  // STAR vs LHC high pT
  //==============================================
  TGraphErrors *cmsRaaVsNpart = (TGraphErrors*)fin->Get("CMS_JpsiRaaVsNpart_HighPt");
  TGraphErrors *cmsRaaVsNpartSys = (TGraphErrors*)fin->Get("CMS_JpsiRaaVsNpart_HighPt_Sys");
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
  hRaaVsNpart->SetMaximum(1.8);
  hRaaVsNpart->DrawCopy();
  gPubNpartErr->Draw("e3sames");
  TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  cmsRaaVsNpartSys->Draw("samesE5");
  cmsRaaVsNpart->Draw("samesPEZ");
  bCmsHighPt->Draw("fsame");

  grStarHighPtSys->Draw("samesE5");
  grStarHighPt->Draw("samesPEZ");
  bStarHighPt->Draw("fsame");

  leg_ncoll_2->Draw();

  TLegend *leg4 = new TLegend(0.16,0.85,0.4,0.95);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(0);
  leg4->SetTextFont(62);
  leg4->SetTextSize(0.035);
  leg4->AddEntry(grStarHighPt,"STAR: Au+Au, #sqrt{s_{NN}} = 200 GeV, |y| < 0.5, p_{T} > 5 GeV/c","P");
  leg4->AddEntry(cmsRaaVsNpart,"CMS: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 2.4, p_{T} > 6.5 GeV/c","P");
  leg4->Draw();
  TPaveText *star_2 = GetPaveText(0.75,0.85,0.68,0.73,26,23);
  star_2->AddText("STAR preliminary");
  star_2->SetTextColor(2);
  star_2->Draw();

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
  star_2->Draw();

  TPaveText *tmodel2 = GetPaveText(0.21,0.38,0.73,0.83);
  tmodel2->SetTextFont(62);
  tmodel2->SetTextAlign(11);
  tmodel2->SetTextSize(0.035);
  tmodel2->AddText("Tsinghua Model");
  tmodel2->AddText("TAMU Model");
  tmodel2->Draw();

  TLegend *leg41 = new TLegend(0.42,0.74,0.6,0.84);
  leg41->SetBorderSize(0);
  leg41->SetFillColor(0);
  leg41->SetTextFont(62);
  leg41->SetTextSize(0.035);
  leg41->AddEntry(gTsuRaaVsNpartRHIC[1],"RHIC","L");
  leg41->AddEntry(gTamuRaaVsNpartRHIC[1],"RHIC","L");
  leg41->Draw();
  TLegend *leg42 = new TLegend(0.55,0.74,0.7,0.84);
  leg42->SetBorderSize(0);
  leg42->SetFillColor(0);
  leg42->SetTextFont(62);
  leg42->SetTextSize(0.035);
  leg42->AddEntry(gTsuRaaVsNpartLHC[1],"LHC","L");
  leg42->AddEntry(gTamuRaaVsNpartLHC[1],"LHC","L");
  leg42->Draw();
}

//-----------------------------------------
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size = 0.04, const Int_t font = 42)
{
  TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}


//-----------------------------------------
void ScaleHistoTitle(TH1 *h, const Double_t xTitleSize = 28, const Double_t xTitleOffset = 0.9, const Double_t xLabelSize = 20,
		     const Double_t yTitleSize = 28, const Double_t yTitleOffset = 0.9, const Double_t yLabelSize = 20,
		     const Int_t font = 42)
{
  if(!h) return;
  h->GetXaxis()->SetTitleFont(font);
  h->GetXaxis()->SetLabelFont(font);
  h->GetYaxis()->SetTitleFont(font);
  h->GetYaxis()->SetLabelFont(font);

  h->GetXaxis()->SetTitleSize(xTitleSize);
  h->GetXaxis()->SetTitleOffset(xTitleOffset);
  h->GetXaxis()->SetLabelSize(xLabelSize);
  h->GetYaxis()->SetTitleSize(yTitleSize);
  h->GetYaxis()->SetTitleOffset(yTitleOffset);
  h->GetYaxis()->SetLabelSize(yLabelSize);
}

//-----------------------------------------
void SetPadMargin(TVirtualPad *pad, const Double_t bottomMargin = 0.12, const Double_t leftMargin = 0.12, const Double_t rightMargin = 0.05, const Double_t topMargin = 0.10)
{
  if(!pad) return;
  pad->SetLeftMargin(leftMargin);
  pad->SetBottomMargin(bottomMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetTopMargin(topMargin);
}

//--------------------------------------------
TLine *GetLine(Double_t xl, Double_t yl, Double_t xh, Double_t yh, Color_t color=2, Width_t width=2, Style_t style=2)
{
  TLine *line = new TLine(xl,yl,xh,yh);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  line->SetLineStyle(style);
  return line;
}
