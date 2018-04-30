void ana_Upsilon()
{
  gStyle->SetOptStat(0);
  //AuAu();
  //ppdata();
  //Xinjie();
  //ptShape();
  abnomalPoint();
}

//================================================
void abnomalPoint(const int savePlot = 1)
{
  const int year = 2016;
  TFile *fin = TFile::Open(Form("output/Run%d_AuAu200.InvMassStudy.root",year-2000));
  THnSparseF *hnInvMass = (THnSparseF*)fin->Get("mhnInvMass_di_mu");
  THnSparseF *hnInvMassMtd = (THnSparseF*)fin->Get("mhnInvMassMtd_di_mu");

  // invariant mass
  const char *pairType[2] = {"LS","UL"};
  const int rebin[5] = {1, 2, 5, 10, 20};
  TH1F *hInvMass[2][5];
  for(int i=0; i<2; i++)
    {
      if(year==2014) hnInvMass->GetAxis(0)->SetRange(i+1, i+1);
      if(year==2016) hnInvMass->GetAxis(0)->SetRange(2-i, 2-i);
      //hnInvMass->GetAxis(2)->SetRangeUser(0,10);
      for(int j=0; j<5; j++)
	{
	  hInvMass[i][j] = (TH1F*)hnInvMass->Projection(1);
	  hInvMass[i][j]->SetName(Form("hInvMass_%s_rebin%d",pairType[i],rebin[j]));
	  hInvMass[i][j]->Rebin(rebin[j]);
	}
      hnInvMass->GetAxis(0)->SetRange(0,-1);
    }
  TString legName[2] = {"LS", "UL"};
  TList *list = new TList();
  for(int j=0; j<5; j++)
    {
      for(int i=0; i<2; i++)
	{
	  list->Add(hInvMass[i][j]);
	}
      c = drawHistos(list,Form("cInvMass_rebin%d",rebin[j]),Form("Invariant mass distribution (BinWidth = %1.0f MeV/c^{2});M_{#mu#mu} [GeV/c^{2}];counts",hInvMass[0][j]->GetBinWidth(1)*1000),true,8,12,false,0.4,1.1,false,true,legName,true,Form("%d",year),0.65,0.75,0.7,0.85);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_InvMass_Rebin%d.pdf",run_type,year,rebin[j]));
      list->Clear();
    }

  // Daughter track kinematics
  const int nComp = 3;
  const int compAxis[nComp] = {2, 2, 1};
  const double compMinMass[nComp] = {8.8, 9, 9};
  const double compMaxMass[nComp] = {8.85, 11, 11};
  TH1F *hPairNVsDay[nComp];
  TH1F *hPairPt[nComp];
  TH1F *hPairEta[nComp];
  TH1F *hPairPhi[nComp];
  TH2F *hMuonPtCor[nComp];
  TH2F *hMuonEtaCor[nComp];
  TH2F *hMuonPhiCor[nComp];
  TH1F *hMuonBL[nComp];
  TH1F *hMuonMod[nComp];
  TH1F *hMuonCell[nComp];
  for(int i=0; i<nComp; i++)
    {
      hnInvMass->GetAxis(0)->SetRange(compAxis[i], compAxis[i]);
      hnInvMass->GetAxis(1)->SetRangeUser(compMinMass[i], compMaxMass[i]);

      hPairNVsDay[i] = (TH1F*)hnInvMass->Projection(11);
      hPairNVsDay[i]->SetName(Form("hPairNVsDay_comp%d",i));
      hPairNVsDay[i]->Sumw2();

      hPairPt[i] = (TH1F*)hnInvMass->Projection(2);
      hPairPt[i]->SetName(Form("hPairPt_comp%d",i));
      hPairPt[i]->Sumw2();

      hPairEta[i] = (TH1F*)hnInvMass->Projection(3);
      hPairEta[i]->SetName(Form("hPairEta_comp%d",i));
      hPairEta[i]->Sumw2();

      hPairPhi[i] = (TH1F*)hnInvMass->Projection(4);
      hPairPhi[i]->SetName(Form("hPairPhi_comp%d",i));
      hPairPhi[i]->Sumw2();

      hMuonPtCor[i] = (TH2F*)hnInvMass->Projection(5, 8);
      hMuonPtCor[i]->SetName(Form("hMuonPtCor_comp%d",i));
      hMuonPtCor[i]->Sumw2();

      hMuonEtaCor[i] = (TH2F*)hnInvMass->Projection(6, 9);
      hMuonEtaCor[i]->SetName(Form("hMuonEtaCor_comp%d",i));
      hMuonEtaCor[i]->Sumw2();

      hMuonPhiCor[i] = (TH2F*)hnInvMass->Projection(7, 10);
      hMuonPhiCor[i]->SetName(Form("hMuonPhiCor_comp%d",i));
      hMuonPhiCor[i]->Sumw2();

      hnInvMass->GetAxis(0)->SetRange(0,-1);
      hnInvMass->GetAxis(1)->SetRange(0,-1);

      // matched hit information
      hnInvMassMtd->GetAxis(0)->SetRange(compAxis[i], compAxis[i]);
      hnInvMassMtd->GetAxis(1)->SetRangeUser(compMinMass[i], compMaxMass[i]);

      hMuonBL[i] = (TH1F*)hnInvMassMtd->Projection(2);
      hMuonBL[i]->SetName(Form("hMuonBL_comp%d",i));
      hMuonBL[i]->Sumw2();

      hMuonMod[i] = (TH1F*)hnInvMassMtd->Projection(3);
      hMuonMod[i]->SetName(Form("hMuonMod_comp%d",i));
      hMuonMod[i]->Sumw2();

      hMuonCell[i] = (TH1F*)hnInvMassMtd->Projection(4);
      hMuonCell[i]->SetName(Form("hMuonCell_comp%d",i));
      hMuonCell[i]->Sumw2();

      hnInvMassMtd->GetAxis(0)->SetRange(0,-1);
      hnInvMassMtd->GetAxis(1)->SetRange(0,-1);
    }

  TString legName2[nComp];
  for(int i=0; i<nComp; i++)
    {
      hPairNVsDay[i]->Rebin(2);
      if(i>0) hPairNVsDay[i]->Scale(hPairNVsDay[0]->Integral()/hPairNVsDay[i]->Integral());
      list->Add(hPairNVsDay[i]);
      legName2[i] = Form("%s: [%1.2f, %1.2f] GeV/c^{2}",pairType[compAxis[i]-1], compMinMass[i], compMaxMass[i]);
    }
  c = drawHistos(list,"cPairNVsDay","# of muon pairs vs. day;Day;counts",false,-100,100,true,0,8,false,true,legName2,true,Form("%d",year),0.2,0.4,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_PairNVsDay.pdf",run_type,year));
  list->Clear();

  for(int i=0; i<nComp; i++)
    {
      if(i>0) hPairPt[i]->Scale(hPairPt[0]->Integral()/hPairPt[i]->Integral());
      list->Add(hPairPt[i]);
    }
  c = drawHistos(list,"cPairPt","Pair p_{T} distribution;p_{T}^{#mu#mu} [GeV/c];counts",false,-100,100,false,0.4,1.1,false,true,legName2,true,Form("%d",year),0.5,0.7,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_PairPt.pdf",run_type,year));
  list->Clear();

  TString legName3[nComp];
  for(int i=0; i<nComp; i++)
    {
      double ratio = (hPairEta[i]->Integral(0,0) + hPairEta[i]->Integral(hPairEta[i]->GetNbinsX()+1, hPairEta[i]->GetNbinsX()+1))*1./hPairEta[i]->Integral(0,-1);
      if(i>0) hPairEta[i]->Scale(hPairEta[0]->Integral()/hPairEta[i]->Integral());
      list->Add(hPairEta[i]);
      legName3[i] = Form("%s: |#eta_{#mu#mu}| > 1 = %3.2f%%",legName2[i].Data(), ratio*100);
    }
  c = drawHistos(list,"cPairEta","Pair #eta distribution;#eta_{#mu#mu};counts",false,-100,100,false,0.4,1.1,false,true,legName3,true,Form("%d",year),0.2,0.4,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_PairEta.pdf",run_type,year));
  list->Clear();

  for(int i=0; i<nComp; i++)
    {
      hPairPhi[i]->Rebin(2);
      if(i>0) hPairPhi[i]->Scale(hPairPhi[0]->Integral()/hPairPhi[i]->Integral());
      list->Add(hPairPhi[i]);
    }
  c = drawHistos(list,"cPairPhi","Pair #varphi distribution;#varphi_{#mu#mu};counts",false,-100,100,true,0,10,false,true,legName2,true,Form("%d",year),0.5,0.7,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_PairPhi.pdf",run_type,year));
  list->Clear();

  c = new TCanvas("cMuonPtCor","cMuonPtCor",900,600);
  c->Divide(2,2);
  for(int i=0; i<nComp; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad, 0.13, 0.1, 0.1, 0.13);
      hMuonPtCor[i]->SetTitle(";p_{T,2} [GeV/c];p_{T,1} [GeV/c]");
      hMuonPtCor[i]->GetXaxis()->SetRangeUser(1, 6);
      hMuonPtCor[i]->GetYaxis()->SetRangeUser(3,10);
      ScaleHistoTitle(hMuonPtCor[i],0.05,1,0.045,0.05,0.8,0.045);
      hMuonPtCor[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%d: %s",year,legName2[i].Data()),0.055);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_MuonPtCor.pdf",run_type,year));

  c = new TCanvas("cMuonEtaCor","cMuonEtaCor",900,600);
  c->Divide(2,2);
  for(int i=0; i<nComp; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad, 0.13, 0.1, 0.1, 0.13);
      hMuonEtaCor[i]->SetTitle(";#eta_{2};#eta_{1}");
      ScaleHistoTitle(hMuonEtaCor[i],0.05,1,0.045,0.05,0.8,0.045);
      hMuonEtaCor[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%d: %s",year,legName2[i].Data()),0.055);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_MuonEtaCor.pdf",run_type,year));

  c = new TCanvas("cMuonPhiCor","cMuonPhiCor",900,600);
  c->Divide(2,2);
  for(int i=0; i<nComp; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad, 0.13, 0.1, 0.1, 0.13);
      hMuonPhiCor[i]->SetTitle(";#varphi_{2};#varphi_{1}");
      ScaleHistoTitle(hMuonPhiCor[i],0.05,1,0.045,0.05,0.8,0.045);
      hMuonPhiCor[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%d: %s",year,legName2[i].Data()),0.055);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_MuonPhiCor.pdf",run_type,year));

  for(int i=0; i<nComp; i++)
    {
      if(i>0) hMuonBL[i]->Scale(hMuonBL[0]->Integral()/hMuonBL[i]->Integral());
      list->Add(hMuonBL[i]);
    }
  c = drawHistos(list,"cMuonBL","Backleg of associated MTD hits for muon daughters;BL;counts",false,-100,100,true,0,15,false,true,legName2,true,Form("%d",year),0.5,0.7,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_MuonBL.pdf",run_type,year));
  list->Clear();

  for(int i=0; i<nComp; i++)
    {
      if(i>0) hMuonMod[i]->Scale(hMuonMod[0]->Integral()/hMuonMod[i]->Integral());
      list->Add(hMuonMod[i]);
    }
  c = drawHistos(list,"cMuonMod","Module of associated MTD hits for muon daughters;Module;counts",false,-100,100,true,0,30,false,true,legName2,true,Form("%d",year),0.2,0.4,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_MuonMod.pdf",run_type,year));
  list->Clear();

  for(int i=0; i<nComp; i++)
    {
      if(i>0) hMuonCell[i]->Scale(hMuonCell[0]->Integral()/hMuonCell[i]->Integral());
      list->Add(hMuonCell[i]);
    }
  c = drawHistos(list,"cMuonCell","Strip number of associated MTD hits for muon daughters;Cell;counts",false,-100,100,true,0,20,false,true,legName2,true,Form("%d",year),0.2,0.4,0.65,0.85);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/Check_%d_MuonCell.pdf",run_type,year));
  list->Clear();
}


//================================================
void Xinjie(const int saveHisto = 0)
{
  TFile *fin = TFile::Open("Rootfiles/160606.all.root","read");
  const char *hName[3] = {"mhUpsInfoUL","mhUpsInfoLSPos","mhUpsInfoLSNeg"};
  THnSparseF *hnInvMass[3] = {0x0};
  TH1F *hInvMass[3] = {0x0};
  
  // same event
  char name[512];
  for(Int_t j=0; j<3; j++) // pair type
    { 
      hnInvMass[j] = (THnSparseF*)fin->Get(hName[j]);
      hnInvMass[j]->GetAxis(3)->SetRange(0,13);
      hnInvMass[j]->GetAxis(4)->SetRange(0,6);
      hnInvMass[j]->GetAxis(5)->SetRange(0,23);
      hnInvMass[j]->GetAxis(6)->SetRange(0,14);

      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("InvMassHft_jpsi_%d",j));
      hInvMass[j]->SetTitle("");
      //hInvMass[j]->Rebin(2);
      hInvMass[j]->Sumw2();
    }
  //hInvMass[1]->Add(hInvMass[2]);

  TH1F *hSeUL  = (TH1F*)hInvMass[0]->Clone("hSeUL");
  TH1F *hSeLS  = (TH1F*)hInvMass[1]->Clone("hSeLS");

  hSeUL->SetMarkerStyle(21);
  hSeUL->SetMarkerColor(2);
  hSeUL->SetLineColor(2);
  hSeUL->GetXaxis()->SetRangeUser(8,12);
  c = draw1D(hSeUL,Form("Dimuon events: invariant mass distribution"),false);
  hSeLS->Draw("sames HIST");

  TFile *fout = TFile::Open("Rootfiles/sQM2016.Upsilon.root","recreate");
  hInvMass[0]->Write("Upsilon_UL_cent0080");
  hInvMass[1]->Write("Upsilon_LS_Pos_cent0080");
  hInvMass[2]->Write("Upsilon_LS_Neg_cent0080");
  fout->Close();
}

//================================================
void AuAu(const int savePlot = 0)
{
  gStyle->SetOptStat(0);
  const int ipt = 0, icent = 0;

  TFile *fin = TFile::Open("output/Pico.Run14.AuAu200.jpsi.root","read");
  TH1F *hStat = (TH1F*)fin->Get("hEventStat");
  printf("+++ check this +++\n");
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3] = {0x0};
  TH1F *hInvMass[3] = {0x0};
  
  // same event
  char name[512];
  for(Int_t j=0; j<3; j++) // pair type
    { 
      hnInvMass[j] = (THnSparseF*)fin->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(3)->SetRangeUser(1.5+0.01,100);
      hnInvMass[j]->GetAxis(4)->SetRangeUser(1.5+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRange(centBins_low[icent],centBins_high[icent]);
      hnInvMass[j]->GetAxis(1)->SetRangeUser(ptBins_low[ipt]+0.01,ptBins_high[ipt]-0.01);

      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("InvMassHft_jpsi_%d",j));
      hInvMass[j]->Rebin(25);
      hInvMass[j]->Sumw2();
    }
  hInvMass[1]->Add(hInvMass[2]);

  TH1F *hSeUL  = (TH1F*)hInvMass[0]->Clone("hSeUL");
  TH1F *hSeLS  = (TH1F*)hInvMass[1]->Clone("hSeLS");
 
  hSeUL->SetMarkerStyle(21);
  hSeUL->SetMarkerColor(2);
  hSeUL->SetLineColor(2);
  hSeUL->GetXaxis()->SetRangeUser(6,14);
  c = draw1D(hSeUL,Form("Dimuon events: %1.1f < p_{T} < %1.1f GeV/c (%s%%)",ptBins_low[ipt],ptBins_high[ipt],cent_Name[icent]),true);
  hSeLS->Draw("sames HIST");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/AuAu200_UpsilonInvMass.pdf",run_type));

  TH1F *hSeDiff  = (TH1F*)hInvMass[0]->Clone("hSeDiff");
  hSeDiff->SetMarkerStyle(20);
  hSeDiff->GetXaxis()->SetRangeUser(6,14);
  hSeDiff->Add(hSeLS,-1); 
  c = draw1D(hSeDiff,Form("Dimuon events: %1.1f < p_{T} < %1.1f GeV/c (%s%%, US-LS)",ptBins_low[ipt],ptBins_high[ipt],cent_Name[icent]),false);
  TLine *line = GetLine(6,0,14,0);
  line->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/AuAu200_UpsilonInvMass_USminusLS.pdf",run_type));

}

//================================================
void ppdata()
{
  TFile *fin = TFile::Open("Rootfiles/Upsilon/ups.root","read");
  TTree *tree = (TTree*)fin->Get("tree");
  tree->Draw("ups_pt>>hups(100,0,10)","abs(ups_y)<0.5");
  TH1F *hups = (TH1F*)gDirectory->Get("hups");
  draw1D(hups);
  TFile *fout = TFile::Open("Rootfiles/Upsion.pp200.root","recreate");
  hups->Write("Upsilon_pp200_midRap");
}


//================================================
void ptShape()
{
  TF1 *fBol = new TF1("Boltzmann","x/(exp(x/[0]+1))",0,10);
  fBol->SetParameter(0,1.11);
  fBol->SetNpx(1000);
  TH1F *hBol = fBol->GetHistogram();
  TCanvas *c = new TCanvas("Boltzmann","Boltzmann",800,600);
  fBol->Draw();

  TFile *fin = TFile::Open("Rootfiles/Upsilon/Upsion.pp200.root","update");
  TH1F *hups = (TH1F*)fin->Get("Upsilon_pp200_midRap");
  hups->Sumw2();
  hups->Scale(1./hups->Integral());
  c = draw1D(hups);
  hBol->Scale(1./hBol->Integral());
  hBol->Draw("sames");
  fBol->Write("Upsilon_pp200_midRap_fBoltzmann");
  hBol->Write("Upsilon_pp200_midRap_hBoltzmann");
  
}

