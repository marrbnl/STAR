TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;

const char *run_config = "";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void ana_InvMass()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  if(year==2013)
    {
      run_type = "Run13_pp500";

      if(iPico)
	f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else
	f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";

      if(iPico)
	f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else
	f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }

  run_cfg_name = Form("%s",run_config);

  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  //Jpsi();
  upsilon();
  //daughters();
  //pt();
  //dimuon();
}

//================================================
void Jpsi(Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.4e\n",hStat->GetBinContent(9));
  printf("single-muon events: %4.4e\n",hStat->GetBinContent(10));
  printf("e-muon      events: %4.4e\n",hStat->GetBinContent(11));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3];
  Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("%s_%s_InvMass",hName[j],trigName[kTrigType]));
      hInvMass[j]->Sumw2();
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }
  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->SetYTitle("Counts");

  // signal
  TH1F *hInvMass_jpsi[2];
  for(Int_t j=0; j<2; j++)
    {
      hInvMass_jpsi[j] = (TH1F*)hInvMass[j]->Rebin(nSpecMBinsJpsi, Form("%s_jpsi",hInvMass[j]->GetName()),specMJpsi);
      scaleHisto(hInvMass_jpsi[j], 1, 1,kTRUE,kFALSE, kTRUE);
      hInvMass_jpsi[j]->GetXaxis()->SetRangeUser(1,3.5);
    }
  c = draw1D(hInvMass_jpsi[0],Form("%s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});dN/dM",trigName[kTrigType],hlt_name[hlt_index]));
  hInvMass_jpsi[1]->Draw("HIST sames");
  leg = new TLegend(0.15,0.63,0.3,0.85);
  if(hlt_index==1) leg = new TLegend(0.25,0.2,0.35,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass_jpsi[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass_jpsi[1],"Like sign (++)+(--)","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_jpsi_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_jpsi_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  // same event background
  Int_t low_bin = hInvMass_jpsi[0]->GetXaxis()->FindFixBin(low_mass+0.001);
  Int_t high_bin = hInvMass_jpsi[0]->GetXaxis()->FindFixBin(high_mass-0.001);
  double errorLS, errorUL;
  double nLS = hInvMass_jpsi[1]->IntegralAndError(low_bin,high_bin,errorLS,"width");
  double nUL = hInvMass_jpsi[0]->IntegralAndError(low_bin,high_bin,errorUL,"width");
  double nSignal = nUL - nLS;
  double nBackground = nLS;
  double errorSignal = TMath::Sqrt(errorLS*errorLS+errorUL*errorUL);
  TPaveText *signif = GetPaveText(0.13,0.3,0.55,0.75);
  signif->SetTextAlign(11);
  signif->SetTextFont(62);
  signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
  signif->AddText(Form("S/B = %1.0f/%1.0f = 1:%3.0f",nSignal,nBackground,nBackground/nSignal));
  signif->AddText(Form("Significance: %3.2f",nSignal/errorSignal));

  const double min = -8e3, max = 2.5e4;
  TH1F *hdiff_jpsi = (TH1F*)hInvMass_jpsi[0]->Clone(Form("InvMass_US_minus_LS_jpsi"));
  hdiff_jpsi->Add(hInvMass_jpsi[1],-1);
  hdiff_jpsi->GetYaxis()->SetRangeUser(min,max);
  hdiff_jpsi->GetXaxis()->SetRangeUser(2.5,3.5);
  //hdiff_jpsi->GetYaxis()->SetNdivisions(505);
  const int min_bin = hdiff_jpsi->FindFixBin(2.8);
  const int max_bin = hdiff_jpsi->FindFixBin(3.3);
  // for(int bin=min_bin; bin<=max_bin; bin++)
  //   {
  //     printf("%3.2f: %3.2f+/-%3.2f - %3.2f+/-%3.2f = %3.2f+/-%3.2f\n",
  // 	     hdiff_jpsi->GetBinCenter(bin),
  // 	     hInvMass_jpsi[0]->GetBinContent(bin)*hInvMass_jpsi[0]->GetBinWidth(bin),hInvMass_jpsi[0]->GetBinError(bin)*hInvMass_jpsi[0]->GetBinWidth(bin),
  // 	     hInvMass_jpsi[1]->GetBinContent(bin)*hInvMass_jpsi[1]->GetBinWidth(bin),hInvMass_jpsi[1]->GetBinError(bin)*hInvMass_jpsi[1]->GetBinWidth(bin),
  // 	     hdiff_jpsi->GetBinContent(bin)*hdiff_jpsi->GetBinWidth(bin),hdiff_jpsi->GetBinError(bin)*hdiff_jpsi->GetBinWidth(bin));
  //   }
  c = draw1D(hdiff_jpsi,Form("%s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});US-LS",trigName[kTrigType],hlt_name[hlt_index]),kFALSE);
  gPad->SetGridy();
  TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85,0.04);
  t1->AddText(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  t1->SetTextFont(62);
  t1->Draw();
  signif->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_US-LS_jpsi_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_US-LS_jpsi_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }
  return;
  
  // Mixed event background
  TFile *fmix = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.MixEvent.%sroot",run_config),"read");
  TH2F *hMixUL2D = (TH2F*)fmix->Get(Form("US_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt1_cut,pt2_cut,trigName[kTrigType]));
  hMixUL2D->Sumw2();
  TH2F *hMixLS2D = (TH2F*)fmix->Get(Form("LS_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt1_cut,pt2_cut,trigName[kTrigType]));
  hMixLS2D->Sumw2();
  TH1F *hMixUL = (TH1F*)hMixUL2D->ProjectionX("MixEvent_UL");
  TH1F *hMixLS = (TH1F*)hMixLS2D->ProjectionX("MixEvent_LS");
  hInvMass[1]->Rebin(10);
  hMixLS->Rebin(10);
  TH1F *hRatio = DivideTH1WithTH1(hInvMass[1],hMixLS);
  //(TH1F*)hInvMass[1]->Clone(Form("%s_ratio",hInvMass[1]->GetName()));
  //draw1D(hInvMass[1]);
  //draw1D(hMixLS);
  //hRatio->Rebin(10);
  //hMixLS->Rebin(10);
  hRatio->GetXaxis()->SetRangeUser(0,8);
  hRatio->GetYaxis()->SetRangeUser(1e-3,2.2e-3);
  hRatio->SetMarkerStyle(25);
  //hRatio->Divide(hMixLS);
  TF1 *func = new TF1(Form("%s_func",hRatio->GetName()),"pol0",2.8,3.2);
  hRatio->Fit(func,"IR0");
  c = draw1D(hRatio,Form("Like-sign: SameEvent/MixEvent (N_{++}+N_{--})"));
  func->SetLineColor(2);
  func->Draw("same");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sLS_SameToMix_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sLS_SameToMix_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  hMixUL->Scale(func->GetParameter(0));
  TH1F *hMixULBkg = (TH1F*)hMixUL->Rebin(nSpecMBinsJpsi, Form("%s_rebin",hMixUL->GetName()),specMJpsi);
  scaleHisto(hMixULBkg, 1, 1,kTRUE,kFALSE, kTRUE);
  hMixULBkg->SetLineColor(4);

  TH1F *hInvMassClone = (TH1F*)hInvMass_jpsi[0]->Clone(Form("%s_clone",hInvMass_jpsi[0]->GetName()));
  c = draw1D(hInvMassClone,Form("%s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});dN/dM",trigName[kTrigType],hlt_name[hlt_index]));
  hInvMass_jpsi[1]->Draw("HIST sames");
  hMixULBkg->Draw("HIST sames");
  leg = new TLegend(0.15,0.63,0.3,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass_jpsi[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass_jpsi[1],"Like sign (++)+(--)","L");
  leg->AddEntry(hMixULBkg,"Mixed: unlike-sign","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_jpsi_MixUL_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_jpsi_MixUL_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  TH1F *hdiff_jpsi_2 = (TH1F*)hInvMass_jpsi[0]->Clone(Form("InvMass_US_minus_LS_jpsi_Mix"));
  hdiff_jpsi_2->Add(hMixULBkg,-1);
  hdiff_jpsi_2->GetYaxis()->SetRangeUser(min,max);
  hdiff_jpsi_2->GetXaxis()->SetRangeUser(2.5,3.5);
  // for(int bin=min_bin; bin<=max_bin; bin++)
  //   {
  //     printf("%3.2f: %3.2f+/-%3.2f - %3.2f+/-%3.2f = %3.2f+/-%3.2f\n",
  // 	     hdiff_jpsi_2->GetBinCenter(bin),
  // 	     hInvMass_jpsi[0]->GetBinContent(bin)*hInvMass_jpsi[0]->GetBinWidth(bin),hInvMass_jpsi[0]->GetBinError(bin)*hInvMass_jpsi[0]->GetBinWidth(bin),
  // 	     hMixULBkg->GetBinContent(bin)*hMixULBkg->GetBinWidth(bin),hMixULBkg->GetBinError(bin)*hMixULBkg->GetBinWidth(bin),
  // 	     hdiff_jpsi_2->GetBinContent(bin)*hdiff_jpsi_2->GetBinWidth(bin),hdiff_jpsi_2->GetBinError(bin)*hdiff_jpsi_2->GetBinWidth(bin));
  //   }
  c = draw1D(hdiff_jpsi_2,Form("Invariant mass of di-muon pairs (mixed-event UL background);M_{#mu#mu} (GeV/c^{2});US-LS"),kFALSE);
  gPad->SetGridy();
  t1->Draw();
  double errorLS;
  double nLS = hMixULBkg->IntegralAndError(low_bin,high_bin,errorLS,"width");
  double nSignal = nUL - nLS;
  double nBackground = nLS;
  double errorSignal = TMath::Sqrt(errorLS*errorLS+errorUL*errorUL);
  TPaveText *signif = GetPaveText(0.13,0.3,0.55,0.75);
  signif->SetTextAlign(11);
  signif->SetTextFont(62);
  signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
  signif->AddText(Form("S/B = %1.0f/%1.0f = 1:%3.0f",nSignal,nBackground,nBackground/nSignal));
  signif->AddText(Form("Significance: %3.2f",nSignal/errorSignal));
  signif->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_US-MixUS_jpsi_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_US-MixUS_jpsi_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }
}



//================================================
void daughters(Int_t save = 0)
{
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH2F *hDaugPtCorr[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hDaugPtCorr[j] = (TH2F*)hnInvMass[j]->Projection(4,5);
      hDaugPtCorr[j]->SetName(Form("%s_DaugPtCorr",hnInvMass[j]->GetName()));
    }
  hDaugPtCorr[1]->Add(hDaugPtCorr[2]);
  const char *hTitle[2] = {"Unlike-sign","Like-sign"};
  const char *pName[2] = {"UL","LS"};
  for(Int_t j=0; j<2; j++)
    {
      c = draw2D(hDaugPtCorr[j],Form("%s: p_{T} correlation of di-muon pairs%s",trigName[kTrigType],hlt_name[hlt_index]));
      TPaveText *t1 = GetPaveText(0.6,0.8,0.3,0.4);
      t1->AddText(hTitle[j]);
      t1->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sDaugPtCorr_%s.pdf",run_type,run_cfg_name.Data(),pName[j]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sDaugPtCorr_%s.png",run_type,run_cfg_name.Data(),pName[j]));
	}
    }

  // pt cut on daughter
  Double_t pt_cut_2 = 1.0;
  TH2F *hInvMassVsPt[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt_cut_2+0.01,100);
      hInvMassVsPt[j] = (TH2F*)hnInvMass[j]->Projection(0,4);
      hInvMassVsPt[j]->SetName(Form("%s_%s_InvMassVsPt",hName[j],trigName[kTrigType]));
      hInvMassVsPt[j]->Sumw2();
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }


  const Int_t nCuts = 4;
  Double_t pt_cuts[nCuts] = {1.5,2.0,2.5,3.0};
  TH1F *hInvMass[3][nCuts];
  for(Int_t j=0; j<3; j++)
    {
      for(Int_t i=0; i<nCuts; i++)
	{
	  hInvMassVsPt[j]->GetXaxis()->SetRangeUser(pt_cuts[i]+0.01,100);
	  hInvMass[j][i] = (TH1F*)hInvMassVsPt[j]->ProjectionY(Form("%s_%s_InvMass_pt1%1.1f_pt2%1.1f",hName[j],trigName[kTrigType],pt_cuts[i],1.2));
	}
    }

  for(Int_t i=0; i<nCuts; i++)
    {
      hInvMass[0][i]->SetMarkerStyle(20);
      hInvMass[0][i]->SetMarkerColor(2);
      hInvMass[0][i]->SetLineColor(2);
      hInvMass[0][i]->GetXaxis()->SetRangeUser(0,4);
      c = draw1D(hInvMass[0][i],Form("%s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
      hInvMass[1][i]->Add(hInvMass[2][i]);
      hInvMass[1][i]->Draw("HIST sames");
      TLegend *leg = new TLegend(0.55,0.65,0.7,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt_cuts[i],pt_cut_2));
      leg->AddEntry(hInvMass[0][i],"Unlike sign","L");
      leg->AddEntry(hInvMass[1][i],"Like sign (++)+(--)","L");
      leg->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts[i]*10,pt_cut_2*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts[i]*10,pt_cut_2*10));
	}
    }


  TList *list = new TList;
  TString legName[nCuts];
  TH1F *hDiff[nCuts];
  for(Int_t i=0; i<nCuts; i++)
    {
      hDiff[i] = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US-LS_pt1%1.1f_pt2%1.1f",pt_cuts[i],pt_cut_2));
      hDiff[i]->Add(hInvMass[1][i],-1);
      hDiff[i]->Rebin(4);
      hDiff[i]->SetMaximum(0.5*hDiff[i]->GetMaximum());
      list->Add(hDiff[i]);

      Int_t low_bin = hInvMass[0][i]->GetXaxis()->FindFixBin(low_mass+0.001);
      Int_t high_bin = hInvMass[0][i]->GetXaxis()->FindFixBin(high_mass-0.001);
      Double_t nBackground = hInvMass[1][i]->Integral(low_bin,high_bin);
      Double_t nSignal = hInvMass[0][i]->Integral(low_bin,high_bin) - nBackground;
      legName[i] = Form("p_{T,1}>%1.1f GeV/c, S/B = %1.0f/%1.0f",pt_cuts[i],nSignal,nBackground);
    }
  c = drawHistos(list,"InvMass_ComparePt1Pt2Cuts",Form("%s: invariant mass of di-muon pairs (US-LS);M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kTRUE,2.0,3.8,kTRUE,-20,300,kFALSE,kTRUE,legName,kTRUE,Form("%1.1f < M_{#mu#mu} < %1.1f GeV/c^{2}",low_mass,high_mass),0.15,0.25,0.6,0.88,kTRUE);
  gPad->SetGridy();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_ComparePt1Cut.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_ComparePt1Cut.png",run_type,run_cfg_name.Data()));
    }
}


//================================================
void pt(Int_t save = 0)
{
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hJpsiPt[3];
  TH1F *hJpsiEta[3];
  TH1F *hJpsiPhi[3];
  Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);

      hJpsiPt[j] = (TH1F*)hnInvMass[j]->Projection(1);
      hJpsiPt[j]->SetName(Form("%s_%s_Jpsi_pt",hName[j],trigName[kTrigType]));
      hJpsiPt[j]->Sumw2();

      hJpsiEta[j] = (TH1F*)hnInvMass[j]->Projection(2);
      hJpsiEta[j]->SetName(Form("%s_%s_Jpsi_eta",hName[j],trigName[kTrigType]));
      hJpsiEta[j]->Sumw2();

      hJpsiPhi[j] = (TH1F*)hnInvMass[j]->Projection(3);
      hJpsiPhi[j]->SetName(Form("%s_%s_Jpsi_phi",hName[j],trigName[kTrigType]));
      hJpsiPhi[j]->Sumw2();

      hnInvMass[j]->GetAxis(0)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }

  TList *list = new TList;
  TString legName[3] = {"Unlike-sign","Like-sign","US-LS"};

  // pt distribution
  hJpsiPt[1]->Add(hJpsiPt[2]);
  TH1F *hDiffPt = (TH1F*)hJpsiPt[0]->Clone("hJpsiPt_US_minus_LS");
  hDiffPt->Add(hJpsiPt[1],-1);

  list->Add(hJpsiPt[0]);
  list->Add(hJpsiPt[1]);
  list->Add(hDiffPt);

  c = drawHistos(list,"Jpsi_pt",Form("%s: p_{T} distribution of J/psi candidates (%1.1f < M_{#mu#mu} < %1.1f GeV/c^{2});p_{T} (GeV/c);counts",trigName[kTrigType],low_mass,high_mass),kFALSE,0,100,kFALSE,-2,2,kFALSE,kTRUE,legName,kTRUE,Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt1_cut,pt2_cut),0.55,0.75,0.55,0.8,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  // eta distribution
  hJpsiEta[1]->Add(hJpsiEta[2]);
  TH1F *hDiffEta = (TH1F*)hJpsiEta[0]->Clone("hJpsiEta_US_minus_LS");
  hDiffEta->Add(hJpsiEta[1],-1);

  list->Clear();
  list->Add(hJpsiEta[0]);
  list->Add(hJpsiEta[1]);
  list->Add(hDiffEta);

  c = drawHistos(list,"Jpsi_eta",Form("%s: #eta distribution of J/psi candidates (%1.1f < M_{#mu#mu} < %1.1f GeV/c^{2});#eta;counts",trigName[kTrigType],low_mass,high_mass),kFALSE,0,100,kFALSE,-2,2,kFALSE,kTRUE,legName,kTRUE,Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt1_cut,pt2_cut),0.15,0.35,0.55,0.8,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_eta_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_eta_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  // phi distribution
  hJpsiPhi[1]->Add(hJpsiPhi[2]);
  TH1F *hDiffPhi = (TH1F*)hJpsiPhi[0]->Clone("hJpsiPhi_US_minus_LS");
  hDiffPhi->Add(hJpsiPhi[1],-1);

  list->Clear();
  list->Add(hJpsiPhi[0]);
  list->Add(hJpsiPhi[1]);
  list->Add(hDiffPhi);

  c = drawHistos(list,"Jpsi_phi",Form("%s: #varphi distribution of J/psi candidates (%1.1f < M_{#mu#mu} < %1.1f GeV/c^{2});#varphi;counts",trigName[kTrigType],low_mass,high_mass),kFALSE,0,100,kFALSE,-2,2,kFALSE,kTRUE,legName,kTRUE,Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt1_cut,pt2_cut),0.55,0.75,0.55,0.8,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_phi_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_phi_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

}


//================================================
void upsilon(Int_t save = 0)
{
  const int type = 0;

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));

  if(type==0) { const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"}; }
  else if(type==1) { const char *hName[3] = {"hUpsilonInfo","hUpsilonBkgLSPos","hUpsilonBkgLSNeg"}; }
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3];
  Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      if(type==0)
	{
	  hnInvMass[j]->GetAxis(1)->SetRangeUser(4,100);
	  hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
	}
      else if(type==1)
	{
	  hnInvMass[j]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[j]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
	}
      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("%s_%s_InvMass",hName[j],trigName[kTrigType]));
      hInvMass[j]->Sumw2();
    }
  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->SetYTitle("Counts");

  cout << hInvMass[0]->Integral(hInvMass[0]->FindFixBin(9.1),hInvMass[0]->FindFixBin(10.9)) << endl;
  cout << hInvMass[1]->Integral(hInvMass[1]->FindFixBin(9.1),hInvMass[1]->FindFixBin(10.9)) << endl;

  TH1F *hDiff = (TH1F*)hInvMass[0]->Clone("InvMass_US_minus_LS");
  hDiff->Add(hInvMass[1],-1);
  hDiff->SetYTitle("US-LS");

  TLegend *leg = new TLegend(0.5,0.63,0.7,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Like sign (++)+(--)","L");

  // Upsilon mass range
  TH1F *hInvMass_upsilon[2];
  for(Int_t i=0; i<2; i++)
    {
      hInvMass_upsilon[i] = (TH1F*)hInvMass[i]->Clone(Form("%s_upsilon",hInvMass[i]->GetName()));
      if(type==0) hInvMass_upsilon[i]->Rebin(100);
      if(type==1) hInvMass_upsilon[i]->Rebin(10);
      hInvMass_upsilon[i]->GetXaxis()->SetRangeUser(8,14);
    }
  c = draw1D(hInvMass_upsilon[0],Form("Invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",hlt_name[hlt_index]),kTRUE,kTRUE);
  hInvMass_upsilon[1]->Draw("HIST sames");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_upsion_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_upsion_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  TH1F *hdiff_upsilon = (TH1F*)hInvMass_upsilon[0]->Clone(Form("%s_upsilon",hDiff->GetName()));
  hdiff_upsilon->Add(hInvMass_upsilon[1],-1);
  c = draw1D(hdiff_upsilon,Form("Invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});UL-LS",hlt_name[hlt_index]),kFALSE,kTRUE);
  gPad->SetGridy();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_upsion_UL-LS_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_upsion_UL-LS_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  
}

//================================================
void dimuon(Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(9));

  const Int_t nSpecMBins = 58;//PRL mass bin
  Double_t specM[nSpecMBins+1] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 
				  0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.17, 0.2, 0.31, 0.4, 0.51, 
				  0.63, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.89, 0.96, 0.99, 1, 1.01, 1.02, 
				  1.03, 1.04, 1.13, 1.25, 1.45, 1.65, 1.87, 2.07, 2.25, 2.47, 2.66, 2.85, 
				  2.99, 3.02, 3.03, 3.05, 3.07, 3.08, 3.09, 3.10, 3.11, 3.12, 3.13, 
				  3.22, 3.4, 3.85, 4.4, 5.5};

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  double pt1_cut[4] = {1.0,1.2,1.2,1.5};
  double pt2_cut[4] = {1.0,1.0,1.2,1.0};
  TH1F *hInvMass[3][4];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      for(int i=0; i<4; i++)
	{
	  hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut[i]+0.01,100);
	  hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut[i]+0.01,100);
	  TH1F *htmp = (TH1F*)hnInvMass[j]->Projection(0);
	  htmp->SetName(Form("%s_InvMass_%d",hName[j],i));
	  hInvMass[j][i] = (TH1F*)htmp->Rebin(nSpecMBins,Form("%s_%s_InvMass_%d",hName[j],trigName[kTrigType],i),specM);
	  scaleHisto(hInvMass[j][i], 1, 1,kTRUE,kFALSE, kTRUE);
	}
    }

  double max = hInvMass[0][0]->GetMaximum();
  for(int i=0; i<4; i++)
    {
      hInvMass[1][i]->Add(hInvMass[2][i]);
      hInvMass[0][i]->SetMarkerStyle(20);
      hInvMass[0][i]->SetMarkerColor(2);
      hInvMass[0][i]->SetLineColor(2);
      hInvMass[0][i]->SetYTitle("Counts");

      TH1F *hDiff = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US_minus_LS_pt1%1.1f_pt2%1.1f",pt1_cut[i],pt2_cut[i]));
      hDiff->Add(hInvMass[1][i],-1);
      hDiff->SetYTitle("US-LS");

      TLegend *leg = new TLegend(0.5,0.63,0.7,0.83);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut[i],pt2_cut[i]));
      leg->AddEntry(hInvMass[0][i],"Unlike sign","PLE");
      leg->AddEntry(hInvMass[1][i],"Like sign (++)+(--)","L");
      //leg->AddEntry(hDiff,"US-LS","L");


      hInvMass[0][i]->SetMaximum(1.1*max);
      hInvMass[0][i]->GetYaxis()->SetNdivisions(505);
      hInvMass[0][i]->SetLineColor(4);
      hInvMass[0][i]->SetMarkerColor(4);
      //hInvMass[0][i]->GetXaxis()->SetRangeUser(0,1.5);
      c = draw1D(hInvMass[0][i],Form("Invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2});counts/bin"),kFALSE,kFALSE);
      hInvMass[1][i]->Draw("HIST sames");
      hDiff->Draw("sames P");
      leg->Draw();

      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sDiMuonContinuum_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut[i]*10,pt2_cut[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sDiMuonContinuum_full_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut[i]*10,pt2_cut[i]*10));
	}
    }

  // BUR Run16/17
  if(0)
    {
      TFile *fout = TFile::Open("Rootfiles/Run13.pp500.DiMuonContinuum.root","recreate");
      hInvMass[0][2]->Write(Form("InvMass_US_pt1%1.0f_pt2%1.0f",pt1_cut[2]*10,pt2_cut[2]*10));
      hInvMass[1][2]->Write(Form("InvMass_LS_pt1%1.0f_pt2%1.0f",pt1_cut[2]*10,pt2_cut[2]*10));
      fout->Close();
    }
}
