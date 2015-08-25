TFile *f;
const char *run_config = "";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void ana_dimuon()
{
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.15); 
  TString fileName;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) fileName = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      fileName = Form("Run13.pp500.jpsi.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) fileName = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      fileName = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");

  run_cfg_name = "";

  dimuon();
}

//================================================
void dimuon(Int_t save = 1)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(9));

  double pt1_cut = 1.5, pt2_cut = 1.0;

  // Acceptance from mixed event
  TFile *fmix = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.MixEvent.root"),"read");
  TH2F *hMixUL = (TH2F*)fmix->Get(Form("US_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s_Rebin",pt1_cut,pt2_cut,trigName[kTrigType]));
  TH2F *hMixLS = (TH2F*)fmix->Get(Form("LS_Geom_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt1_cut,pt2_cut,trigName[kTrigType]));
  TH2F *hMixRatio = DivideTH2WithTH2(hMixUL, hMixLS);
  
  c = draw2D(hMixUL,Form("Mixed event: unlike-sign pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)",pt1_cut,pt2_cut));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/%sInvMass_vs_pt_UL_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/%sInvMass_vs_pt_UL_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  c = draw2D(hMixLS,Form("Mixed event: like-sign pairs with geometic mean (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)",pt1_cut,pt2_cut));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/%sInvMass_vs_pt_LS_Geom_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/%sInvMass_vs_pt_LS_Geom_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  c = draw2D(hMixRatio,Form("Mixed event: US/LS (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)",pt1_cut,pt2_cut),0.04,kFALSE,"colz");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/%sAcceptance_Geom_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/%ssAcceptance_Geom_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10));
    }

  // acceptance correction
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH2F *hInvMassVsPt[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      hInvMassVsPt[j] = (TH2F*)hnInvMass[j]->Projection(1,0);
      hInvMassVsPt[j]->SetName(Form("%s_%s_InvMassVsPt",hName[j],trigName[kTrigType]));
    }

  TH2F *htmp = (TH2F*)hInvMassVsPt[1]->Clone(Form("SameEvent_LS_InvMass_vs_pt_tmp"));
  htmp->Add(hInvMassVsPt[2]);

  TH2F *hSameBkg = new TH2F(Form("SameEvent_LS_InvMass_vs_pt"), htmp->GetTitle(), nSpecMBins, specM, nSpecPtBins, specPt);
  CopyTH2(htmp, hSameBkg);
  hSameBkg->Sumw2();
  hMixRatio->Sumw2();
  TH2F *hSameBkgCorr = MultiplyTH2WithTH2(hSameBkg, hMixRatio);
  hSameBkgCorr->SetName(Form("SameEvent_LS_InvMass_vs_pt_AccCorr"));
  
  TList *list = new TList;
  TString legName[2] = {"Raw","Acceptance corrected"};
  TH1F *hLSBkg = (TH1F*)hSameBkg->ProjectionX("SameEvent_LS_bkg");
  TH1F *hLSBkgCorr = (TH1F*)hSameBkgCorr->ProjectionX("SameEvent_LS_bkg_AccCorr");


  TH1F *hLSBkgRebin = (TH1F*)hLSBkg->Rebin(nSpecMBins2,Form("%s_rebin",hLSBkg->GetName()),specM2);
  TH1F *hLSBkgCorrRebin = (TH1F*)hLSBkgCorr->Rebin(nSpecMBins2,Form("%s_rebin",hLSBkgCorr->GetName()),specM2);
  scaleHisto(hLSBkgRebin, 1, 1,kTRUE,kFALSE, kTRUE);
  scaleHisto(hLSBkgCorrRebin, 1, 1,kTRUE,kFALSE, kTRUE);
  list->Add(hLSBkgRebin);
  list->Add(hLSBkgCorrRebin);
  c = sysCompare(list,Form("LS_pt1_%1.1f_pt2_%1.1f",pt1_cut,pt2_cut),Form("SameEvent: like-sign di-muon pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt1_cut,pt2_cut),"Acceptance factor",kTRUE,0,5,kFALSE,0.1,10,kTRUE,0.7,1.3,kFALSE,kTRUE,legName,kTRUE,"Like-sign",0.2,0.4,0.2,0.4,kTRUE);
  TLine *line = GetLine(0,1,5.5,1,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/SameEvent_LS_AccCorr_pt1_%1.0f_pt2_%1.0f.pdf",run_type,pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/SameEvent_LS_AccCorr_pt1_%1.0f_pt2_%1.0f.png",run_type,pt1_cut*10,pt2_cut*10));
    }

  // di-muon continuum
  scaleHisto(hLSBkg, 1, 1,kTRUE,kFALSE, kTRUE);
  scaleHisto(hLSBkgCorr, 1, 1,kTRUE,kFALSE, kTRUE);
  hInvMassVsPt[0]->Sumw2();
  TH1F *h1 = (TH1F*)hInvMassVsPt[0]->ProjectionX("SameEvent_UL");
  TH1F *hUL = (TH1F*)h1->Rebin(nSpecMBins,"SameEvent_UL_Rebin",specM);
  scaleHisto(hUL, 1, 1,kTRUE,kFALSE, kTRUE);
  
  hUL->SetMarkerStyle(20);
  TH1F *hDiff = (TH1F*)hUL->Clone(Form("InvMass_US_minus_LS_pt1%1.1f_pt2%1.1f",pt1_cut,pt2_cut));
  hDiff->Add(hLSBkgCorr,-1);
  hUL->SetMaximum(1.3*hUL->GetMaximum());
  hUL->SetMinimum(-8000);
  hUL->GetYaxis()->SetNdivisions(505);
  hUL->SetLineColor(4);
  hUL->SetMarkerColor(4);
  hUL->GetXaxis()->SetRangeUser(0,4);
  c = draw1D(hUL,Form("Invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2});dN/dM"),kFALSE,kTRUE);
  hLSBkgCorr->SetLineColor(1);
  hLSBkgCorr->Draw("HIST sames");
  hDiff->SetMarkerColor(2);
  hDiff->SetLineColor(2);
  hDiff->Draw("sames P");
  hLSBkg->SetLineColor(6);
  hLSBkg->Draw("HIST sames");

  TLegend *leg = new TLegend(0.12,0.65,0.35,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hUL,"Unlike sign","PEL");
  leg->AddEntry(hLSBkg,"Like sign","L");
  leg->AddEntry(hLSBkgCorr,"Like sign (accep. corr.)","L");
  leg->AddEntry(hDiff,"US-LS","PEL");
  leg->Draw();

  TLine *line = GetLine(0,0,4,0,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/DiMuonContinuum_pt1_%1.0f_pt2_%1.0f.pdf",run_type,pt1_cut*10,pt2_cut*10));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_dimuon/DiMuonContinuum_pt1_%1.0f_pt2_%1.0f.png",run_type,pt1_cut*10,pt2_cut*10));
    }
}
