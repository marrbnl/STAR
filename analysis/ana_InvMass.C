TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "PrimTrk.ClosePrimVtx.DCA1cm";

//================================================
void ana_InvMass()
{
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  InvMass();
  //upsilon();
  //daughters();
}


//================================================
void daughters(Int_t save = 0)
{
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH2F *hDaugPtCorr[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hDaugPtCorr[j] = (TH2F*)hnInvMass[j]->Projection(4,5);
      hDaugPtCorr[j]->SetName(Form("%s_DaugPtCorr",hnInvMass[j]->GetName()));
    }
  hDaugPtCorr[1]->Add(hDaugPtCorr[2]);
  const char *hTitle[2] = {"Unlike-sign","Like-sign"};
  const char *pName[2] = {"UL","LS"};
  for(Int_t j=0; j<2; j++)
    {
      c = draw2D(hDaugPtCorr[j],Form("Au+Au %s: p_{T} correlation of di-muon pairs%s",trigName[kTrigType],hlt_name[hlt_index]));
      TPaveText *t1 = GetPaveText(0.6,0.8,0.7,0.8);
      t1->AddText(hTitle[j]);
      t1->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.DaugPtCorr_%s.png",run_config,pName[j]));
    }

  // pt cut on daughter
  TH2F *hInvMassVsPt[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(5)->SetRangeUser(1.2+0.01,100);
      hInvMassVsPt[j] = (TH2F*)hnInvMass[j]->Projection(0,4);
      hInvMassVsPt[j]->SetName(Form("%s_%s_InvMassVsPt",hName[j],trigName[kTrigType]));
      hInvMassVsPt[j]->Sumw2();
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }


  Double_t pt_cuts[5] = {1.2,1.5,2.0,2.5,3.0};
  TH1F *hInvMass[3][5];
  for(Int_t j=0; j<3; j++)
    {
      for(Int_t i=0; i<5; i++)
	{
	  hInvMassVsPt[j]->GetXaxis()->SetRangeUser(pt_cuts[i]+0.01,100);
	  hInvMass[j][i] = (TH1F*)hInvMassVsPt[j]->ProjectionY(Form("%s_%s_InvMass_pt1%1.1f_pt2%1.1f",hName[j],trigName[kTrigType],pt_cuts[i],1.2));
	}
    }

  TH1F *hDiff[5];
  for(Int_t i=0; i<5; i++)
    {
      hDiff[i] = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US-LS_pt1%1.1f_pt2%1.1f",pt_cuts[i],1.2));
      hInvMass[1][i]->Add(hInvMass[2][i]);
      hDiff[i]->Add(hInvMass[1][i],-1);
      hDiff[i]->Rebin(4);
      // list->Add(hDiff[i]);
      // legName[i] = Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt1_cut[i],pt2_cut[i]);

      hInvMass[0][i]->SetMarkerStyle(20);
      hInvMass[0][i]->SetMarkerColor(2);
      hInvMass[0][i]->SetLineColor(2);
      hInvMass[0][i]->GetXaxis()->SetRangeUser(2.5,3.5);
      if(i>0) hInvMass[0][i]->SetMaximum(1.3*hInvMass[0][i]->GetMaximum());
      c = draw1D(hInvMass[0][i],Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE);
      hInvMass[1][i]->Draw("HIST sames");
      TLegend *leg = new TLegend(0.55,0.65,0.7,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt_cuts[i],1.2));
      leg->AddEntry(hInvMass[0][i],"Unlike sign","PLE");
      leg->AddEntry(hInvMass[1][i],"Like sign (++)+(--)","L");
      leg->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.jpsi_InvMass_pt1_%1.1f_pt2_%1.1f.png",run_config,pt_cuts[i],1.2));

      hDiff[i]->GetXaxis()->SetRangeUser(2.5,3.5);
      hDiff[i]->SetMarkerStyle(21);
      hDiff[i]->SetMarkerColor(2);
      hDiff[i]->SetLineColor(2);
      hDiff[i]->GetYaxis()->SetRangeUser(-400,500);
      c = draw1D(hDiff[i],Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});US-LS",trigName[kTrigType],hlt_name[hlt_index]),kFALSE);
      gPad->SetGridy();
      TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85,0.04);
      t1->AddText(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt_cuts[i],1.2));
      t1->Draw();
      TLine *line = GetLine(3.097,hDiff[i]->GetMinimum(),3.097,hDiff[i]->GetMaximum()*0.8,1);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.jpsi_US-LS_pt1_%1.1f_pt2_%1.1f.png",run_config,pt_cuts[i],1.2));
    }
  // c = drawHistos(list,"InvMass_ComparePt1Pt2Cuts",Form("Au+Au %s: invariant mass of di-muon pairs%s (US-LS);M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kTRUE,3.0,3.2,kTRUE,-200,300,kFALSE,kTRUE,legName,kTRUE,"",0.15,0.25,0.7,0.88);
  // gPad->SetGridy();
  // if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_ComparePt1Pt2Cut.png",run_config));
}


//================================================
void InvMass(Int_t save = 1)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));
  printf("single-muon events: %4.2e\n",hStat->GetBinContent(8));
  printf("e-muon      events: %4.2e\n",hStat->GetBinContent(9));
  printf("HLT         events: %4.2f%%\n",100*hStat->GetBinContent(10)/hStat->GetBinContent(7));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3][2];
  Double_t pt1_cut = 2.0, pt2_cut = 1.2;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));

      hInvMass[j][0] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j][0]->SetName(Form("%s_%s_InvMass_All",hName[j],trigName[kTrigType]));
      hInvMass[j][0]->Sumw2();

      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      hInvMass[j][1] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j][1]->SetName(Form("%s_%s_InvMass_WithCut",hName[j],trigName[kTrigType]));
      hInvMass[j][1]->Sumw2();

      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }

  for(Int_t i=0; i<2; i++)
    {
      hInvMass[1][i]->Add(hInvMass[2][i]);
      hInvMass[0][i]->SetMarkerStyle(20);
      hInvMass[0][i]->SetMarkerColor(2);
      hInvMass[0][i]->SetLineColor(2);
      hInvMass[0][i]->SetYTitle("Counts");
    }

  TH1F *hDiff = (TH1F*)hInvMass[0][0]->Clone("InvMass_US_minus_LS_all");
  hDiff->Add(hInvMass[1][0],-1);
  hDiff->SetYTitle("US-LS");

  TLegend *leg = new TLegend(0.5,0.63,0.7,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",1.0,1.0));
  leg->AddEntry(hInvMass[0][0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass[1][0],"Like sign (++)+(--)","L");

  // full mass range
  c = draw1D(hInvMass[0][0],Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kTRUE);
  hInvMass[1][0]->Draw("HIST sames");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_full_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  TH1F *hdiff_full = (TH1F*)hDiff->Clone(Form("%s_full",hDiff->GetName()));
  hdiff_full->Rebin(5);
  c = draw1D(hdiff_full,Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
  gPad->SetGridy();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_full_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  // J/psi mass range
  TH1F *hInvMass_jpsi[2];
  for(Int_t j=0; j<2; j++)
    {
      hInvMass_jpsi[j] = (TH1F*)hInvMass[j][1]->Clone(Form("%s_jpsi",hInvMass[j][1]->GetName()));
      hInvMass_jpsi[j]->GetXaxis()->SetRangeUser(2.5,3.5);
    }
  c = draw1D(hInvMass_jpsi[0],Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]));
  hInvMass_jpsi[1]->Draw("HIST sames");
  leg = new TLegend(0.2,0.63,0.4,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass_jpsi[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass_jpsi[1],"Like sign (++)+(--)","L");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_jpsi_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  TH1F *hdiff_jpsi = (TH1F*)hInvMass_jpsi[0]->Clone(Form("%s_jpsi",hDiff->GetName()));
  hdiff_jpsi->Add(hInvMass_jpsi[1],-1);
  hdiff_jpsi->GetXaxis()->SetRangeUser(3.05,3.15);
  c = draw1D(hdiff_jpsi,Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});US-LS",trigName[kTrigType],hlt_name[hlt_index]),kFALSE);
  gPad->SetGridy();
  TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85,0.04);
  t1->AddText(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_jpsi_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));
  hdiff_jpsi->SetName(Form("%s_rebin4",hdiff_jpsi->GetName()));
  hdiff_jpsi->Rebin(4);
  hdiff_jpsi->GetXaxis()->SetRangeUser(2.5,3.5);
  c = draw1D(hdiff_jpsi,Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2});US-LS",trigName[kTrigType],hlt_name[hlt_index]),kFALSE);
  gPad->SetGridy();
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_jpsi_pt1_%1.1f_pt2_%1.1f_rebin4.png",run_config,pt1_cut,pt2_cut));
  
}

//================================================
void upsilon(Int_t save = 1)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3];
  Double_t pt1_cut = 4.5, pt2_cut = 1.2;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("%s_%s_InvMass",hName[j],trigName[kTrigType]));
      hInvMass[j]->Sumw2();
    }
  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->SetYTitle("Counts");

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
      hInvMass_upsilon[i]->Rebin(10);
      hInvMass_upsilon[i]->GetXaxis()->SetRangeUser(8.5,11);
    }
  c = draw1D(hInvMass_upsilon[0],Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kTRUE);
  hInvMass_upsilon[1]->Draw("HIST sames");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_upsion_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  TH1F *hdiff_upsilon = (TH1F*)hInvMass_upsilon[0]->Clone(Form("%s_upsilon",hDiff->GetName()));
  hdiff_upsilon->Add(hInvMass_upsilon[1],-1);
  c = draw1D(hdiff_upsilon,Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kTRUE);
  gPad->SetGridy();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_upsilon_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));
  
}

