TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "GlobalTrk.NoDCA.HLT";

//================================================
void ana_InvMass()
{
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  InvMass();
  //daughters();
}


//================================================
void daughters(Int_t save = 1)
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
}


//================================================
void InvMass(Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));
  printf("single-muon events: %4.2e\n",hStat->GetBinContent(8));
  printf("e-muon      events: %4.2e\n",hStat->GetBinContent(9));
  printf("HLT         events: %4.2f%%\n",100*hStat->GetBinContent(10)/hStat->GetBinContent(7));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3];
  Double_t pt1_cut = 1.5, pt2_cut = 1;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut,100);
      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("%s_%s_InvMass",hName[j],trigName[kTrigType]));
    }

  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->GetXaxis()->SetRangeUser(2.5,3.5);
  hInvMass[0]->SetYTitle("Counts");

  TH1F *hDiff = (TH1F*)hInvMass[0]->Clone("InvMass_US_minus_LS");
  hDiff->Sumw2();
  hDiff->Add(hInvMass[1],-1);
  hDiff->SetYTitle("US-LS");
  c = draw1D(hDiff,Form("Au+Au %s: invariant mass of di-muon pairs%s",trigName[kTrigType],hlt_name[hlt_index]));
  TLine *line = GetLine(2.5,1,3.5,1,1);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_pt1%1.1f_pt2%1.1f.png",run_config,pt1_cut,pt2_cut));


  c = draw1D(hInvMass[0],Form("Au+Au %s: invariant mass of di-muon pairs%s",trigName[kTrigType],hlt_name[hlt_index]));
  hInvMass[1]->Draw("HIST sames");
  TLegend *leg = new TLegend(0.6,0.75,0.75,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hInvMass[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Like sign (++)+(--)","L");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_pt1%1.1f_pt2%1.1f.png",run_config,pt1_cut,pt2_cut));
}
