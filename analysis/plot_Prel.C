const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;

//================================================
void plot_Prel()
{
  Run13_pp500_raw();
}


//================================================
void Run13_pp500_raw(const Int_t save = 0)
{
  TFile *f = TFile::Open(Form("./output/Run13.pp500.jpsi.EventQA.root"),"read");

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));
  printf("single-muon events: %4.2e\n",hStat->GetBinContent(8));
  printf("e-muon      events: %4.2e\n",hStat->GetBinContent(9));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3][2];
  Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));

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
      hInvMass[0][i]->SetYTitle("Counts");
    }

  // J/psi mass range
  TH1F *hInvMass_jpsi[2];
  for(Int_t j=0; j<2; j++)
    {
      hInvMass_jpsi[j] = (TH1F*)hInvMass[j][1]->Clone(Form("%s_jpsi",hInvMass[j][1]->GetName()));
      hInvMass_jpsi[j]->Rebin(2);
      hInvMass_jpsi[j]->SetMarkerColor(4-j*3);
      hInvMass_jpsi[j]->SetLineColor(4-j*3);
      hInvMass_jpsi[j]->GetXaxis()->SetRangeUser(2,4);
    }
  hInvMass_jpsi[0]->GetYaxis()->SetRangeUser(-20,400);
  //c = draw1D(hInvMass_jpsi[0],Form("Invariant mass of di-muon pairs reconstructed using MTD;M_{#mu#mu} (GeV/c^{2})"),kFALSE,kFALSE);
  c = draw1D(hInvMass_jpsi[0],Form(";M_{#mu#mu} (GeV/c^{2})"),kFALSE,kFALSE);
  gPad->SetGridy();
  hInvMass_jpsi[1]->Draw("HIST sames");

  TH1F *hdiff_jpsi = (TH1F*)hInvMass_jpsi[0]->Clone(Form("InvMass_US_minus_LS_jpsi"));
  hdiff_jpsi->Add(hInvMass_jpsi[1],-1);
  hdiff_jpsi->SetMarkerStyle(20);
  hdiff_jpsi->SetMarkerColor(2);
  hdiff_jpsi->SetLineColor(2);
  hdiff_jpsi->Draw("sames P");
  
  leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  //leg->SetHeader("MTD");
  leg->AddEntry(hInvMass_jpsi[0],"Unlike sign","LE");
  leg->AddEntry(hInvMass_jpsi[1],"Like sign (++)+(--)","LE");
  leg->AddEntry(hdiff_jpsi,"US - LS","PLE");
  leg->Draw();

  Int_t low_bin = hInvMass_jpsi[0]->GetXaxis()->FindFixBin(low_mass+0.001);
  Int_t high_bin = hInvMass_jpsi[0]->GetXaxis()->FindFixBin(high_mass-0.001);
  Double_t nBackground = hInvMass_jpsi[1]->Integral(low_bin,high_bin);
  Double_t nSignal = hInvMass_jpsi[0]->Integral(low_bin,high_bin) - nBackground;

  t1 = GetPaveText(0.13,0.4,0.75,0.88,0.04,62);
  t1->AddText("Run13 pp #sqrt{s} = 500 GeV");
  t1->AddText("MTD di-muon trigger, 100M events");
  t1->SetTextAlign(11);
  t1->SetFillColor(0);
  t1->Draw();

  TPaveText *signif = GetPaveText(0.13,0.4,0.5,0.7,0.04,62);
  signif->AddText(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  signif->AddText(Form("%1.1f < M_{#mu#mu} < %1.1f GeV/c^{2}",low_mass,high_mass));
  signif->AddText(Form("S/B = %1.0f/%1.0f = %1.2f:1",nSignal,nBackground,nSignal/nBackground));
  signif->SetTextAlign(11);
  signif->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run13_pp500/Preliminary/Run13_pp500_MTD_Jpsi_Raw.pdf"));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run13_pp500/Preliminary/Run13_pp500_MTD_Jpsi_Raw.png"));
    }
}


