void inv_mass(const TString invFile,
	      const TString outdir)
{
  gStyle->SetOptStat(0);

  // invariant mass plot
  TFile *finv = TFile::Open(invFile.Data(),"read");
  TH1F *hStat = (TH1F*)finv->Get("mhEventStat_dimuon");
  const char *histoName[3] = {"mhInvMvsPtUL_dimuon","mhInvMvsPtLSpos_dimuon","mhInvMvsPtLSneg_dimuon"};
  TH1F *hInvMass[3];
  for(Int_t j=0; j<3; j++)
    {
      TH2F *h2 = (TH2F*)finv->Get(histoName[j]);
      h2->Sumw2();
      hInvMass[j] = (TH1F*)h2->ProjectionX(Form("hInvMass_%d",j));
    }
  hInvMass[1]->Add(hInvMass[2]);

  TCanvas *c = new TCanvas(Form("hInvMass"),Form("hInvMass"),800,600);
  gPad->SetGrid(0,0);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->GetXaxis()->SetRangeUser(0,4.5);
  hInvMass[0]->SetMaximum(1.5*hInvMass[0]->GetMaximum());
  hInvMass[0]->SetTitle(";M_{#mu#mu} (GeV/c^{2});counts");
  hInvMass[0]->Draw("P");
  hInvMass[1]->Draw("sames HIST");

  TPaveText* title = new TPaveText(0.3530151,0.8968531,0.6532663,0.9965035,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->AddText(0.,0.,"Invariant mass of di-muon pairs");
  title->SetTextSize(0.045);
  title->SetTextFont(62);
  title->Draw();

  TLegend *leg = new TLegend(0.15,0.65,0.3,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",1.5,1.0));
  leg->AddEntry(hInvMass[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Like sign (++)+(--)","L");
  leg->Draw();

  Int_t low_bin = hInvMass[0]->GetXaxis()->FindFixBin(3.0+0.001);
  Int_t high_bin = hInvMass[0]->GetXaxis()->FindFixBin(3.2-0.001);
  Double_t nBackground = hInvMass[1]->Integral(low_bin,high_bin);
  Double_t nSignal = hInvMass[0]->Integral(low_bin,high_bin) - nBackground;
  TPaveText *signif = new TPaveText(0.5,0.65,0.8,0.85,"brNDC");
  signif->SetFillStyle(0);
  signif->SetBorderSize(0);
  signif->SetTextSize(0.045);
  signif->SetTextFont(62);
  signif->SetTextAlign(11);
  signif->AddText("Run14 Au+Au 200 GeV/c");
  signif->AddText(Form("%3.2fM events (filtered)",hStat->GetBinContent(2)/1e6));
  signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",3.0,3.2));
  if(nBackground>0) signif->AddText(Form("S/B = %1.0f/%1.0f = %1.2e:1",nSignal,nBackground,nSignal/nBackground));
  else              signif->AddText(Form("S/B = %1.0f/%1.0f",nSignal,nBackground));
  signif->Draw();
  c->SaveAs(Form("%s/DimuonInvMass.pdf",outdir.Data()));
}
