void ana_Run11_eff()
{

  const int rebin = 1;

  TFile *f = TFile::Open("Rootfiles/Run11_eff.root","read");

  // Run with weigth
  TH2F *hMcPtVsEta = (TH2F*)f->Get("mcJpsiPtY");
  TH2F *hRcPtVsEta = (TH2F*)f->Get("hHt2JpsiPE");

  draw2D(hMcPtVsEta);
  draw2D(hRcPtVsEta);

  hMcPtVsEta->GetXaxis()->SetRangeUser(-1+1e-6,1-1e-6);
  TH1F *hMcPt = (TH1F*)hMcPtVsEta->ProjectionY("hMcPt");
  hMcPt->Rebin(rebin);
  hMcPt->SetMarkerStyle(20);
  draw1D(hMcPt,"",kTRUE);

  hRcPtVsEta->GetXaxis()->SetRangeUser(-1+1e-6,1-1e-6);
  TH1F *hRcPt = (TH1F*)hRcPtVsEta->ProjectionY("hRcPt");
  hRcPt->Rebin(rebin);
  hRcPt->SetMarkerStyle(21);
  hRcPt->SetMarkerColor(2);
  hRcPt->SetLineColor(2);
  hRcPt->Draw("sames P");

  TH1F *hRatio = (TH1F*)hRcPt->Clone("hRatio_fromRunning");
  hRatio->Rebin(100);
  hMcPt->Rebin(100);
  hRatio->Divide(hMcPt);
  cEff = draw1D(hRatio,"");

  // Run without weight
  TH2F *hMcPtVsEtaNoWeight = (TH2F*)f->Get("mcJpsiPtY_Or");
  draw2D(hMcPtVsEtaNoWeight);  

  TH2F *hMcPtVsRcNoWeight = (TH2F*)f->Get("hJpsiRcvsMC_Cut1");
  hMcPtVsRcNoWeight->RebinX(rebin);
  hMcPtVsRcNoWeight->RebinY(rebin);
  draw2D(hMcPtVsRcNoWeight);

  hMcPtVsEtaNoWeight->GetXaxis()->SetRangeUser(-1+1e-6,1-1e-6);
  TH1F *hMcPtNoWeight = (TH1F*)hMcPtVsEtaNoWeight->ProjectionY("hMcPtNoWeight");
  hMcPtNoWeight->Rebin(rebin);
  hMcPtNoWeight->SetMarkerStyle(20);
  hMcPtNoWeight->SetMinimum(1);
  draw1D(hMcPtNoWeight,"",kTRUE);

  TH1F *hRcPtNoWeight = (TH1F*)hMcPtVsRcNoWeight->ProjectionX("hRcPtNoWeight");
  hRcPtNoWeight->SetMarkerStyle(21);
  hRcPtNoWeight->SetMarkerColor(2);
  hRcPtNoWeight->SetLineColor(2);
  hRcPtNoWeight->Draw("sames P");

  TH1F *hRatioNoWeight = (TH1F*)hRcPtNoWeight->Clone("hRatioNoWeight");
  hRatioNoWeight->Divide(hMcPtNoWeight);
  cEff->cd();
  hRatioNoWeight->SetMarkerColor(4);
  hRatioNoWeight->Draw("samesP");

  // weight with input histogram
  TH1F *hMcPtWeight = (TH1F*)hMcPtNoWeight->Clone("hMcPtWeight");
  TH2F *hMcPtVsRcWeight = (TH2F*)hMcPtVsRcNoWeight->Clone("hMcPtVsRcWeight");
  for(int ibin=1; ibin<=hMcPtVsRcNoWeight->GetNbinsX(); ibin++)
    {
      double scale = hMcPt->GetBinContent(ibin);
      hMcPtWeight->SetBinContent(ibin,hMcPtWeight->GetBinContent(ibin)*scale);
      hMcPtWeight->SetBinError(ibin,hMcPtWeight->GetBinError(ibin)*scale);
      for(int jbin=1; jbin<=hMcPtVsRcNoWeight->GetNbinsY(); jbin++)
	{
	  hMcPtVsRcWeight->SetBinContent(ibin,jbin,hMcPtVsRcWeight->GetBinContent(ibin,jbin)*scale);
	  hMcPtVsRcWeight->SetBinError(ibin,jbin,hMcPtVsRcWeight->GetBinError(ibin,jbin)*scale);
	}
    }
  TH1F *hRcPtWeight = (TH1F*)hMcPtVsRcWeight->ProjectionY("hRcPtWeight");
  hRcPtWeight->SetMarkerStyle(21);
  hRcPtWeight->SetMarkerColor(2);
  hRcPtWeight->SetLineColor(2);

  draw2D(hMcPtVsRcWeight);
  draw1D(hMcPtWeight,"",kTRUE);
  hRcPtWeight->Draw("sames P");
  
  TH1F *hRatioWeight = (TH1F*)hRcPtWeight->Clone("hRatioWeight");
  hRatioWeight->Divide(hMcPtWeight);
  cEff->cd();
  hRatioWeight->SetMarkerColor(6);
  hRatioWeight->Draw("samesP");

  TH1F *hCheck = (TH1F*)hRatioWeight->Clone("check");
  hCheck->Divide(hRatio);
  draw1D(hCheck);

  // weight with fitted function
  TCanvas *c = new TCanvas("Fit","Fit",800,600);
  SetPadMargin(gPad,0.15,0.15);
  gPad->SetLogy();
  TH1F *h = new TH1F("histogram",";;;",7,0,30);
  h->GetYaxis()->SetRangeUser(1e-7,100);
  h->Draw();

  TFile *fdata = TFile::Open("Rootfiles/Spectrum_in_bin.root","read");
  TGraphErrors	*gr = (TGraphErrors*)fdata->Get("gall");
  gr->SetMarkerColor(1);
  gr->SetLineColor(1);
  gr->GetXaxis()->SetRangeUser(0,30);
  gr->Draw("sames PE");
  TF1 *func = new TF1("func",InvPt,0,30,4);
  func->SetParameters(0.4,-0.4796,4.229,-7.54);
  gr->Fit(func,"RL");

  TH1F *hMcPtFunc = (TH1F*)hMcPtNoWeight->Clone("hMcPtFunc");
  TH2F *hMcPtVsRcFunc = (TH2F*)hMcPtVsRcNoWeight->Clone("hMcPtVsRcFunc");
  for(int ibin=1; ibin<=hMcPtVsRcFunc->GetNbinsX(); ibin++)
    {
      double scale = func->Eval(hMcPtFunc->GetBinCenter(ibin));
      hMcPtFunc->SetBinContent(ibin,hMcPtFunc->GetBinContent(ibin)*scale);
      hMcPtFunc->SetBinError(ibin,hMcPtFunc->GetBinError(ibin)*scale);
      for(int jbin=1; jbin<=hMcPtVsRcNoWeight->GetNbinsY(); jbin++)
	{
	  hMcPtVsRcFunc->SetBinContent(ibin,jbin,hMcPtVsRcFunc->GetBinContent(ibin,jbin)*scale);
	  hMcPtVsRcFunc->SetBinError(ibin,jbin,hMcPtVsRcFunc->GetBinError(ibin,jbin)*scale);
	}
    }
  TH1F *hRcPtFunc = (TH1F*)hMcPtVsRcFunc->ProjectionY("hRcPtFunc");
  hRcPtFunc->SetMarkerStyle(21);
  hRcPtFunc->SetMarkerColor(2);
  hRcPtFunc->SetLineColor(2);
  hMcPtVsRcFunc->GetZaxis()->SetRangeUser(1e-4,1e2);
  draw2D(hMcPtVsRcFunc);
  hMcPtFunc->GetYaxis()->SetRangeUser(1e-4,5e4);
  draw1D(hMcPtFunc,"",kTRUE);
  hRcPtFunc->Draw("sames P");

  TH1F *hRatioFunc = (TH1F*)hRcPtFunc->Clone("hRatioFunc");
  hRatioFunc->Rebin(100);
  hMcPtFunc->Rebin(100);
  hRatioFunc->Divide(hMcPtFunc);
  cEff->cd();
  hRatioFunc->SetMarkerColor(5);
  hRatioFunc->Draw("samesP");
  
  TH1F *hCheck2 = (TH1F*)hRatioFunc->Clone("check2");
  hCheck2->Divide(hRatio);
  draw1D(hCheck2);

}

Double_t InvPt(Double_t *x,Double_t *par)
{
  return par[0]*TMath::Power(TMath::Exp(par[1]*x[0])+x[0]/par[2],par[3]);
}
