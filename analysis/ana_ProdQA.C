TFile *f;
const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;
TString run_cfg_name;

//================================================
void ana_ProdQA()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.17);                
  gStyle->SetStatH(0.16);

  f = TFile::Open("./output/AllDay.MuDst.qa.root","read");
  run_cfg_name = "ProdQA.";
  run_type = "Run14_AuAu200";

  fitSignal();
  //fitBackground();
}



//================================================
void fitSignal(Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("mhEventStat_dimuon");
  for(int bin=1; bin<=hStat->GetNbinsX(); bin++)
    {
      if(hStat->GetBinContent(bin)==0) continue;
      printf("%s: %3.2e\n",hStat->GetXaxis()->GetBinLabel(bin),hStat->GetBinContent(bin));
    }

  const Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  const Int_t nPtBins = 3;
  const Double_t ptBins[nPtBins+1] = {0,2,4,10};
  const char *hName[3] = {"mhInvMvsPtUL","mhInvMvsPtLSpos","mhInvMvsPtLSneg"};

  TH2F *hnInvMass[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (TH2F*)f->Get(Form("%s_dimuon",hName[j]));
      hnInvMass[j]->Sumw2();
      draw2D(hnInvMass[j]);
    }
  hnInvMass[1]->Add(hnInvMass[2]);
 
  TH1F *hInvMass[2][nPtBins];
  TH1F *hSignal[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hnInvMass[j]->GetYaxis()->SetRangeUser(ptBins[i]+0.01,ptBins[i+1]-0.01);
	  hInvMass[j][i] = (TH1F*)hnInvMass[j]->ProjectionX(Form("%s_%s_InvMass_jpsi_%1.1f_%1.1f",hName[j],trigName[kTrigType],ptBins[i],ptBins[i+1]));
	  hInvMass[j][i]->Rebin(2);
	}
      hSignal[i] = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US_minus_LS_jpsi_%1.1f_%1.1f",ptBins[i],ptBins[i+1]));
      hSignal[i]->Add(hInvMass[1][i],-1);
      hSignal[i]->SetYTitle("US-LS");
    }

  double nAll[nPtBins];
  double allErr[nPtBins];
  double nBkg[nPtBins];
  double bkgErr[nPtBins];
  double nSignal[nPtBins];
  double sigErr[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      TH1F *hFit = hInvMass[0][i];
      TF1 *func = new TF1(Form("func_%d",i),"gaus(0)+pol3(3)",2.5,3.5);
      if(i==-1) func->FixParameter(1,3.09);
      else
	{
	  func->SetParLimits(1,3.07,3.09);
	  func->SetParLimits(2,0.02,0.06);
	}
      if(i==0)func->SetRange(2.9,3.5);
      hFit->Fit(func,"IR0");
      hInvMass[0][i]->GetXaxis()->SetRangeUser(2.7,3.5);
      hInvMass[0][i]->SetMarkerStyle(21);
      hInvMass[0][i]->SetMarkerColor(2);
      hInvMass[0][i]->SetLineColor(2);
      hInvMass[0][i]->SetMaximum(1.2*hInvMass[0][i]->GetMaximum());
      c = draw1D(hInvMass[0][i],Form("Invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2});Counts"));
      hInvMass[1][i]->Draw("samesHIST");
      func->SetLineColor(4);
      func->Draw("sames");
      TF1 *func1 = new TF1(Form("signal_%d",i),"gaus",func->GetXmin(),func->GetXmax());
      func1->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
      func1->SetLineColor(2);

      TF1 *funcBkg = new TF1(Form("Bkg_%d",i),"pol3",func->GetXmin(),func->GetXmax());
      funcBkg->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5),func->GetParameter(6));
      funcBkg->SetLineColor(4);
      funcBkg->Draw("sames");

      // Singal counts
      int low_bin = hInvMass[0][i]->FindFixBin(low_mass+0.001);
      int high_bin = hInvMass[0][i]->FindFixBin(high_mass-0.001);
      nAll[i] = hInvMass[0][i]->Integral(low_bin,high_bin);
      allErr[i] = TMath::Sqrt(nAll[i]);
      nBkg[i] = (func->Integral(low_mass,high_mass) - func1->Integral(low_mass,high_mass)) * 1./hInvMass[0][i]->GetBinWidth(1);
      bkgErr[i] = TMath::Sqrt(nBkg[i]);
      nSignal[i] = nAll[i] - nBkg[i];
      sigErr[i] = TMath::Sqrt(nAll[i]+nBkg[i]);
      TPaveText *signif = GetPaveText(0.13,0.3,0.12,0.3);
      signif->SetTextAlign(11);
      signif->SetTextFont(62);
      //signif->AddText(Form("p_{T,1} > 1.5 GeV/c"));
      signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
      signif->AddText(Form("Signal = %3.1f #pm %3.1f",nSignal[i],sigErr[i]));
      signif->AddText(Form("S/B = %1.0f/%1.0f = 1:%3.0f",nSignal[i],nBkg[i],nBkg[i]/nSignal[i]));
      signif->Draw();

      TPaveText *t1 = GetPaveText(0.23,0.33,0.8,0.85,0.045);
      t1->SetTextFont(62);
      t1->AddText(Form("%1.0f < p_{T,J/#psi} < %1.0f GeV/c",ptBins[i],ptBins[i+1]));
      t1->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sFitJpsi_UL_%1.0f_%1.0fGeV.pdf",run_type,run_cfg_name.Data(),ptBins[i],ptBins[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sFitJpsi_UL_%1.0f_%1.0fGeV.png",run_type,run_cfg_name.Data(),ptBins[i],ptBins[i+1]));
	}
    }
  
  TH1F *hJpsiPt = new TH1F("hJpsiPt","p_{T} distribution of J/#psi;p_{T} (GeV/c);Counts",nPtBins,ptBins);
  for(int bin=1; bin<=nPtBins; bin++)
    {
      hJpsiPt->SetBinContent(bin,nSignal[bin-1]);
      hJpsiPt->SetBinError(bin,sigErr[bin-1]);
    }
  hJpsiPt->SetMarkerStyle(21);
  hJpsiPt->SetMarkerColor(4);
  c = draw1D(hJpsiPt);
  TPaveText *t1 = GetPaveText(0.5,0.7,0.65,0.85,0.04);
  t1->SetTextFont(62);
  t1->SetTextAlign(11);
  t1->AddText("Run14 Au+Au @ 200 GeV");
  t1->AddText(Form("%3.0fM dimuon events",hStat->GetBinContent(2)/1e6));
  t1->AddText(Form("%3.0f J/#psi candidates",hJpsiPt->Integral(1,3)));
  t1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt.png",run_type,run_cfg_name.Data()));
    }

  TGraphErrors *gr[3];
  for(int i=0; i<3; i++)
    {
      gr[i] = new TGraphErrors(nPtBins);
      for(int ipoint=0; ipoint<nPtBins; ipoint++)
	{
	  gr[i]->SetPoint(ipoint,hJpsiPt->GetBinCenter(ipoint+1)+0.2*i,1);
	  double err = sigErr[ipoint]/nSignal[ipoint];
	  if(i==2) err = err / sqrt(4.6);
	  if(i==1) err = err / sqrt(2.9);
	  gr[i]->SetPointError(ipoint,0,err);
	}
      gr[i]->SetMarkerStyle(21);
      gr[i]->SetMarkerColor(color[i]);
      gr[i]->SetLineColor(color[i]);
    }
  gr[0]->GetYaxis()->SetRangeUser(0.5,1.5);
  c = drawGraph(gr[0],"Projection of relative statistical errors;J/#psi p_{T} (GeV/c)");
  gPad->SetGridy();
  gr[1]->Draw("same PE");
  gr[2]->Draw("same PE");
  TLine *line = GetLine(0.5,1,7.5,1,1,2,2);
  line->Draw();
  leg = new TLegend(0.5,0.63,0.7,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("Production scenario");
  leg->AddEntry(gr[0],"Current statistics","P");
  leg->AddEntry(gr[1],"50/50 MTD/HFT","P");
  leg->AddEntry(gr[2],"90/10 MTD/HFT","P");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt_err.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt_err.png",run_type,run_cfg_name.Data()));
    }

  // upsilon
  TH1F *hUpsilon[2];
  const double ptCut = 4;
  for(Int_t j=0; j<2; j++)
    {
      hnInvMass[j]->GetYaxis()->SetRangeUser(ptCut,10);
      hUpsilon[j] = (TH1F*)hnInvMass[j]->ProjectionX(Form("%s_%s_InvMass_upsilon",hName[j],trigName[kTrigType]));
      hUpsilon[j]->Rebin(40);
      hUpsilon[j]->GetXaxis()->SetRangeUser(8,14);
    }
  hUpsilon[0]->SetMarkerStyle(20);
  hUpsilon[0]->SetMarkerColor(2);
  hUpsilon[0]->SetLineColor(2);
  hUpsilon[0]->SetYTitle("Counts");
  c = draw1D(hUpsilon[0],Form("Invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2})"),kFALSE,kTRUE);
  hUpsilon[1]->Draw("HIST sames");

  TLegend *leg = new TLegend(0.5,0.63,0.7,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,#mu#mu} > %1.0f GeV/c",ptCut));
  leg->AddEntry(hUpsilon[0],"Unlike sign","PLE");
  leg->AddEntry(hUpsilon[1],"Like sign (++)+(--)","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sUpsilon_invMass_Pt%1.0f.pdf",run_type,run_cfg_name.Data(),ptCut));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sUpsilon_invMass_Pt%1.0f.png",run_type,run_cfg_name.Data(),ptCut));
    }

}


//================================================
void fitBackground(Int_t save = 0)
{
  const Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  const Int_t nPtBins = 3;
  const Double_t ptBins[nPtBins+1] = {0,2,4,10};
  const char *hName[3] = {"mhInvMvsPtUL","mhInvMvsPtLSpos","mhInvMvsPtLSneg"};

  TH2F *hnInvMass[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (TH2F*)f->Get(Form("%s_dimuon",hName[j]));
      hnInvMass[j]->Sumw2();
      draw2D(hnInvMass[j]);
    }
  hnInvMass[1]->Add(hnInvMass[2]);
 
  TH1F *hInvMass[2][nPtBins];
  TH1F *hSignal[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hnInvMass[j]->GetYaxis()->SetRangeUser(ptBins[i]+0.01,ptBins[i+1]-0.01);
	  hInvMass[j][i] = (TH1F*)hnInvMass[j]->ProjectionX(Form("%s_%s_InvMass_jpsi_%1.1f_%1.1f",hName[j],trigName[kTrigType],ptBins[i],ptBins[i+1]));
	  hInvMass[j][i]->Rebin(2);
	}
      hSignal[i] = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US_minus_LS_jpsi_%1.1f_%1.1f",ptBins[i],ptBins[i+1]));
      hSignal[i]->Add(hInvMass[1][i],-1);
      hSignal[i]->SetYTitle("US-LS");
    }
  //return;
  double nAll[nPtBins];
  double allErr[nPtBins];
  double nBkg[nPtBins];
  double bkgErr[nPtBins];
  double nSignal[nPtBins];
  double sigErr[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      TH1F *hFit = hInvMass[1][i];
      TF1 *func = new TF1(Form("func_%d",i),"pol3",2.5,3.5);
      if(i==0)func->SetRange(2.8,3.5);
      hFit->Fit(func,"IR0");
      hInvMass[0][i]->GetXaxis()->SetRangeUser(2.5,3.5);
      hInvMass[0][i]->SetMarkerStyle(21);
      hInvMass[0][i]->SetMarkerColor(2);
      hInvMass[0][i]->SetLineColor(2);
      c = draw1D(hInvMass[0][i],Form("Invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2});Counts"));
      hInvMass[1][i]->Draw("samesHIST");
      func->SetLineColor(4);
      func->Draw("sames");

      // Singal counts
      int low_bin = hInvMass[0][i]->FindFixBin(low_mass+0.001);
      int high_bin = hInvMass[0][i]->FindFixBin(high_mass-0.001);
      nAll[i] = hInvMass[0][i]->Integral(low_bin,high_bin);
      allErr[i] = TMath::Sqrt(nAll[i]);
      nBkg[i] = func->Integral(low_mass,high_mass) * 1./hInvMass[0][i]->GetBinWidth(1);
      bkgErr[i] = TMath::Sqrt(nBkg[i]);
      nSignal[i] = nAll[i] - nBkg[i];
      sigErr[i] = TMath::Sqrt(nAll[i]+nBkg[i]);
      TPaveText *signif = GetPaveText(0.13,0.3,0.12,0.3);
      signif->SetTextAlign(11);
      signif->SetTextFont(62);
      signif->AddText(Form("p_{T,1} > 1.5 GeV/c"));
      signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
      signif->AddText(Form("Signal = %3.1f #pm %3.1f",nSignal[i],sigErr[i]));
      signif->AddText(Form("S/B = %1.0f/%1.0f = 1:%3.0f",nSignal[i],nBkg[i],nBkg[i]/nSignal[i]));
      signif->Draw();
    }

}



