TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;

const char *run_config = "looseCut.";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void ana_PtBin()
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

  ptBins();
}

//================================================
void ptBins(Int_t save = 0)
{
  const Double_t pt1_cut = 1.5, pt2_cut = 1.0;
  const Int_t nPtBins = 3;
  const Double_t ptBins[nPtBins+1] = {0,2,4,10};
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};

  THnSparseF *hnInvMass[3];
  TH1F *hJpsiPt[3];
  TH1F *hInvMass[3][nPtBins];
  TH1F *hSignal[nPtBins];

  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);

      for(Int_t i=0; i<nPtBins; i++)
	{
	  hnInvMass[j]->GetAxis(1)->SetRangeUser(ptBins[i]+0.01,ptBins[i+1]-0.01);
	  TH1F *htmp = (TH1F*)hnInvMass[j]->Projection(0);
	  htmp->SetName(Form("%s_%s_InvMass_jpsi_%1.1f_%1.1f",hName[j],trigName[kTrigType],ptBins[i],ptBins[i+1]));
	  htmp->Sumw2();
	  hInvMass[j][i] = (TH1F*)htmp->Rebin(nSpecMBinsJpsi, Form("%s_jpsi",htmp->GetName()),specMJpsi);
	  scaleHisto(hInvMass[j][i], 1, 1,kTRUE,kFALSE, kTRUE);
	  hnInvMass[j]->GetAxis(1)->SetRange(0,-1);
	}
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }

  for(Int_t i=0; i<nPtBins; i++)
    {
      hInvMass[1][i]->Add(hInvMass[2][i]);
      hSignal[i] = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US_minus_LS_jpsi_%1.1f_%1.1f",ptBins[i],ptBins[i+1]));
      hSignal[i]->Add(hInvMass[1][i],-1);
      hSignal[i]->SetYTitle("US-LS");
    }

  TCanvas *canvas[nPtBins/6+1];
  for(Int_t i=0; i<nPtBins; i++)
    {
      if(i%6==0) 
	{
	  if(nPtBins<5)
	    {
	      canvas[i/6] = new TCanvas(Form("canvas_%d",i/6),Form("canvas_%d",i/6),1100,700);
	      canvas[i/6]->Divide(2,2);
	    }
	  else
	    {
	      canvas[i/6] = new TCanvas(Form("canvas_%d",i/6),Form("canvas_%d",i/6),1100,700);
	      canvas[i/6]->Divide(3,2);
	    }
	}
      canvas[i/6]->cd(i%6+1);
      SetPadMargin(gPad,0.13,0.13);

      hInvMass[0][i]->GetXaxis()->SetRangeUser(2.5,3.5);
      hInvMass[0][i]->SetMaximum(1.4*hInvMass[0][i]->GetMaximum());
      hInvMass[0][i]->SetMinimum(0);
      hInvMass[0][i]->SetMarkerColor(1);
      hInvMass[0][i]->SetLineColor(1);
      ScaleHistoTitle(hInvMass[0][i],0.06,0.9,0.05,0.06,1.2,0.05,62);
      hInvMass[0][i]->SetTitle(";M_{#mu#mu} (GeV/c^{2});dN/dM");

      hInvMass[1][i]->SetMarkerColor(4);
      hInvMass[1][i]->SetLineColor(4);

      hSignal[i]->SetMarkerStyle(21);
      hSignal[i]->SetMarkerColor(2);
      hSignal[i]->SetLineColor(2);

      hInvMass[0][i]->Draw("HIST");
      hInvMass[1][i]->Draw("sames HIST");
      hSignal[i]->Draw("sames");

      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T,J/#psi} < %1.1f GeV/c",ptBins[i],ptBins[i+1]),0.06);
      t1->Draw();

      if(i%6==0)
	{
	  leg = new TLegend(0.15,0.23,0.3,0.45);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.06);
	  leg->AddEntry(hInvMass[0][i],"Unlike sign","L");
	  leg->AddEntry(hInvMass[1][i],"Like sign (++)+(--)","L");
	  leg->AddEntry(hSignal[i],"UL-LS","P");
	  leg->Draw();
	}
      
      Int_t low_bin = hInvMass[0][i]->GetXaxis()->FindFixBin(low_mass+0.001);
      Int_t high_bin = hInvMass[0][i]->GetXaxis()->FindFixBin(high_mass-0.001);
      double errorLS, errorUL;
      double nLS = hInvMass[1][i]->IntegralAndError(low_bin,high_bin,errorLS,"width");
      double nUL = hInvMass[0][i]->IntegralAndError(low_bin,high_bin,errorUL,"width");
      double nSignal = nUL - nLS;
      double nBackground = nLS;
      double errorSignal = TMath::Sqrt(errorLS*errorLS+errorUL*errorUL);

      t1 = GetPaveText(0.15,0.7,0.7,0.85,0.055);
      t1->AddText(Form("S/B = %1.0f/%1.0f = 1:%3.0f",nSignal,nBackground,nBackground/nSignal));
      t1->AddText(Form("Significance: %2.1f",nSignal/errorSignal));
      t1->SetTextAlign(11);
      t1->SetTextFont(62);
      t1->Draw();
      if((i%6==5 || i==nPtBins-1) && save)
	{
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_Jpsi_PtBins_pt1_%1.0f_pt2_%1.0f_Bin%d.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10,i/6));
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_Jpsi_PtBins_pt1_%1.0f_pt2_%1.0f_Bin%d.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10,i/6));
	}
    }

  if(nPtBins==3)
    {
      canvas[0]->cd(4);
      t1 = GetPaveText(0.2,0.5,0.5,0.8,0.07);
      t1->SetTextAlign(11);
      t1->AddText(Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f",pt1_cut,pt2_cut));
      t1->AddText(Form("Mass window: %1.1f < M_{#mu#mu} < %1.1f",low_mass,high_mass));
      t1->SetTextFont(62);
      t1->Draw();
    }
}

