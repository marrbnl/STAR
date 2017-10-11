TFile *f = 0x0;
const int year = YEAR;

//================================================
void ana_LowPtJpsi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.15);

  f = TFile::Open(Form("./output/%s.jpsi.root",run_type),"read");
  run_cfg_name = Form("%s",run_config);

  if(f)
    {
      TH1F *hStat = (TH1F*)f->Get("hEventStat");
      printf("all         events: %4.4e\n",hStat->GetBinContent(1));
      printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
      printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
    }

  lowPt();
}

//================================================
void lowPt(const int savePlot = 1)
{
  const int nPtBins = 2;
  const double ptBinsLow[nPtBins] = {0, 0.2};
  const double ptBinsHigh[nPtBins] = {1, 1};

  // same-event
  TH1F *hInvMass[nPtBins] = {0x0};
  THnSparseF *hnInvMass = (THnSparseF*)f->Get("mhJpsiInfoLowPtWeight_di_mu");
  hnInvMass->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
  hnInvMass->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
  hnInvMass->GetAxis(4)->SetRange(1, 4);
  for(Int_t i=0; i<nPtBins; i++) // pt bins
    {
      hnInvMass->GetAxis(1)->SetRangeUser(ptBinsLow[i]+0.01,ptBinsHigh[i]-0.01);
      hInvMass[i] = (TH1F*)hnInvMass->Projection(0);
      hInvMass[i]->SetName(Form("JpsiInvMass_PtBin%d",i));
      hInvMass[i]->Sumw2();
      hnInvMass->GetAxis(1)->SetRange(0,-1);
    }

  // mix-event
  TFile *fmix = TFile::Open(Form("Output/%s.Mix.%spt%1.1f.pt%1.1f.root",run_type,run_config,pt1_cut,pt2_cut),"read");
  TH1F *hMixInvMass[nPtBins];
  TH3D *hMixMmumuvsPtCen = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
  for(int i=0; i<nPtBins; i++)
    {
      int ybin_min = hMixMmumuvsPtCen->GetYaxis()->FindFixBin(ptBinsLow[i]+1e-4);
      int ybin_max = hMixMmumuvsPtCen->GetYaxis()->FindFixBin(ptBinsHigh[i]-1e-4);
      TH1F *htmp = (TH1F*)hMixMmumuvsPtCen->ProjectionZ(Form("mix_JpsiInvMass_PtBin%d_tmp",i),1,4,ybin_min,ybin_max);
      hMixInvMass[i] = new TH1F(Form("mix_JpsiInvMass_PtBin%d",i),htmp->GetTitle(),200,2,4);
      for(int bin=1; bin<=hMixInvMass[i]->GetNbinsX(); bin++)
	{
	  int jbin = htmp->FindFixBin(hMixInvMass[i]->GetBinCenter(bin));
	  hMixInvMass[i]->SetBinContent(bin,htmp->GetBinContent(jbin));
	  hMixInvMass[i]->SetBinError(bin,htmp->GetBinError(jbin));
	}
    }

  TH1F *hJpsi[nPtBins];
  TF1 *funcJpsi[nPtBins];
  const double min_mass[2] = {2.6, 3.3};
  const double max_mass[2] = {2.9, 3.6};
  for(Int_t i=0; i<nPtBins; i++) // pt bins
    {
      hInvMass[i]->Rebin(4);
      hInvMass[i]->SetMarkerStyle(20);
      hInvMass[i]->GetXaxis()->SetRangeUser(2.5,4);
      hInvMass[i]->GetYaxis()->SetRangeUser(-50,300);
     c = draw1D(hInvMass[i],Form("%s: invariant mass distribution (%1.1f < p_{T} < %1.0f GeV/c, 60-80%%)",run_type,ptBinsLow[i],ptBinsHigh[i]));

      hMixInvMass[i]->Rebin(4);
      hMixInvMass[i]->SetMarkerStyle(24);
      double counts_se = 0, counts_me = 0;
      for(int j=0; j<2; j++)
	{
	  int min_bin = hInvMass[i]->FindFixBin(min_mass[j]);
	  int max_bin = hInvMass[i]->FindFixBin(max_mass[j]);
	  for(int bin=min_bin; bin<=max_bin; bin++)
	    {
	      counts_se += hInvMass[i]->GetBinContent(bin);
	      counts_me += hMixInvMass[i]->GetBinContent(bin);
	    }
	}

      hMixInvMass[i]->Scale(counts_se/counts_me);
      hMixInvMass[i]->Draw("sames");

      hJpsi[i] = (TH1F*)hInvMass[i]->Clone(Form("Jpsi_Pt%d",i));
      hJpsi[i]->Add(hMixInvMass[i], -1);
      hJpsi[i]->SetMarkerStyle(21);
      hJpsi[i]->SetMarkerColor(2);
      hJpsi[i]->SetLineColor(2);
      hJpsi[i]->Draw("sames");

      funcJpsi[i] = new TF1(Form("FuncJpsi_Pt%d",i),"gausn(0)+pol3(3)",2.5,4);
      funcJpsi[i]->SetParameter(1, 3.09);
      funcJpsi[i]->SetParameter(2, 0.05);
      hJpsi[i]->Fit(funcJpsi[i], "IR0Q");
      funcJpsi[i]->SetLineStyle(2);
      funcJpsi[i]->SetLineColor(4);
      funcJpsi[i]->Draw("sames");
      TPaveText *t1 = GetPaveText(0.15,0.3,0.7,0.75,0.04);
      t1->AddText(Form("N_{J/#Psi} = %1.0f #pm %1.0f",funcJpsi[i]->GetParameter(0)/hJpsi[i]->GetBinWidth(1), funcJpsi[i]->GetParError(0)/hJpsi[i]->GetBinWidth(1)));
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/JpsiCounts_%1.1f-%1.0f_cent6080.pdf",run_type,ptBinsLow[i],ptBinsHigh[i]));
    }
}
