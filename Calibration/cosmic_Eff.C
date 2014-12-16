void cosmic_Eff()
{
  gStyle->SetOptStat(0);
  //TFile *f = TFile::Open("Output/cosmic.eff.root","read");
  TFile *f = TFile::Open("Output/test.histos.root","read");
  TH2F *hTrkPtVsBL[2];
  hTrkPtVsBL[0] = (TH2F*)f->Get("nProBlvsPt");
  hTrkPtVsBL[1] = (TH2F*)f->Get("nMatchBlvsPt2"); 

  TH1F *hTrkPt[2][2];
  for(Int_t i=0; i<2; i++)
    {
      //draw2D(hTrkPtVsBL[i]);

      hTrkPt[i][0] = (TH1F*)hTrkPtVsBL[i]->ProjectionX(Form("hTrkPt_%d_Bottom",i),9,23);
      hTrkPt[i][1] = (TH1F*)hTrkPtVsBL[i]->ProjectionX(Form("hTrkPt_%d_Top",i),1,7);
      for(Int_t j=0; j<2; j++)
	{
	  hTrkPt[i][j]->Sumw2();
	  hTrkPt[i][j]->SetMarkerStyle(21);
	  hTrkPt[i][j]->SetMarkerColor(color[i]);
	  hTrkPt[i][j]->SetLineColor(color[i]);
	}
      hTrkPt[i][1]->Add((TH1F*)hTrkPtVsBL[i]->ProjectionX(Form("hTrkPt_%d_Top_2",i),24,30));
    }

  const char *title[2] = {"bottom","top"};
  for(Int_t j=0; j<2; j++)
    {
      c = draw1D(hTrkPt[0][j],Form("p_{T} of tracks matched to backlegs on %s;p_{T} (GeV/c);Counts",title[j]));   
      hTrkPt[1][j]->Draw("samesP");
    }

  TLegend *leg = new TLegend(0.25,0.15,0.45,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  TH1F *hTrkPtClone[2][2];
  TGraphAsymmErrors *hEff[2];
  for(Int_t j=0; j<2; j++)
    {
      hEff[j] = new TGraphAsymmErrors();
      hEff[j]->SetName(Form("hEff_%s",title[j]));
      for(Int_t i=0; i<2; i++)
	{
	  hTrkPtClone[i][j] = (TH1F*)hTrkPt[i][j]->Clone(Form("%s_clone",hTrkPt[i][j]->GetName()));
	}
      hEff[j]->Divide(hTrkPt[1][j],hTrkPt[0][j]);

      hTrkPtClone[1][j]->Divide(hTrkPtClone[0][j]);
      hTrkPtClone[1][j]->SetMarkerStyle(25-j*4);
      hTrkPtClone[1][j]->SetMarkerColor(color[j]);
      hTrkPtClone[1][j]->SetLineColor(color[j]);
      hTrkPtClone[1][j]->GetYaxis()->SetRangeUser(0,1.2);
      leg->AddEntry(hTrkPtClone[1][j],title[j],"PLE");
      
      hEff[j]->SetMarkerStyle(25-j*4);
      hEff[j]->SetMarkerColor(color[j]);
      hEff[j]->SetLineColor(color[j]);
      hEff[j]->SetMaximum(1.1);
    }
  c = draw1D(hTrkPtClone[1][0]);
  hTrkPtClone[1][1]->Draw("sames");
  leg->Draw();

  //cout << hEff[j]->GetErrorYhigh(0) << endl;
  // drawGraph(hEff[0],"Track matching efficiency to MTD;p_{T} (GeV/c);Efficiency");
  // hEff[1]->Draw("samesP");
  // leg->Draw();
}
