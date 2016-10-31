
//================================================
void ana_Projection()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);


  InnerVsOuter();
}

//================================================
void InnerVsOuter(const int savePlot = 1)
{
  TFile *f = TFile::Open("output/Run14.AuAu200.Projection.DcaGeometry.root","read");
  
  const char *name[5] = {"Y","Z","Px","Py","Pz"};
  const int nHistos = 6;
  const double ptCuts[nHistos+1] = {1,1.2,1.5,2,3,5,10};
  const double min[5] = {-1,-1,-0.02,-0.02,-0.02};
  const double max[5] = {1,1,0.02,0.02,0.02};

  TH2F *hTrkDiffVsPt[5];
  TH1F *hTrkDiff[5][nHistos];
  for(int i=0; i<5; i++)
    {
      hTrkDiffVsPt[i] = (TH2F*)f->Get(Form("mhTrk%sDiff",name[i]));
      hTrkDiffVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      c1 = draw2D(hTrkDiffVsPt[i]);
      if(savePlot)
	c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Projection/Trk%sDiffVsPt.pdf",run_type,name[i]));

      TCanvas *c = new TCanvas(Form("%s",name[i]),Form("%s",name[i]),1100,700);
      c->Divide(3,2);
      for(int k=0; k<nHistos; k++)
	{
	  int low_bin = hTrkDiffVsPt[i]->GetXaxis()->FindFixBin(ptCuts[k]+0.001);
	  int hig_bin = hTrkDiffVsPt[i]->GetXaxis()->FindFixBin(ptCuts[k+1]-0.001);
	  hTrkDiff[i][k] = (TH1F*)hTrkDiffVsPt[i]->ProjectionY(Form("hTrkDiff_%d%d",i,k),low_bin,hig_bin);
	  c->cd(k+1);
	  if(i>1) gPad->SetLogy();
	  hTrkDiff[i][k]->Draw();
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{track} < %1.1f GeV/c",ptCuts[k],ptCuts[k+1]),0.055);
	  t1->Draw();
	  double accepted = hTrkDiff[i][k]->Integral(hTrkDiff[i][k]->FindFixBin(min[i]+1e-4),hTrkDiff[i][k]->FindFixBin(max[i]-1e-4));
	  double fraction = accepted/hTrkDiff[i][k]->Integral() * 100;
	  TPaveText *t1 = GetPaveText(0.7,0.75,0.7,0.8,0.05);
	  t1->AddText(Form("f = %2.1f%%",fraction));
	  t1->Draw();
	  double ymax = hTrkDiff[i][k]->GetMaximum() * 0.8;
	  TLine *line1 = GetLine(min[i],0,min[i],ymax);
	  line1->Draw();
	  TLine *line2 = GetLine(max[i],0,max[i],ymax);
	  line2->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Projection/Trk%sDiffInPtBins.pdf",run_type,name[i]));
    }
}
