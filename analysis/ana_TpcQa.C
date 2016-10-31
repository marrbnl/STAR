
TFile *f;
const int year = YEAR;

//================================================
void ana_TpcQa()
{
  gStyle->SetOptStat(0);

  if(year==2014) f = TFile::Open("output/Run14.AuAu200.TpcQa.root","read");

  TpcSector20();
}


//================================================
void TpcSector20(const int savePlot = 1)
{
  // phi vs eta for pT > 1 GeV/c
  const double ptCut = 1;

  THnSparseF *hnTrkInfo[3][2];
  for(int t=0; t<3; t++)
    {
      hnTrkInfo[t][0] = (THnSparseF*)f->Get(Form("mhQaTrkEtaPhi_di_mu%s",gTrgSetupTitle[t+2]));
      hnTrkInfo[t][1] = (THnSparseF*)f->Get(Form("mhQaTrkEtaPhiWithCuts_di_mu%s",gTrgSetupTitle[t+2]));
    }
  const char *cut_title[2] = {"w/o","w"};
  const char *cut_name[2] = {"wo","w"};
  for(int t=0; t<3; t++)
    {
      for(int i=0; i<2; i++)
	{
	  hnTrkInfo[t][i]->GetAxis(0)->SetRangeUser(1,100);
	  TH2F *h2 = (TH2F*)hnTrkInfo[t][i]->Projection(2,1);
	  h2->SetName(Form("TrkEtaPhi_%d_%d",t,i));
	  h2->RebinY(4);
	  c = draw2D(h2,Form("AuAu200%s (%s quality cuts, p_{T} > %1.0f GeV/c);#eta;#varphi",gTrgSetupTitle[t+2],cut_title[i],ptCut),0.04,false);
	  if(savePlot)
	    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TpcQa/TrkPhiEta_%sCuts%s.pdf",run_type,cut_name[i],gTrgSetupTitle[t+2]));
	}
    }

  // track quality cuts
  const int nHistos = 5;
  TH2F *hTrkCuts[4][3][2];
  TH1F *hTrkDis[4][3][2][nHistos];
  const double ptCuts[nHistos+1] = {0.2,0.5,1,1.5,2,10};
  const char *var_name[4] = {"NHitsFit","NHitsDedx","Dca","NSigmaPi"};
  const char *var_title[4] = {"NHitsFit","NHitsDedx","Dca","n#sigma_{#pi}"};
  const char *good_name[2] = {"bad","good"};
  const char *legend[3][2] = { {"prod_low: sector 20","prod_low: other sectors"}, 
			       {"prod_mid: sector 20","prod_mid: other sectors"},
			       {"prod_high: sector 20","prod_high: other sectors"}};
  for(int j=0; j<4; j++)
    {
      TCanvas *c = new TCanvas(Form("%s",var_name[j]),Form("%s",var_name[j]),1100,700);
      c->Divide(3,2);
      TLegend *leg = new TLegend(0.1,0.3,0.7,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06);
      leg->SetHeader(var_title[j]);
      for(int t=0; t<3; t++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hTrkCuts[j][t][i] = (TH2F*)f->Get(Form("mhQaTrk%s_di_mu%s_%s",var_name[j],gTrgSetupTitle[t+2],good_name[i]));
	      hTrkCuts[j][t][i]->Sumw2();
	      for(int k=0; k<nHistos; k++)
		{
		  int low_bin = hTrkCuts[j][t][i]->GetXaxis()->FindFixBin(ptCuts[k]+0.001);
		  int hig_bin = hTrkCuts[j][t][i]->GetXaxis()->FindFixBin(ptCuts[k+1]-0.001);
		  hTrkDis[j][t][i][k] = (TH1F*)hTrkCuts[j][t][i]->ProjectionY(Form("hTrkDis_%d%d%d%d",j,t,i,k),low_bin,hig_bin);
		  hTrkDis[j][t][i][k]->Scale(1./hTrkDis[j][t][i][k]->Integral());
		  hTrkDis[j][t][i][k]->SetLineColor(color[t]);
		  hTrkDis[j][t][i][k]->SetMarkerColor(color[t]);
		  hTrkDis[j][t][i][k]->SetMarkerStyle(20+t);
		  hTrkDis[j][t][i][k]->SetTitle(Form(";%s;prob.",var_title[j]));
		  if(j==1) hTrkDis[j][t][i][k]->GetXaxis()->SetRangeUser(5,50);
		  if(j==0 || j==1) hTrkDis[j][t][i][k]->SetMaximum(1.5*hTrkDis[j][t][i][k]->GetMaximum());
		  if(j==2) hTrkDis[j][t][i][k]->SetMaximum(2*hTrkDis[j][t][i][k]->GetMaximum());
		  if(j==3) hTrkDis[j][t][i][k]->SetMaximum(1.2*hTrkDis[j][t][i][k]->GetMaximum());

		  c->cd(k+1);
		  if(t==0 && i==0)
		    {
		      hTrkDis[j][t][i][k]->Draw("PE");
		      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{track} < %1.1f GeV/c",ptCuts[k],ptCuts[k+1]),0.055);
		      t1->Draw();
		    }
		  else
		    {
		      if(i==1) hTrkDis[j][t][i][k]->Draw("HIST same");
		      else     hTrkDis[j][t][i][k]->Draw("PE   same");
		    }
		  if(k==0) 
		    {
		      if(i==1) leg->AddEntry(hTrkDis[j][t][i][k],legend[t][i],"L");
		      else     leg->AddEntry(hTrkDis[j][t][i][k],legend[t][i],"P");
		    }
		}
	    }
	}
      c->cd(6);
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TpcQa/Trk%s.pdf",run_type,var_name[j]));
    }
}
