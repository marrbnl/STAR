
//================================================
void compare_EventMixing()
{
  gStyle->SetOptStat(1);

  run_type = "Run14_AuAu200";

  compareToShuai();
  //makeShuai();
}

//================================================
void compareToShuai(const int save = 0)
{
  const char *type_name[2] = {"","Mix_"};
  const int type = 1;

  TH1F *hMix[2][nPtBins][nCentBins][2];

  // My mixed events
  TFile *f = TFile::Open("Rootfiles/Pico.Run14.AuAu200.jpsi.root","read");
  //TFile *f = TFile::Open("Rootfiles/Shuai.Pico.Run14.AuAu200.jpsi.mix.check.20List.root","read");
  for(int i=0; i<nPtBins; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hMix[0][i][j][0] = (TH1F*)f->Get(Form("%sInvMass_UL_pt%s_cent%s",type_name[type],pt_Name[i],cent_Name[j]));
	  hMix[1][i][j][0] = (TH1F*)f->Get(Form("%sInvMass_LS_pt%s_cent%s",type_name[type],pt_Name[i],cent_Name[j]));
	}
    }

  //return;
  // Shuai's mixed events
  TFile *fShuai = TFile::Open("Rootfiles/Shuai.Pico.Run14.AuAu200.jpsi.mix.check.20List.root","read");
  for(int i=0; i<nPtBins; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hMix[0][i][j][1] = (TH1F*)fShuai->Get(Form("Shuai_%sInvMass_UL_pt%s_cent%s",type_name[type],pt_Name[i],cent_Name[j]));
	  hMix[1][i][j][1] = (TH1F*)fShuai->Get(Form("Shuai_%sInvMass_LS_pt%s_cent%s",type_name[type],pt_Name[i],cent_Name[j]));
	}
    }
 
  // compare
  const char *name[2] = {"unlike-sign","like-sign"};
  const char *title[2] = {"US","LS"};
  for(int m=0; m<2; m++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  TCanvas *cMix = new TCanvas(Form("%s%s_%s",type_name[type],name[m],cent_Name[j]),Form("%s%s_%s",type_name[type],name[m],cent_Name[j]),1100,650);
	  cMix->Divide(3,2);

	  TCanvas *cRatio = new TCanvas(Form("ratio_%s%s_%s",type_name[type],name[m],cent_Name[j]),Form("ratio_%s%s_%s",type_name[type],name[m],cent_Name[j]),1100,650);
	  cRatio->Divide(3,2);

	  for(int i=0; i<nPtBins; i++)
	    {
	      for(int k=0; k<2; k++)
		{
		  hMix[m][i][j][k]->Sumw2();
		  hMix[m][i][j][k]->Rebin(20);
		  //hMix[m][i][j][k]->Scale(1./hMix[m][i][j][k]->GetBinContent(hMix[m][i][j][k]->FindFixBin(3)));
		  hMix[m][i][j][k]->Scale(1./hMix[m][i][j][k]->Integral());
		  hMix[m][i][j][k]->SetMarkerStyle(20);
		  hMix[m][i][j][k]->SetMarkerColor(color[k]);
		  hMix[m][i][j][k]->SetLineColor(color[k]);
		  hMix[m][i][j][k]->SetTitle("");
		}
	      printf("Rongrong: %s %s %1.0f\n",pt_Name[i],cent_Name[j],hMix[m][i][j][0]->GetEntries());
	      printf("Shuai: %s %s %1.0f\n",pt_Name[i],cent_Name[j],hMix[m][i][j][1]->GetEntries());
	      cMix->cd(i+1);
	      hMix[m][i][j][0]->Draw("P");
	      hMix[m][i][j][1]->Draw("samesP");
	      TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[j]),0.06);
	      t->Draw();
	      if(i==0)
		{
		  leg = new TLegend(0.45,0.63,0.8,0.85);
		  leg->SetBorderSize(0);
		  leg->SetFillColor(0);
		  leg->SetTextSize(0.05);
		  leg->SetHeader(Form("Mix %s",name[m]));
		  leg->AddEntry(hMix[m][i][j][0],"Rongrong","P");
		  leg->AddEntry(hMix[m][i][j][1],"Shuai","P");
		  leg->Draw();
		}

	      TH1F *hRatio = (TH1F*)hMix[m][i][j][1]->Clone(Form("Ratio_%sInvMass_%s_pt%s_cent%s",type_name[type],name[m],pt_Name[i],cent_Name[j]));
	      for(int bin=1; bin<=hRatio->GetNbinsX(); bin++)
		{
		  double x1 = hMix[m][i][j][1]->GetBinContent(bin);
		  double e1 = hMix[m][i][j][1]->GetBinError(bin);
		  double x2 = hMix[m][i][j][0]->GetBinContent(bin);
		  double e2 = hMix[m][i][j][0]->GetBinError(bin);
		  double x,e;
		  if(x2>0 && x1>0)
		    {
		      x = x1/x2;
		      e = x*TMath::Sqrt(e1*e1/x1/x1+e2*e2/x2/x2);
		    }
		  else
		    {
		      x = 0;
		      e = 0;
		    }
		  hRatio->SetBinContent(bin,x);
		  hRatio->SetBinError(bin,e);
		}
	      cRatio->cd(i+1);
	      hRatio->GetYaxis()->SetRangeUser(0.8,1.2);
	      hRatio->GetXaxis()->SetRangeUser(0,4);
	      hRatio->Draw();
	      TLine *line = GetLine(0,1,4,1,1);
	      line->Draw();
	      t->Draw();
	      if(i==0)
		{
		  t = GetPaveText(0.2,0.7,0.75,0.85,0.05);
		  t->SetTextFont(62);
		  t->AddText(Form("Mix %s: Shuai/Rongrong",name[m]));
		  t->Draw();
		}
	    }
	  if(save) 
	    {
	      cMix->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/compare.%s%s_Cent%s.pdf",run_type,type_name[type],title[m],cent_Title[j]));
	      cMix->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/compare.%s%s_Cent%s.png",run_type,type_name[type],title[m],cent_Title[j]));
	    }
	  if(save) 
	    {
	      cRatio->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/compare.%s%s_ratio_Cent%s.pdf",run_type,type_name[type],title[m],cent_Title[j]));
	      cRatio->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/compare.%s%s_ratio_Cent%s.png",run_type,type_name[type],title[m],cent_Title[j]));
	    }
	}
    }
}

//================================================
void makeShuai()
{
  TFile *fShuai = TFile::Open("Rootfiles/5ListsPerJob.root","read");
  TH3D *hMixULMmumuvsPtCen = (TH3D*)fShuai->Get("hMixULMmumuvsPtCen");
  TH3D *hMixLPosMmumuvsPtCen = (TH3D*)fShuai->Get("hMixLPosMmumuvsPtCen");
  TH3D *hMixLNegMmumuvsPtCen = (TH3D*)fShuai->Get("hMixLNegMmumuvsPtCen");
  TH1F *hMixUL[nPtBins][nCentBins];
  TH1F *hMixLS[nPtBins][nCentBins];
  TH1F *hMixLSpos[nPtBins][nCentBins];
  TH1F *hMixLSneg[nPtBins][nCentBins];
  for(int i=0; i<nPtBins; i++)
    {
      int ybin_min = hMixULMmumuvsPtCen->GetYaxis()->FindFixBin(ptBins_low[i]+1e-4);
      int ybin_max = hMixULMmumuvsPtCen->GetYaxis()->FindFixBin(ptBins_high[i]-1e-4);
      for(int j=0; j<nCentBins; j++)
	{
	  hMixUL[i][j] = (TH1F*)hMixULMmumuvsPtCen->ProjectionZ(Form("Shuai_Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[j]),centBins_low[j],centBins_high[j],ybin_min,ybin_max);
	  hMixLSpos[i][j] = (TH1F*)hMixLPosMmumuvsPtCen->ProjectionZ(Form("Shuai_Mix_InvMass_LSpos_pt%s_cent%s",pt_Name[i],cent_Name[j]),centBins_low[j],centBins_high[j],ybin_min,ybin_max);
	  hMixLSneg[i][j] = (TH1F*)hMixLNegMmumuvsPtCen->ProjectionZ(Form("Shuai_Mix_InvMass_LSneg_pt%s_cent%s",pt_Name[i],cent_Name[j]),centBins_low[j],centBins_high[j],ybin_min,ybin_max);
	  hMixLS[i][j] = (TH1F*)hMixLSpos[i][j]->Clone(Form("Shuai_Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[j]));
	  hMixLS[i][j]->Add(hMixLSneg[i][j]);
	}
    }

  TH3D *hULMmumuvsPtCen = (TH3D*)fShuai->Get("hULMmumuvsPtCen");
  TH3D *hLPosMmumuvsPtCen = (TH3D*)fShuai->Get("hLPosMmumuvsPtCen");
  TH3D *hLNegMmumuvsPtCen = (TH3D*)fShuai->Get("hLNegMmumuvsPtCen");
  TH1F *hUL[nPtBins][nCentBins];
  TH1F *hLS[nPtBins][nCentBins];
  TH1F *hLSpos[nPtBins][nCentBins];
  TH1F *hLSneg[nPtBins][nCentBins];
  for(int i=0; i<nPtBins; i++)
    {
      int ybin_min = hULMmumuvsPtCen->GetYaxis()->FindFixBin(ptBins_low[i]+1e-4);
      int ybin_max = hULMmumuvsPtCen->GetYaxis()->FindFixBin(ptBins_high[i]-1e-4);
      for(int j=0; j<nCentBins; j++)
	{
	  hUL[i][j] = (TH1F*)hULMmumuvsPtCen->ProjectionZ(Form("Shuai_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[j]),centBins_low[j],centBins_high[j],ybin_min,ybin_max);
	  hLSpos[i][j] = (TH1F*)hLPosMmumuvsPtCen->ProjectionZ(Form("Shuai_InvMass_LSpos_pt%s_cent%s",pt_Name[i],cent_Name[j]),centBins_low[j],centBins_high[j],ybin_min,ybin_max);
	  hLSneg[i][j] = (TH1F*)hLNegMmumuvsPtCen->ProjectionZ(Form("Shuai_InvMass_LSneg_pt%s_cent%s",pt_Name[i],cent_Name[j]),centBins_low[j],centBins_high[j],ybin_min,ybin_max);
	  hLS[i][j] = (TH1F*)hLSpos[i][j]->Clone(Form("Shuai_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[j]));
	  hLS[i][j]->Add(hLSneg[i][j]);
	}
    }

  TFile *fout = TFile::Open("Rootfiles/Shuai.Pico.Run14.AuAu200.jpsi.mix.check.5List.root","recreate");
  for(int i=0; i<nPtBins; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hUL[i][j]->Write();
	  hLS[i][j]->Write();
	  hLSpos[i][j]->Write();
	  hLSneg[i][j]->Write();
	  hMixUL[i][j]->Write();
	  hMixLS[i][j]->Write();
	  hMixLSpos[i][j]->Write();
	  hMixLSneg[i][j]->Write();
	}
    }
}
