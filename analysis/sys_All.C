TFile *f;
const int year = YEAR;
TString file_name;
const int nSys = 5;
const char *name[nSys] = {"SigExt","TpcTracking","PidCuts","MtdTrigEff","MtdRespEff"};
const TString legName[nSys+1] = {"Total","Signal extraction", "Tpc tracking", "Muon PID", "Trigger efficiency", "MTD matching"};
const int color[nSys+1] = {1, 2, 4, 6, kCyan, kGreen+2};

//================================================
void sys_All()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //collectSys();
  mergeSystematics();
}


//================================================
void mergeSystematics(int savePlot = 1, int saveHisto = 1)
{
  TFile *fout = 0x0;
  if(saveHisto) fout = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type),"update");
  else          fout = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type),"read");

  TList *list = new TList;
  TH1F *hJpsiSysVsPt[nCentBins_pt][nSys+1];
  TH1F *hJpsiSysVsCent[nPtBins_npart][nSys+1];
  for(int s=0; s<nSys; s++)
    {
      for(int i=0; i<nCentBins_pt; i++)
	{
	  hJpsiSysVsPt[i][s+1] = (TH1F*)fout->Get(Form("JpsiSysVsPt_%s_cent%s",name[s],cent_Title_pt[i]));
	  if(s==0)
	    {
	      hJpsiSysVsPt[i][0] = (TH1F*)hJpsiSysVsPt[i][s+1]->Clone(Form("JpsiSysVsPt_All_cent%s",cent_Title_pt[i]));
	      hJpsiSysVsPt[i][0]->Reset();
	    }
	  for(int bin=1; bin<=hJpsiSysVsPt[i][0]->GetNbinsX(); bin++)
	    {
	      double err1 = hJpsiSysVsPt[i][s+1]->GetBinContent(bin);
	      double err2 = hJpsiSysVsPt[i][0]->GetBinContent(bin);
	      hJpsiSysVsPt[i][0]->SetBinContent(bin, sqrt(err1*err1+err2*err2));
	    }
	}

      for(int i=0; i<nPtBins_npart; i++)
	{
	  hJpsiSysVsCent[i][s+1] = (TH1F*)fout->Get(Form("JpsiSysVsCent_%s_Pt%s",name[s],pt_Name_npart[i]));
	  if(s==0)
	    {
	      hJpsiSysVsCent[i][0] = (TH1F*)hJpsiSysVsCent[i][s+1]->Clone(Form("JpsiSysVsCent_All_Pt%s",pt_Name_npart[i]));
	      hJpsiSysVsCent[i][0]->Reset();
	    }
	  for(int bin=1; bin<=hJpsiSysVsCent[i][0]->GetNbinsX(); bin++)
	    {
	      double err1 = hJpsiSysVsCent[i][s+1]->GetBinContent(bin);
	      double err2 = hJpsiSysVsCent[i][0]->GetBinContent(bin);
	      hJpsiSysVsCent[i][0]->SetBinContent(bin, sqrt(err1*err1+err2*err2));
	    }
	}
    }

  // Systematics vs. pt for 0-80%
  TCanvas *c = new TCanvas("SysVsPt_cent0080","SysVsPt_cent0080",800,600);
  TLegend *leg1[2];
  for(int i=0; i<2; i++)
    {
      leg1[i] = new TLegend(0.15+i*0.35, 0.7, 0.5+i*0.35, 0.85);
      leg1[i]->SetBorderSize(0);
      leg1[i]->SetFillColor(0);
      leg1[i]->SetTextSize(0.035);
    }
  for(int s=0; s<nSys+1; s++)
    {
      hJpsiSysVsPt[0][s]->SetLineColor(color[s]);
      hJpsiSysVsPt[0][s]->SetLineWidth(2);
      hJpsiSysVsPt[0][s]->GetYaxis()->SetRangeUser(0,0.3);
      hJpsiSysVsPt[0][s]->GetXaxis()->SetRangeUser(0.5,12);
      hJpsiSysVsPt[0][s]->SetTitle(";p_{T} (GeV/c)");
      if(s>0) hJpsiSysVsPt[0][s]->SetLineStyle(2);
      else    hJpsiSysVsPt[0][s]->SetLineStyle(1);
      if(s==0) hJpsiSysVsPt[0][s]->DrawCopy();
      else     hJpsiSysVsPt[0][s]->DrawCopy("sames");
      leg1[s/3]->AddEntry(hJpsiSysVsPt[0][s],legName[s],"L");
      printf("[i] %s: sys = %2.2f for pT = %2.1f\n",legName[s].Data(),hJpsiSysVsPt[0][s]->GetBinContent(1)*100,hJpsiSysVsPt[0][s]->GetBinCenter(1));
      printf("[i] %s: sys = %2.2f for pT = %2.1f\n",legName[s].Data(),hJpsiSysVsPt[0][s]->GetBinContent(6)*100,hJpsiSysVsPt[0][s]->GetBinCenter(6));
    }
  TPaveText *t1 = GetTitleText(Form("%s: systematic uncertainty of J/psi vs. p_{T} (%s%%)",run_type,cent_Name_pt[0]),0.04);
  t1->Draw();
  for(int i=0; i<2; i++)
    {
      leg1[i]->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/JpsiSysVsPt_0080.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_TotSysVsPt_cent0080.pdf"));
    }

  // Systematics vs. cent
  for(int i=0; i<nPtBins_npart; i++)
    {
      c = new TCanvas(Form("SysVsCent_Pt%s",pt_Name_npart[i]),Form("SysVsCent_Pt%s",pt_Name_npart[i]),800,600);
      for(int s=0; s<nSys+1; s++)
	{
	  hJpsiSysVsCent[i][s]->SetLineColor(color[s]);
	  hJpsiSysVsCent[i][s]->SetLineWidth(2);
	  hJpsiSysVsCent[i][s]->GetYaxis()->SetRangeUser(0,0.4);
	  hJpsiSysVsCent[i][s]->GetXaxis()->SetLabelSize(0.05);
	  if(s>0) hJpsiSysVsCent[i][s]->SetLineStyle(2);
	  else    hJpsiSysVsCent[i][s]->SetLineStyle(1);
	  if(s==0) hJpsiSysVsCent[i][s]->DrawCopy();
	  else     hJpsiSysVsCent[i][s]->DrawCopy("sames");
	}
      t1 = GetTitleText(Form("%s: systematic uncertainty of J/psi (%s)",run_type,pt_Title_npart[i]),0.04);
      t1->Draw();
      for(int j=0; j<2; j++)
	{
	  leg1[j]->Draw();
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/JpsiSysVsCent_Pt%s.pdf",run_type,pt_Name_npart[i]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_TotSysVsCent_Pt%1.0f.pdf",ptBins_low_npart[i]));
	}
    }


  // systematics vs. pt for different centrality
  c = new TCanvas("SysVsPt_centComp","SysVsPt_centComp",800,600);
  TLegend *leg = new TLegend(0.5, 0.6, 0.7, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiSysVsPt[i][0]->SetLineColor(color[i]);
      hJpsiSysVsPt[i][0]->SetLineWidth(2);
      hJpsiSysVsPt[i][0]->GetYaxis()->SetRangeUser(0,0.4);
      hJpsiSysVsPt[i][0]->SetTitle(";p_{T} (GeV/c)");
      hJpsiSysVsPt[i][0]->SetLineStyle(1);
      if(i==2 || i==3) hJpsiSysVsPt[i][0]->SetBinContent(9, 0);
      if(i==4)
	{
	  for(int bin=7; bin<=9; bin++)
	    hJpsiSysVsPt[i][0]->SetBinContent(bin, 0);
	}
      if(i==0) hJpsiSysVsPt[i][0]->DrawCopy();
      else     hJpsiSysVsPt[i][0]->DrawCopy("sames");
      leg->AddEntry(hJpsiSysVsPt[i][0],Form("%s%%",cent_Name_pt[i]),"L");
    }
  t1 = GetTitleText(Form("%s: total systematic uncertainty of J/psi vs. p_{T}",run_type),0.04);
  t1->Draw();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/JpsiSysVsPt_CompCent.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_TotSysVsPt_CompCent.pdf"));
    }

  if(saveHisto)
    {
      fout->cd();
      for(int i=0; i<nCentBins_pt; i++)
	{
	  hJpsiSysVsPt[i][0]->Write("",TObject::kOverwrite);
	}

      for(int i=0; i<nPtBins_npart; i++)
	{
	  hJpsiSysVsCent[i][0]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void collectSys(int saveHisto = 1)
{
// const int nSys = 5;
// const char *name[nSys] = {"SigExt","TpcTracking","PidCuts","MtdTrigEff","MtdRespEff"};

  const int nbins = nPtBins_pt -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low_pt[i+1];
  xbins[nbins] = ptBins_high_pt[nbins];

  TH1F *hJpsiSysVsPt[nCentBins_pt][nSys];
  TH1F *hJpsiSysVsCent[nPtBins_npart][nSys];

  for(int s=0; s<nSys; s++)
    {
      for(int i=0; i<nPtBins_npart; i++)
	{
	  hJpsiSysVsCent[i][s] = new TH1F(Form("JpsiSysVsCent_%s_Pt%s",name[s],pt_Name_npart[i]), "", nCentBins_npart[i], 0, nCentBins_npart[i]);
	  for(int bin=1; bin<=nCentBins_npart[i]; bin++)
	    {
	      hJpsiSysVsCent[i][s]->GetXaxis()->SetBinLabel(bin, Form("%s%%",cent_Name_npart[i*nCentBins_npart[0]+nCentBins_npart[i]-bin]));
	    }
	}
    }

  // signal extraction
  TFile *fSigExt = TFile::Open(Form("Rootfiles/%s.Sys.JpsiYield.root",run_type),"read");
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiSysVsPt[i][0] = (TH1F*)fSigExt->Get(Form("Sys_signalExt_%s",cent_Title_pt[i]));
      for(int bin=1; bin<=nbins; bin++)
	{
	  hJpsiSysVsPt[i][0]->SetBinContent(bin, hJpsiSysVsPt[i][0]->GetBinContent(bin)-1);
	}
    }
  for(int i=0; i<nPtBins_npart; i++)
    {
      TH1F *htmp = (TH1F*)fSigExt->Get(Form("Sys_signalExt_pt%s",pt_Name_npart[i]));
      for(int bin=1; bin<=nCentBins_npart[i]; bin++)
	{
	  hJpsiSysVsCent[i][0]->SetBinContent(bin, htmp->GetBinContent(nCentBins_npart[i]-bin+1)-1);
	}
    }

  // TpcTracking
  TFile *fTpc = TFile::Open(Form("Rootfiles/%s.Sys.TpcTracking.root",run_type),"read");
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiSysVsPt[i][1] = (TH1F*)fTpc->Get(Form("FinalSys_TpcTrackingVsPt_cent%s",cent_Title_pt[i]));
      for(int bin=1; bin<=nbins; bin++)
	{
	  hJpsiSysVsPt[i][1]->SetBinContent(bin, hJpsiSysVsPt[i][1]->GetBinContent(bin)-1);
	}
    }
  for(int i=0; i<nPtBins_npart; i++)
    {
      TH1F *htmp = (TH1F*)fTpc->Get(Form("FinalSys_TpcTrackingVsCent_Pt%s",pt_Name_npart[i]));
      for(int bin=1; bin<=nCentBins_npart[i]; bin++)
	{
	  hJpsiSysVsCent[i][1]->SetBinContent(bin, htmp->GetBinContent(nCentBins_npart[i]-bin+1)-1);
	}
    }

  // MuonPid
  TFile *fDtof = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"read");
  TH1F *hDtofSysVsPt   = (TH1F*)fDtof->Get("Run14_AuAu200_JpsiEffVsPt_Sys_Dtof0.75Eff");
  TH1F *hDtofSysVsCent = (TH1F*)fDtof->Get("Run14_AuAu200_JpsiEffVsCent_Sys_Dtof0.75Eff");

  TFile *fPid = TFile::Open(Form("Rootfiles/%s.Sys.MuonPid.root",run_type),"read");
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiSysVsPt[i][2] = (TH1F*)fPid->Get(Form("FinalSys_MuonPidVsPt_cent%s",cent_Title_pt[i]));
      for(int bin=1; bin<=nbins; bin++)
	{
	  double err1 = hJpsiSysVsPt[i][2]->GetBinContent(bin)-1;
	  double err2 = hDtofSysVsPt->GetBinError(bin);
	  hJpsiSysVsPt[i][2]->SetBinContent(bin, TMath::Sqrt(err1*err1+err2*err2));
	}
    }
  for(int i=0; i<nPtBins_npart; i++)
    {
      TH1F *htmp = (TH1F*)fPid->Get(Form("FinalSys_MuonPidVsCent_Pt%s",pt_Name_npart[i]));
      double err2 = hDtofSysVsCent->GetBinError(i+1);
      for(int bin=1; bin<=nCentBins_npart[i]; bin++)
	{
	  double err1 = htmp->GetBinContent(nCentBins_npart[i]-bin+1)-1;
	  hJpsiSysVsCent[i][2]->SetBinContent(bin, TMath::Sqrt(err1*err1+err2*err2));
	}
    }

  // MtdTrigEff
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"read");
  TH1F *hTrigEffSysVsPt   = (TH1F*)fTrigEff->Get("Run14_AuAu200_JpsiEffVsPt_Sys_TacDiffEff");
  TH1F *hTrigEffSysVsCent = (TH1F*)fTrigEff->Get("Run14_AuAu200_JpsiEffVsCent_Sys_TacDiffEff");
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiSysVsPt[i][3] = (TH1F*)hTrigEffSysVsPt->Clone(Form("JpsiSysVsPt_%s_cent%s",name[3],cent_Title_pt[i]));
      for(int bin=1; bin<=nbins; bin++)
	{
	  hJpsiSysVsPt[i][3]->SetBinContent(bin, hJpsiSysVsPt[i][3]->GetBinError(bin));
	  hJpsiSysVsPt[i][3]->SetBinError(bin, 0);
	}
    }
  for(int i=0; i<nPtBins_npart; i++)
    {
      for(int bin=1; bin<=nCentBins_npart[i]; bin++)
	{
	  hJpsiSysVsCent[i][3]->SetBinContent(bin, hTrigEffSysVsCent->GetBinError(i+1));
	}
    }

  // RespEff
  TFile *fRespEff = TFile::Open(Form("Rootfiles/%s.Sys.MtdRespEff.root",run_type),"read");
  TH1F *hRespEffSysVsPt   = (TH1F*)fRespEff->Get("Jpsi_hMtdRespEffSysAll");
  TH1F *hRespEffSysVsCent = (TH1F*)fRespEff->Get("Jpsi_Npart_hMtdRespEffSysAll");
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiSysVsPt[i][4] = (TH1F*)hRespEffSysVsPt->Clone(Form("JpsiSysVsPt_%s_cent%s",name[4],cent_Title_pt[i]));
      for(int bin=1; bin<=nbins; bin++)
	{
	  hJpsiSysVsPt[i][4]->SetBinContent(bin, hJpsiSysVsPt[i][4]->GetBinError(bin));
	  hJpsiSysVsPt[i][4]->SetBinError(bin, 0);
	}
    }
  for(int i=0; i<nPtBins_npart; i++)
    {
      for(int bin=1; bin<=nCentBins_npart[i]; bin++)
	{
	  hJpsiSysVsCent[i][4]->SetBinContent(bin, hRespEffSysVsCent->GetBinError(i+1));
	}
    }
  
  

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type),"recreate");
      for(int s=0; s<nSys; s++)
	{
	  for(int i=0; i<nCentBins_pt; i++)
	    {
	      hJpsiSysVsPt[i][s]->SetTitle("");
	      hJpsiSysVsPt[i][s]->SetName(Form("JpsiSysVsPt_%s_cent%s",name[s],cent_Title_pt[i]));
	      hJpsiSysVsPt[i][s]->Write("",TObject::kOverwrite);
	    }
	  for(int i=0; i<nPtBins_npart; i++)
	    {
	      hJpsiSysVsCent[i][s]->Write("",TObject::kOverwrite);
	    }
	}
    }
}
