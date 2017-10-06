TFile *f;
const int year = YEAR;
// const int nSys = 12;
// const TString sysName[nSys] = {"default",
// 			       "dcaUp","dcaDown","NHitsUp","NDedxUp",
// 			       "dcaUp_NHitsUp", "dcaUp_NDedxUp", "dcaUp_NHitsUp_NDedxUp",
// 			       "dcaDown_NHitsUp", "dcaDown_NDedxUp", "dcaDown_NHitsUp_NDedxUp",
// 			       "NHitsUp_NDedxUp"};
// const TString outName = "TpcTracking";
// const int nSysMarker = 5;

const int nSys = 27;
const TString sysName[nSys] = {"default",
			       "dzUp","dzDown","dyUp","dyDown","nSigPiUp","nSigPiDown",
			       "dzUp_dyUp","dzUp_dyDown","dzUp_nSigPiUp","dzUp_nSigPiDown","dzUp_dyUp_nSigPiUp","dzUp_dyUp_nSigPiDown","dzUp_dyDown_nSigPiUp","dzUp_dyDown_nSigPiDown",
			       "dzDown_dyUp","dzDown_dyDown","dzDown_nSigPiUp","dzDown_nSigPiDown","dzDown_dyUp_nSigPiUp","dzDown_dyUp_nSigPiDown","dzDown_dyDown_nSigPiUp","dzDown_dyDown_nSigPiDown",
			       "dyUp_nSigPiUp","dyUp_nSigPiDown",
			       "dyDown_nSigPiUp","dyDown_nSigPiDown"};
const TString outName = "MuonPid";
const int nSysMarker = 7;

const int color[30] = {1, 2, 4, 6, 8, 1, 2, 4, 6, 8, 1, 2, 4, 6, 8, 1, 2, 4, 6, 8, 1, 2, 4, 6, 8, 1, 2, 4, 6, 8};
const char* typeName[2] = {"UL", "LS"};

//================================================
void sys_TpcAndPid()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //mergeHistosData();
  //mergeHistosEmb();
  //makeHistos();
  anaSys();
}

//================================================
void anaSys(const int savePlot = 1, const int saveHisto = 1)
{
  const int uncerMode = 0; // 0 - max; 1 - rms
  const double min_mass = 3.0, max_mass = 3.2;
  TFile *fsys = 0;
  if(saveHisto) fsys = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.%s.root",outName.Data()),"update");
  else          fsys = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.%s.root",outName.Data()),"read");

  if(uncerMode==0) { const int nSysUsed = nSysMarker; }
  if(uncerMode==1) { const int nSysUsed = nSys; }

  //----------------------------------------------
  // pT dependence
  //----------------------------------------------
  const int nbins = nPtBins_pt -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low_pt[i+1];
  xbins[nbins] = ptBins_high_pt[nbins];


  TH1F *hInvMassVsPt[nSysUsed][nCentBins_pt][nbins][2];
  TH1F *hJpsiCountVsPt[nSysUsed][nCentBins_pt];
  for(int s=0; s<nSysUsed; s++)
    {
      for(int k=0; k<nCentBins_pt; k++)
	{
	  hJpsiCountVsPt[s][k] = new TH1F(Form("Sys%s_JpsiCountsVsPt_cent%s",sysName[s].Data(),cent_Title_pt[k]),";p_{T} [GeV/c];Counts",nbins,xbins);
	  for(Int_t i=0; i<nbins; i++)
	    {
	      for(Int_t j=0; j<2; j++)
		{
		  hInvMassVsPt[s][k][i][j] = (TH1F*)fsys->Get(Form("Sys%s_InvMassData_%s_Pt%1.0f-%1.0f_cent%s",sysName[s].Data(),typeName[j],xbins[i],xbins[i+1],cent_Title_pt[k]));
		}
	    }
	}
    }

  // get the Jpsi width from embedding
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"read");
  TH1F *hEmbJpsiWidth = (TH1F*)fscan->Get("SmearEmb_JpsiWidth_def");
  TCanvas *cFit[nCentBins_pt];
  TCanvas *cFitSys[nbins];
  for(Int_t i=0; i<nbins; i++)
    {
      cFitSys[i] = new TCanvas(Form("cFitSys_Pt%d",i), Form("cFitSys_Pt%1.0f-%1.0f",xbins[i],xbins[i+1]), 1100, 700);
      if(nSysUsed==12)  cFitSys[i]->Divide(4,3);
      if(nSysUsed==5)   cFitSys[i]->Divide(3,2);
      if(nSysUsed==27)  cFitSys[i]->Divide(4,3);
      if(nSysUsed==7)   cFitSys[i]->Divide(4,2);
    }

  TF1 *funcUL[nSysUsed][nCentBins_pt][nbins];
  for(int s=0; s<nSysUsed; s++)
    {
      for(int k=0; k<nCentBins_pt; k++)
	{
	  if(s==0)
	    {
	      cFit[k] = new TCanvas(Form("cFit_%s",cent_Title_pt[k]),Form("cFit_%s",cent_Title_pt[k]),1100,700);
	      cFit[k]->Divide(3,3);
	    }
	  for(Int_t i=0; i<nbins; i++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hInvMassVsPt[s][k][i][j]->SetMarkerStyle(20+j*4);
		  if(outName=="TpcTracking")
		    {
		      if(i<4) hInvMassVsPt[s][k][i][j]->Rebin(5);
		      else    hInvMassVsPt[s][k][i][j]->Rebin(10);
		    }
		  if(outName=="MuonPid") 
		    {
		      if(i>=4) hInvMassVsPt[s][k][i][j]->Rebin(2);
		    }
		  hInvMassVsPt[s][k][i][j]->GetXaxis()->SetRangeUser(2.6,4);
		}
	      double fit_min = 2.6, fit_max = 4;
	      if(k<4) 
		{
		  if(i==0) {fit_min = 2.8; fit_max = 3.6; }
		  if(i==1) {fit_min = 2.6; fit_max = 4.0; }
		}
	      int nPar = 4;
	      if(i>=4) nPar = 1;
	      TF1 *funcLS = 0x0;
	      if(nPar==4) funcLS = new TF1(Form("%s_func",hInvMassVsPt[s][k][i][1]->GetName()), "pol3", fit_min, fit_max);
	      if(nPar==2) funcLS = new TF1(Form("%s_func",hInvMassVsPt[s][k][i][1]->GetName()), "pol1", fit_min, fit_max);
	      if(nPar==1) funcLS = new TF1(Form("%s_func",hInvMassVsPt[s][k][i][1]->GetName()), "pol0", fit_min, fit_max);
	      hInvMassVsPt[s][k][i][1]->Fit(funcLS,"IR0Q");

	      if(nPar==4) funcUL[s][k][i] = new TF1(Form("%s_func",hInvMassVsPt[s][k][i][0]->GetName()), "gausn(0)+pol3(3)", fit_min, fit_max);
	      if(nPar==2) funcUL[s][k][i] = new TF1(Form("%s_func",hInvMassVsPt[s][k][i][0]->GetName()), "gausn(0)+pol1(3)", fit_min, fit_max);
	      if(nPar==1) funcUL[s][k][i] = new TF1(Form("%s_func",hInvMassVsPt[s][k][i][0]->GetName()), "gausn(0)+pol0(3)", fit_min, fit_max);
	      for(int ipar=0; ipar<3+nPar; ipar++)
		{
		  if(ipar==0) funcUL[s][k][i]->SetParameter(ipar, 100);
		  else if(ipar==1) funcUL[s][k][i]->FixParameter(ipar, 3.096);
		  else if(ipar==2) funcUL[s][k][i]->FixParameter(2, hEmbJpsiWidth->GetBinContent(i+1));
		  else             funcUL[s][k][i]->SetParameter(ipar, funcLS->GetParameter(ipar-3));
		}
	      hInvMassVsPt[s][k][i][0]->Add(hInvMassVsPt[s][k][i][1], -1);
	      hInvMassVsPt[s][k][i][0]->Fit(funcUL[s][k][i],"IR0Q");
	      for(int ipar=0; ipar<nPar; ipar++)
		{
		  funcLS->SetParameter(ipar, funcUL[s][k][i]->GetParameter(ipar+3));
		}
	      double bin_width = hInvMassVsPt[s][k][i][0]->GetBinWidth(1);
	      double nJpsi = funcUL[s][k][i]->GetParameter(0)/bin_width;
	      double nJpsi_err = funcUL[s][k][i]->GetParError(0)/bin_width;
	      hJpsiCountVsPt[s][k]->SetBinContent(i+1, nJpsi);
	      hJpsiCountVsPt[s][k]->SetBinError(i+1, nJpsi_err);
	      hInvMassVsPt[s][k][i][0]->SetMaximum(1.4*hInvMassVsPt[s][k][i][0]->GetMaximum());

	      if(s==0)
		{
		  if( (k==2 || k==3) && i>7) continue;
		  if( (k==4) && i>5) continue;
		  cFit[k]->cd(i+1);
		  hInvMassVsPt[s][k][i][0]->SetMaximum(1.3*hInvMassVsPt[s][k][i][0]->GetMaximum());
		  hInvMassVsPt[s][k][i][0]->Draw("PE");
		  funcLS->SetLineStyle(2);
		  funcLS->SetLineColor(4);
		  funcLS->Draw("sames");
		  funcUL[s][k][i]->SetLineStyle(2);
		  funcUL[s][k][i]->SetLineColor(2);
		  funcUL[s][k][i]->Draw("sames");
		  TPaveText *t1 = GetTitleText(Form("%s%%, %1.0f < p_{T} < %1.0f GeV/c",cent_Name_pt[k],xbins[i],xbins[i+1]),0.08);
		  t1->Draw();
		}

	      if(k==0)
		{
		  cFitSys[i]->cd(s+1);
		  if(s!=0) hInvMassVsPt[s][k][i][0]->SetMaximum(1.3*hInvMassVsPt[s][k][i][0]->GetMaximum());
		  hInvMassVsPt[s][k][i][0]->Draw("PE");
		  funcLS->SetLineStyle(2);
		  funcLS->SetLineColor(4);
		  funcLS->Draw("sames");
		  funcUL[s][k][i]->SetLineStyle(2);
		  funcUL[s][k][i]->SetLineColor(2);
		  funcUL[s][k][i]->Draw("sames");
		  TPaveText *t1 = GetTitleText(sysName[s].Data(),0.07);
		  t1->Draw();
		  if(s==0)
		    {
		      t1 = GetPaveText(0.65,0.75,0.65,0.85,0.06,62);
		      t1->AddText(Form("%s%%",cent_Name_pt[k]));
		      t1->AddText(Form("%1.0f < p_{T} < %1.0f GeV/c",xbins[i],xbins[i+1]));
		      t1->SetTextColor(2);
		      t1->Draw();
		    }
		  t1 = GetPaveText(0.65,0.75,0.4,0.5,0.06,62);
		  t1->AddText(Form("N=%1.0f#pm%1.0f",nJpsi,nJpsi_err));
		  t1->Draw();
		}
	    }
	}
    }

  if(savePlot)
    {
      for(int k=0; k<nCentBins_pt; k++)
	{
	  cFit[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_InvMassVsPt_cent%s.pdf",run_type,outName.Data(),cent_Title_pt[k]));
	}
      for(Int_t i=0; i<nbins; i++)
	{
	  cFitSys[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_InvMassVsSys_Pt%1.0f-%1.0f_cent%s.pdf",run_type,outName.Data(),xbins[i],xbins[i+1],cent_Title_pt[0]));
	}
    }

  TLegend *leg1[2];
  for(int i=0; i<2; i++)
    {
      leg1[i] = new TLegend(0.1+i*0.35, 0.45, 0.45+i*0.35, 0.7);
      leg1[i]->SetBorderSize(0);
      leg1[i]->SetFillColor(0);
      leg1[i]->SetTextSize(0.035);
    }

  TCanvas *cJpsiCountVsPt = new TCanvas("cJpsiCountVsPt", "cJpsiCountVsPt", 1100, 700);
  cJpsiCountVsPt->Divide(3,2);
  for(int k=0; k<nCentBins_pt; k++)
    {
      cJpsiCountVsPt->cd(k+2);
      gPad->SetLogy();
      for(int s=0; s<nSysUsed; s++)
	{
	  hJpsiCountVsPt[s][k]->SetMarkerStyle(20+s);
	  hJpsiCountVsPt[s][k]->SetMarkerColor(color[s]);
	  hJpsiCountVsPt[s][k]->SetMarkerSize(1.3);
	  hJpsiCountVsPt[s][k]->SetLineColor(color[s]);
	  if(k==2 || k==3) hJpsiCountVsPt[s][k]->GetXaxis()->SetRangeUser(0,10);
	  if(k==4) hJpsiCountVsPt[s][k]->GetXaxis()->SetRangeUser(0,6);
	  if(s==0) hJpsiCountVsPt[s][k]->Draw();
	  else     hJpsiCountVsPt[s][k]->Draw("sames");
	  if(k==0 && s<nSysMarker)
	    {
	      leg1[s/4]->AddEntry(hJpsiCountVsPt[s][k], sysName[s].Data(), "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name_pt[k]),0.06);
      t1->Draw();
    }
  cJpsiCountVsPt->cd(1);
  for(int i=0; i<2; i++)
    {
      leg1[i]->Draw();
    }
  TPaveText *system = GetPaveText(0.25,0.45,0.75,0.8,0.06);
  system->AddText(run_type);
  system->Draw();
  if(savePlot)
    cJpsiCountVsPt->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_JpsiCountVsPt.pdf",run_type,outName.Data()));

  // Get efficiency from embedding
  TCanvas *cJpsiEffVsPt = new TCanvas("cJpsiEffVsPt", "cJpsiEffVsPt", 1100, 700);
  cJpsiEffVsPt->Divide(3,2);
  TH1F *hJpsiEffVsPt[nSysUsed][nCentBins_pt];
  for(int k=0; k<nCentBins_pt; k++)
    {
      cJpsiEffVsPt->cd(k+2);
      gPad->SetLogy();
      for(int s=0; s<nSysUsed; s++)
	{
	  hJpsiEffVsPt[s][k] = (TH1F*)fsys->Get(Form("Sys%s_JpsiEffVsPtEmb_cent%s",sysName[s].Data(),cent_Title_pt[k]));
	  hJpsiEffVsPt[s][k]->SetMarkerStyle(20+s);
	  hJpsiEffVsPt[s][k]->SetMarkerColor(color[s]);
	  hJpsiEffVsPt[s][k]->SetLineColor(color[s]);
	  if(s==0) hJpsiEffVsPt[s][k]->Draw("PE");
	  else     hJpsiEffVsPt[s][k]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name_pt[k]),0.06);
      t1->Draw();
    }
  cJpsiEffVsPt->cd(1);
  for(int i=0; i<2; i++)
    {
      leg1[i]->Draw();
    }
  system->Draw();
  if(savePlot)
    cJpsiEffVsPt->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_JpsiEffVsPt.pdf",run_type,outName.Data()));

  // correct J/psi counts to derive uncertainty
  TCanvas *cJpsiSysVsPt = new TCanvas("cJpsiSysVsPt", "cJpsiSysVsPt", 1100, 700);
  cJpsiSysVsPt->Divide(3,2);
  TH1F *hJpsiCountCorrVsPt[nSysUsed][nCentBins_pt];
  TH1F *hJpsiSysVsPt[nSysUsed][nCentBins_pt];
  for(int k=0; k<nCentBins_pt; k++)
    {
      cJpsiSysVsPt->cd(k+2);
      for(int s=0; s<nSysUsed; s++)
	{
	  hJpsiCountCorrVsPt[s][k] = (TH1F*)hJpsiCountVsPt[s][k]->Clone(Form("Sys%s_JpsiCountCorrVsPt_cent%s",sysName[s].Data(),cent_Title_pt[k]));
	  hJpsiCountCorrVsPt[s][k]->Divide(hJpsiEffVsPt[s][k]);

	  hJpsiSysVsPt[s][k] = (TH1F*)hJpsiCountCorrVsPt[s][k]->Clone(Form("Sys%s_JpsiSysVsPt_cent%s",sysName[s].Data(),cent_Title_pt[k]));
	  for(int i=0; i<nbins; i++)
	    {
	      if(hJpsiCountCorrVsPt[0][k]->GetBinContent(i+1)==0) continue;
	      hJpsiSysVsPt[s][k]->SetBinContent(i+1, hJpsiSysVsPt[s][k]->GetBinContent(i+1)/hJpsiCountCorrVsPt[0][k]->GetBinContent(i+1));
	      if(s==0) hJpsiSysVsPt[s][k]->SetBinError(i+1, hJpsiCountCorrVsPt[0][k]->GetBinError(i+1)/hJpsiCountCorrVsPt[0][k]->GetBinContent(i+1));
	      else     hJpsiSysVsPt[s][k]->SetBinError(i+1, 0);
	    }
	  hJpsiSysVsPt[s][k]->GetYaxis()->SetRangeUser(0.7,1.3);
	  hJpsiSysVsPt[s][k]->SetTitle(";p_{T} (GeV/c);Sys. uncert.");
	  if(s==0) hJpsiSysVsPt[s][k]->Draw("PE");
	  else     hJpsiSysVsPt[s][k]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name_pt[k]),0.06);
      t1->Draw();
    }
  cJpsiSysVsPt->cd(1);
  for(int i=0; i<2; i++)
    {
      leg1[i]->Draw();
    }
  system->Draw();
  if(savePlot) cJpsiSysVsPt->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_JpsiSysVsPt.pdf",run_type,outName.Data()));

  // estimate final uncertainty
  TH1F *hEstSysVsPt[nCentBins_pt];
  TH1F *hFinalSysVsPt[nCentBins_pt];
  TH1F *hFinalSysVsPtLower[nCentBins_pt];
  for(int k=0; k<nCentBins_pt; k++)
    {
      hEstSysVsPt[k] = new TH1F(Form("EstSys_%sVsPt_cent%s",outName.Data(),cent_Title_pt[k]),"",nbins,xbins);
      hFinalSysVsPt[k] = new TH1F(Form("FinalSys_%sVsPt_cent%s",outName.Data(),cent_Title_pt[k]),"",nbins,xbins);
      hFinalSysVsPtLower[k] = new TH1F(Form("FinalSysLower_%sVsPt_cent%s",outName.Data(),cent_Title_pt[k]),"",nbins,xbins);
    }
  TH1F *htmp = new TH1F("htmp_sys_distribution","",100,0.5,1.5);
  TF1 *funcUncert = 0x0;
  if(uncerMode==0)
    {
      // maximum
      for(int k=0; k<nCentBins_pt; k++)
	{
	  for(int i=0; i<nbins; i++)
	    {
	      double max = 0;
	      for(int s=0; s<nSysUsed; s++)
		{
		  double err = fabs(hJpsiSysVsPt[s][k]->GetBinContent(i+1)-1);
		  if(err>max) max = err;
		}
	      hEstSysVsPt[k]->SetBinContent(i+1, 1+max);
	      if(k==0)
		{
		  hFinalSysVsPt[k]->SetBinContent(i+1, max+1);
		  if(i>=7) hFinalSysVsPt[k]->SetBinContent(i+1, hFinalSysVsPt[k]->GetBinContent(7));
		}
	      else
		{
		  hFinalSysVsPt[k]->SetBinContent(i+1, hFinalSysVsPt[0]->GetBinContent(i+1));
		}
	      hFinalSysVsPtLower[k]->SetBinContent(i+1, 2-hFinalSysVsPt[k]->GetBinContent(i+1));
	    }
	}
      if(outName=="TpcTracking")
	{
	  TH1F *hFitUncert = (TH1F*)hEstSysVsPt[0]->Clone("hfit_tmp");
	  for(int i=0; i<nbins; i++)
	    {
	      hFitUncert->SetBinError(i+1, hJpsiSysVsPt[0][0]->GetBinError(i+1));
	    }
	  funcUncert = new TF1("funcUncert","pol0",0,15);
	  hFitUncert->Fit(funcUncert,"RQ0");
	  draw1D(hFitUncert,"Fit to estimate uncertainty");
	  funcUncert->Draw("sames");
	  for(int k=0; k<nCentBins_pt; k++)
	    {
	      for(int i=0; i<nbins; i++)
		{
		  hFinalSysVsPt[k]->SetBinContent(i+1, funcUncert->GetParameter(0));
		  hFinalSysVsPtLower[k]->SetBinContent(i+1, 2-hFinalSysVsPt[k]->GetBinContent(i+1));
		}
	    }	
	}
    }
  else if(uncerMode==1)
    {
      // RMS
      for(int k=0; k<nCentBins_pt; k++)
	{
	  for(int i=0; i<nbins; i++)
	    {
	      htmp->Reset();
	      for(int s=0; s<nSysUsed; s++)
		{
		  htmp->Fill(hJpsiSysVsPt[s][k]->GetBinContent(i+1));
		}
	      hEstSysVsPt[k]->SetBinContent(i+1, htmp->GetRMS()+1);
	      if(k==0)
		{
		  hFinalSysVsPt[k]->SetBinContent(i+1, htmp->GetRMS()+1);
		  if(i>=7) hFinalSysVsPt[k]->SetBinContent(i+1, hFinalSysVsPt[k]->GetBinContent(7));
		}
	      else
		{
		  hFinalSysVsPt[k]->SetBinContent(i+1, hFinalSysVsPt[0]->GetBinContent(i+1));
		}
	      hFinalSysVsPtLower[k]->SetBinContent(i+1, 2-hFinalSysVsPt[k]->GetBinContent(i+1));
	    }
	}
    }
  for(int k=0; k<nCentBins_pt; k++)
    {
      cJpsiSysVsPt->cd(k+2);
      if(k==0) hEstSysVsPt[k]->SetLineColor(2);
      else     hEstSysVsPt[k]->SetLineColor(4);
      hEstSysVsPt[k]->SetLineStyle(2);
      hEstSysVsPt[k]->SetLineWidth(2);
      hEstSysVsPt[k]->Draw("sames");

      hFinalSysVsPt[k]->SetLineColor(2);
      hFinalSysVsPt[k]->SetLineWidth(2);
      hFinalSysVsPt[k]->Draw("sames");
      hFinalSysVsPtLower[k]->SetLineColor(2);
      hFinalSysVsPtLower[k]->SetLineWidth(2);
      hFinalSysVsPtLower[k]->Draw("sames");
    }
  cJpsiSysVsPt->cd(1);
  TLegend *leg = new TLegend(0.1, 0.2, 0.45, 0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  if(outName=="TpcTracking") leg->AddEntry(hFinalSysVsPtLower[0],"0-80% fitted");
  if(outName=="MuonPid") leg->AddEntry(hFinalSysVsPtLower[0],"0-80% re-assigned");
  if(uncerMode==0)
    {
      leg->AddEntry(hEstSysVsPt[0],"Max: 0-80% ");
      leg->AddEntry(hEstSysVsPt[1],"Max: individual cent. bin");
    }
  if(uncerMode==1)
    {
      leg->AddEntry(hEstSysVsPt[0],"RMS: 0-80% ");
      leg->AddEntry(hEstSysVsPt[1],"RMS: individual cent. bin");
    }
  leg->Draw();
  if(savePlot) cJpsiSysVsPt->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_FinalSysVsPt.pdf",run_type,outName.Data()));
  if(saveHisto)
    {
      fsys->cd();
      for(int k=0; k<nCentBins_pt; k++)
	{
	  hFinalSysVsPt[k]->Write("",TObject::kOverwrite);
	}
    }

  //----------------------------------------------
  // centrality dependence
  //----------------------------------------------
  const int nCentBins_npart[nPtBins_npart] = {9,8};
  const char *cent_Name_npart[nPtBins_npart*nCentBins_npart[0]] = {"0-80","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80", 
								   "0-80","0-10","10-20","20-30","30-40","40-50","50-60","60-80","60-70"};
  const char *cent_Title_npart[nPtBins_npart*nCentBins_npart[0]] = {"0080","0010","1020","2030","3040","4050","5060","6070","7080",
								    "0080","0010","1020","2030","3040","4050","5060","6080","6070"};
  TH1F *hInvMassVsCent[nSysUsed][nPtBins_npart][nCentBins_npart[0]][2];
  TH1F *hJpsiCountVsCent[nSysUsed][nPtBins_npart];
  for(int s=0; s<nSysUsed; s++)
    {
      for(Int_t i=0; i<nPtBins_npart; i++)
	{
	  hJpsiCountVsCent[s][i] = new TH1F(Form("Sys%s_JpsiCountsVsCent_Pt%s",sysName[s].Data(),pt_Name_npart[i]),";;Counts",nCentBins_npart[i],0,nCentBins_npart[i]);
	  for(int k=1; k<nCentBins_npart[i]; k++)
	    {
	      for(Int_t j=0; j<2; j++)
		{
		  hInvMassVsCent[s][i][k][j] = (TH1F*)fsys->Get(Form("Sys%s_InvMassData_%s_Pt%s_cent%s",sysName[s].Data(),typeName[j],pt_Name_npart[i],cent_Title_npart[k]));
		  if(i==1 && k==nCentBins_npart[i]-1)
		    {
		      TH1F *h1tmp = (TH1F*)fsys->Get(Form("Sys%s_InvMassData_%s_Pt%s_cent%s",sysName[s].Data(),typeName[j],pt_Name_npart[i],cent_Title_npart[k+1]));
		      hInvMassVsCent[s][i][k][j]->Add(h1tmp);
		    }
		  if(k==1) hInvMassVsCent[s][i][0][j] = (TH1F*)hInvMassVsCent[s][i][k][j]->Clone(Form("Sys%s_InvMassData_%s_Pt%s_cent%s",sysName[s].Data(),typeName[j],pt_Name_npart[i],cent_Title_npart[0]));
		  else     hInvMassVsCent[s][i][0][j]->Add(hInvMassVsCent[s][i][k][j]);
		}
	    }
	}
    }


  // get the Jpsi width from embedding
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"read");
  TH1F *hEmbJpsiWidth = (TH1F*)fscan->Get("SmearEmb_JpsiWidthIntegr_def");
  TCanvas *cFitCent[nPtBins_npart];
  TCanvas *cFitSysCent = new TCanvas(Form("cFitSysCent"), Form("cFitSysCent"), 1100, 700);
  if(nSysUsed==12)  cFitSysCent->Divide(4,3);
  if(nSysUsed==5)   cFitSysCent->Divide(3,2);
  if(nSysUsed==27)  cFitSysCent->Divide(4,3);
  if(nSysUsed==7)   cFitSysCent->Divide(4,2);

  TF1 *funcCentUL[nSysUsed][nPtBins_npart][nCentBins_npart[0]];
  for(int s=0; s<nSysUsed; s++)
    {
      for(Int_t i=0; i<nPtBins_npart; i++)
	{
	  if(s==0)
	    {
	      cFitCent[i] = new TCanvas(Form("cFit_%s",pt_Name_npart[i]),Form("cFit_%s",pt_Name_npart[i]),1100,700);
	      cFitCent[i]->Divide(3,3);
	    }
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hInvMassVsCent[s][i][k][j]->SetMarkerStyle(20+j*4);
		  if(outName=="TpcTracking")
		    {
		      if(i==0) hInvMassVsCent[s][i][k][j]->Rebin(5);
		      else     hInvMassVsCent[s][i][k][j]->Rebin(10);
		    }
		  if(outName=="MuonPid") 
		    {
		      if(i!=-0) hInvMassVsCent[s][i][k][j]->Rebin(2);
		    }
		  hInvMassVsCent[s][i][k][j]->GetXaxis()->SetRangeUser(2.6,4);
		}
	      double fit_min = 2.6, fit_max = 4;
	      if(k<4) 
		{
		  if(i==0) {fit_min = 2.6; fit_max = 4; }
		}

	      int nPar = 4;
	      if(i==1) nPar = 1;
	      TF1 *funcLS = 0x0;
	      if(nPar==4) funcLS = new TF1(Form("%s_func",hInvMassVsCent[s][i][k][1]->GetName()), "pol3", fit_min, fit_max);
	      if(nPar==2) funcLS = new TF1(Form("%s_func",hInvMassVsCent[s][i][k][1]->GetName()), "pol1", fit_min, fit_max);
	      if(nPar==1) funcLS = new TF1(Form("%s_func",hInvMassVsCent[s][i][k][1]->GetName()), "pol0", fit_min, fit_max);
	      hInvMassVsCent[s][i][k][1]->Fit(funcLS,"IR0Q");

	      if(nPar==4) funcCentUL[s][i][k] = new TF1(Form("%s_func",hInvMassVsCent[s][i][k][0]->GetName()), "gausn(0)+pol3(3)", fit_min, fit_max);
	      if(nPar==2) funcCentUL[s][i][k] = new TF1(Form("%s_func",hInvMassVsCent[s][i][k][0]->GetName()), "gausn(0)+pol1(3)", fit_min, fit_max);
	      if(nPar==1) funcCentUL[s][i][k] = new TF1(Form("%s_func",hInvMassVsCent[s][i][k][0]->GetName()), "gausn(0)+pol0(3)", fit_min, fit_max);
	      for(int ipar=0; ipar<3+nPar; ipar++)
		{
		  if(ipar==0) funcCentUL[s][i][k]->SetParameter(ipar, 100);
		  else if(ipar==1) funcCentUL[s][i][k]->FixParameter(ipar, 3.096);
		  else if(ipar==2) funcCentUL[s][i][k]->FixParameter(2, hEmbJpsiWidth->GetBinContent(i+1));
		  else             funcCentUL[s][i][k]->SetParameter(ipar, funcLS->GetParameter(ipar-3));
		}
	      hInvMassVsCent[s][i][k][0]->Add(hInvMassVsCent[s][i][k][1], -1);
	      hInvMassVsCent[s][i][k][0]->Fit(funcCentUL[s][i][k],"IR0Q");
	      for(int ipar=0; ipar<nPar; ipar++)
		{
		  funcLS->SetParameter(ipar, funcCentUL[s][i][k]->GetParameter(ipar+3));
		}
	      double bin_width = hInvMassVsCent[s][i][k][0]->GetBinWidth(1);
	      double nJpsi = funcCentUL[s][i][k]->GetParameter(0)/bin_width;
	      double nJpsi_err = funcCentUL[s][i][k]->GetParError(0)/bin_width;
	      hJpsiCountVsCent[s][i]->SetBinContent(k+1, nJpsi);
	      hJpsiCountVsCent[s][i]->SetBinError(k+1, nJpsi_err);
	      hJpsiCountVsCent[s][i]->GetXaxis()->SetBinLabel(k+1, Form("%s%%",cent_Name_npart[k+i*nCentBins_npart[0]]));
	      hInvMassVsCent[s][i][k][0]->SetMaximum(1.4*hInvMassVsCent[s][i][k][0]->GetMaximum());

	      if(s==0)
		{
		  cFitCent[i]->cd(k+1);
		  hInvMassVsCent[s][i][k][0]->SetMaximum(1.3*hInvMassVsCent[s][i][k][0]->GetMaximum());
		  hInvMassVsCent[s][i][k][0]->Draw("PE");
		  funcLS->SetLineStyle(2);
		  funcLS->SetLineColor(4);
		  funcLS->Draw("sames");
		  funcCentUL[s][i][k]->SetLineStyle(2);
		  funcCentUL[s][i][k]->SetLineColor(2);
		  funcCentUL[s][i][k]->Draw("sames");
		  TPaveText *t1 = GetTitleText(Form("%s%%, p_{T} > %1.0f GeV/c",cent_Name_npart[k+i*nCentBins_npart[0]],ptBins_low_npart[i]),0.08);
		  t1->Draw();
		}

	      if(i==0 && k==1)
		{
		  cFitSysCent->cd(s+1);
		  hInvMassVsCent[s][i][k][0]->Draw("PE");
		  funcLS->SetLineStyle(2);
		  funcLS->SetLineColor(4);
		  funcLS->Draw("sames");
		  funcCentUL[s][i][k]->SetLineStyle(2);
		  funcCentUL[s][i][k]->SetLineColor(2);
		  funcCentUL[s][i][k]->Draw("sames");
		  TPaveText *t1 = GetTitleText(sysName[s].Data(),0.07);
		  t1->Draw();
		  if(s==0)
		    {
		      t1 = GetPaveText(0.65,0.75,0.65,0.85,0.06,62);
		      t1->AddText(Form("%s%%",cent_Name_npart[k+i*nCentBins_npart[0]]));
		      t1->AddText(Form("p_{T} > %1.0f GeV/c",ptBins_low_npart[i]));
		      t1->SetTextColor(2);
		      t1->Draw();
		    }
		  t1 = GetPaveText(0.65,0.75,0.4,0.5,0.06,62);
		  t1->AddText(Form("N=%1.0f#pm%1.0f",nJpsi,nJpsi_err));
		  t1->Draw();
		}
	    }
	}
    }
  if(savePlot)
    {
      for(Int_t i=0; i<nPtBins_npart; i++)
	{
	  cFitCent[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_InvMassVsCent_Pt%s.pdf",run_type,outName.Data(),pt_Name_npart[i]));
	}
    }

  const int nLeg = 2;
  TLegend *leg2[nLeg];
  for(int i=0; i<nLeg; i++)
    {
      leg2[i] = new TLegend(0.3+i*0.2, 0.15, 0.5+i*0.2, 0.3);
      leg2[i]->SetBorderSize(0);
      leg2[i]->SetFillColor(0);
      leg2[i]->SetTextSize(0.02);
    }

  TCanvas *cJpsiCountVsCent = new TCanvas("cJpsiCountVsCent", "cJpsiCountVsCent", 1200, 500);
  cJpsiCountVsCent->Divide(2,1);
  for(Int_t i=0; i<nPtBins_npart; i++)
    {
      cJpsiCountVsCent->cd(i+1);
      gPad->SetLogy();
      for(int s=0; s<nSysUsed; s++)
	{
	  hJpsiCountVsCent[s][i]->SetMarkerStyle(s+20);
	  hJpsiCountVsCent[s][i]->SetMarkerColor(color[s]);
	  hJpsiCountVsCent[s][i]->SetMarkerSize(1.3);
	  hJpsiCountVsCent[s][i]->SetLineColor(color[s]);
	  if(s==0) hJpsiCountVsCent[s][i]->Draw();
	  else     hJpsiCountVsCent[s][i]->Draw("sames");
	  if(i==0 && s<nSysMarker)
	    {
	      leg2[s/4]->AddEntry(hJpsiCountVsCent[s][i], sysName[s].Data(), "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s: p_{T} > %1.0f GeV/c",run_type,ptBins_low_npart[i]),0.04);
      t1->Draw();
    }
  if(savePlot)
    cJpsiCountVsCent->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_JpsiCountVsCent.pdf",run_type,outName.Data()));

  // Get efficiency from embedding
  TCanvas *cJpsiEffVsCent = new TCanvas("cJpsiEffVsCent", "cJpsiEffVsCent", 1200, 500);
  cJpsiEffVsCent->Divide(2,1);
  TH1F *hJpsiEffVsCent[nSysUsed][nPtBins_npart];
  for(Int_t i=0; i<nPtBins_npart; i++)
    {
      cJpsiEffVsCent->cd(i+1);
      for(int s=0; s<nSysUsed; s++)
	{
	  hJpsiEffVsCent[s][i] = (TH1F*)fsys->Get(Form("Sys%s_JpsiEffVsCentEmb_Pt%1.0f",sysName[s].Data(),ptBins_low_npart[i]));
	  hJpsiEffVsCent[s][i]->SetMarkerStyle(s+20);
	  hJpsiEffVsCent[s][i]->SetMarkerColor(color[s]);
	  hJpsiEffVsCent[s][i]->SetLineColor(color[s]);
	  if(outName=="TpcTracking") hJpsiEffVsCent[s][i]->GetYaxis()->SetRangeUser(0.06,0.14);
	  if(outName=="MuonPid")     hJpsiEffVsCent[s][i]->GetYaxis()->SetRangeUser(0.005,0.03);
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      hJpsiEffVsCent[s][i]->GetXaxis()->SetBinLabel(k+1, Form("%s%%",cent_Name_npart[k+i*nCentBins_npart[0]]));
	    }
	  if(s==0) hJpsiEffVsCent[s][i]->Draw("PE");
	  else     hJpsiEffVsCent[s][i]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("%s: p_{T} > %1.0f GeV/c",run_type,ptBins_low_npart[i]),0.04);
      t1->Draw();
    }
  if(savePlot)
    cJpsiEffVsCent->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_JpsiEffVsCent.pdf",run_type,outName.Data()));

  // correct J/psi counts to derive uncertainty
  TCanvas *cJpsiSysVsCent = new TCanvas("cJpsiSysVsCent", "cJpsiSysVsCent", 1200, 500);
  cJpsiSysVsCent->Divide(2,1);
  TH1F *hJpsiCountCorrVsCent[nSysUsed][nPtBins_npart];
  TH1F *hJpsiSysVsCent[nSysUsed][nPtBins_npart];
  for(Int_t i=0; i<nPtBins_npart; i++)
    {
      cJpsiSysVsCent->cd(i+1);
      for(int s=0; s<nSysUsed; s++)
	{
	  hJpsiCountCorrVsCent[s][i] = (TH1F*)hJpsiCountVsCent[s][i]->Clone(Form("Sys%s_JpsiCountCorrVsCent_Pt%1.0f",sysName[s].Data(),ptBins_low_npart[i]));
	  hJpsiCountCorrVsCent[s][i]->Divide(hJpsiEffVsCent[s][i]);

	  hJpsiSysVsCent[s][i] = (TH1F*)hJpsiCountCorrVsCent[s][i]->Clone(Form("Sys%s_JpsiSysVsCent_Pt%1.0f",sysName[s].Data(),ptBins_low_npart[i]));
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      if(hJpsiCountCorrVsCent[0][i]->GetBinContent(k+1)==0) continue;
	      hJpsiSysVsCent[s][i]->SetBinContent(k+1, hJpsiSysVsCent[s][i]->GetBinContent(k+1)/hJpsiCountCorrVsCent[0][i]->GetBinContent(k+1));
	      if(s==0) hJpsiSysVsCent[s][i]->SetBinError(k+1, hJpsiCountCorrVsCent[0][i]->GetBinError(k+1)/hJpsiCountCorrVsCent[0][i]->GetBinContent(k+1));
	      else     hJpsiSysVsCent[s][i]->SetBinError(k+1, 0);
	    }
	  hJpsiSysVsCent[s][i]->GetYaxis()->SetRangeUser(0.7,1.3);
	  hJpsiSysVsCent[s][i]->SetTitle(";;Sys. uncert.");
	  if(s==0) hJpsiSysVsCent[s][i]->Draw("PE");
	  else     hJpsiSysVsCent[s][i]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("%s: uncertainty for %s (p_{T} > %1.0f GeV/c)",run_type,outName.Data(),ptBins_low_npart[i]),0.04);
      t1->Draw();
    }
  cJpsiSysVsCent->cd(1);
  for(int i=0; i<nLeg; i++)
    {
      leg2[i]->Draw();
    }
  if(savePlot)
    cJpsiSysVsCent->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_JpsiSysVsCent.pdf",run_type,outName.Data()));

  // estimate final uncertainty
  TH1F *hEstSysVsCent[nPtBins_npart];
  TH1F *hFinalSysVsCent[nPtBins_npart];
  TH1F *hFinalSysVsCentLower[nPtBins_npart];
  for(int i=0; i<nPtBins_npart; i++)
    {
      hEstSysVsCent[i] = new TH1F(Form("EstSys_%sVsCent_Pt%s",outName.Data(),pt_Name_npart[i]),"",nCentBins_npart[i],0,nCentBins_npart[i]);
      hFinalSysVsCent[i] = new TH1F(Form("FinalSys_%sVsCent_Pt%s",outName.Data(),pt_Name_npart[i]),"",nCentBins_npart[i],0,nCentBins_npart[i]);
      hFinalSysVsCentLower[i] = new TH1F(Form("FinalSysLower_%sVsCent_Pt%s",outName.Data(),pt_Name_npart[i]),"",nCentBins_npart[i],0,nCentBins_npart[i]);
    }

  if(uncerMode==0)
    {
      // maximum
      for(int i=0; i<nPtBins_npart; i++)
	{
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      double max = 0;
	      for(int s=0; s<nSysUsed; s++)
		{
		  double err = fabs(hJpsiSysVsCent[s][i]->GetBinContent(k+1)-1);
		  if(err>max) max = err;
		}
	      hEstSysVsCent[i]->SetBinContent(k+1, 1+max);
	      hFinalSysVsCent[i]->SetBinContent(k+1, hEstSysVsCent[i]->GetBinContent(1));
	      hFinalSysVsCentLower[i]->SetBinContent(k+1, 2-hFinalSysVsCent[i]->GetBinContent(k+1));
	    }
	}
    }
  else if(uncerMode==1)
    {
      // RMS
      for(int i=0; i<nPtBins_npart; i++)
	{
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      htmp->Reset();
	      for(int s=0; s<nSysUsed; s++)
		{
		  htmp->Fill(hJpsiSysVsCent[s][i]->GetBinContent(k+1));
		}
	      hEstSysVsCent[i]->SetBinContent(k+1, htmp->GetRMS()+1);
	      hFinalSysVsCent[i]->SetBinContent(k+1, htmp->GetRMS()+1);
	      hFinalSysVsCentLower[i]->SetBinContent(k+1, 2-hFinalSysVsCent[i]->GetBinContent(k+1));
	    }
	}
    }

  if(outName=="TpcTracking")
    {
      for(int i=0; i<nPtBins_npart; i++)
	{
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      hFinalSysVsCent[i]->SetBinContent(k+1, funcUncert->GetParameter(0));
	      hFinalSysVsCentLower[i]->SetBinContent(k+1, 2-hFinalSysVsCent[i]->GetBinContent(k+1));
	    }
	}
    }
  for(Int_t i=0; i<nPtBins_npart; i++)
    {
      cJpsiSysVsCent->cd(i+1);
      hEstSysVsCent[i]->SetLineColor(4);
      hEstSysVsCent[i]->SetLineStyle(2);
      hEstSysVsCent[i]->SetLineWidth(2);
      hEstSysVsCent[i]->Draw("sames");

      hFinalSysVsCent[i]->SetLineColor(2);
      hFinalSysVsCent[i]->SetLineWidth(2);
      hFinalSysVsCent[i]->Draw("sames");
      hFinalSysVsCentLower[i]->SetLineColor(2);
      hFinalSysVsCentLower[i]->SetLineWidth(2);
      hFinalSysVsCentLower[i]->Draw("sames");
    }
  cJpsiSysVsCent->cd(1);
  TLegend *leg = new TLegend(0.5, 0.7, 0.65, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hFinalSysVsCent[0],"0-80%");
  leg->AddEntry(hEstSysVsCent[0],"Individual bins");
  leg->Draw();
  if(savePlot) cJpsiSysVsCent->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_FinalSysVsCent.pdf",run_type,outName.Data()));
  if(saveHisto)
    {
      fsys->cd();
      for(Int_t i=0; i<nPtBins_npart; i++)
	{
	  hFinalSysVsCent[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void makeHistos(const int type = 1, const int mode = 1)
{
  // re-assign global constants
  if(mode==0)
    {
      //----------------------------------------------
      // pT dependence
      const int nPtBins         = nPtBins_pt;
      const double* ptBinsLow  = ptBins_low_pt;
      const double* ptBinsHigh = ptBins_high_pt;
      const char** ptName      = pt_Name_pt;
      const int nCentBins       = nCentBins_pt; 
      const int* centBins_low   = centBins_low_pt;
      const int* centBins_high  = centBins_high_pt;
      const char** cent_Name    = cent_Name_pt;
      const char** cent_Title   = cent_Title_pt;
    }
  else if(mode==1)
    {
      //----------------------------------------------
      // centrality dependence
      const int nPtBins         = nPtBins_npart;
      const double* ptBinsLow  = ptBins_low_npart;
      const double* ptBinsHigh = ptBins_high_npart;
      const char** ptName      = pt_Name_npart;
      const int nCentBins       = nCentBins_npart[0]; 
      const int* centBins_low   = centBins_low_npart;
      const int* centBins_high  = centBins_high_npart;
      const char** cent_Name    = cent_Name_npart;
      const char** cent_Title   = cent_Title_npart;
    }

  if(type==0)
    {
      //Data
      TFile *fdata = TFile::Open(Form("output/Sys.Run14_AuAu200.jpsi.%s.root",outName.Data()),"read");
      TH1F *hInvMassData[nSys][nCentBins][nPtBins][2];
      for(int s=0; s<nSys; s++)
	{
	  printf("+++ Process %s +++\n",sysName[s].Data());
	  for(Int_t j=0; j<2; j++)
	    {
	      THnSparseF *hnInvMass = (THnSparseF*)fdata->Get(Form("mhJpsi%s_%s",typeName[j],sysName[s].Data()));
	      hnInvMass->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
	      hnInvMass->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
	      for(Int_t i=0; i<nPtBins; i++) // pt bins
		{
		  hnInvMass->GetAxis(1)->SetRangeUser(ptBinsLow[i]+0.01,ptBinsHigh[i]-0.01);
		  for(int k=0; k<nCentBins; k++) // centrality bins
		    {
		      hnInvMass->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
		      hInvMassData[s][k][i][j] = (TH1F*)hnInvMass->Projection(0);
		      hInvMassData[s][k][i][j]->SetName(Form("Sys%s_InvMassData_%s_Pt%1.0f-%1.0f_cent%s",sysName[s].Data(),typeName[j],ptBinsLow[i],ptBinsHigh[i],cent_Title[k]));
		      hInvMassData[s][k][i][j]->Sumw2();
		      hnInvMass->GetAxis(4)->SetRange(0,-1);
		    }
		  hnInvMass->GetAxis(1)->SetRange(0,-1);
		}
	      hnInvMass->GetAxis(2)->SetRange(0,-1);
	      hnInvMass->GetAxis(3)->SetRange(0,-1);
	    }
	}
      TFile *fout = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.%s.root",outName.Data()),"update");
      for(int s=0; s<nSys; s++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	      for(Int_t i=0; i<nPtBins; i++)
		{
		  for(int k=0; k<nCentBins; k++)
		    {
		      hInvMassData[s][k][i][j]->SetTitle("");
		      hInvMassData[s][k][i][j]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
    }
  else if(type==1)
    {
      //Embedding
      TFile *fEmbed = TFile::Open(Form("output/Sys.Run14_AuAu200.Embed.jpsi.%s.root",outName.Data()),"read");
      TString embed_name[2] = {"Mc","Tpc"};
      if(outName=="MuonPid") embed_name[1] = "Pid";
      TH1F *hTrkVsPtEmb[nSys][nCentBins][2];
      for(int s=0; s<nSys; s++)
	{
	  printf("+++ Process %s +++\n",sysName[s].Data());

	  for(int i=0; i<2; i++)
	    {
	      THnSparseF *hnTrkEffEmb =  (THnSparseF*)fEmbed->Get(Form("mhTrack%s_%s",embed_name[i].Data(),sysName[s].Data()));
	      hnTrkEffEmb->GetAxis(1)->SetRangeUser(-0.5,0.5); // cut on track eta
	      for(int k=0; k<nCentBins; k++)
		{
		  hnTrkEffEmb->GetAxis(2)->SetRange(centBins_low[k],centBins_high[k]);
		  hTrkVsPtEmb[s][k][i] = (TH1F*)hnTrkEffEmb->Projection(0);
		  hTrkVsPtEmb[s][k][i]->SetName(Form("Sys%s_TrkEffVsPtEmb_%s_cent%s",sysName[s].Data(),embed_name[i].Data(),cent_Title[k]));
		  hTrkVsPtEmb[s][k][i]->SetTitle("");
		  hTrkVsPtEmb[s][k][i]->SetBinContent(hTrkVsPtEmb[s][k][i]->GetNbinsX()+1,0); // reset overflow bin
		  hTrkVsPtEmb[s][k][i]->Sumw2();
		  hnTrkEffEmb->GetAxis(2)->SetRange(0,-1); 
		}
	      hnTrkEffEmb->GetAxis(1)->SetRange(0,-1);
	    }
	}
      TFile *fout = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.%s.root",outName.Data()),"update");
      for(int s=0; s<nSys; s++)
	{
	  for(int i=0; i<2; i++)
	    {
	      for(int k=0; k<nCentBins; k++)
		{
		  hTrkVsPtEmb[s][k][i]->Write("",TObject::kOverwrite);
		}
	    }
	}
      fout->Close();
    }
}


//================================================
void mergeHistosData()
{
  // data
  TH1F *hEventStat[nSys];
  THnSparseF *hnJpsiUL[nSys];
  THnSparseF *hnJpsiLS[nSys];
  for(int i=0; i<nSys; i++)
    {
      printf("[i] Process %s\n",sysName[i].Data());
      TFile *fsys = TFile::Open(Form("output/Sys.Run14_AuAu200.jpsi.%s.root",sysName[i].Data()),"read");
      hEventStat[i] = (TH1F*)fsys->Get("hEventStat");
      hEventStat[i]->SetName(Form("mhEventStat_%s",sysName[i].Data()));

      hnJpsiUL[i] = (THnSparseF*)fsys->Get("mhJpsiInfo_di_mu");
      hnJpsiUL[i]->SetName(Form("mhJpsiUL_%s",sysName[i].Data()));

      hnJpsiLS[i] = (THnSparseF*)fsys->Get("mhBkgLSPos_di_mu");
      THnSparseF *hnTmp = (THnSparseF*)fsys->Get("mhBkgLSNeg_di_mu");
      hnJpsiLS[i]->Add(hnTmp);
      hnJpsiLS[i]->SetName(Form("mhJpsiLS_%s",sysName[i].Data()));
    }
  TFile *fout = TFile::Open(Form("output/Sys.Run14_AuAu200.jpsi.%s.root",outName.Data()),"recreate");
  for(int i=0; i<nSys; i++)
    {
      hEventStat[i]->Write();
      hnJpsiUL[i]->Write();
      hnJpsiLS[i]->Write();
    }
  fout->Close();
}


//================================================
void mergeHistosEmb()
{
  // embedding
  TH1F *hEventStatEmb[nSys];
  THnSparseF *hnTrackMc[nSys];
  THnSparseF *hnTrackTpc[nSys];
  THnSparseF *hnTrackPid[nSys];
  for(int i=0; i<nSys; i++)
    {
      printf("[i] Process %s\n",sysName[i].Data());
      TFile *fsys = TFile::Open(Form("output/Sys.Run14_AuAu200.Embed.Jpsi.%s.root",sysName[i].Data()),"read");
      hEventStatEmb[i] = (TH1F*)fsys->Get("hEventStat");
      hEventStatEmb[i]->SetName(Form("hEventStat_%s",sysName[i].Data()));

      hnTrackMc[i] = (THnSparseF*)fsys->Get("mhMcTrkPtEff_MC_di_mu");
      hnTrackMc[i]->SetName(Form("mhTrackMc_%s",sysName[i].Data()));

      hnTrackTpc[i] = (THnSparseF*)fsys->Get("mhMcTrkPtEff_Tpc_di_mu");
      hnTrackTpc[i]->SetName(Form("mhTrackTpc_%s",sysName[i].Data()));

      hnTrackPid[i] = (THnSparseF*)fsys->Get("mhMcTrkPtEff_MuonPid_di_mu");
      hnTrackPid[i]->SetName(Form("mhTrackPid_%s",sysName[i].Data()));
    }
  TFile *fout = TFile::Open(Form("output/Sys.Run14_AuAu200.Embed.jpsi.%s.root",outName.Data()),"recreate");
  for(int i=0; i<nSys; i++)
    {
      hEventStatEmb[i]->Write();
      hnTrackMc[i]->Write();
      hnTrackTpc[i]->Write();
      hnTrackPid[i]->Write();
    }
  fout->Close();
}
