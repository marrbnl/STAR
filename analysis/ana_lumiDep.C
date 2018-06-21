const int year = YEAR;
const int nCentBins = nCentBins_pt; 
const char** cent_Name = cent_Name_pt;
const char** cent_Title = cent_Title_pt;
const int* centBins_low = centBins_low_pt;
const int* centBins_high = centBins_high_pt;
const TString legName[5] = {"All","prod","prod_low","prod_mid","prod_high"};
//================================================
void ana_lumiDep()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //jpsiWidth();

  //vzDiff();
  //makeVzDiff();

  //makeTrigEff();
  anaTrigEff();
  //MtdTrigEffHadron();

  //muonPID();

  //sigSys();
}


//================================================
void sigSys(int savePlot = 1, int saveHisto = 1)
{
 // re-assign global constants
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;


  // evaluate uncertainty
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  xbins[0] = 0.15;

  TH1F *hSys[nCentBins][gNTrgSetup];
  TH1F *hMax[nCentBins][gNTrgSetup][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hSys[k][j]    = new TH1F(Form("Sys_signalExt_cent%s%s",cent_Title[k],gTrgSetupName[j]),Form("Systematic uncertainty for signal extraction (%s%%)",cent_Name[k]),nbins,xbins);
	  hMax[k][j][0] = new TH1F(Form("Sys_max_%s%s",cent_Title[k],gTrgSetupName[j]),"",nbins,xbins);
	  hMax[k][j][1] = new TH1F(Form("Sys_min_%s%s",cent_Title[k],gTrgSetupName[j]),"",nbins,xbins);
	}
    }

  const char *sys_name[13]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc1","_BkgFunc2","_LargeFit","_SmallFit","_SigFunc","_FixSigDown","_FixSigUp","_LineShape"};
  const char *sys_leg[13] = {"Default","Larger bkg norm range","Smaller bkg norm range","Fit ME/SE w/ pol0","Binning","Res. bkg order-1","Res. bkg order+1","Larger sig fit range","Smaller sig fit range","Crystal-ball","Fix sig. down","Fix sig. up","line shape"};
  const int nSys = 9;

  TFile *fin  = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config),"read");
  TFile *fsys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiYield.root",run_type),"read");

  TH1F *hSignal[nCentBins][gNTrgSetup][nSys];
  double max[nCentBins][gNTrgSetup][nPtBins];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      max[k][j][i] = 0;
	    }
	}
    }
  const char* method = "Fit";
  const int color_sys[nSys] = {1, 2, 3, 4, 6, 7, 8, 1, 2};

  for(int j=0; j<gNTrgSetup; j++)
    {
      TCanvas *c = new TCanvas(Form("Sys_signalExt%s",gTrgSetupName[j]),Form("Sys_signalExt%s",gTrgSetupName[j]),1100,700);
      c->Divide(3,2);
      for(int k=0; k<nCentBins; k++)
	{
	  for(int s=0; s<nSys; s++)
	    {
	      hSignal[k][j][s] = 0x0;
	      if(s==0) hSignal[k][j][s] = (TH1F*)fin->Get(Form("Jpsi_%sYield_cent%s_weight%s%s",method,cent_Title[k],gTrgSetupName[j],sys_name[s]));
	      else     hSignal[k][j][s] = (TH1F*)fsys->Get(Form("Jpsi_%sYield_cent%s_weight%s%s",method,cent_Title[k],gTrgSetupName[j],sys_name[s]));

	      if(!hSignal[k][j][s]) continue;
	      hSignal[k][j][s]->SetMarkerSize(1);
	      hSignal[k][j][s]->SetMarkerStyle(20+s);
	      hSignal[k][j][s]->SetMarkerColor(color_sys[s]);
	  
	      TH1F *htmp = (TH1F*)hSignal[k][j][s]->Clone(Form("%s_tmp",hSignal[k][j][s]->GetName()));
	      htmp->Divide(hSignal[k][j][0]);
	      for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
		{
		  if(s==0) htmp->SetBinError(bin,hSignal[k][j][s]->GetBinError(bin)/hSignal[k][j][s]->GetBinContent(bin));
		  else     htmp->SetBinError(bin,0);
		}
	      htmp->GetYaxis()->SetRangeUser(0.6,1.4);
	      htmp->SetTitle(";p_{T} (GeV/c);Relative difference");
	      if(k==2 || k==3) htmp->GetXaxis()->SetRangeUser(0.15,9);
	      if(k==4) htmp->GetXaxis()->SetRangeUser(0.15,6);
	      c->cd(k+2);
	      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	      ScaleHistoTitle(htmp,0.05,1,0.04,0.05,1,0.04,62);
	      if(s==0) htmp->Draw("P");
	      else htmp->Draw("P sames");
	      TPaveText *t1 = GetTitleText(Form("Signal extraction uncertainty"),0.05);
	      t1->Draw();
	      t1 = GetPaveText(0.25,0.35,0.2,0.3,0.06,62);
	      t1->AddText(Form("%s%%",cent_Name[k]));
	      t1->Draw();

	      for(int i=0; i<nPtBins; i++)
		{
		  if(k==1 && j==3 && i==5 && s==4) continue;
		  if(k==1 && j==2 && i==1 && s==3) continue;
		  double value = fabs(htmp->GetBinContent(i+1)-1);
		  if(value>0.3) continue;
		  if(max[k][j][i]<value) max[k][j][i] = value;
		}
	    }

	  if(k==0)
	    {
	      double xmin = 0, xmax = 0.45, ymin = 0.5, ymax = 0.88;
	      TLegend *leg[2];
	      for(int l=0; l<2; l++)
		{
		  leg[l] = new TLegend(xmin+0.5*l,ymin,xmax+0.5*l,ymax);
		  leg[l]->SetBorderSize(0);
		  leg[l]->SetFillColor(0);
		  leg[l]->SetTextSize(0.04);
		}
	      for(int s=0; s<nSys; s++)
		{
		  leg[s/5]->AddEntry(hSignal[k][j][s],sys_leg[s],"P");
		}
	    }
      
	  TH1F *htmp = (TH1F*)hSys[k][j]->Clone(Form("%s_low",hSys[k][j]->GetName()));
	  for(int bin=1; bin<=hSys[k][j]->GetNbinsX(); bin++)
	    {
	      double sys = max[k][j][bin-1];
	      if(k==0 || k==1)
		{
		  if(bin==9) sys = max[k][j][7];
		}
	      if(k==2 || k==3)
		{
		  if(bin==8) sys = max[k][j][6];
		}
	      if(k==4)
		{
		  if(bin==6) sys = max[k][j][4];
		}
	      hSys[k][j]->SetBinContent(bin,sys+1);
	      htmp->SetBinContent(bin,1-sys);
	      hMax[k][j][0]->SetBinContent(bin,1+max[k][j][bin-1]);
	      hMax[k][j][1]->SetBinContent(bin,1-max[k][j][bin-1]);
	    }
	  hSys[k][j]->SetLineColor(2);
	  hSys[k][j]->SetLineWidth(1);
	  hSys[k][j]->SetLineStyle(2);
	  hSys[k][j]->Draw("samesHIST");
	  htmp->SetLineColor(2);
	  htmp->SetLineWidth(1);
	  htmp->SetLineStyle(2);
	  htmp->Draw("samesHIST");
	  for(int l=0; l<2; l++)
	    {
	      hMax[k][j][l]->SetLineWidth(1);
	      //hMax[k][j][l]->SetLineStyle(2);
	      //hMax[k][j][l]->Draw("samesHIST");
	    }
	}
      c->cd(1);
      leg[0]->Draw();
      leg[1]->Draw();
      c->cd(1);
      TLegend *leg2 =  new TLegend(0.1, 0.3, 0.3, 0.45);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.045);
      leg2->SetHeader(run_type);
      leg2->AddEntry(hSys[j][0],"Assigned uncertainty","L");
      //leg2->AddEntry(hSys[j][0],"Re-assigned uncert.","L");
      leg2->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_sigExtSys%s.pdf",run_type,gTrgSetupTitle[j]));
	}
    }

  TCanvas *c = new TCanvas("cSigSys", "cSigSys", 1100, 700);
  c->Divide(3, 2);
  leg2 = new TLegend(0.2,0.4,0.6,0.85);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.055);
  leg2->SetHeader(run_type);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  for(int bin=1; bin<=hSys[k][j]->GetNbinsX(); bin++)
	    {
	      hSys[k][j]->SetBinContent(bin, hSys[k][j]->GetBinContent(bin)-1);
	    }
	  hSys[k][j]->SetTitle(";p_{T} [GeV/c];sys. uncert.");
	  hSys[k][j]->SetMarkerStyle(20+j);
	  hSys[k][j]->SetMarkerColor(color[j]);
	  hSys[k][j]->SetLineColor(color[j]);
	  hSys[k][j]->SetMarkerSize(1.2);
	  hSys[k][j]->GetYaxis()->SetRangeUser(0, 0.6);
	  if(k==0 || k==1) hSys[k][j]->GetXaxis()->SetRangeUser(0.25,15);
	  if(k==2 || k==3) hSys[k][j]->GetXaxis()->SetRangeUser(0.25,10);
	  if(k==4) hSys[k][j]->GetXaxis()->SetRangeUser(0.25,5);
	  hSys[k][j]->SetLineStyle(1);
	  if(j==0) hSys[k][j]->Draw("HIST");
	  else     hSys[k][j]->Draw("samesHISt");

	  if(k==0) leg2->AddEntry(hSys[k][j], legName[j].Data(), "L");
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(6);
  leg2->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_sigExtSys.pdf",run_type));

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "update");
      for(int k=0; k<nCentBins; k++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      hSys[k][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
  
}


//================================================
void muonPID(const int savePlot = 1, const int saveHisto = 0)
{
  const char* typeName[2] = {"MtdMth", "MuonPid"};
  TFile *fin = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");
  THnSparseF *hnTrkPt[2];
  TH1F *hTrkPt[2][nCentBins][gNTrgSetup];
  TH1F *hTrkPtEff[nCentBins][gNTrgSetup];
  for(int i=0; i<2; i++)
    {
      hnTrkPt[i] = (THnSparseF*)fin->Get(Form("mhMcTrkPtEff_%s_di_mu",typeName[i]));
      hnTrkPt[i]->GetAxis(1)->SetRangeUser(-0.5,0.5);
      for(int k=0; k<nCentBins; k++)
	{
	  hnTrkPt[i]->GetAxis(2)->SetRange(centBins_low[k],centBins_high[k]);
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      if(j>0) hnTrkPt[i]->GetAxis(4)->SetRange(j,j);
	      hTrkPt[i][k][j] = (TH1F*)hnTrkPt[i]->Projection(0);
	      hTrkPt[i][k][j]->SetName(Form("hTrkPt_%s_cent%s%s",typeName[i],cent_Title[k],gTrgSetupName[j]));
	      hTrkPt[i][k][j]->Sumw2();
	      if(k<3) hTrkPt[i][k][j]->Rebin(5);
	      else    hTrkPt[i][k][j]->Rebin(20);
	      if(i==1)
		{
		  hTrkPtEff[k][j] = (TH1F*)hTrkPt[i][k][j]->Clone(Form("hTrkPtEff_%s_cent%s%s",typeName[i],cent_Title[k],gTrgSetupName[j]));
		  hTrkPtEff[k][j]->Divide(hTrkPt[0][k][j]);
		}
		  
	    }
	  hnTrkPt[i]->GetAxis(4)->SetRange(0,-1);
	}
      hnTrkPt[i]->GetAxis(2)->SetRange(0,-1);
    }

  TCanvas *c = new TCanvas("cMuonPid", "cMuonPid", 1100, 700);
  c->Divide(3, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type);
  const char *setupName[5] = {"all","prod","prod_low","prod_mid","prod_high"};
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTrkPtEff[k][j]->SetTitle("");
	  hTrkPtEff[k][j]->SetMarkerStyle(20+j);
	  hTrkPtEff[k][j]->SetMarkerColor(color[j]);
	  hTrkPtEff[k][j]->SetLineColor(color[j]);
	  hTrkPtEff[k][j]->SetMarkerSize(1.2);
	  hTrkPtEff[k][j]->GetYaxis()->SetRangeUser(0.5,1.0);
	  hTrkPtEff[k][j]->GetXaxis()->SetRangeUser(0,10);
	  if(j==0) hTrkPtEff[k][j]->Draw();
	  else     hTrkPtEff[k][j]->Draw("sames");

	  if(k==0) leg->AddEntry(hTrkPtEff[k][j], legName[j].Data(), "PL");
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_MuonPid.pdf",run_type));

  TCanvas *c = new TCanvas("cMuonPidRatio", "cMuonPidRatio", 1100, 700);
  c->Divide(3, 2);
  TH1F *hTrkPtEffRatio[nCentBins][gNTrgSetup-1];
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  hTrkPtEffRatio[k][j] = (TH1F*)hTrkPtEff[k][j+1]->Clone(Form("hTrkPtEffRatio_%s_cent%s%s",typeName[1],cent_Title[k],gTrgSetupName[j+1]));
	  hTrkPtEffRatio[k][j]->Divide(hTrkPtEff[k][0]);
	  hTrkPtEffRatio[k][j]->GetYaxis()->SetRangeUser(0.95,1.05);
	  if(j==0) hTrkPtEffRatio[k][j]->Draw();
	  else     hTrkPtEffRatio[k][j]->Draw("sames");
	  if(k==0)
	    {
	      TF1 *func = new TF1(Form("func_%s",hTrkPtEffRatio[k][j]->GetName()),"pol0",0,15);
	      hTrkPtEffRatio[k][j]->Fit(func,"0R");
	      func->SetLineColor(hTrkPtEffRatio[k][j]->GetLineColor());
	      func->SetLineStyle(2);
	      func->Draw("sames");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_MuonPidRatio.pdf",run_type));
}

//================================================
void MtdTrigEffHadron(const int savePlot = 0, const int saveHisto = 0)
{
  const int nPtBins = 6;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 
  TH2F *hppMuonTacDiffVsTrigUnit[nPtBins];
  TH1F *hppMuonTacDiff[nPtBins];
  TH2F *hppTacDiffVsTrigUnit[nPtBins];
  TH1F *hppTacDiff[nPtBins];
  TF1  *funcppTacDiff[nPtBins];
  TFile *fpp = TFile::Open("Rootfiles/Run15_pp200.MtdTrigEff.root", "read");
  TCanvas *c = new TCanvas("fit_TacDiff_LS", "fit_TacDiff_LS", 1100, 700);
  c->Divide(3,2);
  for(int i=0; i<nPtBins; i++)
    {
      hppMuonTacDiffVsTrigUnit[i] = (TH2F*)fpp->Get(Form("Run15_pp200_hTacDiffVsTrigUnit_Muon_PtBin%d",i+1));
      hppMuonTacDiff[i] = (TH1F*)hppMuonTacDiffVsTrigUnit[i]->ProjectionY(Form("Run15_pp200_hTacDiff_Muon_PtBin%d",i+1));
      hppMuonTacDiff[i]->Rebin(2);
      hppMuonTacDiff[i]->SetMarkerStyle(24);
      hppMuonTacDiff[i]->SetMarkerColor(2);
      hppMuonTacDiff[i]->SetLineColor(2);
      hppMuonTacDiff[i]->Scale(1./hppMuonTacDiff[i]->Integral()*0.5);

      hppTacDiffVsTrigUnit[i] = (TH2F*)fpp->Get(Form("Run15_pp200_hTacDiffVsTrigUnit_LS_PtBin%d",i+1));
      hppTacDiff[i] = (TH1F*)hppTacDiffVsTrigUnit[i]->ProjectionY(Form("Run15_pp200_hTacDiff_LS_PtBin%d",i+1));
      hppTacDiff[i]->Rebin(2);
      hppTacDiff[i]->Scale(1./hppTacDiff[i]->Integral()*0.5);
      if(i<nPtBins-1)
      //if(i<0)
	{
	  funcppTacDiff[i] = new TF1(Form("fit_%s",hppTacDiff[i]->GetName()), CrystalBall, 880, 960, 5);
	  funcppTacDiff[i]->SetParameters(1,925,8,1,0.2);
	}
      else
	{
	  funcppTacDiff[i] = new TF1(Form("fit_%s",hppTacDiff[i]->GetName()), "gaus", 910, 960);
	}

      hppTacDiff[i]->Fit(funcppTacDiff[i], "IR0");
      c->cd(i+1);
      hppTacDiff[i]->SetMarkerStyle(20);
      hppTacDiff[i]->GetXaxis()->SetRangeUser(860, 980);
      hppTacDiff[i]->SetMaximum(1.4*hppTacDiff[i]->GetMaximum());
      hppTacDiff[i]->Draw();
      hppMuonTacDiff[i]->Draw("sames");
      funcppTacDiff[i]->SetLineStyle(2);
      funcppTacDiff[i]->SetLineColor(4);
      funcppTacDiff[i]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("Run15_pp200: %1.1f < p_{T} < %1.1f GeV/c",xPtBins[i],xPtBins[i+1]),0.05);
      t1->Draw();
    }

  // Get the mean and sigma for TacDiff distributions from Au+Au
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root", "update");
  else          fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root", "read");
  const int nLumi = 4;
  const char *lumi_name[nLumi] = {"prod","prod_low","prod_mid","prod_high"};
  const double min_TacDiffCut[nLumi] = {785+1,785+1,788+1,788+1};
  const double max_TacDiffCut[nLumi] = {980, 980, 980, 980};
  TH1F *hTacDiffPtMean[gNTrgSetup];
  TH1F *hTacDiffPtSigma[gNTrgSetup];  
  for(int j=0; j<gNTrgSetup; j++)
    {
      hTacDiffPtMean[j]  = (TH1F*)fin->Get(Form("hTacDiffMeanVsPt%s",gTrgSetupName[j]));
      hTacDiffPtSigma[j] = (TH1F*)fin->Get(Form("hTacDiffSigmaVsPt%s",gTrgSetupName[j]));
    }
  TGraphAsymmErrors* gTacDiffEff[nLumi]; // 0 - bin counting; 1 - fitting;
  double min_TacDiffCut_new[nLumi][nPtBins];
  double pt_arr[2][nPtBins], all_arr[2][nPtBins], all_err_arr[2][nPtBins], acc_arr[2][nPtBins], acc_err_arr[2][nPtBins];
  for(int l=0; l<nLumi; l++)
    {
      printf("+++ %s +++\n",lumi_name[l]);
      // calculate new cut values
      for(int bin=1; bin<=nPtBins; bin++)
	{
	  double pt = hTacDiffPtMean[0]->GetBinCenter(bin);
	  double sigma_auau = hTacDiffPtSigma[l+1]->GetBinContent(bin);
	  double mean_auau = hTacDiffPtMean[l+1]->GetBinContent(bin);
	  double sigma_pp = funcppTacDiff[bin-1]->GetParameter(2);
	  double mean_pp = funcppTacDiff[bin-1]->GetParameter(1);
	  min_TacDiffCut_new[l][bin-1] = mean_pp - (mean_auau - min_TacDiffCut[l]) * sigma_pp/sigma_auau;
	  if(bin==nPtBins) min_TacDiffCut_new[l][bin-1] = min_TacDiffCut_new[l][bin-2];
	  if(l==1 || l==3) printf("[i] %s: auau (%4.2f, %4.2f), pp (%4.2f, %4.2f), old_cut = %4.2f, new_cut = %4.2f\n",gTrgSetupTitle[l+1],mean_auau,sigma_auau,mean_pp,sigma_pp,min_TacDiffCut[l],min_TacDiffCut_new[l][bin-1]);
	}

      // bin counting
      for(int bin=1; bin<=nPtBins; bin++)
	{
	  TH1F *hppTac = hppTacDiff[bin-1];
	  int xbins = hppTac->GetNbinsX(); 
	  pt_arr[0][bin-1] = hppTac->GetBinCenter(bin);
	  all_arr[0][bin-1] = hppTac->IntegralAndError(1, xbins, all_err_arr[0][bin-1]);
	  int low_bin =    hppTac->FindFixBin(min_TacDiffCut_new[l][bin-1]);
	  acc_arr[0][bin-1] = hppTac->IntegralAndError(low_bin+1,
						       hppTac->FindFixBin(max_TacDiffCut[l]-0.1),
						       acc_err_arr[0][bin-1]);
	  double fraction = (hppTac->GetXaxis()->GetBinUpEdge(low_bin)-min_TacDiffCut_new[l][bin-1])/hppTac->GetBinWidth(low_bin);
	  acc_arr[0][bin-1] += hppTac->GetBinContent(low_bin)*fraction;
	  acc_err_arr[0][bin-1] = sqrt(pow(acc_err_arr[0][bin-1],2)+pow(hppTac->GetBinError(low_bin)*fraction,2));
	  if(l==1 || l==3) printf("[i] %s: all = %4.2f, acc = %4.2f, eff = %4.2f%%\n",gTrgSetupTitle[l+1],all_arr[0][bin-1],acc_arr[0][bin-1],acc_arr[0][bin-1]/all_arr[0][bin-1]*100);
	}

      // fitting
      for(int bin=1; bin<=nPtBins; bin++)
	{
	  TF1 *funcppTac    = funcppTacDiff[bin-1];
	  pt_arr[1][bin-1]  = pt_arr[0][bin-1];
	  all_arr[1][bin-1] = funcppTac->Integral(840,max_TacDiffCut[l]);
	  acc_arr[1][bin-1] = funcppTac->Integral(min_TacDiffCut_new[l][bin-1],max_TacDiffCut[l]);
	  all_err_arr[1][bin-1] = all_err_arr[0][bin-1]/all_arr[0][bin-1] * all_arr[1][bin-1];
	  acc_err_arr[1][bin-1] = acc_err_arr[0][bin-1]/acc_arr[0][bin-1] * acc_arr[1][bin-1];	  
	}
      
      gTacDiffEff[l] = GetEfficiencyCurve(nPtBins, pt_arr[1], all_arr[1], all_err_arr[1], acc_arr[1], acc_err_arr[1]);
      gTacDiffEff[l]->SetName(Form("%s_gTacDiffEff_LS_%s",run_type,lumi_name[l]));
    }

  TFile *fRun14 = TFile::Open("Rootfiles/Run14_AuAu200.MtdTrigEff.root","read");
  TGraphAsymmErrors *gRun14 = (TGraphAsymmErrors*)fRun14->Get("Run14_AuAu200_gTacDiffEffFinal_BinCount_prod_Run15_pp200");
  double x, y, x1, y1;
  for(int l=0; l<nLumi; l++)
    {
      int npoints = gTacDiffEff[l]->GetN();
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  gRun14->GetPoint(ipoint, x1, y1);
	  gTacDiffEff[l]->GetPoint(ipoint,x,y);
	  gTacDiffEff[l]->SetPoint(ipoint,x1,y);
	  gTacDiffEff[l]->SetPointEXhigh(ipoint,gRun14->GetErrorXhigh(ipoint));
	  gTacDiffEff[l]->SetPointEXlow(ipoint,gRun14->GetErrorXlow(ipoint));
	}
    }

  TF1 *funcTacEff[nLumi];
  TH1F *hplot = new TH1F("hplot",";p_{T}^{#mu} [GeV/c];Efficiency",100,1,7);
  for(int l=0; l<nLumi; l++)
    {
      funcTacEff[l] = new TF1(Form("%s_gTacDiffEff_LS_%s_func",run_type,lumi_name[l]),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
      funcTacEff[l]->SetParameters(0.9, 3.5, 0.8);
      gTacDiffEff[l]->Fit(funcTacEff[l],"IR0");
      TCanvas *c = new TCanvas(Form("%s_TacDiffEff_%s",run_type,lumi_name[l]),Form("%s_TacDiffEff_%s",run_type,lumi_name[l]),800,600);
      hplot->GetYaxis()->SetRangeUser(0.6,1.1);
      hplot->GetYaxis()->SetNdivisions(505);
      hplot->Draw();
      gTacDiffEff[l]->SetMarkerStyle(21);
      gTacDiffEff[l]->SetMarkerColor(1);
      gTacDiffEff[l]->SetLineColor(1);
      gTacDiffEff[l]->Draw("PEsame");
      funcTacEff[l]->SetLineColor(2);
      funcTacEff[l]->SetLineStyle(2);
      funcTacEff[l]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s: estimated MTD trigger efficiency for LS tracks",run_type));
      t1->Draw();
      TLegend *leg = new TLegend(0.45,0.15,0.65,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(lumi_name[l]);
      leg->AddEntry(gTacDiffEff[l],"Data: bin counting","P");
      leg->AddEntry(funcTacEff[l],"Fit: p0-e^{-p1*(x-p2)}","L");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%s_TacDiffEffWithFit_%s.pdf",run_type,run_type,lumi_name[l]));
    }

  if(saveHisto)
    {
      fin->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  funcppTacDiff[i]->Write("",TObject::kOverwrite);
	}
      for(int l=0; l<nLumi; l++)
	{
	  gTacDiffEff[l]->Write("",TObject::kOverwrite);
	  funcTacEff[l]->Write("",TObject::kOverwrite);
	} 
    }
  
}

//================================================
void anaTrigEff(const int savePlot = 0, const int saveHisto = 0)
{
  //==============================================
  // compare the TacDiff in MT101
  TFile *fStudy = TFile::Open("output/Run14_AuAu200.Study.MtdTrig.root", "read");
  TH2F *h2TacDiff = (TH2F*)fStudy->Get("mhMT101TacDiff_di-muon");
  TH1F *hTacDiffMT101[4];
  TF1 *funcTacDiffMT101[4];
  TCanvas *c = new TCanvas("cTacDiffMT101", "cTacDiffMT101", 800, 600);
  leg = new TLegend(0.15,0.6,0.4,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader(run_type);
  for(int j=0; j<gNTrgSetup-1; j++)
    {
      hTacDiffMT101[j] = (TH1F*)h2TacDiff->ProjectionY(Form("hTacDiffMT101_%s",gTrgSetupName[j+1]), j+1, j+1);
      hTacDiffMT101[j]->Sumw2();
      hTacDiffMT101[j]->SetXTitle("TacSum_{MT101} - TacSum_{VPD}/8 + 1024");
      hTacDiffMT101[j]->SetMarkerStyle(20+j+1);
      hTacDiffMT101[j]->SetMarkerColor(color[j+1]);
      hTacDiffMT101[j]->SetLineColor(color[j+1]);
      hTacDiffMT101[j]->SetMarkerSize(1.2);
      hTacDiffMT101[j]->Scale(1./hTacDiffMT101[j]->GetBinContent(hTacDiffMT101[j]->FindBin(790)));
      hTacDiffMT101[j]->GetYaxis()->SetRangeUser(1e-2,1.5);
      hTacDiffMT101[j]->GetXaxis()->SetRangeUser(760, 820);
      if(j==0) hTacDiffMT101[j]->DrawCopy();
      else     hTacDiffMT101[j]->DrawCopy("sames");
      leg->AddEntry(hTacDiffMT101[j], legName[j+1].Data(), "P");
      funcTacDiffMT101[j] = new TF1(Form("func_%s",hTacDiffMT101[j]->GetName()), "gaus", 790, 805);
      if(j<2) funcTacDiffMT101[j]->SetRange(786, 805);
      if(j==3) funcTacDiffMT101[j]->SetRange(790, 805);
      funcTacDiffMT101[j]->SetParameters(1, 795, 8);
      hTacDiffMT101[j]->Fit(funcTacDiffMT101[j], "IR0Q");
      funcTacDiffMT101[j]->SetLineColor(color[j+1]);
      funcTacDiffMT101[j]->SetLineStyle(2);
      funcTacDiffMT101[j]->Draw("sames");
      cout << hTacDiffMT101[j]->Integral(787,789)/hTacDiffMT101[j]->Integral(787,1024) << endl;
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_MT101TacDiff.pdf",run_type));


  //==============================================
  // matched tracks
  const char* type_name[3] = {"Muon","UL", "LS"};
  const char* type_title[3] = {"Muon = UL - LS","Unlike-sign", "Like-sign"};
  const int nTrigUnit = 28;
  const int nPtBins = 6;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 
  const double fit_min = 789;
  const double fit_max = 806;

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "update");
  else          fin = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "read");

  TCanvas *c = new TCanvas("cTacDiff", "cTacDiff", 1100, 700);
  c->Divide(2, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type);
  TH1F *hTacDiffAll[3][gNTrgSetup];
  for(int i=0; i<3; i++)
    {
      c->cd(i+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffAll[i][j] = (TH1F*)fin->Get(Form("hTacDiffAll_%s%s",type_name[i],gTrgSetupName[j]));
	  hTacDiffAll[i][j]->SetXTitle("#DeltaTacSum");
	  scaleHisto(hTacDiffAll[i][j], 1, 1, true, false, false);
	  hTacDiffAll[i][j]->SetMarkerStyle(20+j);
	  hTacDiffAll[i][j]->SetMarkerColor(color[j]);
	  hTacDiffAll[i][j]->SetLineColor(color[j]);
	  hTacDiffAll[i][j]->SetMarkerSize(1.2);
	  hTacDiffAll[i][j]->Scale(1./hTacDiffAll[i][j]->GetBinContent(hTacDiffAll[i][j]->FindBin(795)));
	  hTacDiffAll[i][j]->GetYaxis()->SetRangeUser(1e-2,1.2);
	  hTacDiffAll[i][j]->GetXaxis()->SetRangeUser(780, 820);
	  if(j==0) hTacDiffAll[i][j]->DrawCopy();
	  else     hTacDiffAll[i][j]->DrawCopy("sames");
	  if(i==0) leg->AddEntry(hTacDiffAll[i][j], legName[j].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(type_title[i],0.05);
      t1->Draw();
    }
  c->cd(4);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_TacDiffDis.pdf",run_type));
  for(int i=0; i<3; i++)
    {
      c->cd(i+1);
      gPad->SetLogy();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_TacDiffDis_Log.pdf",run_type));

  //==============================================
  // matched tracks vs MT101 vs pp shifted muons
  TFile *flumi = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  TF1 *funcppLS = (TF1*)flumi->Get("fit_Run15_pp200_hTacDiff_LS_PtBin1");
  TFile *fMtdTrig = TFile::Open("Rootfiles/Run14_AuAu200.MtdTrigEff.root", "read");
  TH1F *hppShifted = (TH1F*)fMtdTrig->Get("Run15_pp200_hMuonTacDiffCombined_PtBin0");
  hppShifted->Scale(1./hppShifted->GetBinContent(hppShifted->FindBin(795)));
  TCanvas *c = new TCanvas("TacSum_LSvsMT101", "TacSum_LSvsMT101", 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.6,0.35,0.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  leg->SetHeader(run_type);
  for(int j=0; j<gNTrgSetup-1; j++)
    {
      double tacSumCutMin = 786;
      double tacSumCutMax = 837;
      if(j>1) tacSumCutMin = 789;

      c->cd(j+1);
      hTacDiffMT101[j]->SetMarkerStyle(20);
      hTacDiffMT101[j]->SetMarkerColor(1);
      hTacDiffMT101[j]->SetLineColor(1);
      hTacDiffMT101[j]->SetTitle(";#DeltaTacSum;a.u.");
      hTacDiffMT101[j]->DrawCopy();
      funcTacDiffMT101[j]->SetLineColor(kGray);
      funcTacDiffMT101[j]->Draw("sames");
      double tacSumCutMinNew = funcppLS->GetParameter(1) - (funcTacDiffMT101[j]->GetParameter(1)-tacSumCutMin)/funcTacDiffMT101[j]->GetParameter(2)*funcppLS->GetParameter(2);
      double tacSumCutMaxNew = funcppLS->GetParameter(1) + (tacSumCutMax - funcTacDiffMT101[j]->GetParameter(1))/funcTacDiffMT101[j]->GetParameter(2)*funcppLS->GetParameter(2);
      double eff = funcppLS->Integral(tacSumCutMinNew,tacSumCutMaxNew)/funcppLS->Integral(900,1200);
      TPaveText *t1 = GetPaveText(0.65,0.85,0.75,0.85,0.05);
      t1->AddText(Form("#mu = %4.2f, #sigma = %4.2f",funcTacDiffMT101[j]->GetParameter(1),funcTacDiffMT101[j]->GetParameter(2)));
      t1->AddText(Form("eff = %4.3f",eff));
      t1->Draw();

      hppShifted->SetMarkerStyle(25);
      hppShifted->SetMarkerColor(4);
      hppShifted->SetLineColor(4);
      hppShifted->Draw("sames");

      hTacDiffAll[2][j+1]->SetMarkerStyle(24);
      hTacDiffAll[2][j+1]->SetMarkerColor(2);
      hTacDiffAll[2][j+1]->SetLineColor(2);
      hTacDiffAll[2][j+1]->DrawCopy("sames");
      TF1 *functmp = new TF1(Form("fit_%s",hTacDiffAll[2][j+1]->GetName()),"gaus",786,805);
      if(j>1) functmp->SetRange(790,805);
      functmp->SetParameters(1, 795, 8);
      hTacDiffAll[2][j+1]->Fit(functmp, "IR0Q");
      functmp->SetLineStyle(2);
      functmp->SetLineColor(2);
      functmp->Draw("sames");
      tacSumCutMinNew = funcppLS->GetParameter(1) - (functmp->GetParameter(1)-tacSumCutMin)/functmp->GetParameter(2)*funcppLS->GetParameter(2);
      tacSumCutMaxNew = funcppLS->GetParameter(1) + (tacSumCutMax - functmp->GetParameter(1))/functmp->GetParameter(2)*funcppLS->GetParameter(2);
      eff = funcppLS->Integral(tacSumCutMinNew,tacSumCutMaxNew)/funcppLS->Integral(900,1200);
      t1 = GetPaveText(0.65,0.85,0.62,0.72,0.05);
      t1->AddText(Form("#mu = %4.2f, #sigma = %4.2f",functmp->GetParameter(1),functmp->GetParameter(2)));
      t1->AddText(Form("eff = %4.3f",eff));
      t1->SetTextColor(2);
      t1->Draw();
      TPaveText *t1 = GetTitleText(Form("%s: #DeltaTacSum distribution",legName[j+1].Data()),0.055);
      t1->Draw();
      if(j==0)
	{
	  leg->AddEntry(hTacDiffMT101[j], "All hits", "P");
	  leg->AddEntry(hTacDiffAll[2][j+1], "Like-sign pairs", "P");
	  leg->AddEntry(hppShifted, "pp muons (shifted)", "P");
	}
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_CompTacSumDiff.pdf",run_type));
  return;


  //==============================================
  // fit TacSum distributions
  TH1F *hTacDiffPt[nPtBins][gNTrgSetup];
  TF1 *funcTacDiffPt[nPtBins][gNTrgSetup];
  TH1F *hTacDiffPtMean[gNTrgSetup];
  TH1F *hTacDiffPtSigma[gNTrgSetup];  
  for(int j=0; j<gNTrgSetup; j++)
    {
      hTacDiffPtMean[j]  = new TH1F(Form("hTacDiffMeanVsPt%s",gTrgSetupName[j]),"Mean of #DeltaTacSum;p_{T} (GeV/c);<#DeltaTacSum>", nPtBins, xPtBins);
      hTacDiffPtSigma[j] = new TH1F(Form("hTacDiffSigmaVsPt%s",gTrgSetupName[j]),"Width of #DeltaTacSum;p_{T} (GeV/c);#sigma(#DeltaTacSum)", nPtBins, xPtBins);      
    }
  for(int t=0; t<nPtBins; t++)
    {
      TCanvas *c = new TCanvas(Form("cFitTacDiff_%d",t), Form("cFitTacDiff_Pt%d",t), 1100, 700);
      c->Divide(3, 2);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffPt[t][j] = (TH1F*)fin->Get(Form("hTacDiff_Pt%d_%s%s",t+1,type_name[2],gTrgSetupName[j]));
	  if(t==nPtBins-1) hTacDiffPt[t][j]->Rebin(2);
	  scaleHisto(hTacDiffPt[t][j], 1, 1, true, false, false);
	  double fit_min_tmp = fit_min;
	  double fit_max_tmp = fit_max;
	  //if(j==1 || j==2) fit_min_tmp = 786;
	  if(t>=nPtBins-1) fit_max_tmp = 815;
	  funcTacDiffPt[t][j] = new TF1(Form("fit_%s",hTacDiffPt[t][j]->GetName()),"gaus",fit_min_tmp,fit_max_tmp);
	  funcTacDiffPt[t][j]->SetParameters(1, 795, 8);
	  hTacDiffPt[t][j]->Fit(funcTacDiffPt[t][j], "IR0Q");
	  c->cd(j+1);
	  hTacDiffPt[t][j]->GetXaxis()->SetRangeUser(780, 820);
	  hTacDiffPt[t][j]->SetMarkerStyle(20);
	  hTacDiffPt[t][j]->Draw();
	  funcTacDiffPt[t][j]->SetLineColor(2);
	  funcTacDiffPt[t][j]->SetLineStyle(2);
	  funcTacDiffPt[t][j]->Draw("sames");
	  TPaveText *t1 = GetPaveText(0.6, 0.8, 0.7, 0.85, 0.045);
	  t1->AddText(Form("#mu = %4.2f#pm%4.2f",funcTacDiffPt[t][j]->GetParameter(1),funcTacDiffPt[t][j]->GetParError(1)));
	  t1->AddText(Form("#sigma = %4.2f#pm%4.2f",funcTacDiffPt[t][j]->GetParameter(2),funcTacDiffPt[t][j]->GetParError(2)));
	  t1->Draw();
	  hTacDiffPtMean[j]->SetBinContent(t+1, funcTacDiffPt[t][j]->GetParameter(1));
	  hTacDiffPtMean[j]->SetBinError(t+1, funcTacDiffPt[t][j]->GetParError(1));
	  hTacDiffPtMean[j]->SetMarkerStyle(20+j);
	  hTacDiffPtMean[j]->SetMarkerColor(color[j]);
	  hTacDiffPtMean[j]->SetLineColor(color[j]);
	  hTacDiffPtMean[j]->SetMarkerSize(1.2);
	  hTacDiffPtSigma[j]->SetBinContent(t+1, funcTacDiffPt[t][j]->GetParameter(2));
	  hTacDiffPtSigma[j]->SetBinError(t+1, funcTacDiffPt[t][j]->GetParError(2));
	  hTacDiffPtSigma[j]->SetMarkerStyle(20+j);
	  hTacDiffPtSigma[j]->SetMarkerColor(color[j]);
	  hTacDiffPtSigma[j]->SetLineColor(color[j]);
	  hTacDiffPtSigma[j]->SetMarkerSize(1.2);
	  TPaveText *t1 = GetTitleText(Form("%s: %1.1f < p_{T} < %1.1f GeV/c",legName[j].Data(),xPtBins[t],xPtBins[t+1]),0.05);
	  t1->Draw();
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_FitTacDiff_Pt%d.pdf",run_type,t));
    }  

  TCanvas *c = new TCanvas("cTacDiffVsPt", "cTacDiffVsPt", 1100, 500);
  c->Divide(2,1);
  leg = new TLegend(0.5,0.15,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type);
  c->cd(1);
  SetPadMargin(gPad, 0.11, 0.13, 0.1, 0.1);
  for(int j=0; j<gNTrgSetup; j++)
    {
      ScaleHistoTitle(hTacDiffPtMean[j], 0.045, 1.1, 0.04, 0.045, 1.2, 0.04, 62);
      hTacDiffPtMean[j]->SetTitle("");
      hTacDiffPtMean[j]->GetYaxis()->SetRangeUser(792, 802);
      if(j==0) hTacDiffPtMean[j]->Draw();
      else     hTacDiffPtMean[j]->Draw("sames");
      leg->AddEntry(hTacDiffPtMean[j], legName[j].Data(), "P");
    }
  TPaveText *t1 = GetTitleText("Mean of #DeltaTacSum distribution",0.05);
  t1->Draw();
  leg->Draw();
  c->cd(2);
  SetPadMargin(gPad, 0.11, 0.13, 0.1, 0.1);
  for(int j=0; j<gNTrgSetup; j++)
    {
      ScaleHistoTitle(hTacDiffPtSigma[j], 0.045, 1.1, 0.04, 0.045, 1.2, 0.04, 62);
      hTacDiffPtSigma[j]->SetTitle("");
      hTacDiffPtSigma[j]->GetYaxis()->SetRangeUser(4, 9);
      if(j==0) hTacDiffPtSigma[j]->Draw();
      else     hTacDiffPtSigma[j]->Draw("sames");
    }
  TPaveText *t1 = GetTitleText("Width of #DeltaTacSum distribution",0.05);
  t1->Draw();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_TacDiffFitVsPt.pdf",run_type));

  TH1F *hTacDiffTrigUnit[nPtBins][nTrigUnit][gNTrgSetup];
  TF1 *funcTacDiffTrigUnit[nPtBins][nTrigUnit][gNTrgSetup];
  TH1F *hTacDiffMean[nPtBins][gNTrgSetup];
  TH1F *hTacDiffSigma[nPtBins][gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {
      for(int t=0; t<nPtBins-1; t++)
	{
	  hTacDiffMean[t][j]  = new TH1F(Form("hTacDiffMean_Pt%d%s",t,gTrgSetupName[j]), "Mean of #DeltaTacSum;TrigUnit;<#DeltaTacSum>", nTrigUnit, 0, nTrigUnit);
	  hTacDiffSigma[t][j] = new TH1F(Form("hTacDiffSigma_Pt%d%s",t,gTrgSetupName[j]), "Width of #DeltaTacSum;TrigUnit;#sigma(#DeltaTacSum)", nTrigUnit, 0, nTrigUnit);	  
	  TCanvas *c = new TCanvas(Form("cFitTacDiff_%d%d",j,t), Form("cFitTacDiff_Pt%d%s",t,gTrgSetupName[j]), 1100, 700);
	  c->Divide(6, 5);
	  for(int k=0; k<nTrigUnit; k++)
	    {
	      hTacDiffTrigUnit[t][k][j] = (TH1F*)fin->Get(Form("hTacDiff_Pt%d_Unit%d_%s%s",t+1,k+1,type_name[2],gTrgSetupName[j]));
	      scaleHisto(hTacDiffTrigUnit[t][k][j], 1, 1, true, false, false);
	      double fit_min_tmp = 789;
	      double fit_max_tmp = 810;
	      if(j==1 || j==2) fit_min_tmp = 786;
	      if(k==1) fit_max_tmp = 799;
	      if(k==4 || k==5) fit_max_tmp = 803;
	      funcTacDiffTrigUnit[t][k][j] = new TF1(Form("fit_%s",hTacDiffTrigUnit[t][k][j]->GetName()),"gaus",fit_min_tmp,fit_max_tmp);
	      funcTacDiffTrigUnit[t][k][j]->SetParameters(1, 795, 8);
	      hTacDiffTrigUnit[t][k][j]->Fit(funcTacDiffTrigUnit[t][k][j], "IR0Q");
	      c->cd(k+1);
	      hTacDiffTrigUnit[t][k][j]->GetXaxis()->SetRangeUser(780, 820);
	      hTacDiffTrigUnit[t][k][j]->SetMarkerStyle(24);
	      hTacDiffTrigUnit[t][k][j]->SetMarkerSize(1.2);
	      hTacDiffTrigUnit[t][k][j]->Draw();
	      funcTacDiffTrigUnit[t][k][j]->Draw("sames");
	      hTacDiffMean[t][j]->SetBinContent(k+1, funcTacDiffTrigUnit[t][k][j]->GetParameter(1));
	      hTacDiffMean[t][j]->SetBinError(k+1, funcTacDiffTrigUnit[t][k][j]->GetParError(1));
	      hTacDiffMean[t][j]->SetMarkerStyle(20+j);
	      hTacDiffMean[t][j]->SetMarkerColor(color[j]);
	      hTacDiffMean[t][j]->SetLineColor(color[j]);
	      hTacDiffMean[t][j]->SetMarkerSize(1.2);
	      hTacDiffSigma[t][j]->SetBinContent(k+1, funcTacDiffTrigUnit[t][k][j]->GetParameter(2));
	      hTacDiffSigma[t][j]->SetBinError(k+1, funcTacDiffTrigUnit[t][k][j]->GetParError(2));
	      hTacDiffSigma[t][j]->SetMarkerStyle(20+j);
	      hTacDiffSigma[t][j]->SetMarkerColor(color[j]);
	      hTacDiffSigma[t][j]->SetLineColor(color[j]);
	      hTacDiffSigma[t][j]->SetMarkerSize(1.2);
	    }
	}
    }
  
  TCanvas *c = new TCanvas("cTacDiffMean", "cTacDiffMean", 1100, 700);
  c->Divide(3, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type);
  for(int t=0; t<nPtBins-1; t++)
    {
      c->cd(t+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffMean[t][j]->SetTitle("");
	  hTacDiffMean[t][j]->GetYaxis()->SetRangeUser(785, 810);
	  if(j==0) hTacDiffMean[t][j]->Draw();
	  else     hTacDiffMean[t][j]->Draw("sames");
	  if(t==0) leg->AddEntry(hTacDiffMean[t][j], legName[j].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xPtBins[t],xPtBins[t+1]),0.05);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_TacDiffFitMeanVsUnit.pdf",run_type));

  TCanvas *c = new TCanvas("cTacDiffSigma", "cTacDiffSigma", 1100, 700);
  c->Divide(3, 2);
  for(int t=0; t<nPtBins-1; t++)
    {
      c->cd(t+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffSigma[t][j]->SetTitle("");
	  hTacDiffSigma[t][j]->GetYaxis()->SetRangeUser(0, 20);
	  if(j==0) hTacDiffSigma[t][j]->Draw();
	  else     hTacDiffSigma[t][j]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xPtBins[t],xPtBins[t+1]),0.05);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_TacDiffFitSigmaVsUnit.pdf",run_type));

  if(saveHisto)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffPtMean[j]->Write("",TObject::kOverwrite);
	  hTacDiffPtSigma[j]->Write("",TObject::kOverwrite);
	}
    }	  
}

//================================================
void makeTrigEff(const int saveHisto = 1)
{
  const int nBinsTacDiff1 = 35;
  const double xBinsTacDiff1[nBinsTacDiff1+1] = {760,765,770,775,780,782,784,786,787,788,789,790,791,792,793,795,797,799,801,803,805,807,809,811,813,815,817,819,821,823,825,827,829,833,837,841};

  const int nBinsTacDiff2 = 19;
  const double xBinsTacDiff2[nBinsTacDiff2+1] = {760,765,770,775,780,782,784,786,788,789,790,792,795,801,807,813,819,825,833,841};
  const char* type_name[3] = {"Muon","UL", "LS"};
  const int nTrigUnit = 28;
  const int nPtBins = 6;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 


  TFile *fin = TFile::Open("output/Run14_AuAu200.JpsiMuon.root", "read");
  THnSparseF *hn = (THnSparseF*)fin->Get("mhJpsiMuonTrigEff_di_mu");
  hn->GetAxis(0)->SetRangeUser(3.0, 3.2); // cut in Jpsi region
  
  TH1F *hTacDiffAll[3][gNTrgSetup];
  TH1F *hTacDiffPt[2][nPtBins][gNTrgSetup];
  TH1F *hTacDiffTrigUnit[2][nPtBins][nTrigUnit][gNTrgSetup];

  // get UL and LS
  for(int j=0; j<gNTrgSetup; j++)
    {
      if(j>0) hn->GetAxis(7)->SetRange(j, j);
      for(int i=0; i<2; i++)
	{
	  hn->GetAxis(5)->SetRange(i+1,i+1);
	  hn->GetAxis(2)->SetRangeUser(1.3, 10);
	  TH1F *h1tmp = (TH1F*)hn->Projection(1);
	  h1tmp->SetName(Form("hTacDiffAll_%d%d",i,j));
	  //hTacDiffAll[i+1][j] = (TH1F*)h1tmp->Rebin(nBinsTacDiff1, Form("hTacDiffAll_%s%s",type_name[i+1],gTrgSetupName[j]), xBinsTacDiff1);
	  hTacDiffAll[i+1][j] = (TH1F*)h1tmp->Clone(Form("hTacDiffAll_%s%s",type_name[i+1],gTrgSetupName[j]));
	  hTacDiffAll[i+1][j]->Sumw2();
	  for(int t=0; t<nPtBins; t++)
	    {
	      hn->GetAxis(2)->SetRangeUser(xPtBins[t], xPtBins[t+1]);
	      TH1F *h1tmp = (TH1F*)hn->Projection(1);
	      h1tmp->SetName(Form("hTacDiffPt%d%d%d",t,i,j));
	      hTacDiffPt[i][t][j] = (TH1F*)h1tmp->Rebin(nBinsTacDiff1, Form("hTacDiff_Pt%d_%s%s",t+1,type_name[i+1],gTrgSetupName[j]), xBinsTacDiff1);
	      TH2F *h2tmp = (TH2F*)hn->Projection(1, 3);
	      h2tmp->SetName(Form("hTacDiffVsTrigUnit_%d%d%d",j,i,t));
	      h2tmp->Sumw2();
	      for(int k=0; k<nTrigUnit; k++)
		{
		  TH1F *h1tmp = (TH1F*)h2tmp->ProjectionY(Form("hTacDiffTrigUnit_%d%d%d%d",i,j,t,k), k+2, k+2);
		  hTacDiffTrigUnit[i][t][k][j] = (TH1F*)h1tmp->Rebin(nBinsTacDiff1, Form("hTacDiff_Pt%d_Unit%d_%s%s",t+1,k+1,type_name[i+1],gTrgSetupName[j]), xBinsTacDiff1);
		  //hTacDiffTrigUnit[i][t][k][j] = (TH1F*)h2tmp->ProjectionY(Form("hTacDiff_Pt%d_Unit%d_%s%s",t+1,k+1,type_name[i+1],gTrgSetupName[j]), k+2, k+2);
		}
	    }
	  hn->GetAxis(2)->SetRange(0,-1);
	}
      hn->GetAxis(5)->SetRange(0,-1);
    }
  hn->GetAxis(7)->SetRange(0,-1);

  // get pure muon
  for(int j=0; j<gNTrgSetup; j++)
    {
      hTacDiffAll[0][j] = (TH1F*)hTacDiffAll[1][j]->Rebin(nBinsTacDiff2, Form("hTacDiffAll_%s%s",type_name[0],gTrgSetupName[j]), xBinsTacDiff2);
      TH1F *h1Rebin = (TH1F*)hTacDiffAll[2][j]->Rebin(nBinsTacDiff2, Form("%s_rebin",hTacDiffAll[2][j]->GetName()), xBinsTacDiff2);
      hTacDiffAll[0][j]->Add(h1Rebin, -1);
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "update");
      for(int j=0; j<gNTrgSetup; j++)
	{
	  for(int i=0; i<3; i++)
	    {
	      hTacDiffAll[i][j]->SetTitle("");
	      hTacDiffAll[i][j]->Write("",TObject::kOverwrite);
	      if(i==0) continue;
	      for(int t=0; t<nPtBins; t++)
		{
		  hTacDiffPt[i-1][t][j]->SetTitle("");
		  hTacDiffPt[i-1][t][j]->Write("",TObject::kOverwrite);
		  for(int k=0; k<nTrigUnit; k++)
		    {
		      hTacDiffTrigUnit[i-1][t][k][j]->SetTitle("");
		      hTacDiffTrigUnit[i-1][t][k][j]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
      fout->Close();
    }
  
}

//================================================
void vzDiff(const int savePlot = 0)
{
  TH1F *hTpcVz[nCentBins][gNTrgSetup];
  TH1F *hVzDiff[nCentBins][gNTrgSetup];
  TFile *fin = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "read");
  
  TCanvas *c = new TCanvas("cTpcVz", "cTpcVz", 1100, 700);
  c->Divide(3, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type);
  const char *setupName[5] = {"all","prod","prod_low","prod_mid","prod_high"};
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTpcVz[k][j] = (TH1F*)fin->Get(Form("hTpcVz_cent%s%s",cent_Title[k],gTrgSetupName[j]));
	  hTpcVz[k][j]->SetMarkerStyle(20+j);
	  hTpcVz[k][j]->SetMarkerColor(color[j]);
	  hTpcVz[k][j]->SetLineColor(color[j]);
	  hTpcVz[k][j]->SetMarkerSize(1.2);
	  hTpcVz[k][j]->Scale(1./hTpcVz[k][j]->GetBinContent(hTpcVz[k][j]->FindBin(0)));
	  hTpcVz[k][j]->GetYaxis()->SetRangeUser(0,1.2);
	  if(j==0) hTpcVz[k][j]->Draw();
	  else     hTpcVz[k][j]->Draw("sames");

	  if(k==0) leg->AddEntry(hTpcVz[k][j], legName[j].Data(), "L");
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_TpcVz.pdf",run_type));

  TCanvas *c = new TCanvas("cVzDiff", "cVzDiff", 1100, 700);
  c->Divide(3, 2);
  leg2 = new TLegend(0.2,0.3,0.6,0.4);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.055);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hVzDiff[k][j] = (TH1F*)fin->Get(Form("hVzDiff_cent%s%s",cent_Title[k],gTrgSetupName[j]));
	  hVzDiff[k][j]->SetMarkerStyle(20+j);
	  hVzDiff[k][j]->SetMarkerColor(color[j]);
	  hVzDiff[k][j]->SetLineColor(color[j]);
	  hVzDiff[k][j]->SetMarkerSize(1.2);
	  hVzDiff[k][j]->Scale(1./hVzDiff[k][j]->GetBinContent(hVzDiff[k][j]->FindBin(0)));
	  hVzDiff[k][j]->GetXaxis()->SetRangeUser(-20,20);
	  hVzDiff[k][j]->GetYaxis()->SetRangeUser(1e-6, 10);
	  if(j==0) hVzDiff[k][j]->Draw();
	  else     hVzDiff[k][j]->Draw("sames");
	  if(j==0 && (k!=1 && k!=2))
	    {
	      TF1 *func = new TF1(Form("fit_%s",hVzDiff[k][j]->GetName()), "gaus(0)", -20, -10);
	      func->SetParameters(1e-3, 0, 20);
	      hVzDiff[k][j]->Fit(func, "R0");
	      func->SetRange(-20,20);
	      func->SetLineColor(2);
	      func->SetLineStyle(2);
	      func->Draw("sames");
	      if(k==0) leg2->AddEntry(func, "Fit background", "L");
	      TPaveText *t1 = GetPaveText(0.35, 0.55, 0.12, 0.2, 0.055);
	      t1->AddText(Form("f_{bkg}(|#Deltaz|<3) = %4.3f%%",func->Integral(-3,3)/hVzDiff[k][j]->Integral(hVzDiff[k][j]->FindBin(-3), hVzDiff[k][j]->FindBin(3)) * 100));
	      t1->Draw();
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  leg2->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_VzDiff.pdf",run_type));
}

//================================================
void makeVzDiff(const int saveHisto = 0)
{
  TFile *fin[2];
  fin[0] = TFile::Open("output/Run14_AuAu200.MB.VtxEff.prod_low.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.MB.VtxEff.prod_high.root", "read");
  THnSparseF *hn[2];
  for(int i=0; i<2; i++)
    {
      hn[i] = (THnSparseF*)fin[i]->Get("mhMbEvtEff");
      hn[i]->SetName(Form("%s_%d",hn[i]->GetName(),i));
    }
  
  TH1F *hTpcVz[nCentBins][gNTrgSetup];
  TH1F *hVzDiff[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTpcVz[k][j] = new TH1F(Form("hTpcVz_cent%s%s",cent_Title[k],gTrgSetupName[j]), ";vz_{TPC} (cm)", 400, -200, 200);
	  hVzDiff[k][j] = new TH1F(Form("hVzDiff_cent%s%s",cent_Title[k],gTrgSetupName[j]), ";vz_{TPC} (cm)", 200, -50, 50);
	}
    }

  TH2F *hTpcVzVsRun[2][nCentBins];
  TH2F *hVzDiffVsRun[2][nCentBins];
  for(int i=0; i<2; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hn[i]->GetAxis(5)->SetRange(centBins_low[k], centBins_high[k]);
	  hTpcVzVsRun[i][k] = (TH2F*)hn[i]->Projection(1, 0);
	  hTpcVzVsRun[i][k]->SetName(Form("hTpcVzVsRun_cent%s_%d",cent_Title[k],i));

	  hVzDiffVsRun[i][k] = (TH2F*)hn[i]->Projection(2, 0);
	  hVzDiffVsRun[i][k]->SetName(Form("hVzDiffVsRun_cent%s_%d",cent_Title[k],i));
	}
    }
	  
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  TH2F *hVzTmp = 0x0, *hVizDiffTmp = 0x0;
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=1; j<gNTrgSetup; j++)
	{
	  if(j<gNTrgSetup-1) 
	    {
	      hVzTmp = hTpcVzVsRun[0][k];
	      hVizDiffTmp = hVzDiffVsRun[0][k];
	    }
	  else  
	    {
	      hVzTmp = hTpcVzVsRun[1][k];
	      hVizDiffTmp = hVzDiffVsRun[1][k];
	    }
	  ifstream fruns;
	  fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s_2014.list",run_type,trgSetupName[j-1]));
	  int runnumber;
	  while(!fruns.eof())
	    {
	      fruns >> runnumber;
	      int xbin = hVzTmp->GetXaxis()->FindFixBin(runnumber);
	      TH1F *h1 = (TH1F*)hVzTmp->ProjectionY(Form("hVzTmp_cent%s_%d",cent_Title[k],runnumber),xbin,xbin);
	      hTpcVz[k][j]->Add(h1);

	      TH1F *h1 = (TH1F*)hVizDiffTmp->ProjectionY(Form("hVizDiffTmp_cent%s_%d",cent_Title[k],runnumber),xbin,xbin);
	      hVzDiff[k][j]->Add(h1);
	    }
	}
    }
  
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hTpcVz[k][0]->Add(hTpcVz[k][j]);
	  hVzDiff[k][0]->Add(hVzDiff[k][j]);
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "update");
      for(int k=0; k<nCentBins; k++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      hTpcVz[k][j]->Write("",TObject::kOverwrite);
	      hVzDiff[k][j]->Write("",TObject::kOverwrite);
	    }
	}
      fout->Close();
    }
}

//================================================
void jpsiWidth(const int savePlot = 0)
{
  TFile *fsmear = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkResScan.root"),"read");
  TH2F *hMassVsPtEmbed[nCentBins];
  TObjArray embedSlice[nCentBins];
  TH1F *hEmbSmearSigma[nCentBins];
  TF1 *funcSmearSigma[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hMassVsPtEmbed[k] = (TH2F*)fsmear->Get(Form("hRcJpsiMass_%s_final_def",cent_Title[k]));
      hMassVsPtEmbed[k]->FitSlicesY(0, 0, -1, 0, "QNR", &embedSlice[k]);
      hEmbSmearSigma[k] = (TH1F*)((TH1F*)embedSlice[k][2])->Clone(Form("hEmbSmearSigma_%s",cent_Title[k])); 
      funcSmearSigma[k] = new TF1(Form("hFitSmearTmp_%d",k),"pol4",0,15);
      hEmbSmearSigma[k]->Fit(funcSmearSigma[k],"IRQ0");
      funcSmearSigma[k]->SetLineColor(6);
      funcSmearSigma[k]->SetLineStyle(2);
    }     

  TH1F *hJpsiSigma[nCentBins][gNTrgSetup];
  TFile *fin = TFile::Open(Form("Rootfiles/%s.TrkResScan.input.root",run_type), "read");
  TCanvas *c = new TCanvas("cJpsiWidth", "cJpsiWidth", 1100, 700);
  c->Divide(3, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiSigma[k][j] = (TH1F*)fin->Get(Form("Jpsi_FitSigma_cent%s_weight%s",cent_Title[k],gTrgSetupName[j]));
	  hJpsiSigma[k][j]->SetMarkerStyle(20+j);
	  hJpsiSigma[k][j]->SetMarkerColor(color[j]);
	  hJpsiSigma[k][j]->SetLineColor(color[j]);
	  hJpsiSigma[k][j]->SetMarkerSize(1.2);
	  hJpsiSigma[k][j]->GetYaxis()->SetRangeUser(0,0.2);
	  if(k<4) hJpsiSigma[k][j]->GetXaxis()->SetRangeUser(0.5,10);
	  else    hJpsiSigma[k][j]->GetXaxis()->SetRangeUser(0.5,6);
	  if(j==0) hJpsiSigma[k][j]->Draw();
	  else     hJpsiSigma[k][j]->Draw("sames");

	  if(k==0) leg->AddEntry(hJpsiSigma[k][j], legName[j].Data(), "P");
	}
      funcSmearSigma[k]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_JpsiWidth.pdf",run_type));
}
