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
