TFile *f;
const int year = YEAR;
TString file_name;

//================================================
void sys_JpsiYield()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //YieldVsPt();
  YieldVsNpart();
}

//================================================
void YieldVsPt(int savePlot = 1, int saveHisto = 1)
{
 // re-assign global constants
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  // evaluate uncertainty
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  xbins[0] = 0.15;

  TH1F *hSys[nCentBins];
  TH1F *hMax[nCentBins][2];
  for(int i=0; i<nCentBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_signalExt_%s",cent_Title[i]),Form("Systematic uncertainty for signal extraction (%s%%)",cent_Name[i]),nbins,xbins);
      hMax[i][0] = new TH1F(Form("Sys_max_%s",cent_Title[i]),"",nbins,xbins);
      hMax[i][1] = new TH1F(Form("Sys_min_%s",cent_Title[i]),"",nbins,xbins);
    }

  const int nSys = 13;
  const char *sys_name[nSys]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc1","_BkgFunc2","_LargeFit","_SmallFit","_SigFunc","_FixSigDown","_FixSigUp","_LineShape"};
  const char *sys_leg[nSys] = {"Default","Larger bkg norm range","Smaller bkg norm range","Fit ME/SE w/ pol0","Binning","Res. bkg order-1","Res. bkg order+1","Larger sig fit range","Smaller sig fit range","Crystal-ball","Fix sig. down","Fix sig. up","line shape"};

  TString outName    = Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outNameSys = Form("Rootfiles/%s.Sys.JpsiYield.root",run_type);

  TFile *fdata = TFile::Open(Form("%s",outName.Data()),"read");
  TFile *fsys = 0x0;
  if(saveHisto) fsys = TFile::Open(Form("%s",outNameSys.Data()),"update");
  else          fsys = TFile::Open(Form("%s",outNameSys.Data()),"read");

  TH1F *hSignal[nCentBins][nSys+1];

  TCanvas *c = new TCanvas(Form("Sys_signalExt"),Form("Sys_signalExt"),1200,700);
  c->Divide(3,2);
  double max[nCentBins][nPtBins];
  for(int i=0; i<nCentBins; i++)
    {
      for(int k=0; k<nPtBins; k++)
	{
	  max[i][k] = 0;
	}
    }

  const char* method = "Fit";
  const int color_sys[nSys+1] = {1, 2, 3, 4, 6, 7, 8, 1, 2, 4, 6, 7, 8, 1};
  const int nSys_used = 14;
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nSys_used; j++)
	{
	  hSignal[i][j] = 0x0;
	  if(j<nSys)
	    {
	      if(j==0) hSignal[i][j] = (TH1F*)fdata->Get(Form("Jpsi_%sYield_cent%s_weight%s",method,cent_Title[i],sys_name[j]));
	      else hSignal[i][j] = (TH1F*)fsys->Get(Form("Jpsi_%sYield_cent%s_weight%s",method,cent_Title[i],sys_name[j]));
	    }
	  if(j==nSys)
	    {
	      hSignal[i][j] = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_cent%s_weight",cent_Title[i]));
	    }

	  if(!hSignal[i][j]) continue;
	  hSignal[i][j]->SetMarkerSize(1);
	  hSignal[i][j]->SetMarkerStyle(20+j);
	  hSignal[i][j]->SetMarkerColor(color_sys[j]);
	  
	  TH1F *htmp = (TH1F*)hSignal[i][j]->Clone(Form("%s_tmp",hSignal[i][j]->GetName()));
	  htmp->Divide(hSignal[i][0]);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      if(j==0) htmp->SetBinError(bin,hSignal[i][j]->GetBinError(bin)/hSignal[i][j]->GetBinContent(bin));
	      else     htmp->SetBinError(bin,0);
	    }
	  htmp->GetYaxis()->SetRangeUser(0.6,1.4);
	  htmp->SetTitle(";p_{T} (GeV/c);Relative difference");
	  if(i==2 || i==3) htmp->GetXaxis()->SetRangeUser(0.15,9);
	  if(i==4) htmp->GetXaxis()->SetRangeUser(0.15,6);
	  c->cd(i+2);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	  ScaleHistoTitle(htmp,0.05,1,0.04,0.05,1,0.04,62);
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("P sames");
	  TPaveText *t1 = GetTitleText(Form("Signal extraction uncertainty"),0.05);
	  t1->Draw();
	  t1 = GetPaveText(0.25,0.35,0.2,0.3,0.06,62);
	  t1->AddText(Form("%s%%",cent_Name[i]));
	  t1->Draw();
	  if(j<nSys_used)
	    {
	      for(int k=0; k<nPtBins; k++)
		{
		  double value = fabs(htmp->GetBinContent(k+1)-1);
		  if(max[i][k]<value) max[i][k] = value;
		}
	    }
	}

      if(i==0)
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
	  for(int j=0; j<nSys_used; j++)
	    {
	      if(j<nSys) leg[j/7]->AddEntry(hSignal[i][j],sys_leg[j],"P");
	      else leg[j/7]->AddEntry(hSignal[i][j],"Bin-counting","P");
	    }
	}
      
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  double sys = max[i][bin-1];
	  if(i==0 || i==1)
	    {
	      if(bin==9) sys = max[i][7];
	    }
	  if(i==2 || i==3)
	    {
	      if(bin==8) sys = max[i][6];
	    }
	  if(i==4)
	    {
	      if(bin==6) sys = max[i][4];
	    }
	  hSys[i]->SetBinContent(bin,sys+1);
	  htmp->SetBinContent(bin,1-sys);
	  hMax[i][0]->SetBinContent(bin,1+max[i][bin-1]);
	  hMax[i][1]->SetBinContent(bin,1-max[i][bin-1]);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->SetLineWidth(1);
      hSys[i]->SetLineStyle(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->SetLineWidth(1);
      htmp->SetLineStyle(2);
      htmp->Draw("samesHIST");
      for(int j=0; j<2; j++)
	{
	  hMax[i][j]->SetLineWidth(1);
	  //hMax[i][j]->SetLineStyle(2);
	  //hMax[i][j]->Draw("samesHIST");
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
  leg2->AddEntry(hSys[0],"Assigned uncertainty","L");
  //leg2->AddEntry(hSys[0],"Re-assigned uncert.","L");
  leg2->Draw();

  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiXsec/Sys_signalExt.pdf",run_type));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_SysSigExtVsPt.pdf"));
    }

  if(saveHisto)
    {
      fsys->cd();
      for(int i=0; i<nCentBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
  
}



//================================================
void YieldVsNpart(int savePlot = 1, int saveHisto = 1)
{
 // re-assign global constants
  const int nPtBins         = nPtBins_npart;
  const double* ptBins_low  = ptBins_low_npart;
  const double* ptBins_high = ptBins_high_npart;
  const char** pt_Name      = pt_Name_npart;
  const int* nCentBins      = nCentBins_npart; 
  const int* centBins_low   = centBins_low_npart;
  const int* centBins_high  = centBins_high_npart;
  const char** cent_Name    = cent_Name_npart;
  const char** cent_Title   = cent_Title_npart;
  const int kNCent          = nCentBins[0];

  // evaluate uncertainty
  TH1F *hSys[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_signalExt_pt%s",pt_Name[i]),Form("Systematic uncertainty for signal extraction (%s GeV/c)",pt_Name[i]),nCentBins[i],0,nCentBins[i]);
    }

  const int nSys = 14;
  const char *sys_name[nSys]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc1","_BkgFunc2","_LargeFit","_SmallFit","_SigFunc","_FixSigDown","_FixSigUp","_LineShape",""};
  const char *sys_leg[nSys] = {"Default","Larger bkg norm range","Smaller bkg norm range","Fit ME/SE w/ pol0","Binning","Res. bkg order-1","Res. bkg order+1","Larger sig fit range","Smaller sig fit range","Crystal-ball","Fix sig. down","Fix sig. up","Line shape","Bin-counting"};

  TString outName    = Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outNameSys = Form("Rootfiles/%s.Sys.JpsiYield.root",run_type);
  TFile *fdata = TFile::Open(Form("%s",outName.Data()),"read");
  TFile *fsys = 0x0;
  if(saveHisto) fsys = TFile::Open(Form("%s",outNameSys.Data()),"update");
  else          fsys = TFile::Open(Form("%s",outNameSys.Data()),"read");


  const int nSys_used = 14;
  const int color_sys[nSys_used] = {1, 2, 3, 4, 6, 7, 8, 1, 2, 4, 6, 7, 8, 1};  

  TH1F *hSignal[nPtBins][nSys];
  double max[nPtBins][kNCent];
  for(int i=0; i<nPtBins; i++)
    {
      for(int k=0; k<nCentBins[i]; k++)
	{
	  max[i][k] = 0;
	}
    }

  const char* method = "Fit";
  for(int i=0; i<nPtBins; i++)
    {
      TCanvas *c = new TCanvas(Form("Sys_signalExt_pt%s",pt_Name[i]),Form("Sys_signalExt_pt%s",pt_Name[i]),800,600);
      for(int j=0; j<nSys_used; j++)
	{
	  hSignal[i][j] = 0x0;
	  if(j==0) hSignal[i][j] = (TH1F*)fdata->Get(Form("Jpsi_%sYield_pt%s_weight%s",method,pt_Name[i],sys_name[j]));
	  else if(j<nSys-1) hSignal[i][j] = (TH1F*)fsys->Get(Form("Jpsi_%sYield_pt%s_weight%s",method,pt_Name[i],sys_name[j]));
	  else hSignal[i][j] = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_pt%s_weight",pt_Name[i]));
	  hSignal[i][j]->SetMarkerSize(1);
	  hSignal[i][j]->SetMarkerStyle(20+j);
	  hSignal[i][j]->SetMarkerColor(color_sys[j]);
	  
	  TH1F *htmp = (TH1F*)hSignal[i][j]->Clone(Form("%s_tmp",hSignal[i][j]->GetName()));
	  htmp->Divide(hSignal[i][0]);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      if(j==0) htmp->SetBinError(bin,hSignal[i][j]->GetBinError(bin)/hSignal[i][j]->GetBinContent(bin));
	      else     htmp->SetBinError(bin,0);
	    }
	  htmp->GetYaxis()->SetRangeUser(0.6,1.6);
	  htmp->SetTitle(";;Relative difference");
	  ScaleHistoTitle(htmp,0.04,1,0.045,0.04,1,0.035,62);
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("P sames");
	  TPaveText *t1 = GetTitleText(Form("%s: signal extraction uncertainty (%s)",run_type,pt_Title_npart[i]));
	  t1->Draw();
	  for(int k=0; k<nCentBins[i]; k++)
	    {
	      if(i==1 && k>=5 && j==9) continue;
	      if(i==0 && k==0 && j==5) continue;
	      double value = fabs(htmp->GetBinContent(k+1)-1);
	      if(max[i][k]<value) max[i][k] = value;
	    }
	}
        
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  double sys = max[i][bin-1];
	  hSys[i]->SetBinContent(bin,sys+1);
	  htmp->SetBinContent(bin,1-sys);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->SetLineWidth(1);
      hSys[i]->SetLineStyle(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->SetLineWidth(1);
      htmp->SetLineStyle(2);
      htmp->Draw("samesHIST");
      double xmin = 0.15, xmax = 0.45, ymin = 0.6, ymax = 0.88;
      TLegend *leg[2];
      for(int l=0; l<2; l++)
	{
	  leg[l] = new TLegend(xmin+0.35*l,ymin,xmax+0.35*l,ymax);
	  leg[l]->SetBorderSize(0);
	  leg[l]->SetFillColor(0);
	  leg[l]->SetTextSize(0.03);
	}
      for(int j=0; j<nSys_used; j++)
	{
	  leg[j/7]->AddEntry(hSignal[i][j],sys_leg[j],"P");
	  if(j==nSys) leg[j/7]->AddEntry(hSignal[i][j],"Bin-counting","P");
	}
      leg[1]->AddEntry(hSys[i],"Maximum deviation","L");
      leg[0]->Draw();
      leg[1]->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiXsec/Npart.Sys_signalExt_pt%s.pdf",run_type,pt_Name[i]));
	}
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_SysSigExtVsCent_Pt%1.0f.pdf",ptBins_low[i]));
	}
    }

  if(saveHisto)
    {
      fsys->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
  
}
