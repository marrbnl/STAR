const int year = 2014;

TFile *f;

//================================================
void ana_MtdTrigEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //getTrigEff();
  //anaTrigEff();
  sysTrigEff();
 
  //compareTacDiff();
  //lumiDepend();
  //trigElecEff();
  //trigUnitMap();
  //extrapolation();
}

//================================================
void anaTrigEff(const int savePlot = 0, const int saveHisto = 0)
{
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.18);

  if(year==2014)
    {
      const char* run_type ="Run14_AuAu200";
      const int nLumi = 2;
      const char *lumi_name[nLumi] = {"prod_high","prod_low"};
    }
  if(year==2015)
    {
      const char* run_type ="Run15_pAu200";
      const int nLumi = 1;
      const char *lumi_name[nLumi] = {"prod"};
    }

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type),"read");
  
  TGraphAsymmErrors *gTacEff[nLumi];
  TF1 *funcTacEff[nLumi];
  TH1F *hplot = new TH1F("hplot",";p_{T}^{#mu} [GeV/c];Efficiency",100,1,7);
  for(int l=0; l<nLumi; l++)
    {
      gTacEff[l] = (TGraphAsymmErrors*)fin->Get(Form("%s_gTacDiffEffFinal_BinCount_%s_Run15_pp200",run_type,lumi_name[l]));
      funcTacEff[l] = new TF1(Form("%s_gTacDiffEffFinalFit_BinCount_%s_Run15_pp200",run_type,lumi_name[l]),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
      if(run_type=="Run15_pAu200") funcTacEff[l]->FixParameter(0, 1);
      gTacEff[l]->Fit(funcTacEff[l],"IR0");
      TCanvas *c = new TCanvas(Form("%s_TacDiffEff_%s",run_type,lumi_name[l]),Form("%s_TacDiffEff_%s",run_type,lumi_name[l]),800,600);
      if(run_type=="Run14_AuAu200") hplot->GetYaxis()->SetRangeUser(0.6,1.1);
      if(run_type=="Run15_pAu200")  hplot->GetYaxis()->SetRangeUser(0.95,1.02);
      hplot->GetYaxis()->SetNdivisions(505);
      hplot->Draw();
      gTacEff[l]->SetMarkerStyle(21);
      gTacEff[l]->SetMarkerColor(1);
      gTacEff[l]->SetLineColor(1);
      gTacEff[l]->Draw("PEsame");
      funcTacEff[l]->SetLineColor(2);
      funcTacEff[l]->SetLineStyle(2);
      funcTacEff[l]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s: estimated MTD trigger efficiency",run_type));
      t1->Draw();
      TLegend *leg = new TLegend(0.45,0.15,0.65,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(lumi_name[l]);
      leg->AddEntry(gTacEff[l],"Data: bin counting","P");
      leg->AddEntry(funcTacEff[l],"Fit: p0-e^{-p1*(x-p2)}","L");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%s_TacDiffEffWithFit_%s.pdf",run_type,run_type,lumi_name[l]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_FitEff_%s.pdf",lumi_name[l]));
	}
    }
  if(saveHisto)
    {
      fin->cd();
      for(int l=0; l<nLumi; l++)
	{
	  funcTacEff[l]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void sysTrigEff(const int savePlot = 1, const int saveHisto = 1)
{
  if(year==2014)
    {
      const char* run_type ="Run14_AuAu200";
      const int nLumi = 2;
      const char *lumi_name[nLumi] = {"prod_high","prod_low"};
      const double min_TacDiffCut[nLumi] = {789,786};
      const double max_TacDiffCut[nLumi] = {837,837};
    }
  if(year==2015)
    {
      const char* run_type ="Run15_pAu200";
      const int nLumi = 1;
      const char *lumi_name[nLumi] = {"prod"};
      const double min_TacDiffCut[nLumi] = {873};
      const double max_TacDiffCut[nLumi] = {962};
    }

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type),"read");

  TGraphAsymmErrors *gTacEff[nLumi];
  TF1 *funcTacEff[nLumi];
  for(int l=0; l<nLumi; l++)
    {
      gTacEff[l] = (TGraphAsymmErrors*)fin->Get(Form("%s_gTacDiffEffFinal_BinCount_%s_Run15_pp200",run_type,lumi_name[l]));
      funcTacEff[l] = (TF1*)fin->Get(Form("%s_gTacDiffEffFinalFit_BinCount_%s_Run15_pp200",run_type,lumi_name[l]));
    }

  const int nbins = 6;
  TH1F *hplot = new TH1F("hplot",";p_{T} (GeV/c);Efficiency",100,0,10);
  hplot->GetXaxis()->SetRangeUser(1.3,7);
  if(run_type=="Run14_AuAu200") hplot->GetYaxis()->SetRangeUser(0.5,1);
  if(run_type=="Run15_pAu200")  hplot->GetYaxis()->SetRangeUser(0.95,1.02);

  // part I: systematic uncertainty
  TF1 *funcSysSys[nLumi][2];
  const int nSys = 7;
  if(run_type=="Run14_AuAu200")
    {
      const char* sys_name[nSys] = {"Deafult","AuAu reso. up","AuAu reso. down","pp reso. up","pp reso. down","AuAu TacDiff mean","Expected AuAu reso."};
    }
  if(run_type=="Run15_pAu200")
    {
      const char* sys_name[nSys] = {"Deafult","pAu reso. up","pAu reso. down","pp reso. up","pp reso. down","pAu TacDiff mean","Expected pAu reso."};
    }
  const int color[nSys] = {1, 2, 3, 4, 9, 6, kGreen+2};
  const char* data_name[2] = {run_type,"Run15_pp200"};
  TH1F *hMuonTacDiffMean[2];
  TH1F *hMuonTacDiffSigma[2];
  TF1  *funcMuonTacDiffSigma[2];
  for(int i=0; i<2; i++)
    {
      hMuonTacDiffMean[i] = (TH1F*)fin->Get(Form("%s_hMuonTacDiffMeanVsPt",data_name[i]));
      hMuonTacDiffSigma[i] = (TH1F*)fin->Get(Form("%s_hMuonTacDiffSigmaVsPt",data_name[i]));
      funcMuonTacDiffSigma[i] = (TF1*)fin->Get(Form("%s_fitMuonTacDiffSigmaVsPt",data_name[i]));
    }
  TH1F *ppTacDiffSigmaCL = (TH1F*)fin->Get("Run15_pp200_fitMuonTacDiffSigmaVsPtCL");
  TH1F *hppTacDiff[nbins];
  TF1  *hppTacDiffFit[nbins];
  for(int bin=0; bin<nbins; bin++)
    {
      hppTacDiff[bin]   = (TH1F*)fin->Get(Form("Run15_pp200_hMuonTacDiffCombined_PtBin%d",bin+1));
      hppTacDiffFit[bin] = (TF1*)fin->Get(Form("Run15_pp200_funcMuonTacDiffCombined_PtBin%d",bin+1));
    }
 
  TGraphAsymmErrors* gTacDiffEff[nLumi][nSys];
  double pt_arr[nbins], all_arr[nbins], all_err_arr[nbins], acc_arr[nbins], acc_err_arr[nbins];
  double x,y;
  TCanvas *cEff[nLumi];
  for(int l=0; l<nLumi; l++)
    {
      cEff[l] = new TCanvas(Form("cEffSys_%s",lumi_name[l]),Form("cEffSys_%s",lumi_name[l]),800,600);
      hplot->DrawCopy();
      TPaveText *t1 = GetTitleText(Form("%s: MTD trigger efficiency uncertainty estimation",run_type));
      t1->Draw();
      const int nleg = 2;
      TLegend *legSys[nleg];
      for(int ileg=0; ileg<nleg; ileg++)
	{
	  legSys[ileg] = new TLegend(0.25+ileg*0.3,0.15,0.5+ileg*0.3,0.4);
	  legSys[ileg]->SetBorderSize(0);
	  legSys[ileg]->SetFillColor(0);
	  legSys[ileg]->SetTextSize(0.035);
	  if(ileg==0) legSys[ileg]->SetHeader(lumi_name[l]);
	}
      for(int s=0; s<nSys; s++)
	{
	  for(int bin=0; bin<nbins; bin++)
	    {
	      int bin_index = bin;
	      if(bin_index==nbins-1) bin_index = bin-1;

	      double pt     = hMuonTacDiffMean[1]->GetBinCenter(bin_index+1);
	      double mean   = hMuonTacDiffMean[1]->GetBinContent(bin_index+1);
	      double AA_res = funcMuonTacDiffSigma[0]->Eval(pt);
	      double pp_res = funcMuonTacDiffSigma[1]->Eval(pt);
	      
	      if(s==1 || s==2)
		{
		  if(run_type=="Run14_AuAu200")
		    {
		      if(bin==0) 
			{
			  if(s==1) AA_res = ppTacDiffSigmaCL->GetBinContent(bin+1) + ppTacDiffSigmaCL->GetBinError(bin+1);
			  if(s==2) AA_res = funcMuonTacDiffSigma[0]->Eval(5);
			}
		      else       AA_res = funcMuonTacDiffSigma[0]->Eval(pt) + (3-2*s)*funcMuonTacDiffSigma[0]->GetParError(0);
		    }
		  if(run_type=="Run15_pAu200")
		    {
		      AA_res = funcMuonTacDiffSigma[0]->Eval(pt) + (3-2*s)*funcMuonTacDiffSigma[0]->GetParError(0);
		    }
		}
	      if(s==3 || s==4)
		{
		  pp_res = ppTacDiffSigmaCL->GetBinContent(bin+1) + (7-2*s)*ppTacDiffSigmaCL->GetBinError(bin+1);
		}
	      if(s==6)
		{
		  if(run_type=="Run14_AuAu200") AA_res = sqrt(pp_res*pp_res-8.27+3.33);
		  if(run_type=="Run15_pAu200") AA_res = sqrt(pp_res*pp_res-8.27+6.89);
		}
	      
	      double scale = pp_res/AA_res;
	      double minTacDiff = mean - (mean - min_TacDiffCut[l])*scale;
	      if(s==5)
		{
		  minTacDiff += mean - hMuonTacDiffMean[0]->GetBinContent(bin_index+1);
		}
	      //cout << lumi_name[l] << ": " << s << "  " << pt << "  " << mean << "  " << pp_res << "  " << AA_res << ":  " << min_TacDiffCut[l] << " -> " << minTacDiff << endl;

	      TH1F *hppAuTac = hppTacDiff[bin];
	      TF1  *hppAuTacFit = hppTacDiffFit[bin];
	      int xbins = hppAuTac->GetNbinsX(); 
	      gTacEff[l]->GetPoint(bin, x, y);
	      pt_arr[bin] = x;
	      all_arr[bin] = hppAuTac->IntegralAndError(1, xbins, all_err_arr[bin]);
	      int low_bin =  hppAuTac->FindFixBin(minTacDiff);
	      acc_arr[bin] = hppAuTac->IntegralAndError(low_bin+1,
	      						hppAuTac->FindFixBin(max_TacDiffCut[l]-0.1),
	      						acc_err_arr[bin]);
	      double fraction = (hppAuTac->GetXaxis()->GetBinUpEdge(low_bin)-minTacDiff)/hppAuTac->GetBinWidth(low_bin);
	      acc_arr[bin] += hppAuTac->GetBinContent(low_bin)*fraction;
	      acc_err_arr[bin] = sqrt(pow(acc_err_arr[bin],2)+pow(hppAuTac->GetBinError(low_bin)*fraction,2));
	      if(acc_arr[bin]>all_arr[bin]) acc_arr[bin] = all_arr[bin];

	      if(year==2014)
		{
	      // recaculate central value using fitted function to avoid statistical fluctuation
		  all_arr[bin]  = hppAuTacFit->Integral(760, max_TacDiffCut[l]);
		  acc_arr[bin]  = hppAuTacFit->Integral(minTacDiff, max_TacDiffCut[l]);
		  
		}
	    }
	  gTacDiffEff[l][s] = GetEfficiencyCurve(nbins, pt_arr, all_arr, all_err_arr, acc_arr, acc_err_arr);
	  gTacDiffEff[l][s]->SetName(Form("%s_gTacDiffEffSys%d_%s",run_type,s,lumi_name[l]));
	  int npoints = gTacDiffEff[l][s]->GetN();
	  for(int ipoint=0; ipoint<npoints; ipoint++)
	    {
	      gTacDiffEff[l][s]->SetPointEXhigh(ipoint,0);
	      gTacDiffEff[l][s]->SetPointEXlow(ipoint,0);
	      if(s>0)
		{
		  gTacDiffEff[l][s]->SetPointEYhigh(ipoint,0);
		  gTacDiffEff[l][s]->SetPointEYlow(ipoint,0);
		}
	      else
		{
		  gTacDiffEff[l][s]->SetPointEYhigh(ipoint,gTacEff[l]->GetErrorYhigh(ipoint));
		  gTacDiffEff[l][s]->SetPointEYlow(ipoint,gTacEff[l]->GetErrorYlow(ipoint));
		}
	    }
	  gTacDiffEff[l][s]->SetMarkerStyle(20+s);
	  gTacDiffEff[l][s]->SetMarkerColor(color[s]);
	  gTacDiffEff[l][s]->Draw("PEsames");
	  legSys[s/4]->AddEntry(gTacDiffEff[l][s],sys_name[s],"P");
	}
      for(int ileg=0; ileg<nleg; ileg++)
	legSys[ileg]->Draw();
    }

  TGraphErrors* gTacEffSysLimit[nLumi][2];
  double y1;
  for(int l=0; l<nLumi; l++)
    {
      for(int i=0; i<2; i++) gTacEffSysLimit[l][i] = new TGraphErrors(nbins);
      for(int bin=0; bin<nbins; bin++)
	{
	  gTacEff[l]->GetPoint(bin, x, y1);
	  double central = funcTacEff[l]->Eval(x);
	  double error = 0;
	  for(int s=1; s<nSys; s++)
	    {
	      gTacDiffEff[l][s]->GetPoint(bin, x, y);
	      if(fabs(y1-y)>error) error = fabs(y1-y);
	    }
	  gTacEffSysLimit[l][0]->SetPoint(bin, x, central-error);
	  gTacEffSysLimit[l][0]->SetPointError(bin, 0, 0.05);
	  gTacEffSysLimit[l][1]->SetPoint(bin, x, central+error);
	  gTacEffSysLimit[l][1]->SetPointError(bin, 0, 0.05);	  
	}
      
      for(int i=0; i<2; i++)
	{
	  funcSysSys[l][i] = new TF1(Form("%s_Sys1_%d",gTacEff[l]->GetName(),i),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	  gTacEffSysLimit[l][i]->Fit(funcSysSys[l][i],"IR0Q");
	  gTacEffSysLimit[l][i]->SetMarkerStyle(29);
	  gTacEffSysLimit[l][i]->SetMarkerColor(2);
	}
      if(year==2014)
	{
	  cEff[l]->cd();
	  for(int i=0; i<2; i++)
	    {
	      funcSysSys[l][i]->SetLineStyle(2);
	      funcSysSys[l][i]->SetLineColor(4);
	      funcSysSys[l][i]->Draw("sames");
	    }
	}
     if(savePlot)
	cEff[l]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/Sys.%s_TacDiffEffSysSys_%s.pdf",run_type,run_type,lumi_name[l]));
     if(gSaveAN)
       {
	 cEff[l]->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_EffTrig_Sys1_%s.pdf",lumi_name[l]));
       }
    }
 
  // part II: statisitcal error
  TF1 *funcSysStat[nLumi][3];
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
  const int nexpr = 100;
  TGraphAsymmErrors *gEffRndm[nLumi][nexpr];
  double x,y;
  for(int l=0; l<nLumi; l++)
    {
      for(int j = 0; j < nexpr; j++)
	{
	  TGraphAsymmErrors *gBaseEff = gTacEff[l]; // bin counting method 
	  int npoint = gBaseEff->GetN();
	  gEffRndm[l][j] = new TGraphAsymmErrors(npoint);
	  gEffRndm[l][j]->SetName(Form("%s_Rndm%d",gBaseEff->GetName(),j));
	  double *exh = gBaseEff->GetEXhigh();
	  double *exl = gBaseEff->GetEXlow();
	  double *eyh = gBaseEff->GetEYhigh();
	  double *eyl = gBaseEff->GetEYlow();
	  for(int ipoint=0; ipoint<npoint; ipoint++)
	    {
	      gBaseEff->GetPoint(ipoint,x,y);
	      double sigma = eyh[ipoint]>eyl[ipoint] ? eyh[ipoint] : eyl[ipoint];
	      y = rndm->Gaus(y,sigma);
	      gEffRndm[l][j]->SetPoint(ipoint,x,y);
	      gEffRndm[l][j]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
	    }
	}
    }
  const int nbins_rndm = 11;
  double xbins_rndm[nbins_rndm] = {1.4,1.5,1.75,2,2.25,2.5,3,4,5,6,8};
  TH1F *hSpread[nLumi][nbins_rndm];
  TGraphErrors *gSysStat[nLumi][3];
  const char *limit_name[3] = {"center","low","up"};
  for(int l=0; l<nLumi; l++)
    {
      c = new TCanvas(Form("cEffStat_%s",lumi_name[l]),Form("cEffStat_%s",lumi_name[l]),800,600);
      hplot->DrawCopy();
      TPaveText *t1 = GetTitleText(Form("%s: trigger efficiency uncertainty from statistical error",run_type));
      t1->Draw();
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("%s",lumi_name[l]));

      TCanvas *cFit = new TCanvas(Form("cEffFit_%s",lumi_name[l]),Form("cEffFit_%s",lumi_name[l]),1100,700);
      cFit->Divide(4,3);
      for(int k=0; k<nbins_rndm; k++)
	{
	  hSpread[l][k] = new TH1F(Form("hSpread_%s_%d",lumi_name[l],k),Form("hSpread_%s_%d",lumi_name[l],k),600,0,1.2);
	}

      for(int j = 0; j < nexpr; j++)
	{
	  TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[l][j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
	  funcTmp->SetParameter(0,1);
	  funcTmp->SetParLimits(0,0,1);
	  funcTmp->SetParLimits(1,0,5);
	  funcTmp->SetParLimits(2,0,1);
	  gEffRndm[l][j]->Fit(funcTmp,"R0Q");
	  c->cd();
	  funcTmp->SetLineStyle(2);
	  funcTmp->Draw("sames");
	  if(j==0) leg->AddEntry(funcTmp,"Randomization","L");
	  for(int k=0; k<nbins_rndm; k++)
	    {
	      hSpread[l][k]->Fill(funcTmp->Eval(xbins_rndm[k]));
	    }
	}

      gSysStat[l][0] = new TGraphErrors(nbins_rndm);
      gSysStat[l][1] = new TGraphErrors(nbins_rndm);
      gSysStat[l][2] = new TGraphErrors(nbins_rndm);
      for(int k=0; k<nbins_rndm; k++)
	{
	  cFit->cd(k+1);
	  double pt    = xbins_rndm[k];
	  double mean  = hSpread[l][k]->GetMean();
	  double width = hSpread[l][k]->GetRMS();
	  double error = hSpread[l][k]->GetRMSError();
	  hSpread[l][k]->GetXaxis()->SetRangeUser(mean-6*width,mean+8*width);
	  hSpread[l][k]->SetTitle(";Efficiency;");
	  hSpread[l][k]->SetMaximum(1.3*hSpread[l][k]->GetMaximum());
	  hSpread[l][k]->Draw();
	  TPaveText *t1 = GetTitleText(Form("p_{T} = %1.1f GeV/c\n",xbins_rndm[k]),0.06);
	  t1->Draw();

	  TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[l][k]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,1);
	  funcTmp2->SetParameter(0,hSpread[l][k]->GetMaximum());
	  funcTmp2->SetParameter(1,hSpread[l][k]->GetMean());
	  funcTmp2->SetParameter(2,hSpread[l][k]->GetRMS());
	  hSpread[l][k]->Fit(funcTmp2,"R0Q");
	  funcTmp2->SetLineColor(4);
	  funcTmp2->SetLineStyle(2);
	  funcTmp2->Draw("sames");
	  
	  mean = funcTmp2->GetParameter(1);
	  error = funcTmp2->GetParError(2);
	  width = funcTmp2->GetParameter(2);
	  gSysStat[l][0]->SetPoint(k,pt,mean);
	  gSysStat[l][0]->SetPointError(k,0,error);
	  gSysStat[l][1]->SetPoint(k,pt,mean-width);
	  gSysStat[l][1]->SetPointError(k,0,error);
	  gSysStat[l][2]->SetPoint(k,pt,mean+width);
	  gSysStat[l][2]->SetPointError(k,0,error);

	  TPaveText *t1 = GetPaveText(0.6,0.7,0.7,0.85,0.06);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("RMS=%4.3f",hSpread[l][k]->GetRMS()));
	  t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcTmp2->GetParameter(2))));
	  t1->Draw();
	}
      cFit->cd(12);
      TPaveText *t1 = GetPaveText(0.3,0.4,0.3,0.6,0.08);
      t1->AddText(run_type);
      t1->AddText(lumi_name[l]);
      t1->Draw();
      if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/Sys.%s_TacDiffEff_RndmFit_%s.pdf",run_type,run_type,lumi_name[l]));

      c->cd();
      leg->Draw();
      for(int t=0; t<3; t++)
	{
	  gSysStat[l][t]->SetMarkerStyle(20);
	  gSysStat[l][t]->SetMarkerColor(2+2*t);
	  gSysStat[l][t]->SetLineColor(2+2*t);
	  gSysStat[l][t]->Draw("samesP");

	  funcSysStat[l][t] = new TF1(Form("%s_Sys2_%s",gTacEff[l]->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	  funcSysStat[l][t]->SetParameters(0.95,2,0.2);
	  gSysStat[l][t]->Fit(funcSysStat[l][t],"R0Q");
	  funcSysStat[l][t]->SetLineColor(gSysStat[l][t]->GetMarkerColor());
	  if(t>0) funcSysStat[l][t]->SetLineStyle(2);
	  funcSysStat[l][t]->DrawCopy("sames");
	}
      leg->AddEntry(gSysStat[l][0],"Central value","PL");
      leg->AddEntry(gSysStat[l][1],"Lower limit","PL");
      leg->AddEntry(gSysStat[l][2],"Upper limit","PL");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/Sys.%s_TacDiffEff_RndmLimits_%s.pdf",run_type,run_type,lumi_name[l]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_EffTrig_Sys2_%s.pdf",lumi_name[l]));
	}
    }

  // combine the systematic uncertainty
  const char *sys_type[2] = {"down", "up"};
  TGraphErrors *gTacEffSysAll[nLumi][2];
  TF1 *funcSysAll[nLumi][2];
  TCanvas *cFit = new TCanvas("cFit_SysAll","cFit_SysAll",1100,700);
  cFit->Divide(2,2);
  for(int l=0; l<nLumi; l++)
    {
      double *exh = gTacEff[l]->GetEXhigh();
      double *exl = gTacEff[l]->GetEXlow();
      double *eyh = gTacEff[l]->GetEYhigh();
      double *eyl = gTacEff[l]->GetEYlow();
      for(int i=0; i<2; i++)
	{
	  int npoint = gTacEff[l]->GetN();
	  gTacEffSysAll[l][i] = new TGraphErrors(npoint);
	  gTacEffSysAll[l][i]->SetName(Form("%s_AllSys%d",gTacEff[l]->GetName(),i));
	  for(int ipoint=0; ipoint<npoint; ipoint++)
	    {
	      gTacEff[l]->GetPoint(ipoint,x,y);
	      double sys1 = fabs(funcSysSys[l][i]->Eval(x) - funcTacEff[l]->Eval(x));
	      double sys2 = fabs(funcSysStat[l][i+1]->Eval(x) - funcTacEff[l]->Eval(x));
	      double sys_all = sqrt(sys1*sys1+sys2*sys2);
	      double new_y = funcTacEff[l]->Eval(x) + sys_all*(i*2-1);
	      gTacEffSysAll[l][i]->SetPoint(ipoint,x,new_y);
	      gTacEffSysAll[l][i]->SetPointError(ipoint,exl[ipoint],eyl[ipoint]/y*new_y);
	      cout << x << "  " << sys1 << "  " << sys2 << endl;
	    }
	  gTacEffSysAll[l][i]->GetYaxis()->SetRangeUser(0.1,1.05);
	  cFit->cd(l*2+i+1);
	  gPad->SetGridy();
	  funcSysAll[l][i] = new TF1(Form("%s_gTacDiffEffFinal_%s_Sys%s",run_type,lumi_name[l],sys_type[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
	  if(run_type=="Run15_pAu200")funcSysAll[l][i]->SetParameters(1,4.5,0.4);
	  gTacEffSysAll[l][i]->Fit(funcSysAll[l][i],"R0Q");
	  hplot->DrawCopy();
	  gTacEffSysAll[l][i]->SetMarkerStyle(25);
	  gTacEffSysAll[l][i]->Draw("samesPE");
	  funcSysAll[l][i]->SetLineStyle(2);
	  funcSysAll[l][i]->SetLineColor(4);
	  funcSysAll[l][i]->Draw("sames");
	}
    }

  for(int l=0; l<nLumi; l++)
    {
      TCanvas *c = new TCanvas(Form("cEff_%s",lumi_name[l]),Form("cEff_%s",lumi_name[l]),800,600);
      hplot->DrawCopy();
      TPaveText *t1 = GetTitleText(Form("%s: MTD trigger efficiency uncertainty",run_type));
      t1->Draw();
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("%s",lumi_name[l]));
      
      funcTacEff[l]->SetLineColor(2);
      funcTacEff[l]->Draw("sames");
      gTacEff[l]->Draw("samesPE");
      for(int i=0; i<2; i++)
	{
	  funcSysAll[l][i]->Draw("sames");
	}
      leg->AddEntry(gTacEff[l],"Data: bin counting","P");
      leg->AddEntry(funcTacEff[l],"Fit to data","L");
      leg->AddEntry(funcSysAll[l][0],"Sys. Uncert.","L");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/Sys.%s_TacDiffEffWithSys_%s.pdf",run_type,run_type,lumi_name[l]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_EffTrig_SysAll_%s.pdf",lumi_name[l]));
	}
    }

  TF1 *funcEffFinal = new TF1(Form("%s_gTacDiffEffFinal_prod",run_type),"[0]-exp(-1*[1]*(x-[2]))",1.2,10);
  TF1 *funcEffSysFinal[2];
  for(int i=0; i<2; i++) funcEffSysFinal[i] = new TF1(Form("%s_Muon_TacDiffEff_Sys%s",run_type,sys_type[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,10);
  if(run_type=="Run15_pAu200")
    {
      funcEffFinal->SetParameters(funcTacEff[0]->GetParameters());
      for(int i=0; i<2; i++) funcEffSysFinal[i]->SetParameters(funcSysAll[0][i]->GetParameters());
    }
  if(run_type=="Run14_AuAu200")
    {
      // Get the average efficiency
      // MB statistics weight: prod_mid/high (80.75%); prod_/low (19.25%)
      const double weight[nLumi] = {0.8075, 0.1925};
      TF1 *hprodlow = 0x0, *hprodhigh = 0x0, *havg = 0x0;
      for(int i=0; i<3; i++)
	{
	  if(i==0)
	    {
	      hprodhigh = funcTacEff[0];
	      hprodlow = funcTacEff[1];
	      havg = funcEffFinal;
	    }
	  else
	    {
	      hprodhigh = funcSysAll[0][i-1];
	      hprodlow = funcSysAll[1][i-1];
	      havg = funcEffSysFinal[i-1];
	    }
	  TH1F *htmp = new TH1F(Form("htmp_%d",i),"",70,1,8);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      double pt = htmp->GetBinCenter(bin);
	      htmp->SetBinContent(bin, weight[0]*hprodhigh->Eval(pt)+weight[1]*hprodlow->Eval(pt));
	      htmp->SetBinError(bin, 1e-10);
	    }
	  TF1 *funcTmp = new TF1(Form("funcTmp_%d",i),"[0]-exp(-1*[1]*(x-[2]))",1.3,10);
	  funcTmp->SetParameters(hprodlow->GetParameters());
	  htmp->GetYaxis()->SetRangeUser(0.5,1);
	  htmp->SetMarkerStyle(24);
	  htmp->Fit(funcTmp,"IR0Q");
	  havg->SetParameters(funcTmp->GetParameters());
	}
      c = new TCanvas("MtdTacDiffEff_AvgFinal","MtdTacDiffEff_AvgFinal",800,600);
      gPad->SetGridy();
      hplot->GetXaxis()->SetRangeUser(1.3,10);
      hplot->DrawCopy();
      TPaveText *t1 = GetTitleText(Form("%s: MTD trigger efficiency",run_type));
      t1->Draw();
      funcEffFinal->Draw("sames");
      for(int i=0; i<2; i++) 
	{
	  funcEffSysFinal[i]->SetLineColor(4);
	  funcEffSysFinal[i]->SetLineStyle(2);
	  funcEffSysFinal[i]->Draw("sames");
	}
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(funcEffFinal, "Weighted efficiency", "L");
      leg->AddEntry(funcEffSysFinal[0], "Systematic uncertainty", "L");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/Sys.%s_AvgTacDiffEffWithSys.pdf",run_type,run_type));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_EffTrig_SysFinal.pdf"));
	}
    }

  if(saveHisto)
    {
      fin->cd();
      funcEffFinal->Write("",TObject::kOverwrite);
      fin->Close();

      TFile *fout = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"recreate");
      funcEffFinal->SetName(Form("%s_Muon_TacDiffEff",run_type));
      funcEffFinal->Write("",TObject::kOverwrite);
      for(int i=0; i<2; i++)
	{
	  funcEffSysFinal[i]->Write("",TObject::kOverwrite);
	}
      funcTacEff[0]->SetName(Form("%s_Muon_TacDiffEff_prod_high",run_type));
      funcTacEff[0]->Write("",TObject::kOverwrite);
      funcTacEff[1]->SetName(Form("%s_Muon_TacDiffEff_prod_low",run_type));
      funcTacEff[1]->Write("",TObject::kOverwrite);
      fout->Close();
    }
  
}

//================================================
void getTrigEff(const int savePlot = 1, const int saveHisto = 1)
{
  // The MTD trigger efficiency is calculated using Run15 pp 200 GeV data
  // as the baseline

  if(year==2014)
    {
      const int nLumi = 2;
      const char *lumi_name[nLumi] = {"prod_high","prod_low"};
      const double min_TacDiffCut[nLumi] = {788+1,785+1};
      const double max_TacDiffCut[nLumi] = {837,837};
      const char* config = "";
      const TString data_name = Form("Run%d_AuAu200",year-2000);
      const char *data_title  = Form("Run%d Au+Au @ 200 GeV",year-2000);
      const double min_fit = 760;
      const double max_fit = 830;
      const int nBinsTacDiff = 22;
      const double xBinsTacDiff[nBinsTacDiff+1] = {760,765,770,775,780,782,784,786,788,789,793,797,801,805,809,813,817,821,825,829,833,837,841};
    }
  else if(year==2015)
    {
      const int nLumi = 1;
      const char *lumi_name[nLumi] = {"prod"};
      const double min_TacDiffCut[nLumi] = {872+1};
      const double max_TacDiffCut[nLumi] = {962};
      const char* config = "";
      const TString data_name = Form("Run%d_pAu200",year-2000);
      const char *data_title  = Form("Run%d p+Au @ 200 GeV",year-2000);
      const double min_fit = 860;
      const double max_fit = 960;
      const int nBinsTacDiff = 23;
      const double xBinsTacDiff[nBinsTacDiff+1] = {860,865,870,871,872,873,875,880,885,890,895,900,905,910,915,920,925,930,935,940,945,950,955,960};
    }
  else if(year==2016)
    {
      const int nLumi = 1;
      const char *lumi_name[nLumi] = {"prod"};
      const double min_TacDiffCut[nLumi] = {950+1};
      const double max_TacDiffCut[nLumi] = {1005};
      const char* config = "";
      const TString data_name = Form("Run%d_AuAu200",year-2000);
      const char *data_title  = Form("Run%d Au+Au @ 200 GeV",year-2000);
    }

  const int nTrigUnit = 28;
  const int nbins = 7;
  const double lowbins[nbins] = {1.3, 1.3, 1.5, 2.0, 2.5, 3.0, 5.0};
  const double upbins[nbins]  = {10,  1.5, 2.0, 2.5, 3.0, 5.0, 10.0};
  const int nPtBins = nbins -1;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 

  // distributions from pp and p/AuAu events
  // 0 - pp; 1 - pAu/AuAu
  const char* ref_data_name = "Run15_pp200";
  const int nData = 2;
  const char* legName[nData] = {ref_data_name, data_name.Data()};
  TH2F *hMuonTacDiffVsTrigUnit[nData][nbins];
  TH1F *hMuonTacDiffInTrigUnit[nData][nbins][nTrigUnit];
  TH1F *hMuonTacDiffMeanVsTrigUnit[nData][nbins];

  TH2F *hLSBkgTacDiffVsTrigUnit[nData][nbins];
  TH1F *hLSBkgTacDiffMeanVsTrigUnit[nData][nbins];

  TH1F *hMeanPt[2];
  
  TFile *fdata[nData];
  for(int i=0; i<nData; i++)
    {
      if(i==0) fdata[i] = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",legName[i]),"read");
      if(i==1)
	{
	  if(saveHisto) fdata[i] = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",legName[i]),"update");
	  else          fdata[i] = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",legName[i]),"read");
	}
      for(int bin=0; bin<nbins; bin++)
	{
	  hMuonTacDiffVsTrigUnit[i][bin] = (TH2F*)fdata[i]->Get(Form("%s_hTacDiffVsTrigUnit_Muon_PtBin%d",legName[i],bin));
	  hMuonTacDiffMeanVsTrigUnit[i][bin] = (TH1F*)fdata[i]->Get(Form("%s_hTacDiffMeanVsTrigUnit_Muon_PtBin%d",legName[i],bin));
	  for(int k=0; k<nTrigUnit; k++)
	    {
	      hMuonTacDiffInTrigUnit[i][bin][k] = (TH1F*)hMuonTacDiffVsTrigUnit[i][bin]->ProjectionY(Form("%s_hTacDiffVsTrigUnit_Muon_TrigUnit%d_PtBin%d",legName[i],k+1,bin),k+2,k+2);
	      if(hMuonTacDiffInTrigUnit[i][bin][k]->Integral()>=1)
		hMuonTacDiffInTrigUnit[i][bin][k]->Scale(1./hMuonTacDiffInTrigUnit[i][bin][k]->Integral());
	    }

	  hLSBkgTacDiffVsTrigUnit[i][bin] = (TH2F*)fdata[i]->Get(Form("%s_hTacDiffVsTrigUnit_LS_PtBin%d",legName[i],bin));
	  hLSBkgTacDiffMeanVsTrigUnit[i][bin] = (TH1F*)fdata[i]->Get(Form("%s_hTacDiffMeanVsTrigUnit_LS_PtBin%d",legName[i],bin));
	}
      hMeanPt[i] = (TH1F*)fdata[i]->Get(Form("%s_hMeanPt_LS",legName[i]));
    }

  // compare multiplicity distribution per trigger unit
  TH1F *hMultWeight[nData][nbins];
  c = new TCanvas("Multiplicity_vs_TrigUnit","Multiplicity_vs_TrigUnit",1000,600);
  c->Divide(4,2);
  TLegend *leg = new TLegend(0.2,0.3,0.5,0.7);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  for(int bin=0; bin<nbins; bin++)
    {
      c->cd(bin+2);
      SetPadMargin(gPad,0.13,0.13,0.02,0.1);
      for(int i=0; i<nData; i++)
	{
	  hMultWeight[i][bin] = (TH1F*)hLSBkgTacDiffVsTrigUnit[i][bin]->ProjectionX(Form("Data%d_Multiplicity_PtBin%d",i,bin));
	  hMultWeight[i][bin]->Scale(1./hMultWeight[i][bin]->Integral());
	  hMultWeight[i][bin]->SetMarkerSize(1.2);
	  hMultWeight[i][bin]->SetMarkerStyle(20+i*4);
	  hMultWeight[i][bin]->SetMarkerColor(2-i);
	  hMultWeight[i][bin]->SetLineColor(2-i);
	  hMultWeight[i][bin]->SetMaximum(1.5*hMultWeight[i][bin]->GetMaximum());
	  hMultWeight[i][bin]->SetTitle(";TrigUnit;Multiplicity");
	  ScaleHistoTitle(hMultWeight[i][bin], 0.05, 0.9, 0.045, 0.05, 1.2, 0.04);
	  if(i==0) hMultWeight[i][bin]->Draw();
	  else     hMultWeight[i][bin]->Draw("sames");
	  if(bin==0) if(i<nData) leg->AddEntry(hMultWeight[i][bin], Form("%s: LS",legName[i]),"P");
	}
      TPaveText *t1 = GetTitleText(Form("%2.1f < p_{T} < %2.1f GeV/c",lowbins[bin],upbins[bin]),0.07);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sComp_MultVsTrigUnit.pdf",data_name.Data(),config));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_LSMultVsTrigUnit.pdf"));
    }

  // compare line-up
  for(int i=0; i<nData; i++)
    {
      c = new TCanvas(Form("%s_TacDiffMeanVsTrigUnit",legName[i]),Form("%s_TacDiffMeanVsTrigUnit",legName[i]),1100,700);
      c->Divide(3,2);
      for(int bin=1; bin<nbins; bin++)
	{
	  TH1F *h1tmp = (TH1F*)hLSBkgTacDiffMeanVsTrigUnit[i][bin]->Clone(Form("hMean_%d_%d",i,bin));
	  h1tmp->SetMarkerStyle(20);
	  h1tmp->SetMarkerColor(1);
	  h1tmp->SetLineColor(1);
	  h1tmp->SetTitle(";TrigUnit;<#DeltaTacSum>");
	  ScaleHistoTitle(h1tmp, 0.05, 0.9, 0.04, 0.05, 1.2, 0.04);
	  if(legName[i]=="Run15_pp200")   h1tmp->GetYaxis()->SetRangeUser(915, 935);	  
	  if(i==1 && data_name=="Run14_AuAu200") h1tmp->GetYaxis()->SetRangeUser(785, 810);
	  if(i==1 && data_name=="Run15_pAu200") h1tmp->GetYaxis()->SetRangeUser(900, 920);
	  c->cd(bin);
	  SetPadMargin(gPad,0.13,0.13,0.02,0.1);
	  h1tmp->Draw();
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]),0.06);
	  t1->Draw();
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.5,0.75,0.8,0.05);
      t1->AddText(legName[i]);
      t1->SetTextColor(2);
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%s%s_LSBkg_TacDiffMeanVsTrigUnit.pdf",data_name.Data(),config,legName[i]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_LSTacDiffMeanVsTrigUnit_%s.pdf",legName[i]));
	}
    }

  if(data_name=="Run15_pAu200")
    {
      //==============================================
      // Compare TacSum in each trigger unit: pp vs. pAu
      TH1F *hMuonTacDiffShift[nData][nbins][nTrigUnit];
      for(int bin=0; bin<nbins; bin++)
	{
	  c = new TCanvas(Form("ppVspAu_TacDiff_PtBin%d",bin),Form("ppVspAu_TacDiff_PtBin%d",bin),1100,650);
	  c->Divide(6,5);
	  TLegend *leg = new TLegend(0.15,0.3,0.5,0.88);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetHeader(Form("%1.1f < p_{T} < %1.1f GeV/c",lowbins[bin],upbins[bin]));
	  leg->SetTextSize(0.09);
	  for(int j=0; j<nTrigUnit; j++)
	    {
	      c->cd(j+1);
	      SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	      for(int i=0; i<nData; i++)
		{
		  hMuonTacDiffShift[i][bin][j] = new TH1F(Form("%s_MuonTacDiffShift_TrigUnit%d_PtBin%d",legName[i],j+1,bin),"",100,-50,50);
		  int nxbins = hMuonTacDiffShift[i][bin][j]->GetNbinsX();
		  double shift = hMuonTacDiffMeanVsTrigUnit[i][bin]->GetBinContent(j+1);
		  for(int ibin=1; ibin<=nxbins; ibin++)
		    {
		      int jbin = hMuonTacDiffInTrigUnit[i][bin][j]->FindFixBin(hMuonTacDiffShift[i][bin][j]->GetBinCenter(ibin) + shift);
		      hMuonTacDiffShift[i][bin][j]->SetBinContent(ibin, hMuonTacDiffInTrigUnit[i][bin][j]->GetBinContent(jbin));
		      hMuonTacDiffShift[i][bin][j]->SetBinError(ibin, hMuonTacDiffInTrigUnit[i][bin][j]->GetBinError(jbin));
		    }
	
		  TH1F *htmp = (TH1F*)hMuonTacDiffShift[i][bin][j]->Clone(Form("Plot_%s",hMuonTacDiffShift[i][bin][j]->GetName()));
		  htmp->Rebin(5);
		  htmp->SetTitle(";#DeltaTacSum;");
		  ScaleHistoTitle(htmp, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
		  htmp->SetMaximum(2*htmp->GetMaximum());
		  htmp->SetMarkerStyle(20+i*4);
		  int color = i + 1;
		  htmp->SetMarkerColor(color);
		  htmp->SetLineColor(color);
		  if(i==0) htmp->Draw();
		  else     htmp->Draw("samesPE");
		  if(j==0)
		    leg->AddEntry(htmp, legName[i], "P");
		}
	      TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	      t1->Draw();
	    }
	  c->cd(29);
	  leg->Draw();
	  if(savePlot) 
	    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sCompTacDiff_PtBin%d.pdf",data_name.Data(),config,bin));
	}
    }

  //==============================================
  // Calculate efficiency
  // 1) align the trigger units in pp according to pAu/AuAu alignment
  // 2) combine all the trigger units
  // 3) fit the distribution and determine the efficiency

  TH1F *hMuonTacDiffCombined[nbins];
  for(int bin=0; bin<nbins; bin++)
    {
      hMuonTacDiffCombined[bin] = new TH1F(Form("%s_hMuonTacDiffCombined_PtBin%d",legName[0],bin),"",300,700,1000);
      for(int j=0; j<nTrigUnit; j++)
	{
	  int index = bin;
	  if(index == nbins-1) index = nbins-2; 
	  double bkg_mean_AuAu = hLSBkgTacDiffMeanVsTrigUnit[1][index]->GetBinContent(j+1);
	  double bkg_mean_ref  = hLSBkgTacDiffMeanVsTrigUnit[0][index]->GetBinContent(j+1);
	  double shift = bkg_mean_ref - bkg_mean_AuAu;

	  TH1F *htmp = (TH1F*)hMuonTacDiffCombined[bin]->Clone(Form("%s_clone%d",hMuonTacDiffCombined[bin]->GetName(),j));
	  htmp->Reset();
	  int nxbins = htmp->GetNbinsX();
	  for(int ibin=1; ibin<=nxbins; ibin++)
	    {
	      int jbin = hMuonTacDiffInTrigUnit[0][bin][j]->FindFixBin(htmp->GetBinCenter(ibin) + shift);
	      htmp->SetBinContent(ibin, hMuonTacDiffInTrigUnit[0][bin][j]->GetBinContent(jbin));
	      htmp->SetBinError(ibin, hMuonTacDiffInTrigUnit[0][bin][j]->GetBinError(jbin));
	    }
	  double weight = hMultWeight[1][bin]->GetBinContent(j+2);
	  hMuonTacDiffCombined[bin]->Add(htmp, weight);
	}
    }

  // Compare TacDiff distributions
  TH1F *hMuonTacDiffRebin[nData][nbins-1];
  TH1F *h1tmp = 0x0;
  for(int bin=1; bin<nbins; bin++)
    {
      for(int i=0; i<nData; i++)
	{
	  if(i==0) h1tmp = hMuonTacDiffCombined[bin];
	  if(i==1) h1tmp = (TH1F*)hMuonTacDiffVsTrigUnit[i][bin]->ProjectionY(Form("DataJpsiMuon_MtdVpdTacDiff_bin%d",bin));
	  hMuonTacDiffRebin[i][bin-1] = (TH1F*)h1tmp->Rebin(nBinsTacDiff,Form("%s_MtdVpdTacDiff_bin%d_Rebin",legName[i],bin),xBinsTacDiff);
	  scaleHisto(hMuonTacDiffRebin[i][bin-1], 1, 1, true, false, false);
	  hMuonTacDiffRebin[i][bin-1]->SetXTitle("#DeltaTacSum");
	}
    }
  c = new TCanvas("CompTacDiff_PtBins","CompTacDiff_PtBins",1100,650);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.15,0.7,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  double scale_bin = 799;
  if(data_name=="Run15_pAu200") scale_bin = 910;
  for(int bin=1; bin<nbins; bin++)
    {
      c->cd(bin);
      for(int i=0; i<nData; i++)
	{
	  hMuonTacDiffRebin[i][bin-1]->Scale(1./hMuonTacDiffRebin[i][bin-1]->GetBinContent(hMuonTacDiffRebin[i][bin-1]->FindFixBin(scale_bin)));
	  if(bin<nbins-1) hMuonTacDiffRebin[i][bin-1]->SetMaximum(1.6);
	  else  hMuonTacDiffRebin[i][bin-1]->SetMaximum(5.6);
	  hMuonTacDiffRebin[i][bin-1]->SetMarkerStyle(20+i*4);
	  hMuonTacDiffRebin[i][bin-1]->SetMarkerColor(2-i);
	  hMuonTacDiffRebin[i][bin-1]->SetLineColor(2-i);
	  if(i==0) hMuonTacDiffRebin[i][bin-1]->Draw();
	  else     hMuonTacDiffRebin[i][bin-1]->Draw("samesPE");
	  if(bin==1) 
	    {
	      if(i==0) leg->AddEntry(hMuonTacDiffRebin[i][bin-1],Form("%s (shifted)",legName[i]),"P");
	      else     leg->AddEntry(hMuonTacDiffRebin[i][bin-1],Form("%s",legName[i]),"P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%2.1f < p_{T}^{#mu} < %2.1f GeV/c",lowbins[bin],upbins[bin]),0.055);
      t1->Draw();
    }
  c->cd(2);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sCompTacDiffCombined.pdf",data_name.Data(),config));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_TacDiffDis.pdf"));
    }

  // Fit the TacSum distributions
  TF1 *funcMuonTacDiffRebin[nData][nbins-1];
  TH1F *hMuonTacDiffMean[nData];
  TH1F *hMuonTacDiffSigma[nData];
  for(int i=0; i<nData; i++)
    {
      hMuonTacDiffMean[i]  = new TH1F(Form("%s_hMuonTacDiffMeanVsPt",legName[i]),"Mean of #DeltaTac distribution as a function of p_{T};p_{T} (GeV/c);<#DeltaTac>", nPtBins, xPtBins);
      hMuonTacDiffSigma[i] = new TH1F(Form("%s_hMuonTacDiffSigmaVsPt",legName[i]),"Width of #DeltaTac distribution as a function of p_{T};p_{T} (GeV/c);#sigma_{#DeltaTac}", nPtBins, xPtBins);
      c = new TCanvas(Form("%s_FitTacDiff",legName[i]),Form("%s_FitTacDiff",legName[i]), 1100,650);
      c->Divide(3,2);
      for(int bin=1; bin<nbins; bin++)
	{
	  TH1F *hFit = (TH1F*)hMuonTacDiffRebin[i][bin-1]->Clone(Form("Fit_%s",hMuonTacDiffRebin[i][bin-1]->GetName()));
	  hFit->SetMarkerStyle(24);
	  hFit->SetMarkerSize(1);
	  hFit->SetMarkerColor(1);
	  hFit->SetLineColor(1);
	  hFit->Scale(1./hFit->Integral());
	  hFit->SetMaximum(1.4*hFit->GetMaximum());
	  hFit->SetTitle(";#DeltaTac;");
	  funcMuonTacDiffRebin[i][bin-1] = new TF1(Form("%s_funcMuonTacDiff_PtBin%d",legName[i],bin),"gaus",789,815);
	  if(data_name=="Run14_AuAu200")
	    {
	      funcMuonTacDiffRebin[i][bin-1]->SetRange(789, 815);
	      funcMuonTacDiffRebin[i][bin-1]->SetParameters(1, 800, 10);
	    }
	  if(data_name=="Run15_pAu200")
	    {
	      funcMuonTacDiffRebin[i][bin-1]->SetRange(895, 930);
	      if(i==1 && bin==3) funcMuonTacDiffRebin[i][bin-1]->SetRange(900, 940);
	      funcMuonTacDiffRebin[i][bin-1]->SetParameters(1, 910, 10);
	    }
	  hFit->Fit(funcMuonTacDiffRebin[i][bin-1],"IR0Q");
	  c->cd(bin);
	  hFit->Draw();
	  TPaveText *t1 = GetTitleText(Form("%2.1f < p_{T}^{#mu} < %2.1f GeV/c",lowbins[bin],upbins[bin]),0.055);
	  t1->Draw();
	  funcMuonTacDiffRebin[i][bin-1]->SetLineColor(2);
	  funcMuonTacDiffRebin[i][bin-1]->Draw("sames");
	  hMuonTacDiffMean[i]->SetBinContent(bin, funcMuonTacDiffRebin[i][bin-1]->GetParameter(1));
	  hMuonTacDiffMean[i]->SetBinError(bin, funcMuonTacDiffRebin[i][bin-1]->GetParError(1));
	  hMuonTacDiffSigma[i]->SetBinContent(bin, funcMuonTacDiffRebin[i][bin-1]->GetParameter(2));
	  hMuonTacDiffSigma[i]->SetBinError(bin, funcMuonTacDiffRebin[i][bin-1]->GetParError(2));
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85,0.05,62);
      t1->AddText(legName[i]);
      t1->Draw();
      if(savePlot) 
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sFitTacDiff_%s.pdf",data_name.Data(),config,legName[i]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_FitTacDiffDis_%s.pdf",legName[i]));
	}
    }

  TF1 *funcMuonTacDiffSigma[nData];
  TH1F *hCLfuncMuonTacDiffSigma[nData];
  c = new TCanvas("Comp_TacDiffWidth", "Comp_TacDiffWidth", 1100,450);
  c->Divide(2,1);
  for(int i=0; i<nData; i++)
    {
      funcMuonTacDiffSigma[i] = new TF1(Form("%s_fitMuonTacDiffSigmaVsPt",legName[i]),"[0]+[1]*exp([2]/x)", 1.3, 5);
      funcMuonTacDiffSigma[i]->SetParameters(7, 1e-5, 15);
      TH1F *hFit = (TH1F*)hMuonTacDiffSigma[i]->Clone(Form("Fit_%s",hMuonTacDiffSigma[i]->GetName()));
      hFit->SetMarkerStyle(20+i*4);
      hFit->SetMarkerSize(1.5);
      hFit->SetMarkerColor(2-i);
      hFit->SetLineColor(2-i);
      hFit->GetYaxis()->SetRangeUser(5,12);
      hFit->GetXaxis()->SetRangeUser(1.3,4);
      hFit->Fit(funcMuonTacDiffSigma[i],"IR0Q");
      hFit->SetTitle("");
      ScaleHistoTitle(hFit, 0.05, 0.9, 0.045, 0.055, 0.8, 0.04, 62);
      c->cd(i+1);
      hFit->Draw();
      funcMuonTacDiffSigma[i]->SetLineColor(4);
      funcMuonTacDiffSigma[i]->SetLineStyle(2);
      funcMuonTacDiffSigma[i]->Draw("sames");
      TPaveText *t1 = GetTitleText(legName[i],0.05);
      t1->Draw();
      hCLfuncMuonTacDiffSigma[i] = new TH1F(Form("%s_fitMuonTacDiffSigmaVsPtCL",legName[i]),"",nPtBins,xPtBins);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hCLfuncMuonTacDiffSigma[i], 0.68);
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sCompTacDiffSigma.pdf",data_name.Data(),config));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_CompTacDiffSigma.pdf"));
    }

  c = new TCanvas("Comp_TacDiffMean", "Comp_TacDiffMean", 800, 600);
  TLegend *leg = new TLegend(0.5,0.25,0.7,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(int i=0; i<nData; i++)
    {
      TH1F *hplot = (TH1F*)hMuonTacDiffMean[i]->Clone(Form("plot_%s",hMuonTacDiffMean[i]->GetName()));
      hplot->SetMarkerStyle(20+i*4);
      hplot->SetMarkerSize(1.5);
      hplot->SetMarkerColor(2-i);
      hplot->SetLineColor(2-i);
      hplot->SetTitle("");
      hplot->GetXaxis()->SetRangeUser(1.3,4);
      if(data_name=="Run14_AuAu200") hplot->GetYaxis()->SetRangeUser(785, 805);
      if(data_name=="Run15_pAu200") hplot->GetYaxis()->SetRangeUser(900, 915);
      if(i==0) hplot->Draw();
      else     hplot->Draw("sames");
      leg->AddEntry(hplot, legName[i], "P");
    }
  leg->Draw();
  TPaveText *t1 = GetTitleText("Mean of #DeltaTac distribution as a function of p_{T}",0.04);
  t1->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sCompTacDiffMean.pdf",data_name.Data(),config));
 if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_CompTacDiffMean.pdf"));
    }

  // Fit the shifted distribution from pp collisions
  TF1  *hMuonTacDiffCombinedFit[nbins];
  TFitResultPtr hFitPtr[nbins];
  int iFuncType = 1; // 0 - gaussian; 1 - crystal ball
  if(savePlot || saveHisto) iFuncType = 1;
  TCanvas *c = new TCanvas(Form("%s_TacDiff_Combined",legName[0]),Form("%s_TacDiff_Combined",legName[0]),1100,650);
  c->Divide(4,2);
  for(int bin=0; bin<nbins; bin++)
    {
      TH1F *hfit = (TH1F*)hMuonTacDiffCombined[bin]->Clone(Form("%s_fit",hMuonTacDiffCombined[bin]->GetName()));
      hfit->Rebin(5);
      hfit->SetTitle(";#DeltaTacSum;");
      if(iFuncType==0)
	{
	  hMuonTacDiffCombinedFit[bin] = new TF1(Form("%s_funcMuonTacDiffCombined_PtBin%d",legName[0],bin),"gaus",min_fit,max_fit);
	  hMuonTacDiffCombinedFit[bin]->SetParameters(0.2,min_TacDiffCut[0]+10,5);
	}
      else if(iFuncType==1)
	{
	  hMuonTacDiffCombinedFit[bin] = new TF1(Form("%s_funcMuonTacDiffCombined_PtBin%d",legName[0],bin),CrystalBall,min_fit,max_fit,5);
	  hMuonTacDiffCombinedFit[bin]->SetParameters(1,min_TacDiffCut[0]+10,8,1,0.2);
	}
      hFitPtr[bin] = hfit->Fit(hMuonTacDiffCombinedFit[bin],"I0RS");
      c->cd(bin+1);
      SetPadMargin(gPad,0.13,0.1,0.01,0.1);
      ScaleHistoTitle(hfit, 0.05, 1, 0.04, 0.05, 1, 0.04);
      hfit->SetMaximum(1.5*hfit->GetMaximum());
      hfit->SetMarkerStyle(20);
      if(data_name=="Run14_AuAu200") hfit->GetXaxis()->SetRangeUser(730, 880);
      if(data_name=="Run15_pAu200") hfit->GetXaxis()->SetRangeUser(840, 980);
      hfit->Draw();
      TPaveText *t1 = GetTitleText(Form("%2.1f < p_{T} < %2.1f GeV/c",lowbins[bin],upbins[bin]),0.06);
      t1->Draw();
      hMuonTacDiffCombinedFit[bin]->SetLineColor(2);
      hMuonTacDiffCombinedFit[bin]->Draw("sames");
    }   
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sFitTacDiffComb_%s.pdf",data_name.Data(),config,legName[0]));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_FitTacDiffComb.pdf"));
    }

  TGraphAsymmErrors* gTacDiffEffComb[nLumi][2]; // 0 - bin counting; 1 - fitting;
  double pt_arr[nbins-1], all_arr[nbins-1], all_err_arr[nbins-1], acc_arr[nbins-1], acc_err_arr[nbins-1];
  for(int l=0; l<nLumi; l++)
    {
      // bin counting
      for(int bin=1; bin<nbins; bin++)
	{
	  int xbins = hMuonTacDiffCombined[bin]->GetNbinsX(); 
	  pt_arr[bin-1] = (upbins[bin]+lowbins[bin])/2;
	  all_arr[bin-1] = hMuonTacDiffCombined[bin]->IntegralAndError(1, xbins, all_err_arr[bin-1]);
	  acc_arr[bin-1] = hMuonTacDiffCombined[bin]->IntegralAndError(hMuonTacDiffCombined[bin]->FindFixBin(min_TacDiffCut[l]+0.1),
								       hMuonTacDiffCombined[bin]->FindFixBin(max_TacDiffCut[l]-0.1),
								       acc_err_arr[bin-1]);
	  if(acc_arr[bin-1]>all_arr[bin-1]) acc_arr[bin-1] = all_arr[bin-1];
	}
      gTacDiffEffComb[l][0] = GetEfficiencyCurve(nbins-1, pt_arr, all_arr, all_err_arr, acc_arr, acc_err_arr);
      gTacDiffEffComb[l][0]->SetName(Form("%s_gTacDiffEffComb_%s_%s_BinCount",data_name.Data(),lumi_name[l],legName[0]));

      // fitting
      for(int bin=1; bin<nbins; bin++)
	{
	  pt_arr[bin-1]      = (upbins[bin]+lowbins[bin])/2;
	  all_arr[bin-1]     = hMuonTacDiffCombinedFit[bin]->Integral(700,max_TacDiffCut[l]);
	  all_err_arr[bin-1] = hMuonTacDiffCombinedFit[bin]->IntegralError(700,max_TacDiffCut[l],
									   hMuonTacDiffCombinedFit[bin]->GetParameters(),
									   hFitPtr[bin]->GetCovarianceMatrix().GetMatrixArray());
	  acc_arr[bin-1]     = hMuonTacDiffCombinedFit[bin]->Integral(min_TacDiffCut[l],max_TacDiffCut[l]);
	  acc_err_arr[bin-1] = hMuonTacDiffCombinedFit[bin]->IntegralError(min_TacDiffCut[l],max_TacDiffCut[l],
									   hMuonTacDiffCombinedFit[bin]->GetParameters(),
									   hFitPtr[bin]->GetCovarianceMatrix().GetMatrixArray());
	}
      gTacDiffEffComb[l][1] = GetEfficiencyCurve(nbins-1, pt_arr, all_arr, all_err_arr, acc_arr, acc_err_arr);
      gTacDiffEffComb[l][1]->SetName(Form("%s_gTacDiffEffComb_%s_%s_Fitting",data_name.Data(),lumi_name[l],legName[0]));
	  
      for(int m=0; m<2; m++)
	{
	  int npoints = gTacDiffEffComb[l][m]->GetN();
	  double x, y;
	  for(int ipoint=0; ipoint<npoints; ipoint++)
	    {
	      gTacDiffEffComb[l][m]->GetPoint(ipoint,x,y);
	      gTacDiffEffComb[l][m]->SetPoint(ipoint,hMeanPt[1]->GetBinContent(ipoint+1),y);
	      gTacDiffEffComb[l][m]->SetPointEXhigh(ipoint,hMeanPt[1]->GetBinError(ipoint+1));
	      gTacDiffEffComb[l][m]->SetPointEXlow(ipoint,hMeanPt[1]->GetBinError(ipoint+1));
	      if(data_name=="Run15_pAu200" && m==1) 
		gTacDiffEffComb[l][m]->SetPointEYlow(ipoint,0);
	    }
	}
    }

  // Recalculate the cut values according to the resolution difference
  TGraphAsymmErrors* gTacDiffEffFinal[nLumi][2]; // 0 - bin counting; 1 - fitting;
  double min_TacDiffCut_new[nLumi][nPtBins];
  for(int l=0; l<nLumi; l++)
    {
      // calculate new cut values
      for(int bin=1; bin<nbins; bin++)
	{
	  double pt = hMuonTacDiffMean[0]->GetBinCenter(bin);
	  double scale = funcMuonTacDiffSigma[0]->Eval(pt)/funcMuonTacDiffSigma[1]->Eval(pt);
	  double mean = hMuonTacDiffMean[0]->GetBinContent(bin);
	  min_TacDiffCut_new[l][bin-1] = mean - (mean - min_TacDiffCut[l])*scale;
	  if(bin==nbins-1) min_TacDiffCut_new[l][bin-1] = min_TacDiffCut_new[l][bin-2];
	  cout << pt << "  " << scale << ":  " << min_TacDiffCut[l] << " -> " << min_TacDiffCut_new[l][bin-1] << endl;
	}

      // bin counting
      for(int bin=1; bin<nbins; bin++)
	{
	  TH1F *hppTac = hMuonTacDiffCombined[bin];
	  int xbins = hppTac->GetNbinsX(); 
	  pt_arr[bin-1] = (upbins[bin]+lowbins[bin])/2;
	  all_arr[bin-1] = hppTac->IntegralAndError(1, xbins, all_err_arr[bin-1]);
	  int low_bin =    hppTac->FindFixBin(min_TacDiffCut_new[l][bin-1]);
	  acc_arr[bin-1] = hppTac->IntegralAndError(low_bin+1,
						    hppTac->FindFixBin(max_TacDiffCut[l]-0.1),
						    acc_err_arr[bin-1]);
	  double fraction = (hppTac->GetXaxis()->GetBinUpEdge(low_bin)-min_TacDiffCut_new[l][bin-1])/hppTac->GetBinWidth(low_bin);
	  acc_arr[bin-1] += hppTac->GetBinContent(low_bin)*fraction;
	  acc_err_arr[bin-1] = sqrt(pow(acc_err_arr[bin-1],2)+pow(hppTac->GetBinError(low_bin)*fraction,2));
	  if(acc_arr[bin-1]>all_arr[bin-1]) acc_arr[bin-1] = all_arr[bin-1];
	}
      gTacDiffEffFinal[l][0] = GetEfficiencyCurve(nbins-1, pt_arr, all_arr, all_err_arr, acc_arr, acc_err_arr);
      gTacDiffEffFinal[l][0]->SetName(Form("%s_gTacDiffEffFinal_BinCount_%s_%s",data_name.Data(),lumi_name[l],legName[0]));

      // fitting
      for(int bin=1; bin<nbins; bin++)
	{
	  TF1 *funcppTac = hMuonTacDiffCombinedFit[bin];
	  pt_arr[bin-1]      = (upbins[bin]+lowbins[bin])/2;
	  all_arr[bin-1]     = funcppTac->Integral(700,max_TacDiffCut[l]);
	  all_err_arr[bin-1] = funcppTac->IntegralError(700,max_TacDiffCut[l],
							funcppTac->GetParameters(),
							hFitPtr[bin]->GetCovarianceMatrix().GetMatrixArray());
	  acc_arr[bin-1]     = funcppTac->Integral(min_TacDiffCut_new[l][bin-1],max_TacDiffCut[l]);
	  acc_err_arr[bin-1] = funcppTac->IntegralError(min_TacDiffCut_new[l][bin-1],max_TacDiffCut[l],
							funcppTac->GetParameters(),
							hFitPtr[bin]->GetCovarianceMatrix().GetMatrixArray());
	}
      gTacDiffEffFinal[l][1] = GetEfficiencyCurve(nbins-1, pt_arr, all_arr, all_err_arr, acc_arr, acc_err_arr);
      gTacDiffEffFinal[l][1]->SetName(Form("%s_gTacDiffEffFinal_Fitting_%s_%s",data_name.Data(),lumi_name[l],legName[0]));

      for(int m=0; m<2; m++)
	{
	  int npoints = gTacDiffEffFinal[l][m]->GetN();
	  double x, y;
	  for(int ipoint=0; ipoint<npoints; ipoint++)
	    {
	      gTacDiffEffFinal[l][m]->GetPoint(ipoint,x,y);
	      gTacDiffEffFinal[l][m]->SetPoint(ipoint,hMeanPt[1]->GetBinContent(ipoint+1),y);
	      gTacDiffEffFinal[l][m]->SetPointEXhigh(ipoint,hMeanPt[1]->GetBinError(ipoint+1));
	      gTacDiffEffFinal[l][m]->SetPointEXlow(ipoint,hMeanPt[1]->GetBinError(ipoint+1));
	      if(data_name=="Run15_pAu200" && m==1) 
		gTacDiffEffFinal[l][m]->SetPointEYlow(ipoint,0);
	    }
	}
    }

  // compare the effects of new cuts
  const char *method_name[2] = {"Bin counting","Fitting"};
  TH1F *hplot = new TH1F("hplot",";p_{T} (GeV/c);eff",100,1,7);
  if(data_name=="Run14_AuAu200") hplot->GetYaxis()->SetRangeUser(0.5,1);
  if(data_name=="Run15_pAu200") hplot->GetYaxis()->SetRangeUser(0.95,1.02);
  for(int l=0; l<nLumi; l++)
    {
      TCanvas *c = new TCanvas(Form("%s_TacDiffEff_AuAuReso_%s",legName[0],lumi_name[l]),Form("%s_TacDiffEff_AuAuReso_%s",legName[0],lumi_name[l]),800,600);
      hplot->DrawCopy();
      TPaveText *t1 = GetTitleText(Form("%s: MTD online trigger efficiency",data_name.Data()));
      t1->Draw();
      TLegend *leg = new TLegend(0.3,0.15,0.5,0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.03);
      leg->SetHeader(Form("%s",lumi_name[l]));
      for(int m=0; m<2; m++)
	{
	  TGraphAsymmErrors *gtmp = (TGraphAsymmErrors*)gTacDiffEffComb[l][m]->Clone(Form("%s_clone",gTacDiffEffComb[l][m]->GetName()));
	  gtmp->SetMarkerStyle(20+m*4);
	  gtmp->SetMarkerSize(1.2);
	  gtmp->SetMarkerColor(1);
	  gtmp->SetLineColor(1);
	  gtmp->Draw("samesPEZ");
	  leg->AddEntry(gtmp, Form("pp resolution: %s",method_name[m]), "PL");
	}

      for(int m=0; m<2; m++)
	{
	  TGraphAsymmErrors *gtmp = (TGraphAsymmErrors*)gTacDiffEffFinal[l][m]->Clone(Form("%s_clone",gTacDiffEffFinal[l][m]->GetName()));
	  gtmp->SetMarkerStyle(21+m*4);
	  gtmp->SetMarkerSize(1.2);
	  gtmp->SetMarkerColor(2);
	  gtmp->SetLineColor(2);
	  gtmp->Draw("samesPEZ");
	  if(data_name=="Run14_AuAu200") leg->AddEntry(gtmp, Form("AuAu resolution: %s",method_name[m]), "PL");
	  if(data_name=="Run15_pAu200") leg->AddEntry(gtmp, Form("pAu resolution: %s",method_name[m]), "PL");
	}
      leg->Draw();
      if(savePlot) 
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sTacDiffEffAuAuReso_%s.pdf",data_name.Data(),config,lumi_name[l]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_AuAuReso_%s.pdf",lumi_name[l]));
	}
    }
  
  // final efficiency
  for(int l=0; l<nLumi; l++)
    {
      TCanvas *c = new TCanvas(Form("%s_TacDiffEff_%s",data_name.Data(),lumi_name[l]),Form("%s_TacDiffEff_%s",data_name.Data(),lumi_name[l]),800,600);
      TH1F *htmp = new TH1F(Form("htmp_%s",lumi_name[l]),";p_{T} (GeV/c);eff",100,1,7);
      if(data_name=="Run14_AuAu200") htmp->GetYaxis()->SetRangeUser(0.5,1);
      if(data_name=="Run15_pAu200") htmp->GetYaxis()->SetRangeUser(0.95,1.02);
      htmp->Draw();
      TPaveText *t1 = GetTitleText(Form("%s: MTD online trigger efficiency",data_name.Data()));
      t1->Draw();
      TLegend *leg = new TLegend(0.3,0.15,0.5,0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.03);
      leg->SetHeader(Form("%s",lumi_name[l]));
      for(int m=0; m<2; m++)
	{
	  gTacDiffEffFinal[l][m]->SetMarkerStyle(20+4*m);
	  gTacDiffEffFinal[l][m]->SetMarkerSize(1.2);
	  gTacDiffEffFinal[l][m]->SetMarkerColor(m+1);
	  gTacDiffEffFinal[l][m]->SetLineColor(m+1);
	  gTacDiffEffFinal[l][m]->Draw("samesPEZ");
	  leg->AddEntry(gTacDiffEffFinal[l][m], Form("%s: %s",legName[0],method_name[m]), "PL");
	}
      leg->Draw();
      if(savePlot) 
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/%sTacDiffEffFinal_BinCountVsFit_%s.pdf",data_name.Data(),config,lumi_name[l]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTirg_FinalEff_%s.pdf",lumi_name[l]));
	}
    }

  // save histograms
  if(saveHisto) 
    {
      fdata[1]->cd();
      for(int bin=0; bin<nbins; bin++)
	{
	  hMuonTacDiffCombined[bin]->Write("",TObject::kOverwrite);
	  hMuonTacDiffCombinedFit[bin]->Write("",TObject::kOverwrite);
	}
    
      for(int i=0; i<nData; i++)
	{
	  for(int bin=1; bin<nbins; bin++)
	    {
	      hMuonTacDiffRebin[i][bin-1]->Write("",TObject::kOverwrite);
	    }
	  hMuonTacDiffMean[i]->Write("",TObject::kOverwrite);
	  hMuonTacDiffSigma[i]->Write("",TObject::kOverwrite);
	  funcMuonTacDiffSigma[i]->Write("",TObject::kOverwrite);
	  hCLfuncMuonTacDiffSigma[i]->Write("",TObject::kOverwrite);
	}
      for(int l=0; l<nLumi; l++)
	{
	  for(int j=0; j<2; j++)
	    {
	      gTacDiffEffFinal[l][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void compareTacDiff(const int savePlot = 0)
{
  // compare TacDiff distribution in different senarios
  const char* type_name[3] = {"UL","LS","Muon"};
  const int nPt = 4;
  const double lowBinsPt[nPt] = {1.3, 1.5, 2.0, 3.0};
  const double upBinsPt[nPt]  = {1.5, 2.0, 3.0, 5.0};

  /*
  const int nFile = 2;
  const char* fileName[nFile] = {"Run15_pp200.JpsiMuon.root", "Run15_pp200.JpsiMuon.Ranking.root"};
  const char* dataName[nFile] = {"NoRank","Rank"};
  const char* legName[nFile]  = {"All","Ranking>0"};
  const char* saveName = "Rank";
  */
  const int nFile = 2;
  const char* fileName[nFile] = {"Run15_pp200.JpsiMuon.root", "Run15_pp200.JpsiMuon.AuAuPidCut.root"};
  const char* dataName[nFile] = {"ppPID","AuAuPID"};
  const char* legName[nFile]  = {"pp PID cuts","AuAu PID cuts"};
  const char* saveName = "PID";

  TH1F *hTacDiff[nFile][nPt][3]; // 0 - UL; 1 - LS; 2 - muon
  TFile *fdata[nFile];
  THnSparseF *hn[nFile];
  for(int i=0; i<nFile; i++)
    {
      fdata[i] = TFile::Open(Form("output/%s",fileName[i]),"read");
      hn[i] = (THnSparseF*)fdata[i]->Get("mhJpsiMuonTrigEff_di_mu");
      hn[i]->SetName(Form("%s_%s",hn[i]->GetName(),dataName[i]));
      hn[i]->GetAxis(0)->SetRangeUser(3.0+1e-4, 3.2-1e-4); // mass cut
      for(int p=0; p<nPt; p++)
	{
	  hn[i]->GetAxis(2)->SetRangeUser(lowBinsPt[p]+1e-4, upBinsPt[p]-1e-4);
	  for(int j=0; j<2; j++)
	    {
	      hn[i]->GetAxis(5)->SetRange(j+1,j+1);
	      hTacDiff[i][p][j] = (TH1F*)hn[i]->Projection(1);
	      hTacDiff[i][p][j]->SetName(Form("hTacDiff_%s_Pt%d_%s",dataName[i],p,type_name[j]));
	      hTacDiff[i][p][j]->Sumw2();
	    }
	  hn[i]->GetAxis(5)->SetRange(0,-1);
	  hTacDiff[i][p][2] = (TH1F*)hTacDiff[i][p][0]->Clone(Form("hTacDiff_%s_Pt%d_%s",dataName[i],p,type_name[2]));
	  hTacDiff[i][p][2]->Add(hTacDiff[i][p][1],-1);
	}
      hn[i]->GetAxis(2)->SetRange(0,-1);
      hn[i]->GetAxis(0)->SetRange(0,-1);
    }

  //==============================================
  // Compare UL vs. LS
  for(int i=0; i<nFile; i++)
    {
      TCanvas *c = new TCanvas(Form("TacDiff_ULvsLS_%s",dataName[i]),Form("TacDiff_ULvsLS_%s",dataName[i]),1100,700);
      c->Divide(2,2);
      TLegend *leg = new TLegend(0.15,0.65,0.3,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader(legName[i]);
      for(int p=0; p<nPt; p++)
	{
	  c->cd(p+1);
	  gPad->SetLogy();
	  for(int j=0; j<3; j++)
	    {
	      TH1F *htmp = (TH1F*)hTacDiff[i][p][j]->Clone(Form("%s_clone",hTacDiff[i][p][j]->GetName()));
	      if(j<2) htmp->SetMarkerStyle(20+j*4);
	      else    htmp->SetMarkerStyle(21);
	      htmp->SetMarkerColor(TMath::Power(2,j));
	      htmp->SetLineColor(TMath::Power(2,j));
	      htmp->Rebin(4);
	      htmp->GetXaxis()->SetRangeUser(860,980);
	      htmp->SetMaximum(10*htmp->GetMaximum());
	      htmp->SetTitle(";#DeltaTacSum;Counts");
	      ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 0.9, 0.045,62);
	      if(j==0) htmp->Draw();
	      else     htmp->Draw("sames");
	      if(p==0)
		leg->AddEntry(htmp,type_name[j],"P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowBinsPt[p], upBinsPt[p]),0.06);
	  t1->Draw();
	  c->cd(2);
	  leg->Draw();
	}
      c->cd();
      TPaveText *t1 = GetPaveText(0.15,0.25,0.75,0.8,0.045);
      t1->AddText("Run15_pp200");
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15_pp200_TacDiff_ULvsLS_%s.pdf",dataName[i]));
    }

  //==============================================
  // Compare TacDiff distributions
  TCanvas *c = new TCanvas(Form("CompTacDiff"),Form("CompTacDiff"),1100,700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.65,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Run15_pp200");
  for(int p=0; p<nPt; p++)
    {
      for(int i=0; i<nFile; i++)
	{
	  TH1F *htmp = (TH1F*)hTacDiff[i][p][2]->Clone(Form("%s_clone2",hTacDiff[i][p][2]->GetName()));
	  htmp->SetMarkerStyle(21+4*i);
	  htmp->SetMarkerColor(TMath::Power(2,i));
	  htmp->SetLineColor(TMath::Power(2,i));
	  htmp->Rebin(4);
	  htmp->GetXaxis()->SetRangeUser(860,980);
	  htmp->Scale(1./htmp->Integral());
	  htmp->SetTitle(";#DeltaTacSum;");
	  ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 0.9, 0.045,62);
	  c->cd(p+1);
	  gPad->SetLogy();
	  if(p==0) htmp->SetMaximum(20*htmp->GetMaximum());
	  else     htmp->SetMaximum(10*htmp->GetMaximum());
	  if(i==0) htmp->DrawCopy();
	  else     htmp->DrawCopy("sames"); 

	  if(p==0)
	    leg->AddEntry(htmp,legName[i],"P");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowBinsPt[p], upBinsPt[p]),0.06);
      c->cd(p+1);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15_pp200_TacDiff_vs_%s.pdf",saveName));

  //==============================================
  // make ratio
  TCanvas *c = new TCanvas("TacDiffRatio","TacDiffRatio",1100,700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.72,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Run15_pp200");
  for(int p=0; p<nPt; p++)
    {
      TH1F *htmp = (TH1F*)hTacDiff[1][p][2]->Clone(Form("%s_ratio",hTacDiff[1][p][2]->GetName()));
      htmp->Scale(1./htmp->Integral());
      htmp->Rebin(4);
      TH1F *htmp2 = (TH1F*)hTacDiff[0][p][2]->Clone(Form("%s_ratio",hTacDiff[0][p][2]->GetName()));
      htmp2->Rebin(4);
      htmp2->Scale(1./htmp2->Integral());
      htmp->Divide(htmp2);

      htmp->SetMarkerStyle(20);
      htmp->SetMarkerColor(1);
      htmp->SetLineColor(1);
      htmp->GetXaxis()->SetRangeUser(860,980);
      htmp->GetYaxis()->SetRangeUser(0,2);
      htmp->SetTitle(";#DeltaTacSum;Ratio");
      ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 0.9, 0.045,62);
      c->cd(p+1);
      htmp->DrawCopy();
      if(p==0) leg->AddEntry(htmp,Form("%s/%s",legName[1],legName[0]),"P");
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowBinsPt[p], upBinsPt[p]),0.06);
      c->cd(p+1);
      t1->Draw();
      TLine *line = GetLine(860,1,980,1,1);
      line->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15_pp200_TacDiffRatio_%s.pdf",saveName));
  
}

//================================================
void lumiDepend(const int savePlot = 0)
{
  const int nPt = 4;
  const double lowBinsPt[nPt] = {1.3, 1.5, 2.0, 3.0};
  const double upBinsPt[nPt]  = {1.5, 2.0, 3.0, 5.0};
  const int nLumi = 3;
  const double lowBinsLumi[nLumi] = {500, 1100, 1500};
  const double upBinsLumi[nLumi]  = {1000, 1200, 3000};
  const char* type_name[3] = {"UL","LS","Muon"};

  TH1F *hTacDiff[nLumi][nPt][3]; // 0 - UL; 1 - LS; 2 - muon
  TFile *fdata = TFile::Open("output/Run15_pp200.JpsiMuon.root","read");
  THnSparseF *hn = (THnSparseF*)fdata->Get("mhJpsiMuonTrigEff_di_mu");
  hn->GetAxis(0)->SetRangeUser(3.0+1e-4, 3.2-1e-4); // mass cut
  for(int l=0; l<nLumi; l++)
    {
      hn->GetAxis(6)->SetRangeUser(lowBinsLumi[l]+1, upBinsLumi[l]-1);
      for(int p=0; p<nPt; p++)
	{
	  hn->GetAxis(2)->SetRangeUser(lowBinsPt[p]+1e-4, upBinsPt[p]-1e-4);
	  for(int i=0; i<2; i++)
	    {
	      hn->GetAxis(5)->SetRange(i+1,i+1);
	      hTacDiff[l][p][i] = (TH1F*)hn->Projection(1);
	      hTacDiff[l][p][i]->SetName(Form("hTacDiff_Lumi%d_Pt%d_%s",l,p,type_name[i]));
	      hTacDiff[l][p][i]->Sumw2();
	    }
	  hn->GetAxis(5)->SetRange(0,-1);
	  hTacDiff[l][p][2] = (TH1F*)hTacDiff[l][p][0]->Clone(Form("hTacDiff_Lumi%d_Pt%d_%s",l,p,type_name[2]));
	  hTacDiff[l][p][2]->Add(hTacDiff[l][p][1],-1);
	}
      hn->GetAxis(2)->SetRange(0,-1);
    }
  hn->GetAxis(6)->SetRange(0,-1);
  hn->GetAxis(0)->SetRange(0,-1);

  //==============================================
  // Compare UL vs. LS
  for(int l=0; l<nLumi; l++)
    {
      TCanvas *c = new TCanvas(Form("TacDiff_ULvsLS_Lumi%d",l),Form("TacDiff_ULvsLS_Lumi%d",l),1100,700);
      c->Divide(2,2);
      TLegend *leg = new TLegend(0.15,0.65,0.3,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader(Form("%1.1f < BBCx < %1.1f MHz",lowBinsLumi[l]/1e3, upBinsLumi[l]/1e3));
      for(int p=0; p<nPt; p++)
	{
	  c->cd(p+1);
	  gPad->SetLogy();
	  for(int i=0; i<3; i++)
	    {
	      TH1F *htmp = (TH1F*)hTacDiff[l][p][i]->Clone(Form("%s_clone",hTacDiff[l][p][i]->GetName()));
	      if(i<2) htmp->SetMarkerStyle(20+i*4);
	      else    htmp->SetMarkerStyle(21);
	      htmp->SetMarkerColor(TMath::Power(2,i));
	      htmp->SetLineColor(TMath::Power(2,i));
	      htmp->Rebin(4);
	      htmp->GetXaxis()->SetRangeUser(860,980);
	      htmp->SetMaximum(10*htmp->GetMaximum());
	      htmp->SetTitle(";#DeltaTacSum;Counts");
	      ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 0.9, 0.045,62);
	      if(i==0) htmp->Draw();
	      else     htmp->Draw("sames");
	      if(p==0)
		leg->AddEntry(htmp,type_name[i],"P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowBinsPt[p], upBinsPt[p]),0.06);
	  t1->Draw();
	  c->cd(2);
	  leg->Draw();
	}
      c->cd();
      TPaveText *t1 = GetPaveText(0.15,0.25,0.75,0.8,0.045);
      t1->AddText("Run15_pp200");
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15_pp200_TacDiff_ULvsLS_Lumi%d.pdf",l));
    }

  //==============================================
  // Luminosity dependence
  TCanvas *c = new TCanvas(Form("TacDiff_vs_Lumi"),Form("TacDiff_vs_Lumi"),1100,700);
  c->Divide(2,2);
  TCanvas *c1 = new TCanvas(Form("TacDiff_vs_Lumi_2"),Form("TacDiff_vs_Lumi_2"),1100,700);
  c1->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.65,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Run15_pp200");
  for(int p=0; p<nPt; p++)
    {
      for(int l=0; l<nLumi; l++)
	{
	  TH1F *htmp = (TH1F*)hTacDiff[l][p][2]->Clone(Form("%s_clone2",hTacDiff[l][p][2]->GetName()));
	  htmp->SetMarkerStyle(24+l);
	  htmp->SetMarkerColor(TMath::Power(2,l));
	  htmp->SetLineColor(TMath::Power(2,l));
	  htmp->Rebin(4);
	  htmp->GetXaxis()->SetRangeUser(860,980);
	  htmp->Scale(1./htmp->Integral());
	  htmp->SetTitle(";#DeltaTacSum;");
	  ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 0.9, 0.045,62);

	  c1->cd(p+1);
	  if(p==0) htmp->SetMaximum(1.5*htmp->GetMaximum());
	  else     htmp->SetMaximum(1.2*htmp->GetMaximum());
	  if(l==0) htmp->DrawCopy();
	  else     htmp->DrawCopy("sames");

	  c->cd(p+1);
	  gPad->SetLogy();
	  if(p==0) htmp->SetMaximum(20*htmp->GetMaximum());
	  else     htmp->SetMaximum(10*htmp->GetMaximum());
	  if(l==0) htmp->DrawCopy();
	  else     htmp->DrawCopy("sames"); 

	  if(p==0)
	    leg->AddEntry(htmp,Form("%1.1f < BBCx < %1.1f MHz",lowBinsLumi[l]/1e3, upBinsLumi[l]/1e3),"P");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowBinsPt[p], upBinsPt[p]),0.06);
      c->cd(p+1);
      t1->Draw();
      c1->cd(p+1);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  c1->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15_pp200_TacDiff_vs_Lumi.pdf"));
  if(savePlot) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15_pp200_TacDiff_vs_Lumi_Linear.pdf"));
}

//================================================
void trigElecEff(const int savePlot = 0, const int saveHisto = 0)
{
  if(year==2014)
    {
      f = TFile::Open("output/Run14_AuAu200.MB.TrigElecEff.root","read");
    }
  else
    {
      printf("[e] No available input file!\n");
      return;
    }

  THnSparseF *hnTrigEff = (THnSparseF*)f->Get("mhMtdTrigElecEff_mb");
  const int nbins = 11;
  const double xbins[nbins+1] = {0,1,1.5,2,2.5,3,3.5,4,5,6,8,10};
  
  TList *list = new TList;
  // Efficiency vs. dTof
  const int nDtof = 4;
  const double dtof_cut[nDtof] = {1, 0.75, 0.5, 0.25};
  TH1F *hMuonPtDtof[nDtof][3];
  TH1F *hMuonEffDtof[nDtof][3];
  for(int i=0; i<nDtof; i++)
    {
      hnTrigEff->GetAxis(3)->SetRangeUser(-3,dtof_cut[i]);
      hMuonPtDtof[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][0]->SetName(Form("hMuonPtMatch_DtofCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtDtof[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][1]->SetName(Form("hMuonPtQT_DtofCut%d",i));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtDtof[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][2]->SetName(Form("hMuonPtHigh2_DtofCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(3)->SetRange(0,-1);
    }

  TString legName1[nDtof];
  for(int i=0; i<nDtof; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hMuonPtDtof[i][j] = (TH1F*)hMuonPtDtof[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtDtof[i][j]->GetName()),xbins);
	  hMuonEffDtof[i][j] = (TH1F*)hMuonPtDtof[i][j]->Clone(Form("hMuonEffDtof_%d_%d",i,j));
	  hMuonEffDtof[i][j]->Sumw2();
	  hMuonEffDtof[i][j]->Divide(hMuonPtDtof[i][0]);
	  hMuonEffDtof[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffDtof[i][j]); 
	}
      legName1[i] = Form("#Deltatof < %2.2f ns",dtof_cut[i]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_Dtof",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName1,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareDtof.pdf",run_type));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTrigElec_vs_PID.pdf"));
 
  // Efficiency vs. centrality
  // use dtof < 1 ns cut
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;
  TH1F *hMuonPtCent[nCentBins][3];
  TH1F *hMuonEffCent[nCentBins][3];
  for(int i=0; i<nCentBins; i++)
    {
      hnTrigEff->GetAxis(5)->SetRange(centBins_low[i], centBins_high[i]);
      hMuonPtCent[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][0]->SetName(Form("hMuonPtMatch_%s",cent_Title[i]));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtCent[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][1]->SetName(Form("hMuonPtQT_%s",cent_Title[i]));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtCent[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][2]->SetName(Form("hMuonPtHigh2_%s",cent_Title[i]));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(5)->SetRange(0,-1);
    }

  TString legName2[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hMuonPtCent[i][j] = (TH1F*)hMuonPtCent[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtCent[i][j]->GetName()),xbins);
	  hMuonEffCent[i][j] = (TH1F*)hMuonPtCent[i][j]->Clone(Form("hMuonEff_%s_%d",cent_Title[i],j));
	  hMuonEffCent[i][j]->Sumw2();
	  hMuonEffCent[i][j]->Divide(hMuonPtCent[i][0]);
	  hMuonEffCent[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffCent[i][j]); 
	}
      legName2[i] = Form("%s%%",cent_Name[i]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_Cent",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName2,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareCent.pdf",run_type));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTrigElec_vs_cent.pdf"));

  // Efficiency vs. TPC vz
  // use dtof < 1 ns cut and 0-60%
  hnTrigEff->GetAxis(5)->SetRange(5,16);
  const int nTpcVz = 6;
  const double tpcvz_cut[7] = {-100,-30,-5,0,5,30,100};
  TH1F *hMuonPtTpcVz[nTpcVz][3];
  TH1F *hMuonEffTpcVz[nTpcVz][3];
  for(int i=0; i<nTpcVz; i++)
    {
      hnTrigEff->GetAxis(4)->SetRangeUser(tpcvz_cut[i]+1e-4,tpcvz_cut[i+1]-1e-4);
      hMuonPtTpcVz[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][0]->SetName(Form("hMuonPtMatch_TpcVzCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtTpcVz[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][1]->SetName(Form("hMuonPtQT_TpcVzCut%d",i));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtTpcVz[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][2]->SetName(Form("hMuonPtHigh2_TpcVzCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(4)->SetRange(0,-1);

      for(int j=0; j<3; j++)
	{
	  hMuonPtTpcVz[i][j] = (TH1F*)hMuonPtTpcVz[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtTpcVz[i][j]->GetName()),xbins);
	}
    }

  TString legName3[nTpcVz];
  for(int i=3; i<6; i++)
    {
      TH1F *hRef = (TH1F*)hMuonPtTpcVz[i][0]->Clone(Form("%s_clone",hMuonPtTpcVz[i][0]->GetName()));
      hRef->Add(hMuonPtTpcVz[5-i][0]);
      for(int j=0; j<3; j++)
	{
	  hMuonEffTpcVz[i][j] = (TH1F*)hMuonPtTpcVz[i][j]->Clone(Form("hMuonEff_%d_%d",i,j));
	  hMuonEffTpcVz[i][j]->Sumw2();
	  hMuonEffTpcVz[i][j]->Add(hMuonPtTpcVz[5-i][j]);
	  hMuonEffTpcVz[i][j]->Divide(hRef);
	  hMuonEffTpcVz[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffTpcVz[i][j]); 
	}
      legName3[i-3] = Form("%1.0f < |vz| < %1.0f cm",tpcvz_cut[i],tpcvz_cut[i+1]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_TpcVz",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName3,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareTpcVz.pdf",run_type));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTrigElec_vs_vz.pdf"));
  hnTrigEff->GetAxis(5)->SetRange(0,-1);

  // Final efficiency: dtof < 1 ns, 0-60%, |tpcVz| < 100 cm
  TH1F *hTrigElecEff = (TH1F*)hMuonEffCent[0][2]->Clone(Form("%s_TrigElecEff",run_type));
  TF1 *func = new TF1(Form("%s_FitFunc",hTrigElecEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,10);
  func->SetParLimits(0,0,1);
  hTrigElecEff->Fit(func,"R0Q");
  hTrigElecEff->GetYaxis()->SetRangeUser(0.85,1.1);
  c = draw1D(hTrigElecEff,Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type));
  func->SetLineColor(2);
  func->SetLineStyle(2);
  func->Draw("sames");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_Fit.pdf",run_type));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffTrigElecFit.pdf"));
  
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type),"update");
      hTrigElecEff->Write("",TObject::kOverwrite);
      func->Write("",TObject::kOverwrite);
      fout->Close();
    }
}


//================================================
void extrapolation(const int savePlot = 0, const int saveHisto = 0)
{
  const char *trgSetupName[4] = {"AuAu_200_production_2014","AuAu_200_production_low_2014","AuAu_200_production_mid_2014","AuAu_200_production_high_2014"};
  TF1 *funcTrigEff[4][nCentBins];
  TF1 *funcTrigEffCorr[4][nCentBins];
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  funcTrigEff[i][k] = new TF1(Form("MuonTrigEff_cent%s%s",cent_Title[k],gTrgSetupName[i+1]),"[0]-exp(-1*[1]*(x-[2]))",1.2,12);
	  funcTrigEffCorr[i][k] = new TF1(Form("MuonTrigEffCorr_cent%s%s",cent_Title[k],gTrgSetupName[i+1]),"[0]-exp(-1*[1]*(x-[2]))",1.2,12);
	  if(i==0 || i==1)
	    {
	      if(k==0) funcTrigEff[i][k]->SetParameters(0.951456, 2.32223, 6.61693e-14);      
	      if(k==1) funcTrigEff[i][k]->SetParameters(0.948947, 2.55808, 0.0659758);
	      if(k==2) funcTrigEff[i][k]->SetParameters(0.950789, 2.51311, 7.49401e-14);
	      if(k==3) funcTrigEff[i][k]->SetParameters(0.968246, 2.45772, 8.27116e-15);

	      if(k==0) funcTrigEffCorr[i][k]->SetParameters(0.96023,  10.5511, 0.802917);      
	      if(k==1) funcTrigEffCorr[i][k]->SetParameters(0.96918,  9.40215, 0.778563);
	      if(k==2) funcTrigEffCorr[i][k]->SetParameters(0.964368, 3.65518, 0.298387);
	      if(k==3) funcTrigEffCorr[i][k]->SetParameters(0.966162, 2.77122, 3.6e-15);
	    }
	  else
	    {
	      if(k==0) funcTrigEff[i][k]->SetParameters(0.908, 2.09256, 1.03406e-14);      
	      if(k==1) funcTrigEff[i][k]->SetParameters(0.905465, 2.05892, 9.27036e-15);
	      if(k==2) funcTrigEff[i][k]->SetParameters(0.920748, 1.89851, 1.57097e-14);
	      if(k==3) funcTrigEff[i][k]->SetParameters(0.883885, 2.1119, 3.33068e-15);

	      if(k==0) funcTrigEffCorr[i][k]->SetParameters(0.936584, 17.1834, 0.893653);      
	      if(k==1) funcTrigEffCorr[i][k]->SetParameters(0.945798, 11.6541, 0.838267);
	      if(k==2) funcTrigEffCorr[i][k]->SetParameters(0.935776, 5.75038, 0.499548);
	      if(k==3) funcTrigEffCorr[i][k]->SetParameters(0.938118, 13.9462, 0.875862);
	    }
	}
    }

  for(int i=0; i<4; i++)
    {
      TCanvas *c = new TCanvas(Form("MuonTrigEff%s",gTrgSetupName[i+1]),Form("MuonTrigEff%s",gTrgSetupName[i+1]),1100,700);
      c->Divide(2,2);
      for(int k=0; k<nCentBins; k++)
	{
	  c->cd(k+1);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	  gPad->SetGrid(1,1);
	  ScaleHistoTitle(funcTrigEff[i][k]->GetHistogram(),0.06,1,0.05,0.06,0.9,0.05,62);
	  funcTrigEff[i][k]->SetTitle(Form("%s: %s%%;p_{T} (GeV/c);Trigger efficiency",trgSetupName[i],cent_Name[k]));
	  funcTrigEff[i][k]->SetMaximum(1.1);
	  funcTrigEff[i][k]->SetMinimum(0.5);
	  funcTrigEff[i][k]->SetLineColor(1);
	  funcTrigEff[i][k]->Draw();
	  funcTrigEffCorr[i][k]->SetLineColor(4);
	  funcTrigEffCorr[i][k]->Draw("sames");
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.45,0.25,0.65,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(funcTrigEff[i][0],"VPDMB5","L");
      leg->AddEntry(funcTrigEffCorr[i][0],"NoVtx/VPDMB5","L");
      leg->Draw();
  
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEff_CentBins%s.pdf",run_type,gTrgSetupName[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEff_CentBins%s.png",run_type,gTrgSetupName[i+1]));
	}
    }

  TH1F *hMuonTrigEff[4][nCentBins];
  TString legName[nCentBins];
  TList *list = new TList;
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hMuonTrigEff[i][k] = new TH1F(Form("CombinedMuonTrigEff_cent%s%s",cent_Title[k],gTrgSetupName[i+1]),"Muon efficiency",1000,0,10);
	  for(int bin=1; bin<=hMuonTrigEff[i][k]->GetNbinsX(); bin++)
	    {
	      double x = hMuonTrigEff[i][k]->GetXaxis()->GetBinCenter(bin);
	      hMuonTrigEff[i][k]->SetBinContent(bin,funcTrigEff[i][k]->Eval(x)*funcTrigEffCorr[i][k]->Eval(x));
	      hMuonTrigEff[i][k]->SetBinError(bin,0);
	    }
	  list->Add(hMuonTrigEff[i][k]);
	  legName[k] = Form("%s%%",cent_Name[k]);
	}
      c = drawHistos(list,Form("MuonTrigEff_CentBins%s",gTrgSetupName[i+1]),Form("Single muon trigger efficiency (%s);p_{T} (GeV/c);Efficiency",trgSetupName[i]),kTRUE,1.2,10,kTRUE,0.5,1.0,kFALSE,kTRUE,legName,kTRUE,"",0.4,0.6,0.3,0.55,kTRUE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEffCombined_CentBins%s.pdf",run_type,gTrgSetupName[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEffCombined_CentBins%s.png",run_type,gTrgSetupName[i+1]));
	}
      list->Clear();
    }
  
  if(saveHisto)
    {
      char *outName = "Run14.AuAu200.MuonTrigEff.root";
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outName),"recreate");
      for(int i=0; i<4; i++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      funcTrigEff[i][k]->Write("",TObject::kOverwrite);
	      funcTrigEffCorr[i][k]->Write("",TObject::kOverwrite);
	      hMuonTrigEff[i][k]->Write("",TObject::kOverwrite);
	    }
	}
      fout->Close();
    }
}

//================================================
void trigUnitMap(const int savePlot = 1)
{
  const Int_t kTrgID[30][5] = {
    { 7, 14, 6,  15,  8}, //BL1 
    { 7, 14, 6,  15,  8}, //BL2
    { 9, 11, 13, 12, 10}, //BL3
    { 9, 11, 13, 12, 10}, //BL4
    { 9, 11, 13, 12, 10}, //BL5
    { 9, 11, 13, 12, 10}, //BL6
    { 9, 11, 13, 12, 10}, //BL7
    {16, 18, 20, 19, 17}, //BL8
    { 0,  0,  0,  0,  0}, //BL9
    {16, 18, 20, 19, 17}, //BL10
    {16, 18, 20, 19, 17}, //BL11
    { 0, 18, 20, 19,  0}, //BL12
    { 0, 27, 21, 28,  0}, //BL13
    { 0, 27, 21, 28,  0}, //BL14
    { 0, 27, 21, 28,  0}, //BL15
    { 0, 27, 21, 28,  0}, //BL16
    { 0, 27, 21, 28,  0}, //BL17
    { 0, 24, 26, 25,  0}, //BL18
    { 0, 24, 26, 25,  0}, //BL19
    { 0, 24, 26, 25,  0}, //BL20
    {22, 24, 26, 25, 23}, //BL21
    {22, 24, 26, 25, 23}, //BL22
    { 0,  0,  0,  0,  0}, //BL23
    { 1,  3,  5,  4,  2}, //BL24
    { 1,  3,  5,  4,  2}, //BL25
    { 1,  3,  5,  4,  2}, //BL26
    { 1,  3,  5,  4,  2}, //BL27
    { 7, 14,  6, 15,  8}, //BL28
    { 7, 14,  6, 15,  8}, //BL29
    { 7, 14,  6, 15,  8}  //BL30
  };

  TH2F *hMap = new TH2F("hMap","Map for trigger units;Backleg;Module",30,0.5,30.5,5,0.5,5.5);
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  hMap->SetBinContent(i+1, j+1, kTrgID[i][j]);
	}
    }
  c = draw2D(hMap,"",0.04,true,"TEXT");
  for(int i=0; i<30; i++)
    {
      TLine *line = GetLine(0.5+i,0.5,0.5+i,5.5,1);
      line->Draw();
    }
  for(int j=0; j<5; j++)
    {
      TLine *line = GetLine(0.5, 0.5+j, 30.5, 0.5+j,1);
      line->Draw();
    }
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/TrigUnitMap.pdf"));
}
