TFile *f;
const int year = YEAR;
TString run_cfg_name;

//================================================
void study_TpcDistortion()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  ana();
}

//================================================
void ana(const bool savePlot = 0)
{
  gStyle->SetOptFit(0);
  f = TFile::Open(Form("./output/Pico.Run14.AuAu200.study.root"),"read");

  const char* cent_name[4]  = {"0060","0020","2040","4060"};
  const char* cent_title[4] = {"00-60%","00-20%","20-40%","40-60%"};
  const char* lumi_name[4]  = {"all","low","mid","high"};

  const int nTry = 10;
  TH2F *mhJpsiMassVsPtUL[4][4][nTry];
  TH2F *mhJpsiMassVsPtLS[4][4][nTry];
  for(int j=1; j<4; j++)
    {
      for(int k=1; k<4; k++)
	{
	  for(int d=0; d<nTry; d++)
	    {
	      mhJpsiMassVsPtUL[j][k][d] = (TH2F*)f->Get(Form("mhJpsiMassVsPtUL_di_mu_%s_%s_try%d",cent_name[j],lumi_name[k],d));
	      mhJpsiMassVsPtLS[j][k][d] = (TH2F*)f->Get(Form("mhJpsiMassVsPtLS_di_mu_%s_%s_try%d",cent_name[j],lumi_name[k],d));
	      mhJpsiMassVsPtUL[j][k][d]->Sumw2();
	      mhJpsiMassVsPtLS[j][k][d]->Sumw2(); 
	    }
	}
    }

  for(int d=0; d<nTry; d++)
    {
      for(int j=0; j<4; j++)
	{
	  for(int k=0; k<4; k++)
	    {
	      if(j!=0 && k!=0) continue;
	      mhJpsiMassVsPtUL[j][k][d] = (TH2F*)mhJpsiMassVsPtUL[1][1][d]->Clone(Form("mhJpsiMassVsPtUL_di_mu_%s_%s_try%d",cent_name[j],lumi_name[j],d));
	      mhJpsiMassVsPtUL[j][k][d]->Reset();
	      mhJpsiMassVsPtLS[j][k][d] = (TH2F*)mhJpsiMassVsPtLS[1][1][d]->Clone(Form("mhJpsiMassVsPtLS_di_mu_%s_%s_try%d",cent_name[j],lumi_name[j],d));
	      mhJpsiMassVsPtLS[j][k][d]->Reset();
	    }
	}

      for(int j=1; j<4; j++)
	{
	  for(int k=1; k<4; k++)
	    {
	      mhJpsiMassVsPtUL[0][k][d]->Add(mhJpsiMassVsPtUL[j][k][d]);
	      mhJpsiMassVsPtUL[j][0][d]->Add(mhJpsiMassVsPtUL[j][k][d]);
	      mhJpsiMassVsPtUL[0][0][d]->Add(mhJpsiMassVsPtUL[j][k][d]);
	      mhJpsiMassVsPtLS[0][k][d]->Add(mhJpsiMassVsPtLS[j][k][d]);
	      mhJpsiMassVsPtLS[j][0][d]->Add(mhJpsiMassVsPtLS[j][k][d]);
	      mhJpsiMassVsPtLS[0][0][d]->Add(mhJpsiMassVsPtLS[j][k][d]);
	    }
	}
    }

  // =========================
  const int icent = 0, ilumi = 0;
  const int nPtBins = 3;
  const double xPtBins[nPtBins+1] = {0,2,5,10};
  TH1F *hJpsiMean[nPtBins];
  TH1F *hJpsiSigma[nPtBins];
  for(int ipt=0; ipt<nPtBins; ipt++)
    {
      hJpsiMean[ipt] = new TH1F(Form("hJpsiMean_pt%1.0f-%1.0f",xPtBins[ipt],xPtBins[ipt+1]),"Mean of J/#Psi peak;A (%);Mean (GeV/c^{2})",nTry,-0.55,0.45);
      hJpsiSigma[ipt] = new TH1F(Form("hJpsiSigma_pt%1.0f-%1.0f",xPtBins[ipt],xPtBins[ipt+1]),"Sigma of J/#Psi peak;A (%);#sigma (MeV/c^{2})",nTry,-0.55,0.45);

      TCanvas *c = new TCanvas(Form("JpsiMass_%d",ipt),Form("JpsiMass_pt%1.0f-%1.0f",xPtBins[ipt],xPtBins[ipt+1]),1200,700);
      c->Divide(5,2);
      for(int d=0; d<nTry; d++)
	{
	  int lowbin = mhJpsiMassVsPtUL[icent][ilumi][d]->GetXaxis()->FindFixBin(xPtBins[ipt]+1e-4);
	  int upbin  = mhJpsiMassVsPtUL[icent][ilumi][d]->GetXaxis()->FindFixBin(xPtBins[ipt+1]-1e-4);
	  TH1F *hUL = (TH1F*)mhJpsiMassVsPtUL[icent][ilumi][d]->ProjectionY(Form("JpsiMassUL_pt%1.0f-%1.0f_%s_%s_try%d",xPtBins[ipt],xPtBins[ipt+1],cent_name[icent],lumi_name[ilumi],d),lowbin,upbin);
	  hUL->SetMarkerStyle(21);
	  hUL->SetMarkerColor(1);
	  hUL->SetLineColor(1);
	  TH1F *hLS = (TH1F*)mhJpsiMassVsPtLS[icent][ilumi][d]->ProjectionY(Form("JpsiMassLS_pt%1.0f-%1.0f_%s_%s_try%d",xPtBins[ipt],xPtBins[ipt+1],cent_name[icent],lumi_name[ilumi],d),lowbin,upbin);
	  hLS->SetMarkerStyle(24);
	  hLS->SetMarkerColor(1);
	  hLS->SetLineColor(1);
	  
	  c->cd(d+1);
	  gPad->SetRightMargin(0.05);
	  hUL->Rebin(2);
	  hLS->Rebin(2);
	  hUL->Add(hLS,-1);
	  hUL->GetXaxis()->SetRangeUser(2.5,4);
	  hUL->SetTitle("");
	  hUL->SetMaximum(1.5*hUL->GetMaximum());
	  hUL->Draw("PE");

	  TF1 *func = new TF1(Form("Fit_pt%1.0f-%1.0f_try%d",xPtBins[ipt],xPtBins[ipt+1],d),"gausn(0)+pol3(3)",2.5,3.9);
	  func->SetParameter(0,hUL->GetMaximum());
	  func->SetParameter(1,3.1);
	  func->SetParameter(2,0.08);
	  hUL->Fit(func,"IR0Q");
	  func->SetLineColor(6);
	  func->SetLineStyle(1);
	  func->Draw("sames");
	  hJpsiMean[ipt]->SetBinContent(d+1, func->GetParameter(1));
	  hJpsiMean[ipt]->SetBinError(d+1, func->GetParError(1));
	  hJpsiSigma[ipt]->SetBinContent(d+1, func->GetParameter(2)*1000);
	  hJpsiSigma[ipt]->SetBinError(d+1, func->GetParError(2)*1000);

	  TPaveText *t1 = GetPaveText(0.2,0.35,0.8,0.85,0.065);
	  t1->AddText(Form("A = %2.2f%%",(0 + 0.001 * (d-5))*100));
	  t1->SetTextColor(2);
	  t1->Draw();

	  TPaveText *t2 = GetPaveText(0.6,0.75,0.6,0.85,0.06);
	  t2->SetTextAlign(11);
	  t2->AddText(Form("%s",cent_title[icent]));
	  t2->AddText(Form("prod_%s",lumi_name[ilumi]));
	  t2->AddText(Form("N_{J/#Psi} = %1.0f#pm%1.0f",func->GetParameter(0),func->GetParError(0)));
	  t2->AddText(Form("#mu = %4.3f",func->GetParameter(1)));
	  t2->AddText(Form("#sigma = %4.3f",func->GetParameter(2)));
	  t2->SetTextColor(1);
	  t2->Draw();

	  TPaveText *t3 = GetTitleText(Form("%1.0f < p_{T}^{J/#Psi} < %1.0f GeV/c",xPtBins[ipt],xPtBins[ipt+1]),0.065);
	  t3->Draw();
	}
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/study_TpcDistortion/ScanA_FitJpsiPeak_pt%1.0f-%1.0f_%s_%s.pdf",run_type,xPtBins[ipt],xPtBins[ipt+1],cent_name[icent],lumi_name[ilumi]));
    }

  TList *list = new TList();
  TString legName[nPtBins];
  for(int ipt=0; ipt<nPtBins; ipt++)
    {
      hJpsiMean[ipt]->SetMarkerStyle(20+ipt);
      hJpsiMean[ipt]->SetMarkerColor(color[ipt]);
      hJpsiMean[ipt]->SetLineColor(color[ipt]);
      hJpsiMean[ipt]->GetYaxis()->SetTitleOffset(1.2);
      list->Add(hJpsiMean[ipt]);
      legName[ipt] = Form("%1.0f < p_{T}^{J/#Psi} < %1.0f GeV/c",xPtBins[ipt],xPtBins[ipt+1]);
    }
  c = drawHistos(list,Form("JpsiMean_ScanA_%s_%s",cent_name[icent],lumi_name[ilumi]),"",false,0,120,true,3.08,3.14,false,true,legName,true,Form("%s, prod_%s",cent_title[icent],lumi_name[ilumi]),0.5,0.8,0.62,0.88,true,0.04,0.04,false,0.01,false,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/study_TpcDistortion/ScanA_FitJpsiMean_%s_%s.pdf",run_type,cent_name[icent],lumi_name[ilumi]));

  list->Clear();
  for(int ipt=0; ipt<nPtBins; ipt++)
    {
      hJpsiSigma[ipt]->SetMarkerStyle(20+ipt);
      hJpsiSigma[ipt]->SetMarkerColor(color[ipt]);
      hJpsiSigma[ipt]->SetLineColor(color[ipt]);
      hJpsiSigma[ipt]->GetYaxis()->SetTitleOffset(1.2);
      list->Add(hJpsiSigma[ipt]);
    }
  c = drawHistos(list,Form("JpsiSigma_ScanA_%s_%s",cent_name[icent],lumi_name[ilumi]),"",false,0,120,true,20,90,false,true,legName,true,Form("%s, prod_%s",cent_title[icent],lumi_name[ilumi]),0.5,0.8,0.62,0.88,true,0.04,0.04,false,0.01,false,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/study_TpcDistortion/ScanA_FitJpsiSigma_%s_%s.pdf",run_type,cent_name[icent],lumi_name[ilumi]));
}
