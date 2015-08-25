TFile *f;
const char *run_config = "EvtAct.";
const Bool_t iPico = 1;
const int year = 2013;
TString run_cfg_name;
const Double_t low_mass = 2.8;
const Double_t high_mass = 3.3;
const double pt_cut = 1.5;
const int nMultBins = 4;
const int low_tofMult[nMultBins] = {0,0,4,7};
const int high_tofMult[nMultBins] = {30,3,6,30};
TString fileName;

//================================================
void sys_Activity()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.15); 
  run_type = "Run13_pp500";

  //effVsTofMult();
  pTshape();
}

//================================================
void pTshape(const bool saveHisto = 1, const bool saveRoot = 0)
{
  TFile *fWeight = TFile::Open("Rootfiles/HTppSpectrum.root","read");
  TF1 *func = (TF1*)fWeight->Get("function");
  //cout << func->GetExpFormula().Data() << endl;
  TCanvas *c = new TCanvas("func","func",800,600);
  gPad->SetLogy();
  func->SetTitle("z_{T} distribution of J/#psi;z_{T} ;A.U.");
  func->Draw();

  TFile *fHT = TFile::Open("Rootfiles/Spectrum_in_bin.root","read");
  TGraphErrors *gSpectrum[5];
  const char *name_temp [5] = {"gall","gbin1","gbin2","gbin3","gbin4"};
  double x, y;
  double firstx;
  for(int i=0; i<5; i++)
    {
      gSpectrum[i] = (TGraphErrors*)fHT->Get(name_temp[i]);
      gSpectrum[i]->GetPoint(0,x,y);
      firstx = x/mean_pt;
      double scale = func->Eval(firstx)/y;
      for(int ipoint=0; ipoint<gSpectrum[i]->GetN(); ipoint++)
	{
	  gSpectrum[i]->GetPoint(ipoint,x,y);
	  gSpectrum[i]->SetPoint(ipoint,x/mean_pt,y*scale);
	  gSpectrum[i]->SetPointError(ipoint,gSpectrum[i]->GetErrorYlow(ipoint)*scale,
				      gSpectrum[i]->GetErrorYhigh(ipoint)*scale);
	}
      gSpectrum[i]->Draw("sames PE");
    }

  const int nFunc = 3;
  TF1 *myFunc[nFunc];
  const double index[nFunc] = {3.92,3.2,5};
  for(int i=0; i<nFunc; i++)
    {
      myFunc[i] = new TF1(Form("myFunc_%d",i),crossSection,0,20,2);
      myFunc[i]->SetParameter(0,index[i]);
      myFunc[i]->SetParameter(1,1);
      myFunc[i]->SetParameter(1,func->Eval(firstx)/myFunc[i]->Eval(firstx));
      myFunc[i]->SetLineColor(color[i+1]);
      if(i>0) myFunc[i]->Draw("sames");
    }
  leg = new TLegend(0.5,0.5,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(func,"World-data parametrization","L");
  leg->AddEntry(myFunc[1],"Shape variation 1","L");
  leg->AddEntry(myFunc[2],"Shape variation 2","L");
  leg->AddEntry(gSpectrum[0],"Data: TofMult >= 0","P");
  leg->AddEntry(gSpectrum[1],"Data: 1<= TofMult <= 4","P");
  leg->AddEntry(gSpectrum[2],"Data: 5<= TofMult <= 8","P");
  leg->AddEntry(gSpectrum[3],"Data: 9<= TofMult <= 12","P");
  leg->AddEntry(gSpectrum[4],"Data: TofMult >= 13","P");
  leg->Draw();
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.Jpsi_pTShape.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.Jpsi_pTShape.png",run_type));
    }

  f = TFile::Open("./output/Run13.pp500.jpsi.Embed.root","read");
  THnSparseF *hJpsiUS = (THnSparseF*)f->Get("hJpsiInfo_di_mu");
  hJpsiUS->GetAxis(9)->SetRange(1,350);
  TH1F *hJpsiPt[nFunc][2];
  TString name[2] = {"MCinput","MTDreco"};
  for(int j=0; j<2; j++)
    {
      hJpsiUS->GetAxis(7)->SetRange(4+j*2,4+j*2);
      if(j>0)
	{
	  hJpsiUS->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
	  hJpsiUS->GetAxis(4)->SetRangeUser(pt_cut+0.01,100);
	  hJpsiUS->GetAxis(5)->SetRangeUser(1+0.01,100);
	}
      hJpsiPt[0][j] = (TH1F*)hJpsiUS->Projection(1);
      hJpsiPt[0][j]->Sumw2();
      hJpsiPt[0][j]->SetName(Form("JpsiPt_%s_0",name[j].Data()));
      hJpsiPt[0][j]->SetTitle(Form("x_{T} distribution j/psi (%s, pt1 > %1.1f GeV/c)",name[j].Data(),pt_cut)); 
      hJpsiUS->GetAxis(0)->SetRange(0,-1);
      hJpsiUS->GetAxis(4)->SetRange(0,-1);
      hJpsiUS->GetAxis(5)->SetRange(0,-1);
      hJpsiUS->GetAxis(7)->SetRange(0,-1);
    }
  for(int i=1; i<nFunc; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hJpsiPt[i][j] = (TH1F*)hJpsiPt[0][j]->Clone(Form("JpsiPt_%s_%d",name[j].Data(),i));
	}
    }

  TList *list = new TList;
  double nJpsi[nFunc][2];
  double nJpsiError[nFunc][2];
  for(int i=0; i<nFunc; i++)
    {
      list->Clear();
      for(int j=0; j<2; j++)
	{
	  nJpsi[i][j] = 0;
	  nJpsiError[i][j] = 0; 
	  for(int bin=1; bin<=hJpsiPt[i][j]->GetNbinsX(); bin++)
	    {
	      double weight = myFunc[i]->Eval(hJpsiPt[i][j]->GetBinCenter(bin)/mean_pt);
	      hJpsiPt[i][j]->SetBinContent(bin,hJpsiPt[i][j]->GetBinContent(bin)*weight);
	      hJpsiPt[i][j]->SetBinError(bin,hJpsiPt[i][j]->GetBinError(bin)*weight);
	      //if(bin==1) cout << bin << "  " << hJpsiPt[i][j]->GetBinContent(bin) << endl;
	      nJpsi[i][j] += hJpsiPt[i][j]->GetBinContent(bin);
	      nJpsiError[i][j] += TMath::Power(hJpsiPt[i][j]->GetBinError(bin),2);
	    }
	  nJpsiError[i][j] = TMath::Sqrt(nJpsiError[i][j]);
	  hJpsiPt[i][j]->Rebin(2);
	  list->Add(hJpsiPt[i][j]);
	}
      cout << nJpsi[i][1] << "  " << nJpsi[i][0] << endl;
      c = sysCompare(list,Form("JpsiEff_func_%d",i),Form("p_{T} distribution of J/#psi"),"J/#psi efficiency;p_{T} (GeV/c);Efficiency",kTRUE,0,10,kTRUE,0.01,2e5,kFALSE,0.7,1.3,kTRUE,kTRUE,name,kTRUE,"",0.6,0.75,0.7,0.88,kTRUE);
    }

  TH1F *hJpsiEff[nFunc];
  list->Clear();
  TString legName[nFunc];

  double eff[nFunc];
  double eff_err[nFunc];
  for(int i=0; i<nFunc; i++)
    {
      eff[i] = nJpsi[i][1]/nJpsi[i][0];
      eff_err[i] = eff[i] * TMath::Sqrt(nJpsiError[i][0]*nJpsiError[i][0]/nJpsi[i][0]/nJpsi[i][0]+nJpsiError[i][1]*nJpsiError[i][1]/nJpsi[i][1]/nJpsi[i][1]);
      printf("Func %d: eff = %4.3e +/- %4.3e, %3.2f%%\n",i,eff[i],eff_err[i],eff_err[i]/eff[i]*100);
      hJpsiPt[i][0]->Rebin(5);
      hJpsiPt[i][1]->Rebin(5);
      hJpsiEff[i] = (TH1F*)hJpsiPt[i][1]->Clone(Form("JpsiEff_Func%d",i));
      hJpsiEff[i]->Divide(hJpsiPt[i][0]);
      list->Add(hJpsiEff[i]);
      legName[i] = Form("Shape %d",i);
    }

  c = sysCompare(list,Form("Sys_JpsiEff_PtShape"),"J/#psi efficiency;x_{T} (GeV/c);Efficiency","Relative difference to MB",kTRUE,0,10,kFALSE,0.1,3e5,kTRUE,0.5,1.5,kTRUE,kTRUE,legName,kTRUE,"",0.55,0.7,0.2,0.4,kTRUE);
  double rel_eff[nFunc-1];
  double rel_eff_err[nFunc-1];
  for(int i=0; i<nFunc-1; i++)
    {
      rel_eff[i] = eff[i+1]/eff[0];
      rel_eff_err[i] = rel_eff[i] * TMath::Sqrt(eff_err[0]*eff_err[0]/eff[0]/eff[0]+eff_err[i+1]*eff_err[i+1]/eff[i+1]/eff[i+1]);
      printf("TofMult [%d,%d]: rel. eff. = %4.3e +/- %4.3e, %3.2f%%\n",low_tofMult[i+1],high_tofMult[i+1],rel_eff[i],rel_eff_err[i],rel_eff_err[i]/rel_eff[i]*100);
    }
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.JpsiEff_vs_pTShape.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.JpsiEff_vs_pTShape.png",run_type));
    }

  if(saveRoot)
    {
      TH1F *hRelEff = new TH1F("hRelEff","Relative efficiency",nMultBins,0,nMultBins);
      for(int bin=1; bin<=nMultBins; bin++)
	{
	  hRelEff->SetBinContent(bin,1);
	  hRelEff->SetBinError(bin,0.1);
	}
      TFile *fout = TFile::Open("Rootfiles/Sys.Run13.EvtAct.Eff.pTShape.root","recreate");
      hRelEff->Write();
      fout->Close();
    }

}

//================================================
double crossSection(double *x, double *par)
{
  double n = 3.92;
  double b = TMath::Gamma(1.5)*TMath::Gamma(n-1.5)/TMath::Gamma(n-1);
  double a = 2 * b * b * (n-1);
  return a * TMath::Power(1+b*b*x[0]*x[0],par[0]*-1) * 2451./2 * 1e-3 * par[1];
}


//================================================
void effVsTofMult(const bool saveHisto = 0, const bool saveRoot = 0)
{
  f = TFile::Open("./output/Run13.pp500.jpsi.Embed.root","read");
  TList *list = new TList;
  TString legName[nMultBins];

  TH2F *hTofMultVsZdcRate = (TH2F*)f->Get("mhTofMultVsZdcRate_di_mu");
  draw2D(hTofMultVsZdcRate);
  TH1F *hZdcRate[nMultBins];
  for(int i=0; i<nMultBins; i++)
    {
      hTofMultVsZdcRate->GetYaxis()->SetRange(low_tofMult[i]+1,high_tofMult[i]+1);
      hZdcRate[i] = (TH1F*)hTofMultVsZdcRate->ProjectionX(Form("hZdcRate_%d",i));
      hZdcRate[i]->Sumw2();
      hZdcRate[i]->Rebin(20);
      hZdcRate[i]->Scale(1./hZdcRate[i]->Integral());
      list->Add(hZdcRate[i]);
      legName[i] = Form("%d <= TofMult <= %d",low_tofMult[i],high_tofMult[i]);
    }
  c = drawHistos(list,Form("ZdcRate_TofMult"),"ZDC rate",kFALSE,0,10,kFALSE,0.1,3e5,kFALSE,kTRUE,legName,kTRUE,"",0.55,0.7,0.2,0.4,kFALSE);

  //return;
  
  THnSparseF *hJpsiUS = (THnSparseF*)f->Get("hJpsiInfo_di_mu");
  hJpsiUS->GetAxis(9)->SetRange(1,350);
  
  TH1F *hJpsiPt[nMultBins][2];
  TString name[2] = {"MCinput","MTDreco"};
  for(int i=0; i<nMultBins; i++)
    {
      hJpsiUS->GetAxis(8)->SetRange(low_tofMult[i]+1,high_tofMult[i]+1);
      for(int j=0; j<2; j++)
	{
	  hJpsiUS->GetAxis(7)->SetRange(4+j*2,4+j*2);
	  if(j>0)
	    {
	      hJpsiUS->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
	      hJpsiUS->GetAxis(4)->SetRangeUser(pt_cut+0.01,100);
	      hJpsiUS->GetAxis(5)->SetRangeUser(1+0.01,100);
	    }
	  hJpsiPt[i][j] = (TH1F*)hJpsiUS->Projection(1);
	  hJpsiPt[i][j]->Sumw2();
	  hJpsiPt[i][j]->SetName(Form("JpsiPt_%s_TofMult_%d_%d",name[j].Data(),low_tofMult[i],high_tofMult[i]));
	  hJpsiPt[i][j]->SetTitle(Form("p_{T} distribution j/psi (%s, pt1 > %1.1f GeV/c)",name[j].Data(),pt_cut)); 
	  hJpsiUS->GetAxis(0)->SetRange(0,-1);
	  hJpsiUS->GetAxis(4)->SetRange(0,-1);
	  hJpsiUS->GetAxis(5)->SetRange(0,-1);
	  hJpsiUS->GetAxis(7)->SetRange(0,-1);
	}
      hJpsiUS->GetAxis(8)->SetRange(0,-1);
    }

  TFile *fWeight = TFile::Open("Rootfiles/HTppSpectrum.root","read");
  TF1 *func = (TF1*)fWeight->Get("function");
  

  double nJpsi[nMultBins][2];
  double nJpsiError[nMultBins][2];
  for(int i=0; i<nMultBins; i++)
    {
      list->Clear();
      for(int j=0; j<2; j++)
	{
	  for(int bin=1; bin<=hJpsiPt[i][j]->GetNbinsX(); bin++)
	    {
	      double weight = func->Eval(hJpsiPt[i][j]->GetBinCenter(bin)/mean_pt);
	      hJpsiPt[i][j]->SetBinContent(bin,hJpsiPt[i][j]->GetBinContent(bin)*weight);
	      hJpsiPt[i][j]->SetBinError(bin,hJpsiPt[i][j]->GetBinError(bin)*weight);
	      //cout << bin << "  " << hJpsiPt[i][j]->GetBinContent(bin) << endl;
	      nJpsi[i][j] += hJpsiPt[i][j]->GetBinContent(bin);
	      nJpsiError[i][j] += TMath::Power(hJpsiPt[i][j]->GetBinError(bin),2);
	    }
	  nJpsiError[i][j] = TMath::Sqrt(nJpsiError[i][j]);
	  hJpsiPt[i][j]->Rebin(2);
	  list->Add(hJpsiPt[i][j]);
	}
      c = sysCompare(list,Form("JpsiEff_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]),Form("p_{T} distribution of J/#psi"),"J/#psi efficiency;p_{T} (GeV/c);Efficiency",kTRUE,0,10,kTRUE,0.01,2e5,kFALSE,0.7,1.3,kTRUE,kTRUE,name,kTRUE,"",0.6,0.75,0.7,0.88,kTRUE);
      if(saveHisto) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.JpsiEff_TofMult_%d_%d.pdf",run_type,low_tofMult[i],high_tofMult[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.JpsiEff_TofMult_%d_%d.png",run_type,low_tofMult[i],high_tofMult[i]));
	}
    }

  TH1F *hJpsiEff[nMultBins];
  list->Clear();


  double eff[nMultBins];
  double eff_err[nMultBins];
  for(int i=0; i<nMultBins; i++)
    {
      cout << nJpsi[i][1] << "  " << nJpsi[i][0] << endl;
      eff[i] = nJpsi[i][1]/nJpsi[i][0];
      eff_err[i] = eff[i] * TMath::Sqrt(nJpsiError[i][0]*nJpsiError[i][0]/nJpsi[i][0]/nJpsi[i][0]+nJpsiError[i][1]*nJpsiError[i][1]/nJpsi[i][1]/nJpsi[i][1]);
      printf("TofMult [%d,%d]: eff = %4.3e +/- %4.3e, %3.2f%%\n",low_tofMult[i],high_tofMult[i],eff[i],eff_err[i],eff_err[i]/eff[i]*100);
      hJpsiPt[i][0]->Rebin(5);
      hJpsiPt[i][1]->Rebin(5);
      hJpsiEff[i] = (TH1F*)hJpsiPt[i][1]->Clone(Form("JpsiEff_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]));
      hJpsiEff[i]->Divide(hJpsiPt[i][0]);
      list->Add(hJpsiEff[i]);
    }

  c = sysCompare(list,Form("Sys_JpsiEff_TofMult"),"J/#psi efficiency;p_{T} (GeV/c);Efficiency","Relative difference to MB",kTRUE,0,10,kFALSE,0.1,3e5,kTRUE,0.5,1.5,kTRUE,kTRUE,legName,kTRUE,"",0.55,0.7,0.2,0.4,kTRUE);
  double rel_eff[nMultBins-1];
  double rel_eff_err[nMultBins-1];
  for(int i=0; i<nMultBins-1; i++)
    {
      rel_eff[i] = eff[i+1]/eff[0];
      rel_eff_err[i] = rel_eff[i] * TMath::Sqrt(eff_err[0]*eff_err[0]/eff[0]/eff[0]+eff_err[i+1]*eff_err[i+1]/eff[i+1]/eff[i+1]);
      printf("TofMult [%d,%d]: rel. eff. = %4.3e +/- %4.3e, %3.2f%%\n",low_tofMult[i+1],high_tofMult[i+1],rel_eff[i],rel_eff_err[i],rel_eff_err[i]/rel_eff[i]*100);
    }
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.JpsiEff_vs_TofMult.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Activity/Sys.JpsiEff_vs_TofMult.png",run_type));
    }

  if(saveRoot)
    {
      TH1F *hRelEff = new TH1F("hRelEff","Relative efficiency",nMultBins-1,0,nMultBins-1);
      for(int bin=1; bin<=nMultBins; bin++)
	{
	  hRelEff->SetBinContent(bin,rel_eff[bin-1]);
	  hRelEff->SetBinError(bin,rel_eff_err[bin-1]);
	}
      TFile *fout = TFile::Open("Rootfiles/Sys.Run13.EvtAct.Eff.TofMult.root","recreate");
      hRelEff->Write();
      fout->Close();
    }
}
