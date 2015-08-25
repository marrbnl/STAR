TFile *f;
const char *run_config = "EvtAct.";
const Bool_t iPico = 1;
const int year = 2013;
TString run_cfg_name;
const Double_t low_mass = 2.8;
const Double_t high_mass = 3.3;
const double pt_cut = 1.5;
// const int nMultBins = 4;
// const int low_tofMult[nMultBins] = {0,0,4,7};
// const int high_tofMult[nMultBins] = {30,3,6,30};
// const int nMultBins = 4;
// const int low_tofMult[nMultBins] = {0,1,3,6};
// const int high_tofMult[nMultBins] = {30,2,5,30};
const int nMultBins = 4;
const int low_tofMult[nMultBins] = {0,1,4,8};
const int high_tofMult[nMultBins] = {30,3,7,30};
TString fileName;
const bool cut_ranking = 1;
const bool cut_vtxdz = 0;
const double jpsi_pt_cut_low = 0;
const double jpsi_pt_cut_high = 10;

//================================================
void ana_Activity()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.15); 

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) fileName = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      fileName = Form("Run13.pp500.jpsi.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) fileName = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      fileName = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all events: %4.4e\n",hStat->GetBinContent(1));
  printf("dimuon  events: %4.4e\n",hStat->GetBinContent(2));

  run_cfg_name = run_config;

  //ana();
  makeCorr();
  //makeYield();
  //compare();
}


//================================================
void makeYield(const bool saveHisto = 0, const bool saveRoot = 1, const int sys = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhJpsiCutStudy_%s",trigName[kTrigType]));
  hn->GetAxis(1)->SetRangeUser(pt_cut+0.001,100); // pt1 cut
  hn->GetAxis(2)->SetRange(1,1); // unlike-sign
  hn->GetAxis(6)->SetRange(1,350);
  hn->GetAxis(7)->SetRangeUser(jpsi_pt_cut_low+1e-6,jpsi_pt_cut_high-1e-6);
  if(cut_ranking)   hn->GetAxis(3)->SetRange(3,3);
  if(cut_vtxdz)     hn->GetAxis(4)->SetRangeUser(-5.2,9.9);
  

  TH1F *hInvMass[nMultBins];
  TH1F *hJpsiYield = new TH1F("hJpsiYield","Jpsi yield",nMultBins,0,nMultBins);
  double nAll[nMultBins];
  double allErr[nMultBins];
  double nBkg[nMultBins];
  double bkgErr[nMultBins];
  double nSignal[nMultBins];
  double sigErr[nMultBins];
  double check = 0;
  for(int i=0; i<nMultBins; i++)
    {
      hn->GetAxis(5)->SetRange(low_tofMult[i]+1,high_tofMult[i]+1);
      hInvMass[i] = (TH1F*)hn->Projection(0);
      hInvMass[i]->SetName(Form("hInvMass_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]));
      hInvMass[i]->SetMarkerStyle(20);
      hInvMass[i]->SetMarkerColor(2);
      hInvMass[i]->Rebin(4);
      hInvMass[i]->GetXaxis()->SetRangeUser(2,4);

      //TF1 *func = new TF1(Form("func_%d",i),"gaus(0)+expo(3)+pol0(5)",2.5,3.5);
      //if(i!=nMultBins-1)func->SetParameter(5,-30);
      TF1 *func;
      if(sys==0)
	{
	  double variation = 0;
	  func = new TF1(Form("func_%d",i),"gaus(0)+expo(3)+pol0(5)",2.5-variation,3.5+variation);
	  func->SetParLimits(1,3.05,3.15);
	  func->SetParLimits(2,0.04,0.1);
	}
      else if(sys==1)
	{
	  func = new TF1(Form("func_%d",i),"gaus(0)+pol3(3)",2,4);
	  func->SetParLimits(1,3.05,3.15);
	  func->SetParLimits(2,0.04,0.1);
	}
      else if(sys==2)
	{
	  func = new TF1(Form("func_%d",i),"gaus(0)+pol1(3)",2.5,3.5);
	  func->SetParameter(1,3.1);
	  func->SetParameter(2,0.1);
	}
      hInvMass[i]->Fit(func,"IR0Q");
      TF1 *func1 = new TF1(Form("signal_%d",i),"gaus",func->GetXmin(),func->GetXmax());
      func1->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
      c = draw1D(hInvMass[i],Form("Invariant mass of di-muon pairs (%d<=TofMult<=%d);M_{#mu#mu} (GeV/c^{2});Counts",low_tofMult[i],high_tofMult[i]));
      func->SetLineColor(4);
      func->Draw("sames");
      func1->SetLineColor(2);

      // Singal counts
      int low_bin = hInvMass[i]->FindFixBin(low_mass+0.001);
      int high_bin = hInvMass[i]->FindFixBin(high_mass-0.001);
      nAll[i] = hInvMass[i]->Integral(low_bin,high_bin);
      nAll[i] -= hInvMass[i]->GetBinContent(low_bin) * (low_mass-hInvMass[i]->GetXaxis()->GetBinLowEdge(low_bin))/hInvMass[i]->GetBinWidth(low_bin);
      nAll[i] -= hInvMass[i]->GetBinContent(high_bin) * (hInvMass[i]->GetXaxis()->GetBinUpEdge(high_bin)-high_mass)/hInvMass[i]->GetBinWidth(high_bin);
      allErr[i] = TMath::Sqrt(nAll[i]);
      nBkg[i] = (func->Integral(low_mass,high_mass) - func1->Integral(low_mass,high_mass)) * 1./hInvMass[i]->GetBinWidth(1);
      bkgErr[i] = TMath::Sqrt(nBkg[i]);
      nSignal[i] = nAll[i] - nBkg[i];
      sigErr[i] = TMath::Sqrt(nAll[i]+nBkg[i]);
      printf("TofMult: [%d,%d]\n",low_tofMult[i],high_tofMult[i]);
      printf("Unlike-sign: %2.3f +/- %2.3f\n",nAll[i],allErr[i]);
      printf("Background: %2.3f +/- %2.3f\n",nBkg[i],bkgErr[i]);
      printf("Signal: %2.3f +/- %2.3f\n",nSignal[i],sigErr[i]);
      printf("Ratio: %2.3f\n",nSignal[i]/nSignal[0]);
      TPaveText *signif = GetPaveText(0.13,0.3,0.12,0.3);
      signif->SetTextAlign(11);
      signif->SetTextFont(62);
      signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
      signif->AddText(Form("Signal = %3.1f #pm %3.1f",nSignal[i],sigErr[i]));
      signif->AddText(Form("S/B = %1.0f/%1.0f = 1:%3.0f",nSignal[i],nBkg[i],nBkg[i]/nSignal[i]));
      signif->Draw();
      if(i>0) check += nSignal[i];
      hJpsiYield->SetBinContent(i+1,nSignal[i]);
      hJpsiYield->SetBinError(i+1,sigErr[i]);   
      int index = (sys>0? 1 : 0);
      char *name = Form("%s%d.",save_name[index],sys);
      if(saveHisto) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%s%sFit_Jpsi_TofMult_%d_%d.pdf",run_type,name,run_cfg_name.Data(),low_tofMult[i],high_tofMult[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%s%sFit_Jpsi_TofMult_%d_%d.png",run_type,name,run_cfg_name.Data(),low_tofMult[i],high_tofMult[i]));
	}
    }
  printf("Check: %3.2f =? %3.2f\n",nSignal[0],check);

  if(saveRoot)
    {
      TFile *fout;
      if(jpsi_pt_cut_low==0 && jpsi_pt_cut_high==10 && !cut_ranking && !cut_vtxdz)
	fout = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.yield.default.root"),"recreate");
      else
	fout = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.yield.Rank%d.VtxDz%d.pt%1.0f-%1.0f.root",cut_ranking,cut_vtxdz,jpsi_pt_cut_low,jpsi_pt_cut_high),"recreate");
      hJpsiYield->Write();
    }
}


//================================================
void makeCorr(bool saveHisto = 0, const bool saveRoot = 0, const bool sys = 0)
{
  if(sys) saveHisto = 0;

// Event Multiplicity
  TFile *fMB = 0;
  if(!cut_vtxdz) fMB = TFile::Open("output/Run13.pp500.jpsi.VPDMB.root","read");
  else           fMB = TFile::Open("output/Run13.pp500.jpsi.VPDMB.VtxDz.root","read");
  THnSparseF *hnMult = (THnSparseF*)fMB->Get("mhEventMult_qa_mb");
  hnMult->GetAxis(0)->SetRange(1,350);
  if(cut_ranking)   hnMult->GetAxis(4)->SetRange(3,3);
  if(cut_vtxdz)     hnMult->GetAxis(5)->SetRangeUser(-5.2,9.9);
  TH1F *hTofMult[nMultBins];
  double meanTofMult[nMultBins];
  double eventFraction[nMultBins];
  cout << "=== Raw ===" << endl;
  for(int i=0; i<nMultBins; i++)
    {
      hnMult->GetAxis(2)->SetRange(low_tofMult[i]+1,high_tofMult[i]+1);
      hTofMult[i] = (TH1F*)hnMult->Projection(2);
      hTofMult[i]->SetName(Form("Raw_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]));
      hTofMult[i]->Sumw2();
      meanTofMult[i] = hTofMult[i]->GetMean();
      eventFraction[i] = hTofMult[i]->Integral(1,30)/hTofMult[0]->Integral(1,30);
      printf("TofMult: [%d,%d]\n",low_tofMult[i],high_tofMult[i]);
      printf("<TofMult> = %3.2f +/- %3.2e\n",meanTofMult[i],hTofMult[i]->GetMeanError());
      printf("Event fraction: %3.2f%%\n",eventFraction[i]*100);
    }
  c = draw1D(hTofMult[0],"Raw TofMult distributin in VPDMB events;TofMult",kTRUE,kFALSE);
  for(int i=1; i<nMultBins-1; i++)
    {
      double max = hTofMult[0]->GetBinContent(high_tofMult[i]+1);
      TLine *line = GetLine(high_tofMult[i]+1,0,high_tofMult[i]+1,max);
      line->Draw();
    }
  TPaveText *t1 = GetPaveText(0.5,0.7,0.65,0.85);
  for(int i=0; i<nMultBins; i++)
    {
      t1->AddText(Form("<TofMult> = %3.2f, frac = %3.1f%%",meanTofMult[i],eventFraction[i]*100));
    }
  t1->SetTextAlign(11);
  t1->Draw();
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sTofMult_Raw.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sTofMult_Raw.png",run_type,run_cfg_name.Data()));
    }

  // Correction factor
  TList *list = new TList;
  TString legName[2] = {"VPDMB","BBCMB"};
  TFile *fCorr = TFile::Open("Rootfiles/VPDMB_Corr.root","read");

  const char *mbName[2] = {"VPDMB","BBCMB"};
  TH1F *hTofMultMB[2];
  for(int i=0; i<2; i++)
    {
      TFile *fMB = TFile::Open(Form("./output/Run13.pp500.jpsi.%s.root",mbName[i]),"read");
      THnSparseF *hn = (THnSparseF*)fMB->Get("mhEventMult_qa_mb");
      hn->GetAxis(0)->SetRange(1,350);
      if(cut_ranking)   hn->GetAxis(4)->SetRange(3,3);
      hTofMultMB[i] = (TH1F*)hn->Projection(2);
      hTofMultMB[i]->SetName(Form("hTofMult_%s",mbName[i]));
      hTofMultMB[i]->Sumw2();
      hTofMultMB[i]->Scale(1./hTofMultMB[i]->GetBinContent(1));
      hTofMultMB[i]->SetLineColor(2);
      list->Add(hTofMultMB[i]);
    }
  c = drawHistos(list,"VPD_vs_BBC","TofMult distribution in different triggers;TofMult;Probability",kTRUE,0,30,kFALSE,1e-6,0.05,kTRUE,kTRUE,legName,kTRUE,"Trigger",0.55,0.75,0.6,0.8,kTRUE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sTofMult_VPD_vs_BBC.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sTofMult_VPD_vs_BBC.png",run_type,run_cfg_name.Data()));
    }

  TH1F *hRatio = (TH1F*)hTofMultMB[0]->Clone("VPDMB_to_BBCMB");
  hRatio->Divide(hTofMultMB[1]);
  TF1 *func;
  if(!sys) func = new TF1("func_VPDMB_to_BBCMB","pol3",10,30);
  else func = new TF1("func_VPDMB_to_BBCMB","pol0",10,30);
  hRatio->Fit(func,"IR0");
  hRatio->SetMarkerStyle(21);
  hRatio->SetMarkerColor(4);
  c = draw1D(hRatio,"Ratio = VPDMB/BBCMB;TofMult;Ratio");
  func->Draw("sames");
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sFit_TofMult_VPD_to_BBC.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sFit_TofMult_VPD_to_BBC.png",run_type,run_cfg_name.Data()));
    }

  TH1F *hBBCcorr = (TH1F*)fCorr->Get("BBC_efficiency");
  draw1D(hBBCcorr,"Efficiency of BBCMB trigger estimated using PYTHIA;TofMult;Efficiency");
  TH1F *hVPDcorr = (TH1F*)hBBCcorr->Clone("VPD_efficiency");
  for(int bin=1; bin<=hVPDcorr->GetNbinsX(); bin++)
    {
      double scale = hRatio->GetBinContent(bin);
      if(bin>10)
	scale = func->Eval(hVPDcorr->GetBinCenter(bin));
      hVPDcorr->SetBinContent(bin,scale*hVPDcorr->GetBinContent(bin));
      hVPDcorr->SetBinError(bin,scale*hVPDcorr->GetBinError(bin));
    }
  hVPDcorr->GetYaxis()->SetRangeUser(0.7,1.8);
  //TF1 *func2 = new TF1("func_BBCMB_Corr","pol3",20,30);
  //hVPDcorr->Fit(func2,"IR0");
  c = draw1D(hVPDcorr,"Relative efficiency of VPDMB trigger;TofMult;Efficiency");
  //func2->Draw("sames");
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sVPD_eff_vs_TofMult.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sVPD_eff_vs_TofMult.png",run_type,run_cfg_name.Data()));
    }

  cout << "=== Corrected ===" << endl;
  TH1F *hTofMultCorr[nMultBins];
  double meanTofMultCorr[nMultBins];
  double eventFractionCorr[nMultBins];
  for(int i=0; i<nMultBins; i++)
    {
      hTofMultCorr[i] = (TH1F*)hTofMult[i]->Clone(Form("Corrected_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]));
      for(int bin=1; bin<=hTofMultCorr[i]->GetNbinsX(); bin++)
	{
	  double scale = hVPDcorr->GetBinContent(hVPDcorr->FindBin(hTofMultCorr[i]->GetBinCenter(bin)));
	  if(scale==0) continue;
	  hTofMultCorr[i]->SetBinContent(bin,hTofMultCorr[i]->GetBinContent(bin)/scale);
	  hTofMultCorr[i]->SetBinError(bin,hTofMultCorr[i]->GetBinError(bin)/scale);
	}
      meanTofMultCorr[i] = hTofMultCorr[i]->GetMean();
      eventFractionCorr[i] = hTofMultCorr[i]->Integral(1,30)/hTofMultCorr[0]->Integral(1,30);
      printf("TofMult: [%d,%d]\n",low_tofMult[i],high_tofMult[i]);
      printf("<TofMult> = %3.2f +/- %3.2e\n",meanTofMultCorr[i],hTofMultCorr[i]->GetMeanError());
      printf("Event fraction: %3.2f%%\n",eventFractionCorr[i]*100);
    }
  c = draw1D(hTofMultCorr[0],"Corrected TofMult distributin in VPDMB events;TofMult",kTRUE,kFALSE);
  for(int i=1; i<nMultBins-1; i++)
    {
      double max = hTofMultCorr[0]->GetBinContent(high_tofMult[i]+1);
      TLine *line = GetLine(high_tofMult[i]+1,0,high_tofMult[i]+1,max);
      line->Draw();
    }
  TPaveText *t1 = GetPaveText(0.5,0.7,0.65,0.85);
  for(int i=0; i<nMultBins; i++)
    {
      t1->AddText(Form("<TofMult> = %3.2f, frac = %3.1f%%",meanTofMultCorr[i],eventFractionCorr[i]*100));
    }
  t1->SetTextAlign(11);
  t1->Draw();

  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sTofMult_Corr.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sTofMult_Corr.png",run_type,run_cfg_name.Data()));
    }

  if(saveRoot)
    {
      TFile *fout;
      if(!cut_ranking && !cut_vtxdz) fout = TFile::Open("Rootfiles/Run13.pp500.VPDMB.corr.root","recreate");
      else             fout = TFile::Open(Form("Rootfiles/Run13.pp500.VPDMB.corr.Rank%d.VtxDz%d.root",cut_ranking,cut_vtxdz),"recreate");
      for(int i=0; i<nMultBins; i++)
	{
	  hTofMult[i]->Write();
	}
      hTofMultMB[0]->Write();
      hTofMultMB[1]->Write();
      hRatio->Write();
      hBBCcorr->Write();
      hVPDcorr->Write();
      for(int i=0; i<nMultBins; i++)
	{
	  hTofMultCorr[i]->Write();
	}
    }
  return;
}

//================================================
void ana(bool save = 1, bool saveRoot = 1)
{
  TString file_name;
  TString corr_name;

  if(jpsi_pt_cut_low==0 && jpsi_pt_cut_high==10 && !cut_ranking && !cut_vtxdz)
    file_name = "Pico.Run13.pp500.jpsi.yield.default.root";
  else
    file_name = Form("Pico.Run13.pp500.jpsi.yield.Rank%d.VtxDz%d.pt%1.0f-%1.0f.root",cut_ranking,cut_vtxdz,jpsi_pt_cut_low,jpsi_pt_cut_high);

  if(!cut_ranking && !cut_vtxdz) 
    corr_name = "Run13.pp500.VPDMB.corr.root";
  else            
    corr_name =Form("Run13.pp500.VPDMB.corr.Rank%d.VtxDz%d.root",cut_ranking,cut_vtxdz);

  TString out_name = file_name;
  out_name.ReplaceAll("yield","EvtAct");
  TFile *fin = TFile::Open(Form("Rootfiles/%s",file_name.Data()),"read");
  TFile *fMult = TFile::Open(Form("Rootfiles/%s",corr_name.Data()),"read");

  cout << out_name.Data() << endl;

  // Jpsi yield
  TH1F *hJpsiYield = (TH1F*)fin->Get("hJpsiYield");
  double nSignal[nMultBins];
  double sigErr[nMultBins];
  for(int i=0; i<nMultBins; i++)
    {
      nSignal[i] = hJpsiYield->GetBinContent(i+1);
      sigErr[i] = hJpsiYield->GetBinError(i+1);
    }

  // event activity
  TH1F *hTofMultCorr[nMultBins];
  double meanTofMult[nMultBins];
  double eventFraction[nMultBins];
  for(int i=0; i<nMultBins; i++)
    {
      hTofMultCorr[i] = (TH1F*)fMult->Get(Form("Corrected_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]));
      meanTofMult[i] = hTofMultCorr[i]->GetMean();
      eventFraction[i] = hTofMultCorr[i]->Integral(1,30)/hTofMultCorr[0]->Integral(1,30);
      printf("TofMult: [%d,%d]\n",low_tofMult[i],high_tofMult[i]);
      printf("<TofMult> = %3.2f +/- %3.2e\n",meanTofMult[i],hTofMultCorr[i]->GetMeanError());
      printf("Event fraction: %3.2f%%\n",eventFraction[i]*100);
    }
  
  TH1F *hYield[2];
  for(int i=0; i<2; i++)
    {
      hYield[i] = new TH1F(Form("hYield_%d",i),"J/psi yield",nMultBins-1,0,nMultBins-1);
    }
  for(int bin=1; bin<=nMultBins-1; bin++)
    {
      hYield[0]->SetBinContent(bin,nSignal[0]);
      hYield[0]->SetBinError(bin,sigErr[0]);
      hYield[1]->SetBinContent(bin,nSignal[bin]);
      hYield[1]->SetBinError(bin,sigErr[bin]);
      printf("Bin %d: %3.2f +/- %3.2f\n",bin,nSignal[bin]/nSignal[0],
	     nSignal[bin]/nSignal[0]*TMath::Sqrt(sigErr[0]*sigErr[0]/nSignal[0]/nSignal[0]+sigErr[bin]*sigErr[bin]/nSignal[bin]/nSignal[bin]));
    }

  TCanvas *c = new TCanvas("Yield_vs_activity","Yield_vs_activity",800,600);
  SetPadMargin(gPad,0.15,0.15);
  TH1F *h = new TH1F("histogram",";#frac{TofMult}{<TofMult>};#frac{N_{J/#psi}}{<N_{J/#psi}>}",7,0,4);
  h->GetYaxis()->SetRangeUser(0,6);
  ScaleHistoTitle(h,0.045,1.2,0.04,0.045,1.1,0.04,62);
  h->DrawCopy();

  TGraphAsymmErrors *ratio = new TGraphAsymmErrors();
  ratio->Divide(hYield[1],hYield[0],"cl=0.683 b(1,1) mode");
  int mode = 1;
  for(int i=0; i<ratio->GetN(); i++)
    {
      double x,y;
      ratio->GetPoint(i,x,y);

      printf("check error for bin %d: %3.2f +/- %3.2e/%3.2e\n",i+1,y,ratio->GetErrorYlow(i),ratio->GetErrorYhigh(i));

      double newx = meanTofMult[i+1]/meanTofMult[0];
      double exl = 0;
      double exh = 0;
      double newy = y/eventFraction[i+1];
      double eyl, eyh;
      if(mode==0)
	{
	  eyl = ratio->GetErrorYlow(i)/eventFraction[i+1];
	  eyh = ratio->GetErrorYhigh(i)/eventFraction[i+1];
	}
      else
	{
	  double error =  nSignal[i+1]/nSignal[0]*TMath::Sqrt(sigErr[0]*sigErr[0]/nSignal[0]/nSignal[0]+sigErr[i+1]*sigErr[i+1]/nSignal[i+1]/nSignal[i+1]);
	  eyl = error/eventFraction[i+1];
	  eyh = error/eventFraction[i+1];
	}
      ratio->SetPoint(i,newx,newy);
      ratio->SetPointError(i,exl,exh,eyl,eyh);
    }
  ratio->SetTitle();
  ratio->SetMarkerStyle(21);
  ratio->SetMarkerColor(2);
  ratio->SetLineColor(2);
  ratio->GetHistogram()->GetXaxis()->SetRangeUser(0,10);
  ratio->Draw("sames PE");
  TGraphAsymmErrors *sys = (TGraphAsymmErrors*)ratio->Clone("gSys");
  for(int i=0; i<sys->GetN(); i++)
    {
      double x,y;
      sys->GetPoint(i,x,y);
      double exl = 0.2;
      double exh = 0.2;
      double eyl = y * TMath::Sqrt(10*10+5*5+1.5*1.5)/100;
      double eyh = eyl;
      sys->SetPoint(i,x,y);
      sys->SetPointError(i,exl,exh,eyl,eyh);
    }
  sys->SetFillStyle(0);
  sys->SetLineColor(ratio->GetLineColor());
  sys->Draw("sameE5");
  TLine *line = GetLine(0,0,4,4,1);
  line->Draw();
 if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct_Run13.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct_Run13.png",run_type,run_cfg_name.Data()));
    }
  
  // Compare to Qian
  TCanvas *c = new TCanvas("Yield_vs_activity_3","Yield_vs_activity_3",800,600);
  SetPadMargin(gPad,0.15,0.15);
  h->DrawCopy();
  ratio->Draw("sames PE");
  sys->Draw("sameE5");
  TFile *fHT = TFile::Open("Rootfiles/Spectrum_in_bin_Rank.root");
  TGraphErrors *gr = (TGraphErrors*)fHT->Get("event_Act");
  gr->Draw("sames P");

  TGraphErrors *grSys = (TGraphErrors*)gr->Clone("grSys");
  for(int i=0; i<grSys->GetN(); i++)
    {
      double x,y;
      grSys->GetPoint(i,x,y);
      double exl = 0.2;
      double exh = 0.2;
      double eyl = y * TMath::Sqrt(1.3*1.3*2)/100;
      double eyh = eyl;
      grSys->SetPoint(i,x,y);
      grSys->SetPointError(i,exl,eyl);
    }
  grSys->SetFillStyle(0);
  grSys->SetLineColor(gr->GetLineColor());
  grSys->Draw("sameE5");


  leg = new TLegend(0.2,0.65,0.36,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(ratio,"MTD: p_{T,J/#psi} > 0 GeV/c","PL");
  leg->AddEntry(gr,"BEMC: p_{T,J/#psi} > 4 GeV/c","PL");
  leg->Draw();
  line->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct.png",run_type,run_cfg_name.Data()));
    }

  TCanvas *c = new TCanvas("Yield_vs_activity_2","Yield_vs_activity_2",800,600);
  SetPadMargin(gPad,0.15,0.15);
  h->SetMaximum(10);
  h->DrawCopy();
  ratio->Draw("sames PE");
  sys->Draw("sameE5");
  TGraphErrors *graph[3];
  graph[0] = (TGraphErrors*)fHT->Get("event_Act");
  graph[1] = (TGraphErrors*)fHT->Get("evAct_lowpT");
  graph[2] = (TGraphErrors*)fHT->Get("evAct_higpT");
  for(int i=1; i<3; i++)
    {
      graph[i]->Draw("sames P");
    }
  leg = new TLegend(0.2,0.65,0.36,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(ratio,"MTD: p_{T,J/#psi} > 0 GeV/c","PL");
  leg->AddEntry(graph[1],"BEMC: 4 < p_{T,J/#psi} < 8 GeV/c","PL");
  leg->AddEntry(graph[2],"BEMC: p_{T,J/#psi} > 8 GeV/c","PL");
  leg->Draw();
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct_2.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct_2.png",run_type,run_cfg_name.Data()));
    }

  if(saveRoot)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s",out_name.Data()),"recreate");
      ratio->Write("Run13_Jpsi_vs_evtAct",TObject::kOverwrite);
      sys->Write("Run13_Jpsi_vs_evtAct_sys",TObject::kOverwrite);
      gr->Write("Run11_Jpsi_vs_evtAct",TObject::kOverwrite);
      grSys->Write("Run11_Jpsi_vs_evtAct_sys",TObject::kOverwrite);
    }
}

//================================================
void compare(bool save = 1)
{
  TFile *fin[2];
  fin[0] = TFile::Open("Rootfiles/Pico.Run13.pp500.jpsi.EvtAct.default.root","read");
  fin[1] = TFile::Open("Rootfiles/Pico.Run13.pp500.jpsi.EvtAct.Rank1.VtxDz1.pt0-10.root","read");

  TGraphAsymmErrors *gr[2];
  for(int i=0; i<2; i++)
    {
      gr[i] = (TGraphAsymmErrors*)fin[i]->Get("Run13_Jpsi_vs_evtAct");
      gr[i]->SetName(Form("%s_%d",gr[i]->GetName(),i));
      gr[i]->SetMarkerColor(color[i]);
      gr[i]->SetLineColor(color[i]);
    }

  TCanvas *c = new TCanvas("Yield_vs_activity","Yield_vs_activity",800,600);
  SetPadMargin(gPad,0.15,0.15);
  TH1F *h = new TH1F("histogram",";#frac{TofMult}{<TofMult>};#frac{N_{J/#psi}}{<N_{J/#psi}>}",7,0,6);
  h->GetYaxis()->SetRangeUser(0,20);
  ScaleHistoTitle(h,0.045,1.2,0.04,0.045,1.1,0.04,62);
  h->Draw();
  gr[0]->Draw("sames PE");
  gr[1]->Draw("sames PE");


  leg = new TLegend(0.2,0.65,0.36,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("MTD: p_{T,J/#psi} > 0 GeV/c");
  leg->AddEntry(gr[0],"No cut on ranking","PL");
  leg->AddEntry(gr[1],"Ranking > 0","PL");
  leg->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct_CompareRanking_3.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Activity/%sJpsi_vs_EvtAct_CompareRanking_3.png",run_type,run_cfg_name.Data()));
    }

}
