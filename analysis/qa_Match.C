TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "FitDz.";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void qa_Match()
{						
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.15); 

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  //Track();
  //DeltaZ();
  //CompareFitDz();
  //DeltaY();
  //qualityCuts();
  //yzDistribution();
  //eLoss();
  magneticField();
}

//================================================
void CompareFitDz(const Int_t savePlot = 1)
{
  TString filename = f->GetName();
  filename.ReplaceAll("output","Rootfiles");
  filename.ReplaceAll(Form(".%sroot",run_config),".FitDz.root");
  TFile *fin = TFile::Open(filename.Data(),"read");
  const char *type_name[2] = {"2Gaus","3Gaus"};

  TList *list = new TList;
  TString legName[2] = {"Two-Gaussian fit","Three-Gaussian fit"};

  list->Clear();
  TH1F *hChi2[2];
  for(int i=0; i<2; i++)
    {
      hChi2[i] = (TH1F*)fin->Get(Form("hChi_%s",type_name[i]));
      list->Add(hChi2[i]);
    }
  c = drawHistos(list,"chi2","Chi2/NDF of #Deltaz fit;p_{T} [GeV/c];Chi2/ndf",kFALSE,2.0,3.8,kTRUE,0.1,300,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_FitDzChi2_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_FitDzChi2_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TH1F *hMean[2];
  for(int i=0; i<2; i++)
    {
      hMean[i] = (TH1F*)fin->Get(Form("hMean_%s",type_name[i]));
      list->Add(hMean[i]);
    }
  c = drawHistos(list,"mean","Mean of fitted #Deltaz distribution",kFALSE,2.0,3.8,kFALSE,0.1,300,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_FitDzMean_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_FitDzMean_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  TH1F *hSigma[2];
  for(int i=0; i<2; i++)
    {
      hSigma[i] = (TH1F*)fin->Get(Form("hSigma_%s",type_name[i]));
      hSigma[i]->SetMarkerStyle(20);
      hSigma[i]->SetMarkerColor(color[i]);
      hSigma[i]->SetLineColor(color[i]);
    }
  hSigma[0]->SetMaximum(16);
  c = draw1D(hSigma[0],"Sigma of fitted #Deltaz distribution");
  hSigma[1]->Draw("sames");
  TF1 *fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)",1,20);
  fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.444217);
  fResDzVsPt->SetLineColor(4);
  fResDzVsPt->Draw("sames");
  TLegend *leg = new TLegend(0.55,0.65,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hSigma[0],legName[0].Data(),"PL");
  leg->AddEntry(hSigma[1],legName[1].Data(),"PL");
  leg->AddEntry(fResDzVsPt,"Run13 embedding","L");
  leg->Draw();
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_FitDzSigma_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_FitDzSigma_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TH1F *hPurity[2];
  for(int i=0; i<2; i++)
    {
      hPurity[i] = (TH1F*)fin->Get(Form("hPurity_%s",type_name[i]));
      list->Add(hPurity[i]);
    }
  c = drawHistos(list,"hPurity","Purity of selected muon sample by #Deltaz cut",kFALSE,2.0,3.8,kFALSE,0.1,300,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.2,0.4,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_DzPurity_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_DzPurity_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TH1F *hEff[2];
  for(int i=0; i<2; i++)
    {
      hEff[i] = (TH1F*)fin->Get(Form("hEff_%s",type_name[i]));
      list->Add(hEff[i]);
    }
  c = drawHistos(list,"hEff","Efficiency of #Deltaz cut",kFALSE,2.0,3.8,kFALSE,0.1,300,kFALSE,kTRUE,legName,kTRUE,"",0.2,0.4,0.2,0.4,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_DzCutEff_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sCompare_DzCutEff_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }
}

//================================================
void DeltaZ(const Int_t savePlot = 0, const Int_t saveHisto = 0)
{
  TH2F *hDzVsPt = (TH2F*)f->Get(Form("mhDzVsPt_%s",trigName[kTrigType]));
  c = draw2D(hDzVsPt,Form("%s: #Deltaz of matched %s track-hit pairs%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sDeltaZ_vs_pt_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sDeltaZ_vs_pt_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  // Fitting
  const int type = 0;
  const char *type_name[2] = {"2Gaus","3Gaus"};

  TH1F *hMthDz = (TH1F*)hDzVsPt->ProjectionY(Form("hTrkDz_%s_proj",trigName[kTrigType]));
  hMthDz->SetTitle(Form("%s: #Deltaz of matched track-hit pairs;#Deltaz (cm)",trigName[kTrigType]));
  Double_t range = 80;
  TF1 *func;
  if(type==0)
    {
      func = new TF1("func","gaus(0)+gaus(3)",-1*range,range);
      func->SetParameter(1,0);
      func->SetParameter(4,0);
      func->SetParameter(2,10);
      func->SetParameter(5,50);
    }
  else if(type==1)
    {
      func = new TF1("func","gaus(0)+gaus(3)+gaus(6)",-1*range,range);
      func->SetParameter(1,0);
      func->SetParameter(4,0);
      func->SetParameter(7,0);
      func->SetParameter(2,10);
      func->SetParameter(5,20);
      func->SetParameter(8,100);
    }
  c = FitDeltaZ(hMthDz,func,range,20.,type);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_%s_All_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_%s_All_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  //return;
  // pt dependence
  const int nTrkPtBin = 26;
  const double trkPtBins[nTrkPtBin+1] = {1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0};
  double trkPtbins[nTrkPtBin+2];
  trkPtbins[0] = 0;
  for(int i=1; i<nTrkPtBin+2; i++) trkPtbins[i] = trkPtBins[i-1];
  TH1F *hMean = new TH1F(Form("hMean_%s",type_name[type]),Form("Mean of fitted #Deltaz distribution;p_{T} [GeV/c];mean [cm]"),nTrkPtBin,trkPtbins);
  TH1F *hSigma = new TH1F(Form("hSigma_%s",type_name[type]),Form("Sigma of fitted #Deltaz distribution;p_{T} [GeV/c];#sigma [cm]"),nTrkPtBin,trkPtbins);
  TH1F *hPurity = new TH1F(Form("hPurity_%s",type_name[type]),Form("Purity of selected muon sample by #Deltaz cut;p_{T} [GeV/c];purity"),nTrkPtBin,trkPtbins);
  TH1F *hEff = new TH1F(Form("hEff_%s",type_name[type]),Form("Efficiency of #Deltaz cut;p_{T} [GeV/c];Efficiency"),nTrkPtBin,trkPtbins);
  TH1F *hChi2 = new TH1F(Form("hChi_%s",type_name[type]),Form("Chi2/NDF of #Deltaz fit;p_{T} [GeV/c];Chi2/ndf"),nTrkPtBin,trkPtbins);
  TF1 *fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)");
  fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.444217);
  const int nCanvas = nTrkPtBin/6+1;
  TCanvas *canvas[nCanvas];
  for(int i=0; i<nTrkPtBin; i++)
    {
      if(i%6==0) 
	{
	  canvas[i/6] = new TCanvas(Form("Fit_dz_%d",i/6),Form("Fit_dz_%d",i/6),1100,650);
	  canvas[i/6]->Divide(3,2);
	}
      int low_bin  = hDzVsPt->GetXaxis()->FindFixBin(trkPtBins[i]+1e-4);
      int high_bin = hDzVsPt->GetXaxis()->FindFixBin(trkPtBins[i+1]-1e-4);
      
      TH1F *htmp = (TH1F*)hDzVsPt->ProjectionY(Form("hTrkDz_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],trigName[kTrigType]),low_bin,high_bin);
      htmp->SetTitle(Form("%s: #Deltaz of matched track-hit pairs (%1.1f < p_{T} < %1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],trkPtBins[i],trkPtBins[i+1]));

      TF1 *func;
      if(type==0)
	{
	  func = new TF1(Form("func_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],type_name[type]),"gaus(0)+gaus(3)",-1*range,range);
	  func->SetParameter(1,0);
	  func->SetParameter(4,0);
	  func->SetParameter(2,10);
	  func->SetParameter(5,50);
	}
      if(type==1)
	{
	  func = new TF1(Form("func_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],type_name[type]),"gaus(0)+gaus(3)+gaus(6)",-1*range,range);
	  func->SetParameter(1,0);
	  func->SetParameter(4,0);
	  func->SetParameter(7,0);
	  func->SetParameter(2,10);
	  func->SetParameter(5,20);
	  func->SetParameter(8,60);
	  if(i<6)
	    {
	      func->SetParameter(2,8);
	      func->SetParameter(5,25);
	      func->SetParameter(8,40);
	    }
	  if(i==1)
	    {
	      func->SetParameter(2,10);
	      func->SetParameter(5,25);
	      func->SetParameter(8,40);
	    }
	  if(i==0)
	    {
	      func->SetParLimits(2,7,9);
	      func->SetParameter(5,20);
	      func->SetParameter(8,54);
	    }
	}
      htmp->Fit(func,"IR0");
      hChi2->SetBinContent(i+2, func->GetChisquare()/func->GetNDF());
      hChi2->SetBinError(i+2, 0);
      htmp->GetYaxis()->SetNdivisions(505);
      canvas[i/6]->cd(i%6+1);
      htmp->SetTitle(";#Deltaz [cm]");
      htmp->GetXaxis()->SetRangeUser(-100,100);
      htmp->SetMarkerStyle(20);
      htmp->Draw("HIST");
      func->SetLineColor(2);
      func->Draw("sames");

      TF1 *func1, *func2, *func3;
      func1 = new TF1("func1","gaus",-1*range,range);
      func1->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
      if(type==0)
	{
	  func2 = new TF1("func2","gaus",-1*range,range);
	  func2->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
	  func2->SetLineColor(4);
	  func2->Draw("sames");
	}
      else if(type==1)
	{
	  func2 = new TF1("func2","gaus(0)+gaus(3)",-1*range,range);
	  func2->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5),func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
	  func2->SetLineColor(kGreen+2);
	  func2->Draw("sames");

	  func3 = new TF1("func3","gaus",-1*range,range);
	  func3->SetParameters(func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
	  func3->SetLineColor(4);
	  func3->Draw("sames");
	}


      TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",trkPtBins[i],trkPtBins[i+1]),0.06);
      t->Draw();
      if(savePlot && (i%6==5 || i==nTrkPtBin-1))
	{
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_%s_Bin%d_%s.png",run_type,run_cfg_name.Data(),type_name[type],i/6,trigName[kTrigType]));
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_%s_Bin%d_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],i/6,trigName[kTrigType]));
	}

      hMean->SetBinContent(i+2,func->GetParameter(1));
      hMean->SetBinError(i+2, func->GetParError(1));
      hSigma->SetBinContent(i+2,func->GetParameter(2));
      hSigma->SetBinError(i+2, func->GetParError(2));

      // pt dependent dz cut
      double pt = hPurity->GetBinCenter(i+1);
      double dz_sigma = fResDzVsPt->Eval(pt);
      double cut = 999.;
      if(pt<3) cut = 2 * dz_sigma;
      else     cut = 2.5 * dz_sigma;
      
      double all = htmp->Integral(htmp->FindFixBin(-1*cut),htmp->FindFixBin(cut));
      double all_err = TMath::Sqrt(all);
      double bkg = (func->Integral(-1*cut,cut)-func1->Integral(-1*cut,cut)) * 1./htmp->GetBinWidth(1);
      double bkg_err = TMath::Sqrt(bkg);
      double signal = all - bkg;
      double signal_err = TMath::Sqrt(all_err*all_err+bkg_err*bkg_err);
      double purity = 1 - bkg/all;
      double purity_err = TMath::Sqrt(bkg_err*bkg_err/all/all+bkg*bkg/all/all/all/all*all_err*all_err);
      hPurity->SetBinContent(i+2,purity);
      hPurity->SetBinError(i+2,purity_err);
      
      double out = 
	func1->Integral(func1->GetParameter(1)-5*func1->GetParameter(2),-1*cut)* 1./htmp->GetBinWidth(1) +
	func1->Integral(cut, func1->GetParameter(1)+5*func1->GetParameter(2))* 1./htmp->GetBinWidth(1);
      double eff = signal/(signal+out);
      hEff->SetBinContent(i+2,eff);
      hEff->SetBinError(i+2,1e-10);

    }
  hChi2->SetMarkerStyle(21);
  hChi2->SetMaximum(1.2*hChi2->GetMaximum());
  c = draw1D(hChi2);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_Chi2_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_Chi2_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hMean->SetMarkerStyle(21);
  hMean->GetYaxis()->SetRangeUser(-0.5,0.5);
  c = draw1D(hMean);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_Mean_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_Mean_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hSigma->SetMarkerStyle(21);
  hSigma->GetYaxis()->SetRangeUser(0,15);
  c = draw1D(hSigma);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_Sigma_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sFitDz_Sigma_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hPurity->SetMarkerStyle(21);
  hPurity->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hPurity);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sDzCut_Purity_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sDzCut_Purity_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hEff->SetMarkerStyle(21);
  hEff->GetYaxis()->SetRangeUser(0,1.1);
  c = draw1D(hEff);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sDzCut_Efficiency_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%sDzCut_Efficiency_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }

  if(saveHisto)
    {
      TString filename = f->GetName();
      filename.ReplaceAll("output","Rootfiles");
      filename.ReplaceAll(Form(".%sroot",run_config),".FitDz.root");
      cout << filename.Data() << endl;
      TFile *fout = TFile::Open(filename.Data(),"update");
      hChi2->Write("",TObject::kOverwrite);
      hMean->Write("",TObject::kOverwrite);
      hSigma->Write("",TObject::kOverwrite);
      hPurity->Write("",TObject::kOverwrite);
      hEff->Write("",TObject::kOverwrite);
      fout->Close();
    }
}

//================================================
TCanvas *FitDeltaZ(TH1 *histo, TF1 *func, const Double_t range1 = 50., Double_t range2 = 20., int type = 0)
{
  histo->Fit(func,"IR0");
  histo->GetYaxis()->SetNdivisions(505);
  c = draw1D(histo,"",kFALSE,kFALSE);
  func->SetLineColor(2);
  func->Draw("sames");

  TF1 *func1, *func2;
  if(type==0)
    {
      func1 = new TF1("func1","gaus",-1*range1,range1);
      func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
      func1->SetLineColor(4);
      func1->Draw("sames");
    }
  else if(type==1)
    {
      func1 = new TF1("func1","gaus(0)+gaus(3)",-1*range1,range1);
      func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5),func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
      func1->SetLineColor(kGreen+2);
      func1->Draw("sames");

      func2 = new TF1("func2","gaus",-1*range1,range1);
      func2->SetParameters(func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
      func2->SetLineColor(4);
      func2->Draw("sames");
    }
  TPaveText *t1 = GetPaveText(0.2,0.3,0.65,0.78,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f] cm",range1,range1));
  t1->AddText(Form("S/B ~ %1.3f:1",((func->Integral(-1*range1,range1))-(func1->Integral(-1*range1,range1)))/(func1->Integral(-1*range1,range1))));
  t1->Draw();
  t1 = GetPaveText(0.2,0.3,0.47,0.6,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f] cm",range2,range2));
  t1->AddText(Form("S/B ~ %1.3f:1",((func->Integral(-1*range2,range2))-(func1->Integral(-1*range2,range2)))/(func1->Integral(-1*range2,range2))));
  t1->Draw();
  return c;
}

//================================================
void DeltaY(const Int_t save = 1)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("hTrkDzDy_%s",trigName[kTrigType]));
  TH2F *hTrkDyVsPt = (TH2F*)hn->Projection(2,0);
  c = draw2D(hTrkDyVsPt,Form("Au+Au %s: #Deltay of matched %s track-hit pairs%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.DeltaY_vs_pt_%s.png",run_config,trigName[kTrigType]));

  c = new TCanvas("hDy_TrkPtBin","hDy_TrkPtBin",1200,650);
  c->Divide(2,2);
  Double_t pt_cuts[5] = {1,2,4,10,20};
  for(Int_t i=0; i<4; i++)
    {
      hTrkDyVsPt->GetXaxis()->SetRangeUser(pt_cuts[i]+0.1,pt_cuts[i+1]-0.1);
      TH1F *hDy = (TH1F*)hTrkDyVsPt->ProjectionY(Form("hTrkDy_pt%1.0f-%1.0f_%s",pt_cuts[i],pt_cuts[i+1]));
      c->cd(i+1);
      hDy->Draw();
      TPaveText *t1 = GetTitleText(Form("Au+Au %s: #Deltay of matched track-hit pairs",trigName[kTrigType]),0.06);
      t1->Draw();
      t1 = GetPaveText(0.15,0.35,0.7,0.75,0.06);
      t1->AddText(Form("%1.0f < p_{T} < %1.0f GeV/c",pt_cuts[i],pt_cuts[i+1]));
      t1->SetTextColor(2);
      t1->Draw();
    }
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.DeltaY_InPtBin_%s.png",run_config,trigName[kTrigType]));
}



//================================================
void qualityCuts(const Int_t save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("hTrkHitDz_qa_%s",trigName[kTrigType]));
  Double_t pt_cut = 2;
  hn->GetAxis(0)->SetRangeUser(pt_cut+0.1,20);

  TH1F *hdz[3];
  Double_t nsigma_low[3] = {-5,-1,0};
  Double_t nsigma_hi[3]  = {5, 3, 3};
  Double_t range = 50, range2 = 20;
  for(Int_t i=0; i<3; i++)
    {
      hn->GetAxis(4)->SetRangeUser(nsigma_low[i]+0.1,nsigma_hi[i]-0.1);
      hdz[i] = (TH1F*)hn->Projection(5);
      hdz[i]->SetName(Form("hdz_nsigma_cut%d",i));
      if(i==0) hdz[i]->SetTitle(Form("Au+Au %s: #Deltaz distribution w/o n#sigma_{#pi} cut (p_{T} > %1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],pt_cut));
      else     hdz[i]->SetTitle(Form("Au+Au %s: #Deltaz distribution w/ %1.0f<n#sigma_{#pi}<%1.0f (p_{T} > %1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],nsigma_low[i],nsigma_hi[i],pt_cut));

      TF1 *func = new TF1(Form("func_%d",i),"gaus(0)+gaus(3)",-1*range,range);
      if(i==0) func->SetParameters(1000,0,15,1000,0,60);
      if(i==1) func->SetParameters(1000,0,15,1000,0,60);
      if(i==2) func->SetParameters(1000,0,15,1000,0,60);
      c = FitDeltaZ(hdz[i],func,range,range2);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.FitDz_nsigmaCut%d_Pt%1.0fGeV_%s.png",run_config,i,pt_cut,trigName[kTrigType]));
    }
}




//================================================
void yzDistribution(const Int_t save = 0)
{
  const char *title[2] = {"z","y"};
  const char *trkname[2] = {"projected tracks","MTD hits"};
  const char *name[2] = {"track","hit"};

  // hit map
  TH2F *hMtdHitMap = (TH2F*)f->Get(Form("hMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMtdHitMap,Form("Au+Au %s: channel vs backleg of good MTD hits%s",trigName[kTrigType],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.MtdGoodHitMap_%s.png",run_config,trigName[kTrigType]));

  // matched hit map
  TH2F *hMthMtdHitMap = (TH2F*)f->Get(Form("hMthMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMthMtdHitMap,Form("Au+Au %s: channel vs backleg of matched MTD hits%s",trigName[kTrigType],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.MthMtdHitMap_%s.png",run_config,trigName[kTrigType]));
  
  THnSparseF *hYZ[2];
  hYZ[0] = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hYZ[1] = (THnSparseF*)f->Get(Form("hHitYZ_qa_%s",trigName[kTrigType]));
  TH2F *hYBL[2], *hZBL[2];
  TH1F *hYAll[2], *hZAll[2];
  for(Int_t i=0; i<2; i++)
    {
      hYBL[i] = (TH2F*)hYZ[i]->Projection(0,2);
      hYBL[i]->SetName(Form("%s_YBL",trkname[i]));
      hYBL[i]->GetYaxis()->SetRangeUser(-50,50);

      hZBL[i] = (TH2F*)hYZ[i]->Projection(1,2);
      hZBL[i]->SetName(Form("%s_ZBL",trkname[i]));

      c1 = draw2D(hYBL[i],Form("Au+Au %s: local y distribution of %s per module;backleg*5+module",trigName[kTrigType],trkname[i]));
      c2 = draw2D(hZBL[i],Form("Au+Au %s: global z distribution of %s per module;backleg*5+module",trigName[kTrigType],trkname[i]));
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_y_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));

      if(i==1)
	{
	  TH2F *h2tmp = (TH2F*)hZBL[i]->Clone(Form("%s_zoomin",hZBL[i]->GetName()));
	  h2tmp->GetXaxis()->SetRangeUser(0,25);
	  c2 = draw2D(h2tmp,Form("Au+Au %s: global z distribution of %s per module;backleg*5+module",trigName[kTrigType],trkname[i]),0.04,kFALSE);
	  if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_vs_Module_zoomin_%s.png",run_config,name[i],trigName[kTrigType]));
	}


      hYAll[i] = (TH1F*)hYBL[i]->ProjectionY(Form("%s_Y_All",trkname[i]));
      hZAll[i] = (TH1F*)hZBL[i]->ProjectionY(Form("%s_Z_All",trkname[i]));     
      hYAll[i]->Scale(1./hYAll[i]->Integral());
      hZAll[i]->Scale(1./hZAll[i]->Integral());

      c1 = draw1D(hYAll[i],Form("Au+Au %s: local y distribution of %s",trigName[kTrigType],trkname[i]),kFALSE,kFALSE);
      c2 = draw1D(hZAll[i],Form("Au+Au %s: global z distribution of %s",trigName[kTrigType],trkname[i]),kFALSE,kFALSE);
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_y_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_%s.png",run_config,name[i],trigName[kTrigType]));
    }

  TH2F *hTrkYZAll = (TH2F*)hYZ[0]->Projection(0,1);
  hTrkYZAll->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYZAll,Form("Au+Au %s: local y vs global z of projected tracks at MTD radius;global z (cm);local y (cm)",trigName[kTrigType]),0.04,kFALSE);
  gPad->SetRightMargin(0.13);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.track_y_vs_z_%s.png",run_config,trigName[kTrigType]));

  return;

  THnSparseF *hTrkYZ = (THnSparseF*)f->Get(Form("hTrkYZ_qa_%s",trigName[kTrigType]));
  const char *trkYZ_name[4] = {"start point","outer magnet","MTD (one module)","MTD (multiple modules)"};
  for(Int_t i=0; i<4; i++)
    {
      hTrkYZ->GetAxis(3)->SetRange(i+1,i+1);
      TH2F *h2 = (TH2F*)hTrkYZ->Projection(2,1);
      h2->SetName(Form("hTrkYZ_%s",trkYZ_name[i]));
      c = draw2D(h2,Form("Au+Au %s: #varphi vs global z of %s tracks at %s",trigName[kTrigType],trk_name[trk_index],trkYZ_name[i]),0.04,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.track_phi_vs_z_%s_%s.png",run_config,trkYZ_name[i],trigName[kTrigType]));

      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("hTrkY_%s",trkYZ_name[i]));
      c = draw1D(h1,Form("Au+Au %s: #varphi distribution of %s tracks at %s",trigName[kTrigType],trk_name[trk_index],trkYZ_name[i]),kFALSE,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.track_phi_%s_%s.png",run_config,trkYZ_name[i],trigName[kTrigType]));
    }


  TH1F *hHitZ[5];
  TH2F *hHitZvsBL[5];
  TList *list1 = new TList;
  TString legName1[5];
  for(Int_t i=0; i<5; i++)
    {
      hHitZvsBL[i] = new TH2F(Form("hHitZvsBL_Mod%d",i+1),Form("Au+Au %s: global z distribution of MTD hits in module %d;backleg;z (cm)",trigName[kTrigType],i+1),30,0,30,hZBL[1]->GetNbinsY(),hZBL[1]->GetYaxis()->GetXmin(), hZBL[1]->GetYaxis()->GetXmax());
      for(Int_t j=0; j<30; j++)
	{
	  TH1F *htmp = (TH1F*)hZBL[1]->ProjectionY(Form("hHitZ_BL%d_Mod%d",j+1,i+1),j*5+i+1,j*5+i+1);
	  if(j==0) hHitZ[i] = (TH1F*)htmp->Clone(Form("hHitZ_Mod%d_AllBL",i+1));
	  else     hHitZ[i]->Add(htmp);
	  for(Int_t ibin=1; ibin<=hHitZvsBL[i]->GetNbinsY(); ibin++)
	    {
	      hHitZvsBL[i]->SetBinContent(j+1,ibin,htmp->GetBinContent(ibin));
	      hHitZvsBL[i]->SetBinError(j+1,ibin,htmp->GetBinError(ibin));
	    }
	}
      hHitZ[i]->Sumw2();
      hHitZ[i]->Scale(1./hHitZ[i]->Integral());
      list1->Add(hHitZ[i]);
      legName1[i] = Form("Module %d",i+1);
      c = draw2D(hHitZvsBL[i],"",0.04,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.hit_z_vs_BL_Mod%d_%s.png",run_config,i+1,trigName[kTrigType]));
    }
  c = drawHistos(list1,list1->At(0)->GetName(),Form("Au+Au %s: z distribution of MTD hits in each module;z (cm)",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-6,0.018,kFALSE,kTRUE,legName1,kTRUE,"",0.6,0.8,0.68,0.88,kFALSE);
  for(Int_t i=0; i<6; i++)
    {
      TLine *line = GetLine(-217.5+i*87, 0, -217.5+i*87, 0.0125, 2, 2, 2);
      line->Draw();
    }
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.hit_z_per_Module_%s.png",run_config,trigName[kTrigType]));

      
}


//================================================
void Track(const Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");

  // track multiplicity
  TH1F *hNTrk       = (TH1F*)f->Get(Form("mhTrkN_%s",trigName[kTrigType]));
  TH1F *hMthMtdHitN = (TH1F*)f->Get(Form("mhMthTrkN_%s",trigName[kTrigType]));
  scaleHisto( hNTrk,    hStat->GetBinContent(9), 1);
  scaleHisto( hMthMtdHitN, hStat->GetBinContent(9), 1);

  TList *list = new TList;
  list->Add(hNTrk);
  list->Add(hMthMtdHitN);
  TString legName[2];
  legName[0] = Form("Good tracks <N> = %2.2f",hNTrk->GetMean());
  legName[1] = Form("Matched tracks <N> = %2.2f",hMthMtdHitN->GetMean());
  c = drawHistos(list,"Track_multiplicity",Form("%s: multiplicity of %s tracks%s;N_{trk};1/N_{evt} dN_{trk}",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kFALSE,0,100,kTRUE,1e-8,10,kTRUE,kTRUE,legName,kTRUE,run_type,0.5,0.7,0.6,0.8);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%s.NMthTrk.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%s.NMthTrk.pdf",run_type,run_cfg_name.Data()));
    }
  c = drawHistos(list,"Track_multiplicity_zoomin",Form("%s: multiplicity of %s tracks%s;N_{trk};1/N_{evt} dN_{trk}",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kTRUE,0,10,kTRUE,1e-8,10,kTRUE,kTRUE,legName,kTRUE,run_type,0.2,0.4,0.2,0.4);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%s.NMthTrk_zoomin.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Match/%s.NMthTrk_zoomin.pdf",run_type,run_cfg_name.Data()));
    }
}


//================================================
void eLoss(const Int_t save = 0, const int saveAN = 1)
{
  TCanvas *c = new TCanvas("energy_loss","energy_loss",800,600);
  TF1 *fEloss = new TF1("f2","[0]*exp(-pow([1]/x,[2]))",1.,20);
  fEloss->SetParameters(1.38147e+00,6.08655e-02,5.03337e-01);
  fEloss->SetMinimum(0);
  fEloss->SetTitle(";p_{track} (GeV/c);#DeltaE (GeV)");
  fEloss->Draw();
 
  TF1 *fEmc = new TF1("fEmc","pol0",0.,100.);
  fEmc->SetParameter(0,0.215);
  fEmc->SetLineColor(2);
  fEmc->Draw("sames");

  TF1 *fCoil = new TF1("fCoil","pol0",0.,100.);
  fCoil->SetParameter(0,0.176);
  fCoil->SetLineColor(4);
  fCoil->Draw("sames");

  TF1 *fMag = new TF1("fMag","[0]*exp(-pow([1]/x,[2]))-[3]",0.,100);
  fMag->SetParameters(1.38147e+00,6.08655e-02,5.03337e-01,fEmc->GetParameter(0),fCoil->GetParameter(0));
  fMag->SetLineColor(6);
  fMag->Draw("sames");

  TLegend *leg = new TLegend(0.6,0.3,0.8,0.6);
  leg->SetHeader("Energy loss");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(fEloss,"Sum","L");
  leg->AddEntry(fEmc,"EMC","L");
  leg->AddEntry(fCoil,"Coil","L");
  leg->AddEntry(fMag,"Magnet","L");
  leg->Draw();

  TPaveText *t1 = GetTitleText("Energy loss of tracks propagated to MTD",0.05);
  t1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/EnergyLoss.png"));
      c->SaveAs("/Users/admin/Dropbox/STAR/MTD/figures/EnergyLoss.pdf");
    }

  if(saveAN)
    c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch2_EnergyLoss.pdf"));
}

//================================================
void magneticField(const Int_t save = 0, const int saveAN = 1)
{
  f = TFile::Open("output/jpsi.test.histos.root");

  TH2F *hMag = (TH2F*)f->Get("hMagneticMap");
  hMag->RebinX(2);
  hMag->RebinY(2);
  hMag->SetZTitle("B (T)");
  hMag->SetTitleOffset(1.2,"Z");
  TCanvas *c = new TCanvas("MagneticMap","MagneticMap",880,800);
  TPaveText *t1 = GetTitleText(hMag->GetTitle());
  hMag->SetTitle("");
  hMag->GetYaxis()->SetTitleOffset(1.4);
  hMag->Draw("colz");
  t1->Draw();
  SetPadMargin(gPad,0.12,0.12,0.15,0.12);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/MagneticFieldMap.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/MagneticFieldMap.pdf",run_type));
    }
  if(saveAN)
    c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch2_MagneticField.pdf"));
}
