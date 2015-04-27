const char *run_config = "muon.";
const Bool_t iPico = 1;
const int year = 2013;
TString run_cfg_name;

TFile *f;

//================================================
void calib_dTof()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2); 

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

  makeCalib();
}


//================================================
void makeCalib(const int save = 0, const int saveHisto = 0)
{
  TH2F *hDtofVsMod = (TH2F*)f->Get(Form("mhDeltaTof_%s",trigName[kTrigType]));
  hDtofVsMod->GetYaxis()->SetRangeUser(-10,10);
  c = draw2D(hDtofVsMod,Form("%s: #Deltatof of tracks matched to MTD;p_{T} (GeV/c);#Deltatof (ns)",trigName[kTrigType]));
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/calib_dTof/Dtof_vs_Mod.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/calib_dTof/Dtof_vs_Mod.pdf",run_type));
    }

  TSpectrum *peakfinder = new TSpectrum();
  TH1F *hDtof[30][5];
  TF1  *func[30][5];
  TH1F *hMean = new TH1F("hDtof_Mean",Form("%s: peak position of #Deltatof distribution;module;#Deltatof_{peak} (ns)",run_type),150,0,150);
  for(int i=0; i<30; i++)
    {
      TCanvas *cTof = 0x0;
      for(int j=0; j<5; j++)
	{
	  int bl = i+1, mod = j+1;
	  hDtof[i][j] = (TH1F*)hDtofVsMod->ProjectionY(Form("hDtof_BL%d_Mod%d",bl,mod),i*5+j+1,i*5+j+1);
	  if(hDtof[i][j]->GetEntries()==0) continue;
	  if(j==0) 
	    {
	      cTof = new TCanvas(Form("hDtof_BL%d",bl),Form("hDtof_BL%d",bl),1100,750);
	      cTof->Divide(3,2);
	    }
	  cTof->cd(j+1);

	  TH1F *htmp = (TH1F*)hDtof[i][j]->Clone(Form("%s_clone",hDtof[i][j]->GetName()));
	  Int_t nfound = peakfinder->Search(htmp,1,"nodrawnew");
	  TList *functions = htmp->GetListOfFunctions();
	  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
	  Double_t *xpeaks = pm->GetX();

	  func[i][j] = new TF1(Form("Fit_BL%d_Mod%d",bl,mod),"gaus",xpeaks[0]-0.3,xpeaks[0]+0.2);
	  if((bl==6||bl==27) && mod==1) func[i][j]->SetRange(xpeaks[0]-0.3,xpeaks[0]+0.4);
	  func[i][j]->SetParameters(hDtof[i][j]->GetMaximum(),xpeaks[0],0.2);
	  hDtof[i][j]->Fit(func[i][j],"IR0");
	  hDtof[i][j]->GetXaxis()->SetRangeUser(-6,8);
	  hDtof[i][j]->Draw();
	  double mean = func[i][j]->GetParameter(1);
	  TLine *line = GetLine(mean,0,mean,func[i][j]->GetParameter(0)*1.02);
	  line->Draw();
	  TPaveText *t = GetTitleText(Form("(BL,Mod) = (%d,%d)",bl,mod),0.05);
	  t->Draw();

	  hMean->SetBinContent(i*5+j+1,mean);
	  hMean->SetBinError(i*5+j+1,func[i][j]->GetParError(1));
	}
      if(saveHisto && cTof) 
	{
	  cTof->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/calib_dTof/Fit_Dtof_BL%d.png",run_type,bl));
	  cTof->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/calib_dTof/Fit_Dtof_BL%d.pdf",run_type,bl));
	}
    }
  hMean->GetYaxis()->SetRangeUser(-0.6,0.5);
  hMean->SetMarkerStyle(20);
  c = draw1D(hMean);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/calib_dTof/DtofMean_vs_Mod.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/calib_dTof/DtofMean_vs_Mod.pdf",run_type));
    }

  if(save) 
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.CalibDtof.offline.root",run_type),"recreate");
      hDtofVsMod->Write();
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      if(hDtof[i][j]->GetEntries()==0) continue;
	      hDtof[i][j]->Write();
	      func[i][j]->Write();
	    }
	}
      fout->Close();

      // database
      std::ofstream outfile;
      outfile.open(Form("Rootfiles/%s.CalibDtof.offline.dat",run_type),std::fstream::out);

      for(Int_t i=0; i<30; i++)
      	{
      	  for(Int_t j=0; j<5; j++)
      	    {
	      double mean = 0;
	      if(hDtof[i][j]->GetEntries()>0) mean = -1*func[i][j]->GetParameter(1);
	      outfile << std::setprecision(4);
      	      outfile << i+1 << "  " << j+1 << "  " << mean << "\n";
      	    }
      	}
      outfile.close();
    }
}



/*
//================================================
double CrystalBall(double *x, double *par)
{
  double t = (x[0]-par[2])/par[3];
  if(par[0]<0) t = -t;

  double absAlpha = TMath::Abs(par[0]);
  double cb;
  if(t>-absAlpha)
    {
      cb = par[4]*exp(-0.5*t*t);
    }
  else
    {
      double aa = TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
      double bb = par[1]/absAlpha-absAlpha;
      cb = par[4]*aa/TMath::Power(bb-t,par[1]);
    }

  double f = cb + par[5]*exp(-0.5*TMath::Power((x[0]-par[6])/par[7],2));
  return cb;
}

//================================================
void makeCalib(const int save = 1)
{
  TH2F *hDtofVsMod = (TH2F*)f->Get(Form("mhDeltaTof_%s",trigName[kTrigType]));
  hDtofVsMod->GetYaxis()->SetRangeUser(-10,10);
  c = draw2D(hDtofVsMod,Form("%s: #Deltatof of tracks matched to MTD;p_{T} (GeV/c);#Deltatof (ns)",trigName[kTrigType]));

  TH1F *hDtof[30][5];
  TF1  *func[30][5];
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  int bl = i+1, mod = j+1;
	  hDtof[i][j] = (TH1F*)hDtofVsMod->ProjectionY(Form("hDtof_BL%d_Mod%d",bl,mod),i*5+j+1,i*5+j+1);
	  if(hDtof[i][j]->GetEntries()==0) continue;
	  if(j==0) 
	    {
	      TCanvas *cTof = new TCanvas(Form("hDtof_BL%d",bl),Form("hDtof_BL%d",bl),1100,750);
	      cTof->Divide(3,2);
	    }
	  cTof->cd(j+1);
	  func[i][j] = new TF1(Form("Fit_BL%d_Mod%d",bl,mod),CrystalBall,-1,5,5);
	  //func->SetParameters(-1, 2, -0.5, 0.5, 1e4, 10, -2.5, 2);
	  func[i][j]->SetParameters(-0.5, 1, -0.5, 0.2, hDtof[i][j]->GetMaximum());
	  hDtof[i][j]->Fit(func[i][j],"IR0N");
	  // hDtof[i][j]->GetXaxis()->SetRangeUser(-3,5);
	  // hDtof[i][j]->Draw();
	  // func->Draw("sames");
	  // TF1 *func1 = new TF1(Form("Gaus1_BL%d_Mod%d",bl,mod),"gaus");
	  // double mean1 = func->GetParameter(6);
	  // double sigma1 = func->GetParameter(7);
	  // func1->SetParameters(func->GetParameter(5),mean1,sigma1);
	  // func1->SetRange(mean1-3*sigma1,mean1+3*sigma1);
	  // func1->SetLineColor(4);
	  // func1->Draw("sames");

	  // TF1 *func2 = new TF1(Form("Gaus2_BL%d_Mod%d",bl,mod),"gaus");
	  // double mean2 = func->GetParameter(2);
	  // double sigma2 = func->GetParameter(3);
	  // func2->SetParameters(func->GetParameter(4),mean2,sigma2);
	  // func2->SetRange(mean2-3*sigma2,mean2+3*sigma2);
	  // func2->SetLineColor(2);
	  // func2->Draw("sames");
	}
    }

  if(save) 
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.CalibDtof.offline.root",run_type),"recreate");
      hDtofVsMod->Write();
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      hDtof[i][j]->Write();
	      func[i][j]->Write();
	    }
	}
      fout->Close();
    }
}


*/
