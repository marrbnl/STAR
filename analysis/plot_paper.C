const char *run_type = "Run14_AuAu200";
TString run_cfg_name = "paper";
const double luminosity = 14.2;
const double ppInelastic = 42.; // mb
const double ppInelasticErr = 3.; // mb

//================================================
void plot_paper()
{  
  //rawSignal();
  //efficiency();
  //xsec();
  //nPart();
  ppRef();
  //makeModel();
}

//================================================
void makeModel(const bool savePlot = 1, const bool saveHisto = 1)
{
  gStyle->SetOptStat(0);
  const int nCentBins = 6;
  const char *cent_Name[nCentBins] = {"0-80","0-20","20-40","40-60","60-80","0-60"};
  const char *cent_Title[nCentBins] = {"0080","0020","2040","4060","6080","0060"};
  const int nPtBins = 2;
  const double ptBins_low[nPtBins]  = {0,5};

  const int nModel = 3;
  const char* model_name[4] = {"TAMU", "Tsinghua", "SHM", "Co-mover"};
  const int color[4] = {1, 2, 4, 6};
  TGraphErrors *gRaaVsPt[nModel][nCentBins];
  TGraphErrors *gRaaVsNpart[nModel][nPtBins];
  TGraphErrors *gV2VsPt[nModel];
  for(int i=0; i<nModel; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  gRaaVsPt[i][j] = 0x0;
	}
      for(int j=0; j<nPtBins; j++)
	{
	  gRaaVsNpart[i][j] = 0x0;
	}
      gV2VsPt[i] = 0x0;
    }
  
  const int gNpoints = 100;
  double npart[nPtBins][gNpoints], raa_npart[nPtBins][gNpoints], raa_npart_min[nPtBins][gNpoints], raa_npart_max[nPtBins][gNpoints];
  double pt[gNpoints], raa_pt[nCentBins][gNpoints], raa_pt_min[nCentBins][gNpoints], raa_pt_max[nCentBins][gNpoints];
  double v2_pt[gNpoints], v2_pt_min[gNpoints], v2_pt_max[gNpoints];
  char tmp_str[256];
  double tmp_d;

  //==============================================
  // TAMU model
  // PRC82(2010)064905, NPA943(2015)147
  ifstream file;
  file.open("Rootfiles/Paper/models/TAMU_AuAu200GeV_Jpsi.txt");
  cout << file.is_open() << endl;
  for(int i=0; i<29; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gRaaVsNpart[0][0] = new TGraphErrors(29);
  gRaaVsNpart[0][0]->SetName(Form("gRaaVsNpart_Pt%1.0f_%s",ptBins_low[0],model_name[0]));
  for(int i=0; i<29; i++)
    {
      file >> npart[0][i] >> raa_npart_min[0][i] >> raa_npart_max[0][i];
      gRaaVsNpart[0][0]->SetPoint(i, npart[0][i], (raa_npart_max[0][i]+raa_npart_min[0][i])/2);
      gRaaVsNpart[0][0]->SetPointError(i, 0, fabs((raa_npart_max[0][i]-raa_npart_min[0][i])/2));
    }
  for(int i=0; i<5; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gRaaVsNpart[0][1] = new TGraphErrors(29);
  gRaaVsNpart[0][1]->SetName(Form("gRaaVsNpart_Pt%1.0f_%s",ptBins_low[1],model_name[0]));
  for(int i=0; i<29; i++)
    {
      file >> npart[1][i] >> raa_npart_min[1][i] >> raa_npart_max[1][i];
      gRaaVsNpart[0][1]->SetPoint(i, npart[1][i], (raa_npart_max[1][i]+raa_npart_min[1][i])/2);
      gRaaVsNpart[0][1]->SetPointError(i, 0, fabs((raa_npart_max[1][i]-raa_npart_min[1][i])/2));
    }
  for(int i=0; i<7; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }

  for(int k=0; k<nCentBins; k++)
    {
      gRaaVsPt[0][k] = new TGraphErrors(76);
      gRaaVsPt[0][k]->SetName(Form("gRaaVsPt_cent%s_%s",cent_Title[k],model_name[0]));
    }
  for(int i=0; i<76; i++)
    {
      file >> pt[i] 
	   >> raa_pt_min[1][i] >> raa_pt_max[1][i]
	   >> raa_pt_min[2][i] >> raa_pt_max[2][i] 
	   >> raa_pt_min[3][i] >> raa_pt_max[3][i]
	   >> raa_pt_min[4][i] >> raa_pt_max[4][i] 
	   >> raa_pt_min[0][i] >> raa_pt_max[0][i] 
	   >> raa_pt_min[5][i] >> raa_pt_max[5][i];
      for(int k=0; k<nCentBins; k++)
	{
	  gRaaVsPt[0][k]->SetPoint(i, pt[i], (raa_pt_max[k][i]+raa_pt_min[k][i])/2);
	  gRaaVsPt[0][k]->SetPointError(i, 0, fabs((raa_pt_max[k][i]-raa_pt_min[k][i])/2));
	}
    }
  for(int i=0; i<7; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gV2VsPt[0] = new TGraphErrors(50);
  gV2VsPt[0]->SetName(Form("gV2VsPt_cent0080_%s",model_name[0]));
  for(int i=0; i<50; i++)
    {
      file >> pt[i] >> v2_pt_min[i] >> v2_pt_max[i];
      gV2VsPt[0]->SetPoint(i, pt[i], (v2_pt_max[i]+v2_pt_min[i])/2);
      gV2VsPt[0]->SetPointError(i, 0, fabs((v2_pt_max[i]-v2_pt_min[i])/2));
    }
  file.close();

  //==============================================
  // Tsinghua model
  file.open("Rootfiles/Paper/models/Tsinghua_AuAu200GeV_Jpsi.txt");
  cout << file.is_open() << endl;
  for(int i=0; i<7; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gRaaVsNpart[1][0] = new TGraphErrors(21);
  gRaaVsNpart[1][0]->SetName(Form("gRaaVsNpart_Pt%1.0f_%s",ptBins_low[0],model_name[1]));
  for(int i=0; i<21; i++)
    {
      file >> npart[0][i] >> tmp_d >> tmp_d >> tmp_d >> raa_npart[0][i];
      gRaaVsNpart[1][0]->SetPoint(i, npart[0][i], raa_npart[0][i]);
      gRaaVsNpart[1][0]->SetPointError(i, 0, 0);
    }
  for(int i=0; i<5; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gRaaVsNpart[1][1] = new TGraphErrors(21);
  gRaaVsNpart[1][1]->SetName(Form("gRaaVsNpart_Pt%1.0f_%s",ptBins_low[1],model_name[1]));
  for(int i=0; i<21; i++)
    {
      file >> npart[1][i] >> tmp_d >> raa_npart[1][i];
      gRaaVsNpart[1][1]->SetPoint(i, npart[1][i], raa_npart[1][i]);
      gRaaVsNpart[1][1]->SetPointError(i, 0, 0);
    }


  for(int k=0; k<nCentBins; k++)
    {
      if(k==0 || k==1 || k==4) gRaaVsPt[1][k] = new TGraphErrors(14);
      else gRaaVsPt[1][k] = new TGraphErrors(15);
      gRaaVsPt[1][k]->SetName(Form("gRaaVsPt_cent%s_%s",cent_Title[k],model_name[1]));

      for(int i=0; i<4; i++)
	{
	  file.getline(tmp_str,256);
	  cout << k << "  " << tmp_str << endl;
	}
      for(int i=0; i<gRaaVsPt[1][k]->GetN(); i++)
	{
	  if(k==4) file >> pt[i] >> raa_pt[k][i];
	  else     file >> pt[i] >> tmp_d >> tmp_d >> tmp_d >> raa_pt[k][i];
	  gRaaVsPt[1][k]->SetPoint(i, pt[i], raa_pt[k][i]);
	  gRaaVsPt[1][k]->SetPointError(i, 0, 0);
	}
    }
  for(int i=0; i<4; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gV2VsPt[1] = new TGraphErrors(15);
  gV2VsPt[1]->SetName(Form("gV2VsPt_cent0080_%s",model_name[1]));
  for(int i=0; i<14; i++)
    {
      file >> pt[i] >> tmp_d >> v2_pt[i];
      gV2VsPt[1]->SetPoint(i, pt[i], v2_pt[i]);
      gV2VsPt[1]->SetPointError(i, 0, 0);
    }
  file.close();

  //==============================================
  // SHM model
  file.open("Rootfiles/Paper/models/SHM_AuAu200GeV_Jpsi.txt");
  cout << file.is_open() << endl;
  for(int i=0; i<19; i++)
    {
      file.getline(tmp_str,256);
      cout << tmp_str << endl;
    }
  gRaaVsNpart[2][0] = new TGraphErrors(17);
  gRaaVsNpart[2][0]->SetName(Form("gRaaVsNpart_Pt%1.0f_%s",ptBins_low[0],model_name[2]));
  for(int i=0; i<17; i++)
    {
      file >> npart[0][i] >> tmp_d >> raa_npart[0][i];
      gRaaVsNpart[2][0]->SetPoint(i, npart[0][i], raa_npart[0][i]);
      gRaaVsNpart[2][0]->SetPointError(i, 0, 0);
    }
  file.close();


  //==============================================
  // compare results from models
  TH1F *hplotVsNpart = new TH1F("hplotVsNpart",";N_{part};R_{AA}",100,0,400);
  TCanvas *cRaaVsNpart = new TCanvas("cRaaVsNpart", "cRaaVsNpart", 1100, 500);
  cRaaVsNpart->Divide(2,1);
  double x,y;
  for(int j=0; j<nPtBins; j++)
    {
      cRaaVsNpart->cd(j+1);
      hplotVsNpart->DrawCopy();
      for(int i=0; i<nModel; i++)
	{
	  if(!gRaaVsNpart[i][j]) continue;
	  gRaaVsNpart[i][j]->SetFillStyle(3004+i);
	  gRaaVsNpart[i][j]->SetFillColor(color[i]);
	  gRaaVsNpart[i][j]->SetLineColor(color[i]);
	  gRaaVsNpart[i][j]->SetMarkerColor(color[i]);
	  if(i==0) gRaaVsNpart[i][j]->Draw("samesE4");
	  else gRaaVsNpart[i][j]->Draw("samesL");
	  int npoints = gRaaVsNpart[i][j]->GetN();
	  TGraph *gRaaUp = new TGraph(npoints);
	  TGraph *gRaaDown = new TGraph(npoints);
	  for(int ipoint=0; ipoint<npoints; ipoint++)
	    {
	      gRaaVsNpart[i][j]->GetPoint(ipoint, x, y);
	      double error = gRaaVsNpart[i][j]->GetErrorY(ipoint);
	      gRaaUp->SetPoint(ipoint, x, y+error);
	      gRaaDown->SetPoint(ipoint, x, y-error);
	    }
	  gRaaUp->SetLineColor(color[i]);
	  gRaaDown->SetLineColor(color[i]);
	  gRaaUp->Draw("samesL");
	  gRaaDown->Draw("samesL");
	}
      TPaveText *t1 = GetTitleText(Form("J/#Psi nuclear modification factor (p_{T} > %1.0f GeV/c)",ptBins_low[j]));
      t1->Draw();
      if(j==0) 
	{
	  TLegend *leg = new TLegend(0.5,0.85-nModel*0.05,0.75,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  for(int i=0; i<nModel; i++)
	    leg->AddEntry(gRaaVsNpart[i][j],model_name[i],"F");
	  leg->Draw();
	}
    }

  TH1F *hplotVsPt = new TH1F("hplotVsPt",";p_{T} (GeV/c);R_{AA}",100,0,16);
  TCanvas *cRaaVsPt = new TCanvas("cRaaVsPt", "cRaaVsPt", 1100, 700);
  cRaaVsPt->Divide(3,2);
  double x,y;
  for(int k=0; k<nCentBins; k++)
    {
      cRaaVsPt->cd(k+1);
      hplotVsPt->DrawCopy();
      for(int i=0; i<nModel; i++)
	{
	  if(!gRaaVsPt[i][k]) continue;
	  gRaaVsPt[i][k]->SetFillStyle(3004+i);
	  gRaaVsPt[i][k]->SetFillColor(color[i]);
	  gRaaVsPt[i][k]->SetLineColor(color[i]);
	  gRaaVsPt[i][k]->SetMarkerColor(color[i]);
	  if(i==0) gRaaVsPt[i][k]->Draw("samesE4");
	  else gRaaVsPt[i][k]->Draw("samesL");
	  int npoints =  gRaaVsPt[i][k]->GetN();
	  TGraph *gRaaUp = new TGraph(npoints);
	  TGraph *gRaaDown = new TGraph(npoints);
	  for(int ipoint=0; ipoint<npoints; ipoint++)
	    {
	      gRaaVsPt[i][k]->GetPoint(ipoint, x, y);
	      double error =  gRaaVsPt[i][k]->GetErrorY(ipoint);
	      gRaaUp->SetPoint(ipoint, x, y+error);
	      gRaaDown->SetPoint(ipoint, x, y-error);
	    }
	  gRaaUp->SetLineColor(color[i]);
	  gRaaDown->SetLineColor(color[i]);
	  gRaaUp->Draw("samesL");
	  gRaaDown->Draw("samesL");
	}
      TPaveText *t1 = GetTitleText(Form("J/#Psi nuclear modification factor (%s%%)",cent_Name[k]));
      t1->Draw();
      if(k==0) 
	{
	  TLegend *leg = new TLegend(0.5,0.85-nModel*0.05,0.75,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  for(int i=0; i<nModel; i++)
	    leg->AddEntry(gRaaVsPt[i][k],model_name[i],"F");
	  leg->Draw();
	}
    }  

  TH1F *hplotV2 = new TH1F("hplotV2",";p_{T} (GeV/c);v_{2}",100,0,10);
  TCanvas *cV2VsPt = new TCanvas("cV2VsPt", "cV2VsPt", 600, 400);
  hplotV2->GetYaxis()->SetRangeUser(0,0.1);
  hplotV2->DrawCopy();
  for(int i=0; i<nModel; i++)
    {
      if(!gV2VsPt[i]) continue;
      gV2VsPt[i]->SetFillStyle(3004+i);
      gV2VsPt[i]->SetFillColor(color[i]);
      gV2VsPt[i]->SetLineColor(color[i]);
      gV2VsPt[i]->SetMarkerColor(color[i]);
      if(i==0) gV2VsPt[i]->Draw("samesE4");
      else gV2VsPt[i]->Draw("samesL");
      int npoints = gV2VsPt[i]->GetN();
      TGraph *gRaaUp = new TGraph(npoints);
      TGraph *gRaaDown = new TGraph(npoints);
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  gV2VsPt[i]->GetPoint(ipoint, x, y);
	  double error = gV2VsPt[i]->GetErrorY(ipoint);
	  gRaaUp->SetPoint(ipoint, x, y+error);
	  gRaaDown->SetPoint(ipoint, x, y-error);
	}
      gRaaUp->SetLineColor(color[i]);
      gRaaDown->SetLineColor(color[i]);
      gRaaUp->Draw("samesL");
      gRaaDown->Draw("samesL");
    }
  TPaveText *t1 = GetTitleText(Form("J/#Psi elliptic flow (0-80%%)"),0.06);
  t1->Draw();
  TLegend *leg = new TLegend(0.5,0.85-nModel*0.05,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(int i=0; i<nModel; i++)
    leg->AddEntry(gV2VsPt[i],model_name[i],"F");
  leg->Draw();

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Paper/models/AuAu200.models.root","recreate");
      for(int i=0; i<nModel; i++)
	{
	  for(int j=0; j<nPtBins; j++)
	    {
	      if(!gRaaVsNpart[i][j]) continue;
	      gRaaVsNpart[i][j]->Write("", TObject::kOverwrite);
	    }
	  for(int k=0; k<nCentBins; k++)
	    {
	      if(!gRaaVsPt[i][k]) continue;
	      gRaaVsPt[i][k]->Write("", TObject::kOverwrite);
	    }
	  if(!gV2VsPt[i]) continue;
	  gV2VsPt[i]->Write("", TObject::kOverwrite);
	} 
    }
}

//================================================
void efficiency(const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptStat(0);

  // Online timing cut efficiency
  TFile *fTrigEff = TFile::Open("Rootfiles/Run14_AuAu200.Sys.MtdTrigEff.root","read");
  TF1 *fucnTrig = (TF1*)fTrigEff->Get("Run14_AuAu200_Muon_TacDiffEff");
  TF1 *fucnTrigUp = (TF1*)fTrigEff->Get("Run14_AuAu200_Muon_TacDiffEff_Sysup");
  TF1 *fucnTrigDown = (TF1*)fTrigEff->Get("Run14_AuAu200_Muon_TacDiffEff_Sysdown");
  TGraphAsymmErrors *gTrigEffSys = new TGraphAsymmErrors(1000);
  for(int ipoint=0; ipoint<gTrigEffSys->GetN(); ipoint++)
    {
      double x = ipoint * 0.01 + 1.3;
      double y = fucnTrig->Eval(x);
      double yh = fucnTrigUp->Eval(x) - y;
      double yl = y - fucnTrigDown->Eval(x);
      gTrigEffSys->SetPoint(ipoint, x, y);
      gTrigEffSys->SetPointError(ipoint, 0.005, 0.005, yl, yh);
    }
  cout << "[i] Trigger eff: "<< fucnTrig->Eval(1.3) << " at 1.3 GeV/c" << endl;
  cout << "[i] Trigger eff: "<< fucnTrig->Eval(10) << " at 10 GeV/c" << endl;

  // PID efficiency
  TFile *fpid = TFile::Open("Rootfiles/Run14_AuAu200.EmbTrkEff.root","read");
  TH1F *hTrkPtMtdMth = (TH1F*)fpid->Get("McTrkPt_MtdMth_cent0080");
  TH1F *hTrkPtPid = (TH1F*)fpid->Get("McTrkPt_MuonPid_cent0080");
  hTrkPtMtdMth->Rebin(2);
  hTrkPtPid->Rebin(2);
  TH1F *hPidEff = (TH1F*)hTrkPtPid->Clone("hPidEff");
  hPidEff->Divide(hTrkPtMtdMth);
  hPidEff->GetXaxis()->SetRangeUser(1.0, 10);
  c = draw1D(hPidEff,"Muon PID efficiency");
  TF1 *funcPidLowPt = new TF1("funcPidLowPt","[0]-exp(-1*[1]*(x-[2]))",1.1,3.0);
  funcPidLowPt->SetParameters(0.795, 2.59, 0.26);
  hPidEff->Fit(funcPidLowPt, "IR0");

  TF1 *funcPidHighPt = new TF1("funcPidHighPt","[0]-exp(-1*[1]*(x-[2]))",3.0,10.0);
  //TF1 *funcPidHighPt = new TF1("funcPidHighPt","pol0",3.0,10.0);
  funcPidHighPt->SetParameters(0.845, 2.59, 0.26);
  hPidEff->Fit(funcPidHighPt, "IR0");
  funcPidLowPt->Draw("sames");
  funcPidHighPt->Draw("sames");

  TFile *fdtof = TFile::Open("Rootfiles/Run14_AuAu200.DtofEff.root","read");
  TF1 *fucnDtof = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_FitFunc");
  TF1 *fucnDtofUp = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_Sysup");
  TF1 *fucnDtofDown = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_Sysdown");
  TGraphAsymmErrors *gPidEffSys[2];
  gPidEffSys[0] = new TGraphAsymmErrors(170);
  gPidEffSys[1] = new TGraphAsymmErrors(730);
  for(int ipoint=0; ipoint<1000; ipoint++)
    {
      double x = ipoint * 0.01 + 1.3;
      double yh = sqrt(0.035*0.035+pow(fabs(fucnDtofUp->Eval(x)/fucnDtof->Eval(x)-1),2))*y;
      double yl = sqrt(0.035*0.035+pow(fabs(fucnDtofDown->Eval(x)/fucnDtof->Eval(x)-1),2))*y;
      if(x<3)
	{
	  double y = funcPidLowPt->Eval(x);
	  gPidEffSys[0]->SetPoint(ipoint, x, y);
	  gPidEffSys[0]->SetPointError(ipoint, 0.005, 0.005, yl, yh);
	}
      else
	{
	  double y = funcPidHighPt->Eval(x);
	  gPidEffSys[1]->SetPoint(ipoint-170, x, y);
	  gPidEffSys[1]->SetPointError(ipoint-170, 0.005, 0.005, yl, yh);
	}
    }

  TCanvas *cEff = new TCanvas("MtdEff", "MtdEff", 600, 700);
  TPad *pads[2];
  pads[0] = GetSinglePad("pad_0", 0.01, 0.98, 0.54, 0.99);
  pads[0]->Draw();
  pads[0]->cd();
  SetPadMargin(gPad,0.01,0.12,0.01,0.05);
  TH1F *hplot = new TH1F("hplot",";p_{T,#mu} [GeV/c];Efficiency", 1100, 0, 11);
  hplot->GetYaxis()->SetRangeUser(0.45, 1.01);
  hplot->GetXaxis()->SetRangeUser(1.2, 6);
  ScaleHistoTitle(hplot,0.06,1.1,0.06,22,1.6,20,63);
  hplot->GetXaxis()->SetLabelSize(0);
  hplot->DrawCopy();
  gTrigEffSys->SetFillStyle(1001);
  gTrigEffSys->SetFillColor(kRed-10);
  gTrigEffSys->SetLineColor(kRed);
  gTrigEffSys->SetLineWidth(3);
  gTrigEffSys->Draw("samesE2L");
  TLegend *leg = new TLegend(0.5,0.2,0.75,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(22);;
  leg->AddEntry(gTrigEffSys,"MTD trigger efficiency","FL");
  leg->Draw();
  cEff->Update();

  cEff->cd();
  pads[1] = GetSinglePad("pad_1", 0.01, 0.98, 0.01, 0.54);
  pads[1]->Draw();
  pads[1]->cd();
  SetPadMargin(gPad,0.18,0.12,0.01,0.01);
  ScaleHistoTitle(hplot,22,2.5,20,22,1.6,20,63);
  hplot->DrawCopy();
  for(int i=0; i<2; i++)
    {
      gPidEffSys[i]->SetFillStyle(1001);
      gPidEffSys[i]->SetFillColor(kGray);
      gPidEffSys[i]->SetLineColor(kBlack);
      gPidEffSys[i]->SetLineWidth(3);
      gPidEffSys[i]->Draw("samesE2L");
    }
  TLegend *leg = new TLegend(0.5,0.3,0.75,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(22);;
  leg->AddEntry(gPidEffSys[0],"Muon PID efficiency","FL");
  leg->Draw();
  cEff->cd();
  cEff->Update();

  if(savePlot)
    {
      cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/MtdEfficiency.pdf",run_type,run_cfg_name.Data()));      
      cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/MtdEfficiency.png",run_type,run_cfg_name.Data()));   
      cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/MtdEfficiency.eps",run_type,run_cfg_name.Data()));
    }
  if(gSavePaper)
    {
      cEff->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_MtdEfficiency.eps",gPaperVersion));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type),"update");
      gTrigEffSys->Write("MTD_TrigEff", TObject::kOverwrite);
      gPidEffSys[0]->Write("MTD_PidEff_LowPt", TObject::kOverwrite);
      gPidEffSys[1]->Write("MTD_PidEff_HighPt", TObject::kOverwrite);
    }
    
}

//================================================
void nPart(const bool savePlot = 0, const bool saveHisto = 0)
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

  gStyle->SetOptStat(0);

  // Ncoll uncertainty
  const double ncoll[kNCent] =    {13.09, 29.66, 62.18, 120.06, 216.45, 366.44, 593.67, 941.24}; 
  const double ncollErr[kNCent] = {5.88, 10.90, 16.85, 23.36, 29.20, 32.33, 30.18, 26.27};
  const double npart[kNCent] =    {14.39, 27.38, 47.79, 76.88, 116.74, 169.04, 237.27, 325.48};
  const double npartErr[kNCent] = {5.2, 7.7, 9.5, 10.7, 11.1, 10.5, 8.5, 3.6};

  const double ncoll2[7] =    {21.57, 62.18, 120.06, 216.45, 366.44, 593.67, 941.24}; 
  const double npart2[7] =    {21.04, 47.79, 76.88, 116.74, 169.04, 237.27, 325.48}; 
  
  double npart_y[kNCent];
  double ncoll_relerr[kNCent];
  for(int k=0; k<kNCent; k++)
    {
      npart_y[k] = 1;
      ncoll_relerr[k] = ncollErr[k]/ncoll[k];
    }

  TGraphErrors *gPubNpartErr = new TGraphErrors(kNCent, npart, npart_y, npartErr, ncoll_relerr);
  gPubNpartErr->SetFillColor(kGreen-10);
  gPubNpartErr->SetLineColor(kGreen-10);

  // plotting template
  TH1F *hRaaVsNpart = new TH1F("hRaaVsNpart",";N_{part};R_{AA}",100,-100,370);
  ScaleHistoTitle(hRaaVsNpart,28,1,24,28,1,24,63);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,2);
  hRaaVsNpart->GetYaxis()->CenterTitle();

 // MTD results
  TFile *fout = 0x0;
  if(saveHisto) fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type),"update");
  else fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type),"read");
  TGraphAsymmErrors *hpp = (TGraphAsymmErrors*)fout->Get("hpp200JpsiVsCentFinalSys");

  TFile *fdata = TFile::Open(Form("Rootfiles/%s.JpsiXsec.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  TFile *fsys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type), "read");
  TH1F *hdata[nPtBins], *hsys[nPtBins];
  TH1F *hRaa[nPtBins];
  TGraphErrors *raaVsNpart[nPtBins], *raaVsNpartSys[nPtBins];
  TBox *globalSys[nPtBins];
  double x,y;
  for(int i=0; i<nPtBins; i++)
    {
      hdata[i] = (TH1F*)fdata->Get(Form("Jpsi_InvYieldVsCent_pt%s",pt_Name[i]));
      hsys[i] = (TH1F*)fsys->Get(Form("JpsiSysVsCent_All_Pt%s",pt_Name[i]));
      int npoints = nCentBins[i];
      raaVsNpart[i] = new TGraphErrors(npoints);
      raaVsNpart[i]->SetName(Form("Jpsi_RaaVsNpart_pt%s",pt_Name[i]));
      raaVsNpartSys[i] = new TGraphErrors(npoints);
      raaVsNpartSys[i]->SetName(Form("Jpsi_RaaVsNpartSys_pt%s",pt_Name[i]));

      hRaa[i] = (TH1F*)hdata[i]->Clone(Form("hRaaVsNpart_pt%s",pt_Name[i]));
      hRaa[i]->Reset();

      hpp->GetPoint(i, x, y);
      double pp_yield = y;
      double pp_err_h = hpp->GetErrorYlow(i);
      double pp_err_l = hpp->GetErrorYhigh(i);
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  double x = npart[ipoint];
	  if(i==1) x = npart2[ipoint];
	  double NCOLL;
	  if(i==0) NCOLL = ncoll[ipoint];
	  if(i==1) NCOLL = ncoll2[ipoint];
	  double auau_yield = hdata[i]->GetBinContent(ipoint+1) * 2 * pi;
	  double auau_err   = hdata[i]->GetBinError(ipoint+1) * 2 * pi;
	  double raa = ppInelastic/NCOLL * 1e6 * auau_yield/pp_yield;
	  double raa_err =  auau_err / auau_yield * raa;
	  double raa_sys = hsys[i]->GetBinContent(ipoint+1) * raa;
	  raaVsNpart[i]->SetPoint(ipoint,x,raa);
	  raaVsNpart[i]->SetPointError(ipoint,0,raa_err);
	  raaVsNpartSys[i]->SetPoint(ipoint,x,raa);
	  raaVsNpartSys[i]->SetPointError(ipoint,5,raa_sys);
	  hRaa[i]->SetBinContent(ipoint+1, raa);
	  hRaa[i]->SetBinError(ipoint+1, raa_err);
	  cout << auau_yield << "  " << pp_yield << "  " << raa << "  " << raa_err << "  " << raa_sys << endl;
	}
      double gSys_h = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err_h*pp_err_h/pp_yield/pp_yield);
      double gSys_l = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err_l*pp_err_l/pp_yield/pp_yield);
      globalSys[i] = new TBox(360,1-gSys_l,365,1+gSys_h);
      globalSys[i]->SetLineWidth(2.);
      globalSys[i]->SetFillStyle(1001);
    }

  //==============================================
  // LHC results
  //==============================================

  // ALICE Low pT
  // http://hepdata.cedar.ac.uk/view/ins1263062
  // Phys. Lett. B 734 (2014) 314-327

  // CMS high pT
  // ARXIV:1610.00613

  TGraphErrors *grLhcRaaVsCent[2];
  TGraphErrors *grLhcRaaVsCentSys[2];

  const int nAlice_lowpT = 3;
  double aliceLowPt_npart[nAlice_lowpT] = {356.0, 191.5, 37.9};
  double aliceLowPt_npart_stat[nAlice_lowpT] = {0, 0, 0};
  double aliceLowPt_npart_sys[nAlice_lowpT] = {6,6,6};
  double aliceLowPt_raa[nAlice_lowpT] = {0.73, 0.70, 0.79};
  double aliceLowPt_raa_stat[nAlice_lowpT] = {0.09, 0.08, 0.15};
  double aliceLowPt_raa_sys[nAlice_lowpT] = {0.06, 0.05, 0.09};
  grLhcRaaVsCent[0] = new TGraphErrors(nAlice_lowpT, aliceLowPt_npart, aliceLowPt_raa, aliceLowPt_npart_stat, aliceLowPt_raa_stat);
  grLhcRaaVsCentSys[0]  = new TGraphErrors(nAlice_lowpT, aliceLowPt_npart, aliceLowPt_raa, aliceLowPt_npart_sys, aliceLowPt_raa_sys);

  const int nCms_highpT = 6;
  double cms_npart[nCms_highpT] = {355.4, 261.4, 187.2, 130.0, 86.3, 22.1};
  double cms_npart_err[nCms_highpT] = {0.,0.,0.,0.,0.,0.};
  double cms_npart_sys[nCms_highpT] = {6,6,6,6,6,6};
  double cms_raa[nCms_highpT] = {0.28, 0.33, 0.44, 0.50, 0.57, 0.77};
  double cms_raa_err[nCms_highpT] = {0.011, 0.013, 0.016, 0.033, 0.033, 0.049};
  double cms_raa_sys[nCms_highpT] = {0.026, 0.034, 0.049, 0.065, 0.081, 0.13};
  grLhcRaaVsCent[1] = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_err, cms_raa_err);
  grLhcRaaVsCentSys[1] = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_sys, cms_raa_sys);

  //==============================================
  // theory
  //==============================================
  // tsinghua group

  char tmp[256];
  double dtmp;
  ifstream ifin;
  ifin.open("Rootfiles/Published/Jpsi_Raa_200/model/data_and_figures_from_Yunpeng/data_pt_5-10_GeV.txt");
  ifin.getline(tmp,256);
  ifin.getline(tmp,256);
  double raa_5_tsu[11];
  double npart_tsu[11];
  for(int i=0; i<11; i++)
    {
      ifin >> npart_tsu[i] >> raa_5_tsu[i];
    }
  ifin.close();
  TGraph *gTsuRaaVsNpartRHIC[2];
  gTsuRaaVsNpartRHIC[0] = new TGraph("Rootfiles/Published/Zebo/QWG2013/Jpsi_Cent/Jpsi_Zhuang_RHIC_mid.dat");
  gTsuRaaVsNpartRHIC[0]->SetName("gTsuRaaVsNpartRHIC_LowPt");
  gTsuRaaVsNpartRHIC[1] = new TGraph(11,npart_tsu,raa_5_tsu);
  gTsuRaaVsNpartRHIC[1]->SetName("gTsuRaaVsNpartRHIC_HighPt");

  TGraphErrors *gTsuRaaVsNpartLHC[2];
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/Mid_Y_RAA_Npart.dat");
  double tsu_lhc_lowPt_npart[26];
  double tsu_lhc_lowPt_npart_err[26];
  double tsu_lhc_lowPt_raa[26];
  double tsu_lhc_lowPt_raa_err[26];
  double tsu_lhc_lowPt_raa_low[26];
  double tsu_lhc_lowPt_raa_high[26];
  for(int i=0; i<11; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<26; i++)
    {
      ifin >> tsu_lhc_lowPt_npart[i] >> dtmp >> dtmp >> tsu_lhc_lowPt_raa_low[i] >> dtmp >> tsu_lhc_lowPt_raa_high[i];
      tsu_lhc_lowPt_npart_err[i] = 0;
      tsu_lhc_lowPt_raa[i] = (tsu_lhc_lowPt_raa_low[i] + tsu_lhc_lowPt_raa_high[i])/2;
      tsu_lhc_lowPt_raa_err[i] = (tsu_lhc_lowPt_raa_high[i] - tsu_lhc_lowPt_raa_low[i])/2;
    }
  ifin.close();
  gTsuRaaVsNpartLHC[0] = new TGraphErrors(26, tsu_lhc_lowPt_npart, tsu_lhc_lowPt_raa, tsu_lhc_lowPt_npart_err, tsu_lhc_lowPt_raa_err);
  gTsuRaaVsNpartLHC[0]->SetName("gTsuRaaVsNpartLHC_LowPt");

  //ifin.open("Rootfiles/Published/Model/Jpsi_Raa/Mid_Y_RAA_Npart_highpt.dat");
  ifin.open("Rootfiles/Published/Model/Jpsi_Raa/cms_raa_npart");
  double tsu_lhc_highPt_npart[24];
  double tsu_lhc_highPt_npart_err[24];
  double tsu_lhc_highPt_raa[24];
  double tsu_lhc_highPt_raa_err[24];
  for(int i=0; i<4; i++)
    {
      ifin.getline(tmp,256);
    }
  for(int i=0; i<24; i++)
    {
      ifin >> tsu_lhc_highPt_npart[i] >> tsu_lhc_highPt_raa[i];
      tsu_lhc_highPt_npart_err[i] = 0;
      tsu_lhc_highPt_raa_err[i] = 0;
      tsu_lhc_highPt_raa[i] = tsu_lhc_highPt_raa[i]*1.2;
    }
  ifin.close();
  gTsuRaaVsNpartLHC[1] = new TGraphErrors(24, tsu_lhc_highPt_npart, tsu_lhc_highPt_raa, tsu_lhc_highPt_npart_err, tsu_lhc_highPt_raa_err);
  gTsuRaaVsNpartLHC[0]->SetName("gTsuRaaVsNpartLHC_HighPt");
  
  // tamu group
  TGraph *gTamuRaaVsNpartRHIC[2], *gTamuRaaVsNpartLHC[2];
  ifin.open("Rootfiles/Published/Jpsi_Raa_200/model/data_from_Xingbo/raa_centra_rhic_pt4.5to10_tozebo_110513.dat");
  double raa_5_tamu[15];
  double npart_tamu[15];
  for(int i=0; i<15; i++)
    {
      ifin >> npart_tamu[i] >> raa_5_tamu[i];
    }
  ifin.close();
  gTamuRaaVsNpartRHIC[1] = new TGraph(15,npart_tamu,raa_5_tamu);
  gTamuRaaVsNpartRHIC[1]->SetName("gTamuRaaVsNpartRHIC_HighPt");

  TFile *ftamu = TFile::Open("Rootfiles/2016sQM/raa200theory1.root","read");
  TGraph *gtmp = (TGraph*)ftamu->Get("gr1");
  gTamuRaaVsNpartRHIC[0] = new TGraph(gtmp->GetN(), gtmp->GetX(), gtmp->GetY());
  gTamuRaaVsNpartRHIC[0]->SetName("gTamuRaaVsNpartRHIC_LowPt");

  gTamuRaaVsNpartLHC[0] = new TGraph("Rootfiles/Published/Zebo/QWG2013/Jpsi_Cent/Jpsi_Rapp_mid.dat");
  gTamuRaaVsNpartLHC[1] = new TGraph("Rootfiles/Published/Zebo/QWG2013/Jpsi_Cent/Jpsi_Rapp_LHC_mid_highPt.dat");

  
  //==============================================
  // Final figures
  //==============================================
  int marker_color[2] = {kRed, kBlue};
  double x_min = -10, x_max = 375;
  for(int i=0; i<nPtBins; i++)
    {
      TCanvas *c = new TCanvas(Form("Jpsi_raa_vs_npart_Pt%s",pt_Name[i]),Form("Jpsi_raa_vs_npart_Pt%s",pt_Name[i]),800,600);
      SetPadMargin(gPad,0.13,0.13,0.05,0.02);
      hRaaVsNpart->GetXaxis()->SetRangeUser(x_min, x_max);
      hRaaVsNpart->SetMaximum(2.0);
      hRaaVsNpart->DrawCopy();
      gPubNpartErr->Draw("e3sames");
      TLine *line = GetLine(x_min,1,x_max,1,1);
      line->Draw();

      gTamuRaaVsNpartLHC[i]->SetFillStyle(1001);
      gTamuRaaVsNpartLHC[i]->SetFillColor(kCyan-9);
      gTamuRaaVsNpartLHC[i]->SetLineColor(kCyan-9);
      gTamuRaaVsNpartLHC[i]->SetLineWidth(1);
      gTamuRaaVsNpartLHC[i]->SetLineStyle(1);

      gTsuRaaVsNpartLHC[i]->SetFillStyle(1001);
      gTsuRaaVsNpartLHC[i]->SetFillColor(kViolet-9);
      gTsuRaaVsNpartLHC[i]->SetLineColor(kViolet-9);
      gTsuRaaVsNpartLHC[i]->SetLineWidth(1);
      gTsuRaaVsNpartLHC[i]->SetLineStyle(1);

      if(i==1)
	{
	  gTamuRaaVsNpartLHC[1]->SetLineColor(marker_color[1]);
	  gTamuRaaVsNpartLHC[1]->SetLineWidth(2);
	  gTamuRaaVsNpartLHC[1]->SetLineStyle(7);
	  gTsuRaaVsNpartLHC[1]->SetLineWidth(2);
	  gTsuRaaVsNpartLHC[1]->SetLineColor(marker_color[1]);
	  gTsuRaaVsNpartLHC[1]->SetLineStyle(1);
	}

      gTamuRaaVsNpartRHIC[i]->SetLineColor(2);
      gTamuRaaVsNpartRHIC[i]->SetLineStyle(7);
      gTamuRaaVsNpartRHIC[i]->SetLineWidth(2);

      gTsuRaaVsNpartRHIC[i]->SetLineColor(2);
      gTsuRaaVsNpartRHIC[i]->SetLineStyle(1);
      gTsuRaaVsNpartRHIC[i]->SetLineWidth(2);

      if(i==0)
	{
	  gTamuRaaVsNpartLHC[i]->Draw("fsames");
	  gTsuRaaVsNpartLHC[i]->Draw("samesE4");
	}
      else
	{
	  gTamuRaaVsNpartLHC[i]->Draw("lsames");
	  gTsuRaaVsNpartLHC[i]->Draw("lsames");
	}
      gTamuRaaVsNpartRHIC[i]->Draw("lsames");
      gTsuRaaVsNpartRHIC[i]->Draw("lsames");


  
      TPaveText *tmodel2 = GetPaveText(0.213,0.35,0.64,0.74);
      tmodel2->SetTextFont(62);
      tmodel2->SetTextAlign(11);
      tmodel2->SetTextSize(0.035);
      tmodel2->AddText("Tsinghua Model");
      tmodel2->AddText("TAMU Model");
      tmodel2->Draw();

      TLegend *leg41 = new TLegend(0.42,0.65,0.6,0.75);
      leg41->SetBorderSize(0);
      leg41->SetFillColor(0);
      leg41->SetTextFont(62);
      leg41->SetTextSize(0.035);
      leg41->AddEntry(gTsuRaaVsNpartRHIC[i],"RHIC","L");
      leg41->AddEntry(gTamuRaaVsNpartRHIC[i],"RHIC","L");
      leg41->Draw();
      TLegend *leg42 = new TLegend(0.55,0.65,0.7,0.75);
      leg42->SetBorderSize(0);
      leg42->SetFillColor(0);
      leg42->SetTextFont(62);
      leg42->SetTextSize(0.035);
      if(i==0)
	{
	  leg42->AddEntry(gTsuRaaVsNpartLHC[i],"LHC","F");
	  leg42->AddEntry(gTamuRaaVsNpartLHC[i],"LHC","F");
	}
      else
	{
	  leg42->AddEntry(gTsuRaaVsNpartLHC[i],"LHC","L");
	  leg42->AddEntry(gTamuRaaVsNpartLHC[i],"LHC","L");
	}
      leg42->Draw();

      // LHC
      grLhcRaaVsCent[i]->SetMarkerStyle(21);
      grLhcRaaVsCent[i]->SetMarkerSize(1.5);
      grLhcRaaVsCent[i]->SetLineColor(marker_color[1]);
      grLhcRaaVsCent[i]->SetMarkerColor(marker_color[1]);
      grLhcRaaVsCentSys[i]->SetMarkerColor(marker_color[1]);
      grLhcRaaVsCentSys[i]->SetLineColor(marker_color[1]);
      grLhcRaaVsCentSys[i]->SetFillStyle(0);
      TBox *bLhc;
      if(i==0) bLhc = new TBox(365,1-0.13,370,1+0.13);
      else     bLhc = new TBox(365,1-0.078,370,1+0.078);
      bLhc->SetLineColor(marker_color[1]);
      bLhc->SetFillColor(marker_color[1]);
      bLhc->SetFillStyle(1001);
      grLhcRaaVsCentSys[i]->Draw("samesE5");
      grLhcRaaVsCent[i]->Draw("samesPEZ");
      bLhc->Draw("fsame");

      // STAR
      TGraphErrors *grStar = (TGraphErrors*)raaVsNpart[i]->Clone(Form("grStar_Pt%s",pt_Name[i]));
      TGraphErrors *grStarSys = (TGraphErrors*)raaVsNpartSys[i]->Clone(Form("grStar_Pt%s_Sys",pt_Name[i]));
      TBox *bStar = (TBox*)globalSys[i]->Clone(Form("bStar_Pt%s",pt_Name[i]));
      grStarSys->SetMarkerSize(0);
      grStarSys->SetFillStyle(0);
      grStarSys->SetMarkerColor(marker_color[0]);
      grStarSys->SetLineColor(marker_color[0]);
      grStarSys->Draw("samesE5");
      grStar->SetMarkerStyle(29);
      grStar->SetMarkerSize(2.5);
      grStar->SetMarkerColor(marker_color[0]);
      grStar->SetLineColor(marker_color[0]);
      grStar->Draw("samesPEZ");
      bStar->SetLineColor(marker_color[0]);
      bStar->SetFillColor(marker_color[0]);
      bStar->SetFillStyle(1001);
      bStar->Draw("fsame");

      TLegend *leg_ncoll_2 = new TLegend(0.15,0.17,0.35,0.22);
      leg_ncoll_2->SetBorderSize(0);
      leg_ncoll_2->SetFillColor(0);
      leg_ncoll_2->SetTextFont(62);
      leg_ncoll_2->SetTextSize(0.035);;
      leg_ncoll_2->AddEntry(gPubNpartErr,"STAR N_{coll} uncertainty","f");
      leg_ncoll_2->Draw();

      TLegend *leg3 = new TLegend(0.16,0.76,0.4,0.96);
      leg3->SetBorderSize(0);
      leg3->SetFillColor(0);
      leg3->SetTextFont(62);
      leg3->SetTextSize(0.04);
      if(i==0)
	{
	  leg3->SetHeader("p_{T,J/#psi} > 0 GeV/c");
	  leg3->AddEntry(grStar,"STAR: Au+Au, #sqrt{s_{NN}} = 200 GeV, J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
	  leg3->AddEntry(grLhcRaaVsCent[i],"ALICE: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, J/#psi#rightarrowe^{+}e^{-}, |y| < 0.8","P");
	}
      else
	{
	  leg3->SetHeader("J/#psi#rightarrow#mu^{+}#mu^{-}");
	  leg3->AddEntry(grStar,"STAR: Au+Au, #sqrt{s_{NN}} = 200 GeV, |y| < 0.5, p_{T} > 5 GeV/c","P");
	  leg3->AddEntry(grLhcRaaVsCent[i],"CMS: Pb+Pb, #sqrt{s_{NN}} = 2.76 TeV, |y| < 1.2, p_{T} > 6.5 GeV/c","P");
	}
      leg3->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiRaaVsNpart_Pt%s_RHICvsLHC.pdf",run_type,run_cfg_name.Data(),pt_Name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiRaaVsNpart_Pt%s_RHICvsLHC.png",run_type,run_cfg_name.Data(),pt_Name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiRaaVsNpart_Pt%s_RHICvsLHC.eps",run_type,run_cfg_name.Data(),pt_Name[i]));
	}
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_RaaVsNpart_Pt%1.0f.pdf",ptBins_low_npart[i]));
	}
      if(gSavePaper)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_RaaVsNpart_Pt%1.0f.pdf",gPaperVersion,ptBins_low_npart[i]));
	}
    }

  if(saveHisto)
    {
      fout->cd();
      for(int i=0; i<2; i++)
	{
	  raaVsNpart[i]->Write(Form("MTD_Run14AuAu_JpsiRaaVsNpart_pt%dGeV",i*5), TObject::kOverwrite);
	  raaVsNpartSys[i]->Write(Form("MTD_Run14AuAu_JpsiRaaVsNpart_pt%dGeV_Sys",i*5), TObject::kOverwrite);
	}
      grLhcRaaVsCent[0]->Write("ALICE_JpsiRaaVsNpart_LowPt", TObject::kOverwrite);
      grLhcRaaVsCentSys[0]->Write("ALICE_JpsiRaaVsNpart_LowPt_Sys", TObject::kOverwrite);
      grLhcRaaVsCent[1]->Write("CMS_JpsiRaaVsNpart_HighPt", TObject::kOverwrite);
      grLhcRaaVsCentSys[1]->Write("CMS_JpsiRaaVsNpart_HighPt_Sys", TObject::kOverwrite);

      gTsuRaaVsNpartRHIC[0]->Write("Tsinghua_RHIC_JpsiRaaVsNpart_LowPt", TObject::kOverwrite);
      gTsuRaaVsNpartRHIC[1]->Write("Tsinghua_RHIC_JpsiRaaVsNpart_HighPt", TObject::kOverwrite);
      gTsuRaaVsNpartLHC[0]->Write("Tsinghua_LHC_JpsiRaaVsNpart_LowPt", TObject::kOverwrite);
      gTsuRaaVsNpartLHC[1]->Write("Tsinghua_LHC_JpsiRaaVsNpart_HighPt", TObject::kOverwrite);
      gTamuRaaVsNpartRHIC[0]->Write("TAMU_RHIC_JpsiRaaVsNpart_LowPt", TObject::kOverwrite);
      gTamuRaaVsNpartRHIC[1]->Write("TAMU_RHIC_JpsiRaaVsNpart_HighPt", TObject::kOverwrite); 
      gTamuRaaVsNpartLHC[0]->Write("TAMU_LHC_JpsiRaaVsNpart_LowPt", TObject::kOverwrite);
      gTamuRaaVsNpartLHC[1]->Write("TAMU_LHC_JpsiRaaVsNpart_HighPt", TObject::kOverwrite);
    }
}

//================================================
void xsec(const bool savePlot = 0, const bool saveHisto = 0)
{
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  gStyle->SetOptStat(0);
  const double ncoll[nCentBins] = {291.9, 766.47, 290.87, 91.33, 21.57};
  const double ncollErr[nCentBins] = {20.46, 28.56, 30.47, 20.01, 8.04};
  const int marker_style[nCentBins] = {kFullCircle, kFullStar, kFullSquare, kFullCross, kFullDiamond};
  const double marker_size[nCentBins] = {1.5,2,1.5,2,2};
  const int marker_color[nCentBins] = {1,2,4,6,kGreen+2};
  const double scale_factor[nCentBins] = {10,1,0.2,0.1,0.01};
  const double x_max[nCentBins] = {15, 15, 12, 12, 7.5};

  TFile *fout = 0x0;
  if(saveHisto) fout = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","update");
  else fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type),"read");
  TGraphAsymmErrors *hJpsipp = (TGraphAsymmErrors*)fout->Get("hpp200JpsiVsPtFinalSys");


  TFile *fSys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type),"read");
  TH1F *hAuAuJpsiSys[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hAuAuJpsiSys[k]= (TH1F*)fSys->Get(Form("JpsiSysVsPt_All_cent%s",cent_Title[k]));
    }

  TFile *fdata = TFile::Open(Form("Rootfiles/%s.JpsiXsec.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiInvYield[nCentBins];
  TGraphAsymmErrors *hJpsiXsec[nCentBins];
  TGraphAsymmErrors *hJpsiXsecSys[nCentBins];
  TGraphAsymmErrors *hJpsiRaa[nCentBins];
  TGraphAsymmErrors *hJpsiRaaSys[nCentBins];
  TGraphAsymmErrors *hJpsiRaaSys2[nCentBins];
  TH1F *hJpsiPtPos = (TH1F*)fdata->Get("hJpsiPtPos_cent0080");
  double x, y;
  double x1, y1;
  const double x_err = 0.25;
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)fdata->Get(Form("Jpsi_InvYieldVsPt_cent%s",cent_Title[k]));
      int npoints = hJpsiInvYield[k]->GetNbinsX();
      hJpsiXsec[k] = new TGraphAsymmErrors(npoints);
      hJpsiXsec[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiXsecSys[k] = new TGraphAsymmErrors(npoints);
      hJpsiXsecSys[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));
      for(int i=0; i<npoints; i++)
	{
	  if( (k==2 || k==3) && i==npoints-1) continue;
	  if(k==4 && i>5) continue;
	  double pt = hJpsiPtPos->GetBinContent(i+1);
	  x = hJpsiInvYield[k]->GetBinCenter(i+1);
	  y = hJpsiInvYield[k]->GetBinContent(i+1)*x/pt;
	  double min_pt = hJpsiInvYield[k]->GetXaxis()->GetBinLowEdge(i+1);
	  double max_pt = hJpsiInvYield[k]->GetXaxis()->GetBinUpEdge(i+1);
	  double stat_rel = hJpsiInvYield[k]->GetBinError(i+1)/hJpsiInvYield[k]->GetBinContent(i+1);
	  hJpsiXsec[k]->SetPoint(i,pt,y*scale_factor[k]);
	  hJpsiXsec[k]->SetPointError(i,pt-min_pt,max_pt-pt,stat_rel*y*scale_factor[k],stat_rel*y*scale_factor[k]);
	  hJpsiXsecSys[k]->SetPoint(i,pt,y*scale_factor[k]);
	  hJpsiXsecSys[k]->SetPointError(i,x_err,x_err,hAuAuJpsiSys[k]->GetBinContent(i+1)*y*scale_factor[k],hAuAuJpsiSys[k]->GetBinContent(i+1)*y*scale_factor[k]);
	}

      hJpsiRaa[k] = new TGraphAsymmErrors(npoints);
      hJpsiRaa[k]->SetName(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
      hJpsiRaaSys[k] = new TGraphAsymmErrors(npoints);
      hJpsiRaaSys[k]->SetName(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
      hJpsiRaaSys2[k] = new TGraphAsymmErrors(npoints);
      hJpsiRaaSys2[k]->SetName(Form("Graph_Jpsi_Raa_cent%s_sys_pp",cent_Title[k]));
      for(int i=0; i<npoints; i++)
	{
	  if( (k==2 || k==3) && i==npoints-1) continue;
	  if(k==4 && i>5) continue;
	  double bin_center = hJpsiInvYield[k]->GetBinCenter(i+1);
	  hJpsiXsec[k]->GetPoint(i,x,y);
	  double AuAu_val = y/scale_factor[k]*x/bin_center;
	  double AuAu_err = hJpsiXsec[k]->GetErrorY(i)/y*AuAu_val;
	  double AuAu_sys = hAuAuJpsiSys[k]->GetBinContent(i+1) * AuAu_val;

	  hJpsipp->GetPoint(i, x1, y1);
	  double pp_val = y1;
	  double pp_err_h = hJpsipp->GetErrorYlow(i);
	  double pp_err_l = hJpsipp->GetErrorYhigh(i);
	  double prefix = ppInelastic/ncoll[k] * 1e6;
	  double val = prefix * AuAu_val / pp_val;
	  double err = prefix * AuAu_err / pp_val;
	  double sys = prefix * AuAu_sys / pp_val;
	  hJpsiRaa[k]->SetPoint(i,bin_center,val);
	  hJpsiRaa[k]->SetPointError(i,hJpsiInvYield[k]->GetBinWidth(i+1)/2,hJpsiInvYield[k]->GetBinWidth(i+1)/2,err,err);
	  //printf("%s: raa -> %2.2f +/- %2.2f%%, AuAu -> %2.2e +/- %2.2f%%, pp -> %2.2e +/-%2.2f%%\n",cent_Title[k],val,err/val*100,AuAu_val,AuAu_err/AuAu_val*100,pp_val,pp_err_h/pp_val*100);
	  //printf("%s: raa -> %2.2f +/- %2.2f%%, AuAu -> %2.2e +/- %2.2f%%, pp -> %2.2e +/-%2.2f%%\n",pt_Name[i+1],val,sys2/val*100,AuAu_val,AuAu_sys/AuAu_val*100,pp_val,pp_err/pp_val*100);

	  hJpsiRaaSys[k]->SetPoint(i,bin_center,val);
	  hJpsiRaaSys[k]->SetPointError(i,x_err*x_max[k]/x_max[0],x_err*x_max[k]/x_max[0],sys,sys);

	  hJpsiRaaSys2[k]->SetPoint(i,bin_center,val);
	  hJpsiRaaSys2[k]->SetPointError(i,x_err*x_max[k]/x_max[0],x_err*x_max[k]/x_max[0],pp_err_l/pp_val*val,pp_err_h/pp_val*val);
	}
      hJpsiRaa[k]->SetMarkerStyle(29);
      hJpsiRaa[k]->SetMarkerColor(2);
      hJpsiRaa[k]->SetLineColor(2);
      hJpsiRaa[k]->SetMarkerSize(2.5);
      hJpsiRaaSys[k]->SetFillStyle(0);
      hJpsiRaaSys[k]->SetLineColor(hJpsiRaa[k]->GetLineColor());
      hJpsiRaaSys2[k]->SetLineColor(kGray);
      hJpsiRaaSys2[k]->SetFillColor(kGray);
      hJpsiRaaSys2[k]->SetFillStyle(1001);
    }


  TCanvas *c1 = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",800,700);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} [GeV/c];B_{ll}d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{-2}]",15,0,15);
  hAuAu->GetYaxis()->SetRangeUser(1e-12,1e-2);
  ScaleHistoTitle(hAuAu,0.05,1,0.04,0.05,1.2,0.04,62);
  gPad->SetLogy();
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hAuAu->Draw();
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiXsec[k]->SetMarkerStyle(marker_style[k]);
      hJpsiXsec[k]->SetMarkerColor(marker_color[k]);
      hJpsiXsec[k]->SetLineColor(marker_color[k]);
      hJpsiXsec[k]->SetMarkerSize(marker_size[k]+0.5);
      hJpsiXsecSys[k]->SetFillStyle(0);
      hJpsiXsecSys[k]->SetLineColor(hJpsiXsec[k]->GetLineColor());
      hJpsiXsecSys[k]->Draw("samesE5");
      hJpsiXsec[k]->Draw("samesPEZ");
    }
  TLegend *leg = new TLegend(0.7,0.65,0.9,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hJpsiXsec[0],"0-80%#times10","P");
  leg->AddEntry(hJpsiXsec[1],"0-20%","P");
  leg->AddEntry(hJpsiXsec[2],"20-40%/5","P");
  leg->AddEntry(hJpsiXsec[3],"40-60%/10","P");
  leg->AddEntry(hJpsiXsec[4],"60-80%/100","P");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.15,0.3,0.85,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV, #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  t1->Draw();
  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.eps",run_type,run_cfg_name.Data()));
    }
  if(gSaveAN)
    {
      c1->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_InvYieldVsPt.pdf"));
    }
  if(gSavePaper)
    {
      c1->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_JpsiInvYieldVsPt.pdf",gPaperVersion));
    }

  //==============================================
  // RAA vs. pT
  //==============================================

  // CMS: Eur. Phys. J. C 77 (2017) 252
  // Prompt Jpsi, |y| < 2.4, 0-100%
  const int ncms = 5;
  double cms_pt[ncms] = {7.5, 9.0, 10.25, 12, 14.5};
  double cms_pt_err[ncms] = {1, 0.5, 0.75, 1, 1.5};
  double cms_pt_sys[ncms] = {0.25, 0.25, 0.25, 0.25, 0.25};
  double cms_raa[ncms] = {0.37, 0.337, 0.374, 0.371, 0.446};
  double cms_raa_err[ncms] = {0.0097, 0.0136, 0.0156, 0.0178, 0.0238};
  double cms_raa_sys[ncms] = {0.0407, 0.0282, 0.0288, 0.0356, 0.0406};
  TGraphErrors *gCmsRaaVsPt = new TGraphErrors(ncms, cms_pt, cms_raa, cms_pt_err, cms_raa_err);
  gCmsRaaVsPt->SetName("CMS_PromptJpsiRaaVsPt_0100");
  TGraphErrors *gCmsRaaVsPtSys = new TGraphErrors(ncms, cms_pt, cms_raa, cms_pt_sys, cms_raa_sys);
  gCmsRaaVsPt->SetName("CMS_PromptJpsiRaaVsPt_0100_sys");
  double cms_gsys = 0.0748;
  
  // ALICE JHEP 1507,051 (2015) (http://hepdata.cedar.ac.uk/view/ins1364887)
  // Inclusive Jpsi, |y| < 0.8, 0-40%
  const int nalice = 2;
  double alice_pt[nalice] = {1.56, 3.33};
  double alice_pt_err_low[nalice] = {1.56, 0.83};
  double alice_pt_err_high[nalice] = {0.94, 2.67};
  double alice_pt_sys_low[nalice] = {0.25, 0.25};
  double alice_pt_sys_high[nalice] = {0.25, 0.25};
  double alice_raa[nalice] = {0.82, 0.58};
  double alice_raa_err_low[nalice] = {0.11, 0.06};
  double alice_raa_err_high[nalice] = {0.11, 0.06};
  double alice_raa_sys_low[nalice] = {0.10,0.08};
  double alice_raa_sys_high[nalice] = {0.10,0.08};
  TGraphAsymmErrors *gAliceRaaVsPt = new TGraphAsymmErrors(nalice, alice_pt, alice_raa, alice_pt_err_low, alice_pt_err_high, alice_raa_err_low, alice_raa_err_high);
  gAliceRaaVsPt->SetName("ALICE_InclusiveJpsiRaaVsPt_040");
  TGraphAsymmErrors *gAliceRaaVsPtSys  = new TGraphAsymmErrors(nalice, alice_pt, alice_raa, alice_pt_sys_low, alice_pt_sys_high, alice_raa_sys_low, alice_raa_sys_high);
  gAliceRaaVsPtSys->SetName("ALICE_InclusiveJpsiRaaVsPt_040_sys");
  double alice_gsys = 0.12;
  
  // PHENIX 
  // four centrality bins: 0-20%, 20-40%, 40-60%, 60-92%
  const char* phenix_cent[4] = {"020","2040","4060","6092"};
  const int nphenix = 5;
  double phenix_pt[nphenix] = {0.5, 1.5, 2.5, 3.5, 4.5};
  double phenix_pt_err[nphenix] = {0.5, 0.5, 0.5, 0.5, 0.5};
  double phenix_pt_sys[nphenix] = {0.25, 0.25, 0.25, 0.25, 0.25};
  double phenix_raa[4][nphenix] = { {0.365, 0.379, 0.318, 0.134, 0.636},
				    {0.487, 0.554, 0.560, 0.649, 0.977},
				    {0.738, 0.566, 0.540, 0.694, 1.610},
				    {1.130, 0.586, 0.431, 0.818, 0.792}};
  double phenix_raa_err[4][nphenix] = { {0.057, 0.054, 0.067, 0.107, 0.339},
				     {0.073, 0.074, 0.101, 0.203, 0.486},
				     {0.111, 0.091, 0.125, 0.258, 0.848},
				     {0.237, 0.151, 0.183, 0.485, 0.833}};
  double phenix_raa_sys[4][nphenix] = { {0.035, 0.036, 0.031, 0.013, 0.061},
				     {0.072, 0.081, 0.082, 0.095, 0.144},
				     {0.107, 0.082, 0.078, 0.101, 0.234},
				     {0.163, 0.084, 0.062, 0.118, 0.114}};
  double phenix_gsys[4] = {0.1, 0.1, 0.13, 0.28};
  TGraphErrors *gPhenixRaaVsPt[4];
  TGraphErrors *gPhenixRaaVsPtSys[4];
  for(int i=0; i<4; i++)
    {
      gPhenixRaaVsPt[i] = new TGraphErrors(nphenix, phenix_pt, phenix_raa[i], phenix_pt_err, phenix_raa_err[i]);
      gPhenixRaaVsPt[i]->SetName(Form("PHENIX_InclusiveJpsiRaaVsPt_%s",phenix_cent[i]));
      gPhenixRaaVsPtSys[i] = new TGraphErrors(nphenix, phenix_pt, phenix_raa[i], phenix_pt_sys, phenix_raa_sys[i]);
      gPhenixRaaVsPtSys[i]->SetName(Form("PHENIX_InclusiveJpsiRaaVsPt_%s_sys",phenix_cent[i]));
      for(int ipoint=0; ipoint<nphenix; ipoint++)
	{
	  gPhenixRaaVsPtSys[i]->SetPointError(ipoint, phenix_pt_sys[ipoint]*x_max[i+1]/x_max[0], phenix_raa_sys[i][ipoint]);
	}
    }

  // STAR
  // High pt: Phys. Lett. B 722 (2013) 55
  // Low pt: Phys. Rev. C 90 (2014) 24906
  TFile *fpub = TFile::Open(Form("Rootfiles/%s/Publication.Jpsi.200GeV.root",run_cfg_name.Data()),"read");
  const char* star_cent[4] = {"0020","2040","4060","0060"};
  TGraphAsymmErrors *gRaaLowPt[4];
  TGraphAsymmErrors *gRaaLowPtSys[4];
  TGraphAsymmErrors *gRaaHighPt[4];
  TGraphAsymmErrors *gRaaHighPtSys[4];
  double x,y;
  for(int k=0; k<4; k++)
    {
      gRaaLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",star_cent[k]));
      gRaaLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",star_cent[k]));
      gRaaHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",star_cent[k]));
      gRaaHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",star_cent[k]));
    }

  //==============================================
  // compare all the previous Raa results from RHIC
  const int lowpt_color = 9;
  const int highpt_color = 12;
  const int phenix_color = kGreen + 2;
  TCanvas *c = new TCanvas("Raa_Jpsi_vs_pub","Raa_Jpsi_vs_pub",1100,700);
  c->Divide(2,2);

  TH1F *hRaa = new TH1F("Raa_Jpsi",";p_{T} (GeV/c);R_{AA}",10,0,15);
  hRaa->GetYaxis()->SetRangeUser(0.05,1.95);
  hRaa->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hRaa,22,1.9,18,22,1.9,18,63);
  for(int k=0; k<4; k++)
    {
      c->cd(k+1);
      SetPadMargin(gPad, 0.14, 0.13, 0.01, 0.01);
      hRaa->GetXaxis()->SetRangeUser(0, x_max[k+1]);
      hRaa->DrawCopy();
      TLine *line = GetLine(0,1,x_max[k+1],1,1);
      line->Draw();
      if(k<3)
	{
	  gRaaLowPt[k]->SetMarkerStyle(kFullCircle);
	  gRaaLowPt[k]->SetMarkerSize(2);
	  gRaaLowPt[k]->SetMarkerColor(lowpt_color);
	  gRaaLowPt[k]->SetLineColor(lowpt_color);
	  gRaaLowPtSys[k]->SetMarkerColor(lowpt_color);
	  gRaaLowPtSys[k]->SetLineColor(lowpt_color);

	  gRaaHighPt[k]->SetMarkerStyle(kOpenCircle);
	  gRaaHighPt[k]->SetMarkerSize(2);
	  gRaaHighPt[k]->SetMarkerColor(highpt_color);
	  gRaaHighPt[k]->SetLineColor(highpt_color);
	  gRaaHighPtSys[k]->SetMarkerColor(highpt_color);
	  gRaaHighPtSys[k]->SetLineColor(highpt_color);

	  gRaaLowPtSys[k]->Draw("sameE5");
	  gRaaHighPtSys[k]->Draw("sameE5");
	  gRaaLowPt[k]->Draw("sames PEZ");
	  gRaaHighPt[k]->Draw("sames PEZ");
	}
      gPhenixRaaVsPt[k]->SetMarkerStyle(kFullCross);
      gPhenixRaaVsPt[k]->SetMarkerSize(2.2);
      gPhenixRaaVsPt[k]->SetMarkerColor(phenix_color);
      gPhenixRaaVsPt[k]->SetLineColor(phenix_color);
      gPhenixRaaVsPtSys[k]->SetMarkerColor(phenix_color);
      gPhenixRaaVsPtSys[k]->SetLineColor(phenix_color);
      gPhenixRaaVsPtSys[k]->SetFillStyle(0);
      gPhenixRaaVsPtSys[k]->Draw("sameE5");
      gPhenixRaaVsPt[k]->Draw("sames PEZ");

      hJpsiRaaSys2[k+1]->Draw("sameE5");
      hJpsiRaaSys[k+1]->Draw("sameE5");
      hJpsiRaa[k+1]->Draw("samesPEZ");

      // centrality label
      TPaveText *t1 = GetPaveText(0.8,0.9,0.9,0.95);
      t1->SetTextFont(63);
      t1->SetTextSize(22);
      t1->AddText(Form("%s%%",cent_Name[k+1]));
      t1->Draw();
      if(k==3)
	{
	  t1 = GetPaveText(0.8,0.9,0.8,0.85);
	  t1->SetTextFont(63);
	  t1->SetTextSize(22);
	  t1->SetTextColor(gPhenixRaaVsPt[k]->GetMarkerColor());
	  t1->AddText("60-92%");
	  t1->Draw();
	}

      // Global systematics
      TBox *box_phenix = new TBox(x_max[k+1]-0.5,1-phenix_gsys[k],x_max[k+1]-0.7,1+phenix_gsys[k]);
      box_phenix->SetFillStyle(1001);
      box_phenix->SetFillColor(gPhenixRaaVsPt[k]->GetMarkerColor());
      box_phenix->Draw();

      double gerr = sqrt(pow(ncollErr[k]/ncoll[k],2)+pow(ppInelasticErr/ppInelastic,2));
      TBox *box_star = new TBox(x_max[k+1]-0.3,1-gerr,x_max[k+1]-0.5,1+gerr);
      box_star->SetFillStyle(1001);
      box_star->SetFillColor(kRed);
      box_star->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.2,0.65,0.45,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(16);
  leg->AddEntry(hJpsiRaa[0],"STAR: J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg->AddEntry(gRaaLowPt[0],"STAR: J/#psi#rightarrowe^{+}e^{-}, |y| < 1 (MB)","P");
  leg->AddEntry(gRaaHighPt[0],"STAR: J/#psi#rightarrowe^{+}e^{-}, |y| < 1 (HT)","P");
  leg->AddEntry(gPhenixRaaVsPt[0],"PHENIX: J/#psi#rightarrowe^{+}e^{-}, |y| < 0.35","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.eps",run_type,run_cfg_name.Data()));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_RaaVsPt_CompPub.pdf"));
    }


  //==============================================
  // final plot

  // model calculation
  ifstream ifin;
  char tmp[256];
  double dtmp;

  // Tsinghua group
  ifin.open("Rootfiles/Published/Jpsi_Raa_200/model/data_and_figures_from_Yunpeng/data_Raa_Pt.dat");
  double pt_thu[11];
  double v2_init[4][11];
  double v2_gen[4][11];
  double v2_all[4][11];
  ifin.getline(tmp,256);
  cout << tmp << endl;
  ifin.getline(tmp,256);
  cout << tmp << endl;
  for(int i=0; i<11; i++)
    {
      ifin >> pt_thu[i] 
	   >> v2_init[1][i] >> v2_gen[1][i] >> v2_all[1][i]
	   >> v2_init[2][i] >> v2_gen[2][i] >> v2_all[2][i]
	   >> v2_init[3][i] >> v2_gen[3][i] >> v2_all[3][i]
	   >> v2_init[0][i] >> v2_gen[0][i] >> v2_all[0][i];
    }
  TGraph *grTHU[nCentBins]; 
  for(int k=0; k<nCentBins; k++)
    {
      grTHU[k] = new TGraph(11,pt_thu,v2_all[k]);
      grTHU[k]->SetName(Form("THU_v2_all_cent%s",cent_Title[k]));
      grTHU[k]->SetLineColor(1);
      grTHU[k]->SetLineWidth(2);
      grTHU[k]->SetLineStyle(1);
    }
  ifin.close();

  // TAMU group
  TGraph *grTAMU[nCentBins];
  const char *fname[nCentBins] = {"060","020","2040","4060"};
  double pt_tamu[51];
  double v2_tamu[51];
  for(int k=0; k<nCentBins; k++)
    {
      ifin.open(Form("Rootfiles/Published/Jpsi_Raa_200/model/data_from_Xingbo/raa_pt_rhic_%s_tozebo_110507.dat",fname[k]));
      for(int i=0; i<51; i++)
	{
	  ifin >> pt_tamu[i] >> v2_tamu[i];
	}
      grTAMU[k] = new TGraph(51,pt_tamu,v2_tamu);
      grTAMU[k]->SetName(Form("TAMU_v2_all_cent%s",cent_Title[k]));
      grTAMU[k]->SetLineColor(1);
      grTAMU[k]->SetLineWidth(2);
      grTAMU[k]->SetLineStyle(7);
      ifin.close();
    }


  TCanvas *c = new TCanvas("Raa_Jpsi_dimuon","Raa_Jpsi_dimuon",1200,700);
  TBox *box_ncoll[nCentBins];
  TPaveText *centLabel[nCentBins];
  TPad *pads[nCentBins];
  const char* alphabet[nCentBins] = {"a", "b", "c", "d", "e"};
  for(int k=0; k<nCentBins; k++)
    {
      c->cd();
      if(k==0) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.05, 0.38, 0.51, 0.99);
      if(k==1) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.38, 0.68, 0.51, 0.99);
      if(k==2) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.68, 0.98, 0.51, 0.99);
      if(k==3) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.05, 0.38, 0.01, 0.49);
      if(k==4) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.38, 0.68, 0.01, 0.49);
      pads[k]->Draw();
      pads[k]->cd();
      SetPadMargin(gPad, 0.15, 0.13, 0.005, 0.005);

      TH1F *hRaa = new TH1F(Form("Raa_Jpsi_cent%d",k),";p_{T} [GeV/c];",150,0,15);
      hRaa->GetYaxis()->SetRangeUser(0.01,1.95);
      hRaa->GetYaxis()->CenterTitle();
      ScaleHistoTitle(hRaa,22,2.1,20,24,1.9,20,63);
      if(k!=0 && k!=3) 
	{
	  hRaa->GetYaxis()->SetLabelSize(0);
	  gPad->SetLeftMargin(0.02);
	}
      if(k==2 || k==3)  hRaa->GetXaxis()->SetRangeUser(0, 11);
      if(k==4)  hRaa->GetXaxis()->SetRangeUser(0, 8.1);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();

      if(k==0)
	{
	  gCmsRaaVsPt->SetMarkerStyle(21);
	  gCmsRaaVsPt->SetMarkerSize(1.8);
	  gCmsRaaVsPt->SetMarkerColor(4);
	  gCmsRaaVsPt->SetLineColor(4);
	  gCmsRaaVsPtSys->SetFillStyle(0);
	  gCmsRaaVsPtSys->SetLineColor(gCmsRaaVsPt->GetLineColor());
	  gCmsRaaVsPtSys->Draw("sameE5");
	  gCmsRaaVsPt->Draw("samesPEZ");

	  gAliceRaaVsPt->SetMarkerStyle(25);
	  gAliceRaaVsPt->SetMarkerSize(1.8);
	  gAliceRaaVsPt->SetMarkerColor(1);
	  gAliceRaaVsPt->SetLineColor(1);
	  gAliceRaaVsPtSys->SetFillStyle(0);
	  gAliceRaaVsPtSys->SetLineColor(gAliceRaaVsPt->GetLineColor());
	  gAliceRaaVsPtSys->Draw("sameE5");
	  gAliceRaaVsPt->Draw("samesPEZ");
	}

      if(k==1)
	{
	  gPhenixRaaVsPtSys[0]->Draw("sameE5");
	  gPhenixRaaVsPt[0]->Draw("sames PEZ");
	}

      hJpsiRaa[k]->SetMarkerStyle(29);
      hJpsiRaa[k]->SetMarkerColor(2);
      hJpsiRaa[k]->SetLineColor(2);
      hJpsiRaa[k]->SetMarkerSize(2.5);
      hJpsiRaaSys[k]->SetFillStyle(0);
      hJpsiRaaSys[k]->SetLineColor(hJpsiRaa[k]->GetLineColor());
      hJpsiRaaSys2[k]->SetLineColor(kGray);
      hJpsiRaaSys2[k]->SetFillColor(kGray);
      hJpsiRaaSys2[k]->SetFillStyle(1001);

      if(k>0 && k<nCentBins-1)
	{
	  grTAMU[k]->Draw("samesC");
	  grTHU[k]->Draw("samesC");
	}

      hJpsiRaaSys2[k]->Draw("sameE5");
      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      if(k==0) centLabel[k] = GetPaveText(0.17,0.22,0.9,0.95);
      else if(k==3) centLabel[k] = GetPaveText(0.27,0.32,0.9,0.95);
      else centLabel[k] = GetPaveText(0.19,0.24,0.9,0.95);
      centLabel[k]->SetTextFont(63);
      centLabel[k]->SetTextSize(22);
      if(k==0) centLabel[k]->AddText(Form("(%s)",alphabet[k]));
      else     centLabel[k]->AddText(Form("(%s) %s%%",alphabet[k],cent_Name[k]));
      centLabel[k]->Draw();

      // Global systematics
      double gerr = sqrt(pow(ncollErr[k]/ncoll[k],2)+pow(ppInelasticErr/ppInelastic,2));
      double x_pos = 14.25;
      if(k==2 || k==3) x_pos = 10.5;
      if(k==4) x_pos = 7.5;
      box_ncoll[k] = new TBox(x_pos,1-gerr,x_pos+0.25,1+gerr);
      box_ncoll[k]->SetFillStyle(1001);
      box_ncoll[k]->SetFillColor(kRed-7);
      box_ncoll[k]->Draw();

      if(k==0)
	{
	  TBox *box_alice = new TBox(x_pos-0.25,1-alice_gsys,x_pos,1+alice_gsys);
	  box_alice->SetFillStyle(1001);
	  box_alice->SetFillColor(gAliceRaaVsPt->GetLineColor());
	  box_alice->Draw();

	  TBox *box_cms = new TBox(x_pos-0.5,1-cms_gsys,x_pos-0.25,1+cms_gsys);
	  box_cms->SetFillStyle(1001);
	  box_cms->SetFillColor(gCmsRaaVsPt->GetLineColor());
	  box_cms->Draw();
	}

      if(k==1)
	{
	  TBox *box_phenix = new TBox(x_pos-0.25,1-phenix_gsys[0],x_pos,1+phenix_gsys[0]);
	  box_phenix->SetFillStyle(1001);
	  box_phenix->SetFillColor(gPhenixRaaVsPt[0]->GetMarkerColor());
	  box_phenix->Draw();
	}
    }
  c->cd();
  pads[0]->cd();
  TLegend *leg = new TLegend(0.25,0.7,0.5,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(18);
  leg->AddEntry(hJpsiRaa[0],"STAR: 0-80%, |y| < 0.5","P");
  leg->AddEntry(gCmsRaaVsPt,"CMS: 0-100%, |y| < 2.4","P");
  leg->AddEntry(gAliceRaaVsPt,"ALICE: 0-40%, |y| < 0.8","P");
  leg->Draw();

  pads[1]->cd();
  TLegend *leg = new TLegend(0.1,0.8,0.5,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(18);
  leg->AddEntry(gPhenixRaaVsPt[0],"PHENIX: 0-20%, |y| < 0.35","P");
  leg->Draw();

  c->cd();
  TPad *pad = GetSinglePad(Form("pad_cent6"), 0.68, 0.98, 0.01, 0.49);
  pad->Draw();
  pad->cd();
  TLegend *leg30 = new TLegend(0.25,0.45,0.5,0.87);
  leg30->SetBorderSize(0);
  leg30->SetFillColor(0);
  leg30->SetTextFont(63);
  leg30->SetTextSize(20);
  leg30->SetHeader("Au+Au @ 200 GeV");
  leg30->AddEntry(hJpsiRaa[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg30->AddEntry(hJpsiRaaSys[0],"Au+Au uncertainty","F");
  leg30->AddEntry(hJpsiRaaSys2[0],"p+p uncertainty","F");
  leg30->Draw();

  TLegend *leg32 = new TLegend(0.25,0.25,0.5,0.45);
  leg32->SetBorderSize(0);
  leg32->SetFillColor(0);
  leg32->SetTextFont(63);
  leg32->SetTextSize(20);
  leg32->AddEntry(grTHU[0],"TM I (Y.-P. Liu #it{et al.})","L");
  leg32->AddEntry(grTAMU[0],"TM II (X. Zhao #it{et al.})","L");
  leg32->Draw();

  c->cd();
  
  TPad *pad = GetSinglePad(Form("pad_new"), 0.01, 0.05, 0.51, 0.99);
  pad->Draw();
  pad->cd();
  TLatex latex;
  latex.SetTextFont(63);
  latex.SetTextSize(26);
  latex.SetTextAngle(90);  //align at top
  latex.DrawLatex(.85,.5,"R_{AA}");
  c->cd();
  TPad *pad = GetSinglePad(Form("pad_new"), 0.01, 0.05, 0.01, 0.49);
  pad->Draw();
  pad->cd();
  TLatex latex;
  latex.SetTextFont(63);
  latex.SetTextSize(26);
  latex.SetTextAngle(90);  //align at top
  latex.DrawLatex(.85,.5,"R_{AA}");

  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.eps",run_type,run_cfg_name.Data()));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_RaaVsPt.pdf"));
    }
  if(gSavePaper)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_JpsiRaaVsPt.pdf",gPaperVersion));
    }

  if(saveHisto)
    {
      fout->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiXsec[k]->Write("", TObject::kOverwrite);
	  hJpsiXsecSys[k]->Write("", TObject::kOverwrite);
	  hJpsiRaa[k]->Write("", TObject::kOverwrite);
	  hJpsiRaaSys[k]->Write("", TObject::kOverwrite);
	  hJpsiRaaSys2[k]->Write("", TObject::kOverwrite);
	}
    }
}

//================================================
void rawSignal(const int icent = 0, const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(kFALSE);
  TFile *fin = TFile::Open(Form("Rootfiles/Run14_AuAu200.Jpsi.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  const int nPtBins = 2;
  const double ptBins_low[nPtBins]  = {0,5};
  const double ptBins_high[nPtBins] = {15,15};
  const char *pt_Name[nPtBins] = {"0-15","5-15"};

  TH1F *hSeUL[nPtBins];
  TH1F *hSeLS[nPtBins];
  TH1F *hMeUL[nPtBins];
  TH1F *hMeLS[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hSeUL[i] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s_weight",pt_Name[i],cent_Name_pt[icent]));
      hSeLS[i] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s_weight",pt_Name[i],cent_Name_pt[icent]));
      hMeUL[i] = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]));
      hMeLS[i] = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]));
    }

 // get the Jpsi width from embedding
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"read");
  TH1F *hEmbJpsiWidth[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hEmbJpsiWidth[i] = (TH1F*)fscan->Get(Form("SmearEmb_JpsiWidthIntegr_Pt%s_def",pt_Name[i]));
    }

  // mixed event scaling
  TH1F *hMixScale[nPtBins];
  TF1 *funcScale[nPtBins];
  TCanvas *cScaling = new TCanvas("mix_scale","mix_scale",1100,500);
  cScaling->Divide(2,1);
  for(int i=0; i<nPtBins; i++)
    {
      hMixScale[i] = (TH1F*)hSeLS[i]->Clone(Form("%s_MixScale",hSeLS[i]->GetName()));
      hMixScale[i]->Rebin(2);
      TH1F *htmp = (TH1F*)hMeLS[i]->Clone(Form("%s_clone",hMeLS[i]->GetName()));
      htmp->Rebin(int(hMixScale[i]->GetBinWidth(1)/hMeLS[i]->GetBinWidth(1)));
      hMixScale[i]->Divide(htmp);

      // fitting method
      funcScale[i] = new TF1(Form("Fit_%s",hMixScale[i]->GetName()),"pol1",2.7,3.8);
      hMixScale[i]->Fit(funcScale[i],"IR0Q");
	      
      // plotting
      cScaling->cd(i+1);
      hMixScale[i]->SetTitle("");
      hMixScale[i]->SetMarkerStyle(21);
      hMixScale[i]->GetXaxis()->SetRangeUser(2.5,4);
      hMixScale[i]->SetMaximum(1.5*hMixScale[i]->GetMaximum());
      hMixScale[i]->SetMinimum(0.5*hMixScale[i]->GetMinimum());
      hMixScale[i]->Draw();
      funcScale[i]->SetLineColor(2);
      funcScale[i]->Draw("sames");
      TPaveText *text = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],cent_Name_pt[icent]),0.05);
      text->Draw();
    }
 
  // raw signal
  TGaxis::SetMaxDigits(3); 
  TH1F *hSignal[2];
  TH1F *hhSeUL, *hhSeLS, *hhMixBkg, *hhSignal;
  double g_bin_width[nPtBins] = {0.05, 0.05};
  for(int i=0; i<nPtBins; i++)
    {
      TCanvas *c = new TCanvas(Form("Jpsi_All_pt%s",pt_Name[i]),Form("Jpsi_All_pt%s",pt_Name[i]),800,650);
      SetPadMargin(gPad,0.12,0.12,0.05,0.05);
      hhSeUL = (TH1F*)hSeUL[i]->Clone(Form("hhSeUL_pt%s",pt_Name[i]));
      hhSeLS = (TH1F*)hSeLS[i]->Clone(Form("hhSeLS_pt%s",pt_Name[i]));

      hhMixBkg = (TH1F*)hMeUL[i]->Clone(Form("hhMixBkg_pt%s",pt_Name[i]));
      int low_bin = hhMixBkg->FindFixBin(2.3);
      int up_bin  = hhMixBkg->FindFixBin(4.2);
      for(int ibin=low_bin; ibin<=up_bin; ibin++)
	{
	  double mass = hhMixBkg->GetBinCenter(ibin);
	  double scale = funcScale[i]->Eval(mass);
	  hhMixBkg->SetBinContent(ibin, scale * hhMixBkg->GetBinContent(ibin));
	  hhMixBkg->SetBinError(ibin, scale * hhMixBkg->GetBinError(ibin));
	}
      hhSeUL->Rebin(g_bin_width[i]/hhSeUL->GetBinWidth(1));
      hhSeLS->Rebin(g_bin_width[i]/hhSeLS->GetBinWidth(1));
      hhMixBkg->Rebin(g_bin_width[i]/hhMixBkg->GetBinWidth(1));

      hhSignal = (TH1F*)hhSeUL->Clone(Form("hhSignal_pt%s",pt_Name[i]));
      hhSignal->Add(hhMixBkg, -1);
      hSignal[i] = (TH1F*)hhSignal->Clone(Form("hhSignal_pt%s_clone",pt_Name[i]));
      

      // plot    
      ScaleHistoTitle(hhSeUL,0.05,1,0.04,0.045,1.1,0.04,62);
      hhSeUL->GetXaxis()->SetRangeUser(2.6,3.8);
      TH1F *hplot = (TH1F*)hhSeUL->Clone(Form("hplot_%d",i));
      if(i==0) 
	{
	  hplot->Scale(0.15);
	  hhSeLS->Scale(0.15);
	  hhMixBkg->Scale(0.15);
	  hplot->SetMaximum(4*hhSignal->GetMaximum());
	  hplot->SetMinimum(-5e3);
	}
      if(i==1) 
	{
	  hplot->SetMaximum(1.45*hplot->GetMaximum());
	  hplot->SetMinimum(-1e2);
	}
      hplot->SetMarkerStyle(20);
      hplot->GetYaxis()->SetNdivisions(505);
      hplot->SetYTitle(Form("Counts/(%1.0f MeV/c^{2})",g_bin_width[i]*1e3));
      hplot->SetXTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
      hplot->Draw();
      hhSeLS->SetMarkerStyle(24);
      hhSeLS->Draw("samesP");
      hhMixBkg->SetLineColor(4);
      hhMixBkg->GetXaxis()->SetRangeUser(2.6,3.8);
      hhMixBkg->Draw("samesHIST");
      hhSignal->SetMarkerColor(2);
      hhSignal->SetLineColor(2);
      hhSignal->SetMarkerStyle(21);
      hhSignal->Draw("sames");

      // fitting
      TF1 *funcSignal = new TF1(Form("Fit_JpsiSignal_pt%s",pt_Name[i]),Form("gausn(0)+pol3(3)"),2.45,3.8);
      funcSignal->SetParameter(1,3.09);
      funcSignal->FixParameter(2,hEmbJpsiWidth[i]->GetBinContent(0));
      cout << hEmbJpsiWidth[i]->GetBinContent(0) << endl;
      TFitResultPtr ptr = hhSignal->Fit(funcSignal,"IR0");
      funcSignal->SetLineColor(1);
      funcSignal->SetLineStyle(1);
      funcSignal->Draw("same");

      TLegend *leg = new TLegend(0.15,0.3,0.3,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(63);
      leg->SetTextSize(18);
      if(i==0)
	{
	  leg->AddEntry(hplot,"Unlike-sign pairs (#times 0.15)","P");
	  leg->AddEntry(hhSeLS,"Like-sign pairs (#times 0.15)","PL");
	  leg->AddEntry(hhMixBkg,"Mixed-event (#times 0.15)","L");
	}
      else
	{
	  leg->AddEntry(hplot,"Unlike-sign pairs","P");
	  leg->AddEntry(hhSeLS,"Like-sign pairs","PL");
	  leg->AddEntry(hhMixBkg,"Mixed-event","L");
	}
      leg->Draw();

      double nJpsiAll = funcSignal->GetParameter(0)/hhSignal->GetBinWidth(1);
      double nJpsiAll_err = funcSignal->GetParError(0)/hhSignal->GetBinWidth(1);
      double min_mass = 3.0, max_mass = 3.2;
      double nAll = hSeUL[i]->Integral(hSeUL[i]->FindFixBin(min_mass),hSeUL[i]->FindFixBin(max_mass));
      double nJpsiMass = funcSignal->Integral(min_mass, max_mass)/funcSignal->Integral(2.9, 3.3)*nJpsiAll;


      TPaveText *t1 = GetPaveText(0.58,0.75,0.68,0.92,0.038,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("J/#psi #rightarrow #mu^{+} + #mu^{-}"));
      t1->AddText(Form("|y_{J/#psi}| < 0.5, p_{T,J/#psi} > %1.0f GeV/c",ptBins_low[i]));
      t1->AddText(Form("N_{J/#psi} = %2.0f #pm %2.0f",nJpsiAll, nJpsiAll_err));
      t1->AddText(Form("S/B = 1:%2.1f",(nAll-nJpsiMass)/nJpsiMass));
      t1->Draw();

      TPaveText *t1 = GetPaveText(0.15,0.2,0.82,0.92,0.038,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("Au+Au @ 200 GeV"));
      t1->AddText(Form("%s%%",cent_Name_pt[icent]));
      t1->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiSig_pt%s_cent%s.pdf",run_type,run_cfg_name.Data(),pt_Name[i],cent_Title_pt[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiSig_pt%s_cent%s.png",run_type,run_cfg_name.Data(),pt_Name[i],cent_Title_pt[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiSig_pt%s_cent%s.eps",run_type,run_cfg_name.Data(),pt_Name[i],cent_Title_pt[icent]));
	}
      if(gSavePaper)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_JpsiSig_pt%1.0f.pdf",gPaperVersion,ptBins_low[i]));
	}
    }
  
  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","update");
      for(int i=0; i<nPtBins; i++)
	{
	  hSeUL[i]->Write(Form("InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]), TObject::kOverwrite);
	  hSeLS[i]->Write(Form("InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]), TObject::kOverwrite);
	  hMeUL[i]->Write(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]), TObject::kOverwrite);
	  hMeLS[i]->Write(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]), TObject::kOverwrite);
	  funcScale[i]->Write(Form("funcScale_LS_MixToSame_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]), TObject::kOverwrite);
	  hSignal[i]->Write(Form("InvMass_Signal_pt%s_cent%s",pt_Name[i],cent_Name_pt[icent]), TObject::kOverwrite);
	}
    }
}
//================================================
void ppRef(const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  xbins[0] = 0.15;

  TFile *f = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.pp200GeV.root","read");

  // pp data
  TH1F *hpp = new TH1F("pp200_Jpsi",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hpp->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hpp,"",kTRUE);
  const int npp = 11;
  double xpp[npp] = {2.25, 2.75, 3.25, 3.75, 4.5, 5.5, 6.5, 7.5, 9, 11, 13};
  double exlpp[npp];
  double sxlpp[npp];
  for(int i=0; i<npp; i++) { exlpp[i] = 0; sxlpp[i] = 0.2; }
  double ypp[npp] = {0.68, 0.318, 0.187, 0.1032, 0.0334, 0.00905, 0.00154, 0.00084, 2.005e-4, 4.55e-5, 9.71e-6};
  double eylpp[npp] = {0.14, 0.050, 0.023, 0.0118, 0.0034, 0.00117, 0.00042, 0.00017, 2.42e-5, 7.2e-6, 2.41e-6};
  double sylpp[npp] = {0.06, 0.028, 0.018, 0.0089, 0.0028, 0.00077, 0.00013, 0.00025, 5.21e-5, 1.18e-5, 2.62e-6};
  double syhpp[npp] = {0.03, 0.007, 0.009, 0.0018, 0.0004, 0.00011, 0.00003, 0.00013, 1.6e-6, 4e-7, 1.6e-7};

  double xsec_5_pub = 0;
  for(int i=0; i<6; i++)
    {
      double pt = xpp[i+5];
      double dpt = 1;
      if(i>2) dpt = 2;
      xsec_5_pub += ypp[i+5] * 2 * pi * pt * dpt;
    }
  printf("[i] published pT > 5, xsec = %4.2e \n",xsec_5_pub);
  

  TGraphAsymmErrors *gHighPtPP = new TGraphAsymmErrors(npp, xpp, ypp, exlpp, exlpp, eylpp, eylpp);
  gHighPtPP->SetName("Jpsi_xsec_pp200_highPt");
  gHighPtPP->SetMarkerStyle(20);
  gHighPtPP->SetMarkerSize(1.5);
  gHighPtPP->Draw("sames PE");

  TGraphAsymmErrors *gHighPtPPSys = new TGraphAsymmErrors(npp, xpp, ypp, sxlpp, sxlpp, sylpp, syhpp);
  gHighPtPPSys->SetName("Jpsi_xsec_pp200_highPt_systematics");
  gHighPtPPSys->SetMarkerStyle(20);
  gHighPtPPSys->SetMarkerSize(0);
  gHighPtPPSys->SetFillStyle(0);
  gHighPtPPSys->Draw("sameE5");

  // PHENIX measurements
  TFile *fjpsi = TFile::Open(Form("Rootfiles/2016sQM/jpsi_xsec_pp200_run12.root",run_cfg_name.Data()),"read");
  TGraphAsymmErrors *gPhenixSys = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix_Systematics");
  gPhenixSys->SetMarkerSize(0);
  gPhenixSys->SetFillStyle(0);
  gPhenixSys->SetLineColor(4);
  gPhenixSys->Draw("sameE5");

  TGraphAsymmErrors *gPhenix = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix");
  gPhenix->SetMarkerStyle(24);
  gPhenix->SetMarkerSize(1.5);
  gPhenix->SetMarkerColor(4);
  gPhenix->SetLineColor(4);
  gPhenix->Draw("sames PE");

  leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  leg->AddEntry(gHighPtPP,"STAR 2009 HT |y|<1","PL");
  leg->AddEntry(gPhenix,"PHENIX |y|<0.35","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/ppRef_Compare.pdf",run_type,run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Comp.pdf"));
    }

  const int nResults = 2;
  TGraphAsymmErrors *hJpsiRebin[nResults];
  TGraphAsymmErrors *hJpsiRebinSys[nResults];
  const char *name[nResults] = {"STAR_2009","PHENIX"};
  const int markercolor[nResults] = {1, 4};
  const int markerstyle[nResults] = {20, 24};


  TF1 *funcJpsiXsec[nResults];
  TGraphAsymmErrors *graph = 0x0, *gSys = 0x0;
  double pT, y1, xh1, xl1, yh1, yl1, syh1, syl1;
  for(int i=0; i<nResults; i++)
    {
      hJpsiRebin[i] = new TGraphAsymmErrors(nbins);
      hJpsiRebin[i]->SetName(Form("Jpsi_xsec_%s",name[i]));
      hJpsiRebinSys[i] = new TGraphAsymmErrors(nbins);
      hJpsiRebinSys[i]->SetName(Form("Jpsi_xsec_%s_Sys",name[i]));
      if(i==0)
	{
	  graph = gHighPtPP;
	  gSys = gHighPtPPSys;
	}
      if(i==1)
	{
	  graph = gPhenix;
	  gSys = gPhenixSys;
	}
      cout << name[i] << endl;
      double fraction2 = 1;
      if(i==0)
	{
	  TCanvas *c = new TCanvas(Form("fit_xsec_%d",i),Form("fit_xsec_%s",name[i]),800,600);
	  hpp->DrawCopy("");
	  gPad->SetLogy();
	  graph->GetXaxis()->SetRangeUser(0,15);
	  graph->GetYaxis()->SetRangeUser(1e-6,10);
	  funcJpsiXsec[i] = new TF1(Form("Func_Jpsi_xsec_%s",name[i]),"exp([0]+[1]*x+[2]*x*x)",4,15);
	  graph->Fit(funcJpsiXsec[i], "IR0Q");
	  graph->Draw("samesPE");
	  funcJpsiXsec[i]->Draw("sames");
	  TF1 *funcTmp = new TF1(Form("%s_yield",funcJpsiXsec[i]->GetName()),Form("%s*x",funcJpsiXsec[i]->GetExpFormula().Data()));
	  funcTmp->SetParameters(funcJpsiXsec[i]->GetParameters());
	  fraction2 = funcTmp->Integral(10,14)/funcTmp->Integral(10,15);
	  TPaveText *t1 = GetTitleText(name[i],0.05);
	  t1->Draw();
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Fit_%s.pdf",run_type,run_cfg_name.Data(),name[i]));
	  if(gSaveAN)
	    {
	      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_FitSTAR.pdf"));
	    }
	  printf("[i] Fraction of [10,14]/[10,15] = %4.2f\n",fraction2);
	}
      double fraction1 = 1;
      if(i==1)
	{
	  TCanvas *c = new TCanvas(Form("fit_xsec_%d",i),Form("fit_xsec_%s",name[i]),800,600);
	  hpp->GetXaxis()->SetRangeUser(0,2);
	  hpp->GetYaxis()->SetRangeUser(0.1,7);
	  hpp->DrawCopy("");
	  graph->GetXaxis()->SetRangeUser(0,1);
	  graph->GetYaxis()->SetRangeUser(1e-1,10);
	  //funcJpsiXsec[i] = new TF1(Form("Func_Jpsi_xsec_%s",name[i]),"exp([0]+[1]*x+[2]*x*x)",0,2);
	  funcJpsiXsec[i] = new TF1(Form("Func_Jpsi_xsec_%s",name[i]),"[0]*(([1]-1)*([1]-2))/(2*3.14*[1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]))**(-1*[1])",0,2);
	  funcJpsiXsec[i]->SetParameters(1, 3, 1, 3.1);
	  graph->Fit(funcJpsiXsec[i], "IR0");
	  graph->Draw("samesPE");
	  funcJpsiXsec[i]->Draw("sames");
	  TF1 *funcTmp = new TF1(Form("%s_yield",funcJpsiXsec[i]->GetName()),Form("%s*x",funcJpsiXsec[i]->GetExpFormula().Data()));
	  funcTmp->SetParameters(funcJpsiXsec[i]->GetParameters());
	  fraction1 = funcTmp->Integral(0.15,0.25)/funcTmp->Integral(0,0.25);
	  TPaveText *t1 = GetTitleText(name[i],0.05);
	  t1->Draw();
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Fit_%s.pdf",run_type,run_cfg_name.Data(),name[i]));
	  if(gSaveAN)
	    {
	      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_FitPHENIX.pdf"));
	    }
	  printf("[i] Fraction of [0.15,0.25]/[0,0.25] = %4.2f\n",fraction1);
	}

      for(int bin=1; bin<=nbins; bin++)
	{
	  double x = 0, y = 0, xh = 0, xl = 0, yh = 0, yl = 0, syh = 0, syl = 0;
	  double dpT = 0;
	  x  = (xbins[bin]+xbins[bin-1])/2.;
	  xh = (xbins[bin]-xbins[bin-1])/2.;
	  xl = (xbins[bin]-xbins[bin-1])/2.;
	  for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
	    {
	      graph->GetPoint(ipoint, pT, y1);
	      if(bin==1)
		{
		  if(pT>xbins[bin]) continue;
		}
	      else
		{
		  if(pT<xbins[bin-1] || pT>xbins[bin]) continue;
		}
	      xh1 = graph->GetErrorXhigh(ipoint);
	      xl1 = graph->GetErrorXlow(ipoint);
	      yh1 = graph->GetErrorYhigh(ipoint);
	      yl1 = graph->GetErrorYlow(ipoint);
	      syh1 = gSys->GetErrorYhigh(ipoint);
	      syl1 = gSys->GetErrorYlow(ipoint);
	      if(i<1)
		{
		  // STAR 2009
		  if(pT<4) dpT = 0.5;
		  else if(pT<8) dpT = 1;
		  else dpT = 2;
		}
	      else
		{
		  // Phenix
		  if(pT<5) dpT = 0.25;
		  else     dpT = 1;
		}
	      double yield = y1 * dpT * pT;
	      if(bin==1 && ipoint==0) 
		{
		  yield = yield * fraction1;
		}
	      y += yield;
	      yh += pow(yh1/y1 * yield, 2);
	      yl += pow(yl1/y1 * yield, 2);
	      syh += syh1/y1 * yield;
	      syl += syl1/y1 * yield;
	    }
	  syh /= y;
	  syl /= y;
	  y = y / (xbins[bin]-xbins[bin-1]) / x;
	  yh = sqrt(yh) / (xbins[bin]-xbins[bin-1]) / x;
	  yl = sqrt(yl) / (xbins[bin]-xbins[bin-1]) / x;
	  syh *= y;
	  syl *= y;
	  if(y==0) continue;
	  if(i==1 && x>8) continue;
	  if(xbins[bin]>14) 
	    {
	      y/ = fraction2;
	      yl/ = fraction2;
	      yh/ = fraction2;
	      syl/ = fraction2;
	      syh/ = fraction2;
	    }
	  hJpsiRebin[i]->SetPoint(bin-1, x, y);
	  hJpsiRebin[i]->SetPointError(bin-1, 0, 0, yl, yh);
	  hJpsiRebinSys[i]->SetPoint(bin-1, x, y);
	  hJpsiRebinSys[i]->SetPointError(bin-1, xl, xh, syl, syh);
	  printf("[i] Rebin: pT = %4.3f, y = %4.3e, stat = (-%4.3f%%, +%4.3f%%), sys = (-%4.3f%%, +%4.3f%%)\n",x,y,yl/y*100,yh/y*100,syl/y*100,syh/y*100);
	}

      hJpsiRebinSys[i]->SetMarkerSize(0);
      hJpsiRebinSys[i]->SetFillStyle(0);
      hJpsiRebinSys[i]->SetLineColor(markercolor[i]);

      hJpsiRebin[i]->SetMarkerColor(markercolor[i]);
      hJpsiRebin[i]->SetMarkerStyle(markerstyle[i]);      
      hJpsiRebin[i]->SetMarkerSize(1.5);
      hJpsiRebin[i]->SetLineColor(markercolor[i]);
    }

  TH1F *hppRebin = new TH1F("pp200_Jpsi_Rebin",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hppRebin->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hppRebin,"",kTRUE);
  for(int i=0; i<nResults; i++)
    {
      hJpsiRebinSys[i]->Draw("sameE5");
      hJpsiRebin[i]->Draw("samePE");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/ppRef_Compare_rebin.pdf",run_type,run_cfg_name.Data()));

  // take the average
  TH1F *hPPJpsiFinal = new TH1F("hpp200JpsiVsPtFinal","",nbins,xbins);
  TGraphAsymmErrors *hPPJpsiFinalSys = new TGraphAsymmErrors(nbins);
  hPPJpsiFinalSys->SetName("hpp200JpsiVsPtFinalSys");

  double xsec_0 = 0, xsec_0_errh = 0, xsec_0_errl = 0;
  double xsec_5 = 0, xsec_5_errh = 0, xsec_5_errl = 0;
  printf("+++ Final avergae +++\n");
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      x = 0, y = 0, xh = 0, xl = 0, yh = 0, yl = 0, syh = 0, syl = 0;
      x  = (xbins[ipoint+1]+xbins[ipoint])/2.;
      xh = (xbins[ipoint+1]-xbins[ipoint])/2.;
      xl = (xbins[ipoint+1]-xbins[ipoint])/2.;

      double total_weight = 0;
      for(int i=0; i<nResults; i++)
	{
	  hJpsiRebin[i]->GetPoint(ipoint, pT, y1);
	  if(y1<=0) continue;
	  yh1 = hJpsiRebin[i]->GetErrorYhigh(ipoint);
	  yl1 = hJpsiRebin[i]->GetErrorYlow(ipoint);
	  syh1 = hJpsiRebinSys[i]->GetErrorYhigh(ipoint);
	  syl1 = hJpsiRebinSys[i]->GetErrorYlow(ipoint);
	  double weight = 1./(pow((yh1+yl1)/2,2) + pow((syh1+syl1)/2,2));
	  total_weight += weight;
	  y += y1 * weight;
	  yh += 1./pow(yh1,2);
	  yl += 1./pow(yl1,2);
	  syh += syh1 * weight;
	  syl += syl1 * weight;
	}
      y /= total_weight;
      yh = 1./sqrt(yh);
      yl = 1./sqrt(yl);
      syh /= total_weight;
      syl /= total_weight;
      yh = sqrt(yh*yh + syh*syh);
      yl = sqrt(yl*yl + syl*syl);
      hPPJpsiFinal->SetBinContent(ipoint+1, y);
      hPPJpsiFinal->SetBinError(ipoint+1, (yh+yl)/2);
      hPPJpsiFinalSys->SetPoint(ipoint, x, y);
      hPPJpsiFinalSys->SetPointError(ipoint, xl, xh, yl, yh);

      double yield = y * 2 * pi * x * (xbins[ipoint+1]-xbins[ipoint]);
      double yield_errh = yh/y * yield;
      double yield_errl = yl/y * yield;
      if(x>0) 
	{
	  xsec_0 += yield;
	  xsec_0_errh += pow(yield_errh, 2);
	  xsec_0_errl += pow(yield_errl, 2);
	}
      if(x>5) 
	{
	  xsec_5 += yield;
	  xsec_5_errh += pow(yield_errh, 2);
	  xsec_5_errl += pow(yield_errl, 2);
	}
      printf("[i] Rebin: pT = %4.3f, y = %4.3e, stat = (-%4.3f%%, +%4.3f%%), sys = (-%4.3f%%, +%4.3f%%)\n",x,y,yl/y*100,yh/y*100,syl/y*100,syh/y*100);
    }
  hPPJpsiFinalSys->SetMarkerStyle(21);
  hPPJpsiFinalSys->SetMarkerSize(1.5);
  hPPJpsiFinalSys->SetMarkerColor(kGreen+2);
  hPPJpsiFinalSys->SetLineColor(kGreen+2);
  hPPJpsiFinalSys->Draw("samesPEZ");

  leg = new TLegend(0.5,0.55,0.8,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hPPJpsiFinalSys,"Combined","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Final.pdf",run_type,run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Comb.pdf"));
    }

  xsec_0_errh = sqrt(xsec_0_errh); xsec_0_errl = sqrt(xsec_0_errl);
  xsec_5_errh = sqrt(xsec_5_errh); xsec_5_errl = sqrt(xsec_5_errl);

  printf("[i] pT > 0, xsec = %4.2f + %4.3f - %4.3f\n",xsec_0,xsec_0_errh,xsec_0_errl);
  printf("[i] pT > 5, xsec = %4.2e + %4.3f - %4.3f\n",xsec_5,xsec_5_errh,xsec_5_errl);
  TH1F *hJpsiIntXsec = new TH1F("hpp200JpsiVsCentFinal",";p_{T} (GeV/c);Bd#sigma/dy [nb]",2,0,2);
  hJpsiIntXsec->SetBinContent(1, xsec_0);
  hJpsiIntXsec->SetBinError(1, (xsec_0_errh+xsec_0_errl)/2);
  hJpsiIntXsec->SetBinContent(2, xsec_5);
  hJpsiIntXsec->SetBinError(2, (xsec_5_errh+xsec_5_errl)/2);

  TGraphAsymmErrors *hJpsiIntXsecSys = new TGraphAsymmErrors(2);
  hJpsiIntXsecSys->SetName("hpp200JpsiVsCentFinalSys");
  hJpsiIntXsecSys->SetPoint(0, 0, xsec_0);
  hJpsiIntXsecSys->SetPointError(0, 0, 0, xsec_0_errl, xsec_0_errh);
  hJpsiIntXsecSys->SetPoint(1, 5, xsec_5);
  hJpsiIntXsecSys->SetPointError(1, 0, 0, xsec_5_errl, xsec_5_errh);
  
  
  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","update");
      fout->cd();
      hJpsiIntXsec->Write("",TObject::kOverwrite);
      hJpsiIntXsecSys->Write("",TObject::kOverwrite);
      hPPJpsiFinal->Write("",TObject::kOverwrite);
      hPPJpsiFinalSys->Write("",TObject::kOverwrite);
    }
}


void scaleGraph(TGraphAsymmErrors *gr, double scale)
{
  int npoints = gr->GetN();
  Double_t x,y;
  for(int i=0; i<npoints; i++)
    {
      gr->GetPoint(i,x,y);
      double err_xlow = gr->GetErrorXlow(i);
      double err_xhigh = gr->GetErrorXhigh(i);
      double err_ylow = gr->GetErrorYlow(i);
      double err_yhigh = gr->GetErrorYhigh(i);
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,err_xlow,err_xhigh,err_ylow*scale,err_yhigh*scale);
    }
}

void drawBoxes(Int_t n, Double_t* x, Double_t dx, Double_t* y, Double_t* dy, Int_t lineWidth, Int_t lineColor){
  for(int i=0;i<n;i++){
    TBox *box = new TBox(x[i]-dx,y[i]-dy[i],x[i]+dx,y[i]+dy[i]);
    box->SetLineWidth(lineWidth);
    box->SetLineColor(lineColor);
    box->SetFillStyle(0);
    box->Draw("lsame");
  }
}

