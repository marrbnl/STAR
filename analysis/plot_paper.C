TString run_type = "Run14_AuAu200";
TString run_cfg_name = "paper";
const double luminosity = 14.2;
const double ppInelastic = 42.; // mb
//const double ppInelasticErr = 3; // mb
const double ppInelasticErr = 4.2; // 10% global error due to PHENIX
const char* model_name[4] = {"TAMU", "Tsinghua", "SHM", "Co-mover"};

//================================================
void plot_paper()
{  
  //rawSignal();
  //efficiency();
  //xsec();
  //xsec_v2();
  //nPart();
  //v2();
  ppRef();
  //ppRef2();
  //makeModel();
  //publication();
}

//================================================
void v2(const bool savePlot = 1, const bool saveHisto = 0)
{
  gStyle->SetOptStat(0);

  //A. Get Experimental Data

  //A-1. STAR Run10 Jpsi->ee Final (0-80%) 
  //PRL 111, 052301 (2013)
  double pt[4] = {1.01, 2.91, 4.98, 7.16};
  double pt_err[4] = {0, 0, 0, 0};
  double pt_sys[4] = {0.1, 0.1, 0.1, 0.1};
  double v2[4] = {0.05, 0.01, 0.05, 0.08};
  double v2_err[4] = {0.03, 0.04, 0.04, 0.05};
  double v2_sys[4] = {0.02, 0.01, 0.02, 0.01};
  double nonflow[4] = {0.01, 0.02, 0.04, 0.1};

  TGraphErrors *gr1 = new TGraphErrors(4, pt, v2, pt_err, v2_err);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(1);
  gr1->SetMarkerSize(1.5);
  gr1->SetLineColor(1);
  gr1->SetLineWidth(1.);

  TBox *box[4];
  for(int i=0; i<4; i++)
    {
      box[i] = new TBox(pt[i]-0.1, v2[i], pt[i]+0.1, v2[i]-nonflow[i]);
      box[i]->SetLineColor(kGray);
      box[i]->SetFillColor(kGray);
      box[i]->SetLineWidth(2.);
      box[i]->SetFillStyle(1);
    }
  double boxtshift = 4.7;
  double boxvshift = 0.198;

  TBox *boxLabel = new TBox(0.0+boxtshift, 0.00+boxvshift, 0.3+boxtshift, 0.02+boxvshift);
  boxLabel->SetLineColor(1);
  boxLabel->SetFillColor(1);
  boxLabel->SetLineWidth(1);
  boxLabel->SetFillStyle(0);

  TGraph* bracketUp[4];
  TGraph* bracketDown[4];
  for(int i=0; i<4; i++)
    {
      double bracketUpX[4] = {pt[i]-0.15, pt[i]-0.1, pt[i]+0.1, pt[i]+0.15};
      double v2Up = v2[i]+v2_sys[i];
      double bracketUpY[4] = {v2Up-0.004, v2Up, v2Up, v2Up-0.004};
      bracketUp[i] = new TGraph(4, bracketUpX, bracketUpY);
      bracketUp[i]->SetLineColor(1);
      bracketUp[i]->SetLineWidth(2.);
      double bracketDownX[4] = {pt[i]-0.15, pt[i]-0.1, pt[i]+0.1, pt[i]+0.15};
      double v2Down = v2[i]-v2_sys[i];
      double bracketDownY[4] = {v2Down+(float)0.004, v2Down, v2Down, v2Down+(float)0.004};
      bracketDown[i] = new TGraph(4, bracketDownX, bracketDownY);
      bracketDown[i]->SetLineColor(1);
      bracketDown[i]->SetLineWidth(2.);
    }

  TGraphErrors *gr1_sys = new TGraphErrors(4, pt, v2, pt_sys, v2_sys);
  gr1_sys->SetFillStyle(0);
  gr1_sys->SetMarkerSize(0);
  gr1_sys->SetMarkerColor(gr1->GetMarkerColor());
  gr1_sys->SetLineColor(gr1->GetMarkerColor());

  TGraphErrors *gNonFlow = new TGraphErrors(4);
  for(int i=0; i<4; i++)
    {
      gNonFlow->SetPoint(i, pt[i], nonflow[i]/2);
      gNonFlow->SetPointError(i, 0.1, nonflow[i]/2);
    }
  gNonFlow->SetFillStyle(1001);
  gNonFlow->SetFillColor(kGray);
  gNonFlow->SetLineColor(kGray);


  //A-2. STAR Run14 Jpsi->ee Prelim(Will be final soon)
  const int coloruu = 2;
  double pt_mtd[3] = {1.03, 3.2, 6.49};
  double ptError_mtd[3] = {0,0,0};
  double v2_mtd[3] = {-0.0316427, 0.0539167, 0.0878805};
  double v2StatErr_mtd[3] = {0.105994 , 0.0898011, 0.0901608};
  double v2SysErr_mtd[3] = {0.0411273, 0.0502348, 0.0168495};
  double ptSysErr_mtd[3] = {0.1,0.1,0.1};

  TGraphErrors *gr2 = new TGraphErrors(3, pt_mtd, v2_mtd, ptError_mtd, v2StatErr_mtd);
  gr2->SetMarkerStyle(29);
  gr2->SetMarkerColor(coloruu);
  gr2->SetMarkerSize(2.5);
  gr2->SetLineColor(coloruu);
  gr2->SetLineWidth(1);

  TGraphErrors *gr2_sys = new TGraphErrors(3, pt_mtd, v2_mtd, ptSysErr_mtd, v2SysErr_mtd);
  gr2_sys->SetMarkerColor(coloruu);
  gr2_sys->SetMarkerSize(0);
  gr2_sys->SetLineColor(coloruu);
  gr2_sys->SetFillStyle(0);
  gr2_sys->SetFillColor(coloruu);

  TGraph* bracketuuUp[4];
  TGraph* bracketuuDown[4];

  for(int i=0; i<3; i++)
    {
      double bracketUpX[4] = {pt_mtd[i]-0.15, pt_mtd[i]-0.1, pt_mtd[i]+0.1, pt_mtd[i]+0.15};
      double v2Up = v2_mtd[i]+v2SysErr_mtd[i];
      double bracketUpY[4] = {v2Up-0.004, v2Up, v2Up, v2Up-0.004};
      bracketuuUp[i] = new TGraph(4, bracketUpX, bracketUpY);
      bracketuuUp[i]->SetLineColor(coloruu);
      bracketuuUp[i]->SetLineWidth(2.);

      double bracketDownX[4] = {pt_mtd[i]-0.15, pt_mtd[i]-0.1, pt_mtd[i]+0.1, pt_mtd[i]+0.15};
      double v2Down = v2_mtd[i]-v2SysErr_mtd[i];
      double bracketDownY[4] = {v2Down+0.004, v2Down, v2Down, v2Down+0.004};
      bracketuuDown[i] = new TGraph(4, bracketDownX, bracketDownY);
      bracketuuDown[i]->SetLineColor(coloruu);
      bracketuuDown[i]->SetLineWidth(2.);
    }

  // B. Theory
  TFile *fmodel = TFile::Open("Rootfiles/Paper/models/AuAu200.models.root","read");
  TGraphErrors *gV2VsPtModel[2];
  gV2VsPtModel[0] = (TGraphErrors*)fmodel->Get("gV2VsPt_cent0080_TAMU");
  gV2VsPtModel[0]->SetFillStyle(3244);
  gV2VsPtModel[0]->SetFillColor(kGray+2);
  gV2VsPtModel[1] = (TGraphErrors*)fmodel->Get("gV2VsPt_cent0080_Tsinghua");
  gV2VsPtModel[1]->SetLineStyle(7);
  gV2VsPtModel[1]->SetLineWidth(2);
  gV2VsPtModel[1]->SetLineColor(kViolet);

  //C. Plot cosmetics
  TPaveText *t1 = GetPaveText(0.25,0.4,0.9,0.95,0.045,62);
  t1->AddText("Au+Au @ 200 GeV 0-80 %");

  TLegend* leg1 = new TLegend(0.15, 0.7, 0.3, 0.88);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.04);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr2,"J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg1->AddEntry(gr1,"J/#psi#rightarrowe^{+}e^{-}, |y| < 1 (PRL 111 (2013) 52301)","P");
  leg1->AddEntry(gNonFlow,"Non-Flow estimation","F");

  TLegend* legt = new TLegend(0.55,0.2,0.8,0.35);
  legt->SetTextFont(62);
  legt->SetTextSize(0.04);
  legt->SetFillStyle(0);
  legt->SetBorderSize(0);
  legt->AddEntry(gV2VsPtModel[1],"TM I: Tsinghua","L");
  legt->AddEntry(gV2VsPtModel[0],"TM II: TAMU","F");


  //D. Plot Jpsi v2 vs pt
  TH1F *hv2Draw = new TH1F("hv2Draw",";p_{T} (GeV/c);v_{2}               ",10,0,8);
  hv2Draw->GetYaxis()->SetRangeUser(-0.20-1e-4,0.35+1e-4);
  hv2Draw->SetLineColor(1);
  ScaleHistoTitle(hv2Draw,28,1,24,28,0.8,24,63);

  TCanvas *c1 = new TCanvas("v2_combined","v2_combined",700,500);
  c1->cd();
  SetPadMargin(gPad,0.13,0.1,0.05,0.02);
  hv2Draw->Draw();
  gNonFlow->Draw("samesE3");

  // theory 
  gV2VsPtModel[0]->Draw("samesE4");
  gV2VsPtModel[1]->Draw("samesL");
  
  // STAR data

  gr1_sys->Draw("samesE5");
  gr1->Draw("psamez");
  for(int i=0; i<4; i++)
    {
      //box[i]->Draw("fsame");
    }

  for(int i=0; i<4; i++)
    {
      //bracketUp[i]->Draw("lsame");
      //bracketDown[i]->Draw("lsame");
    }


  for(int i=0; i<3; i++)
    {
      //bracketuuUp[i]->Draw("lsame");
      //bracketuuDown[i]->Draw("lsame");
    }
  gr2_sys->Draw("samesE5");
  gr2->Draw("samesPEz");


  // legend etc.
  //boxLabel->Draw("fsame");
  t1->Draw("same");
  leg1->Draw("same");
  legt->Draw("same");

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14_model.pdf",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14_model.png",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14_model.jpg",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14_model.eps",run_type.Data(),run_cfg_name.Data()));
    }

  if(gSavePaper)
    {
      c1->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_JpsiV2.pdf",gPaperVersion));
    }


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
  gV2VsPt[1] = new TGraphErrors(14);
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

  // nsigmaPi, dy and dz efficiency
  TFile *femb = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root","read");
  TH1F *hTrkPt[4];
  hTrkPt[0] = (TH1F*)femb->Get("mhTrkPtDis_Tpc_di_mu");
  hTrkPt[1] = (TH1F*)femb->Get("mhTrkPtDis_NsigmaPi_di_mu");
  hTrkPt[2] = (TH1F*)femb->Get("mhTrkPtDis_MtdMth_di_mu");
  hTrkPt[3] = (TH1F*)femb->Get("mhTrkPtDis_DyDz_di_mu");
  for(int i=0; i<4; i++)
    {
      hTrkPt[i]->Sumw2();
    }
  hTrkPt[1]->Divide(hTrkPt[0]);
  hTrkPt[3]->Divide(hTrkPt[2]);
  
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
  leg->SetTextSize(20);;
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

  fucnDtof->Draw("sames"); // dtof efficiency
  hTrkPt[1]->SetMarkerStyle(20);
  hTrkPt[1]->SetMarkerColor(4);
  hTrkPt[1]->SetLineColor(4);
  hTrkPt[1]->Draw("samesPE"); // nSigmaPi efficiency
  hTrkPt[3]->SetMarkerStyle(25);
  hTrkPt[3]->SetMarkerColor(6);
  hTrkPt[3]->SetLineColor(6);
  hTrkPt[3]->Draw("samesPE"); // dy,dz efficiency
  for(int i=0; i<2; i++)
    {
      gPidEffSys[i]->SetFillStyle(1001);
      gPidEffSys[i]->SetFillColor(kGray);
      gPidEffSys[i]->SetLineColor(kBlack);
      gPidEffSys[i]->SetLineWidth(3);
      gPidEffSys[i]->Draw("samesE2L");
    }
  TLegend *leg = new TLegend(0.45,0.22,0.7,0.55);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);;
  leg->AddEntry(hTrkPt[1],"n#sigma_{#pi} cut","P");
  leg->AddEntry(hTrkPt[3],"#Deltay & #Deltaz cut","P");
  leg->AddEntry(fucnDtof,"#DeltaT_{tof} cut (after other PID cuts)","L");
  leg->AddEntry(gPidEffSys[0],"Total muon PID efficiency","FL");
  leg->Draw();
  cEff->cd();
  cEff->Update();

  if(savePlot)
    {
      cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/MtdEfficiency.pdf",run_type.Data(),run_cfg_name.Data()));      
      cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/MtdEfficiency.png",run_type.Data(),run_cfg_name.Data()));   
      cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/MtdEfficiency.eps",run_type.Data(),run_cfg_name.Data()));
    }
  if(gSavePaper)
    {
      cEff->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_MtdEfficiency.eps",gPaperVersion));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type.Data()),"update");
      gTrigEffSys->Write("MTD_TrigEff", TObject::kOverwrite);
      gPidEffSys[0]->Write("MTD_PidEff_LowPt", TObject::kOverwrite);
      gPidEffSys[1]->Write("MTD_PidEff_HighPt", TObject::kOverwrite);
    }
    
}

//================================================
void nPart(const bool savePlot = 1, const bool saveHisto = 1)
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
  const double ncollErr2[7] = {8.04, 16.85, 23.36, 29.20, 32.33, 30.18, 26.27};
  const double npart2[7] =    {21.04, 47.79, 76.88, 116.74, 169.04, 237.27, 325.48}; 
  
  double npart_y[kNCent];
  double ncoll_relerr[kNCent];
  for(int k=0; k<kNCent; k++)
    {
      npart_y[k] = 1;
      ncoll_relerr[k] = ncollErr[k]/ncoll[k];
    }

  TGraphErrors *gPubNpartErr = new TGraphErrors(kNCent, npart, npart_y, npartErr, ncoll_relerr);
  gPubNpartErr->SetFillColor(kBlue-10);
  gPubNpartErr->SetLineColor(kBlue-10);

  // plotting template
  TH1F *hRaaVsNpart = new TH1F("hRaaVsNpart",";N_{part};J/#psi R_{AA}",100,-100,370);
  ScaleHistoTitle(hRaaVsNpart,28,1,24,28,1,24,63);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,2);
  hRaaVsNpart->GetYaxis()->CenterTitle();

  // MTD results
  TFile *fout = 0x0;
  if(saveHisto) fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type.Data()),"update");
  else fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type.Data()),"read");
  TGraphAsymmErrors *hpp = (TGraphAsymmErrors*)fout->Get("hpp200JpsiVsCentFinalSys");

  TFile *fdata = TFile::Open(Form("Rootfiles/%s.JpsiXsec.pt%1.1f.pt%1.1f.root",run_type.Data(),pt1_cut,pt2_cut),"read");
  TFile *fsys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type.Data()), "read");
  TH1F *hdata[nPtBins], *hsys[nPtBins];
  TH1F *hRaa[nPtBins];
  TGraphErrors *raaVsNpart[nPtBins], *raaVsNpartSys[nPtBins];
  TBox *globalSys[nPtBins];
  double x,y;
  for(int i=0; i<nPtBins; i++)
    {
      printf("+++ pt %s +++\n",pt_Name[i]);
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
	  cout << NCOLL << ": " << auau_yield << "  " << pp_yield << "  " << raa << "  " << raa_err << "  " << raa_sys << endl;
	}
      double gSys_h = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err_h*pp_err_h/pp_yield/pp_yield);
      double gSys_l = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err_l*pp_err_l/pp_yield/pp_yield);
      globalSys[i] = new TBox(365,1-gSys_l,370,1+gSys_h);
      globalSys[i]->SetLineWidth(2.);
      globalSys[i]->SetFillStyle(1001);
      cout << "Global err = " << gSys_h << "  " << gSys_l << endl;
    }

  //==============================================
  // LHC results
  //==============================================

  // ALICE Low pT
  // http://hepdata.cedar.ac.uk/view/ins1263062
  // Phys. Lett. B 734 (2014) 314-327

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

  /*
  // CMS high pT
  // JHEP05(2012)063
  const int nCms_highpT = 6;
  double cms_npart[nCms_highpT] = {355.4, 261.4, 187.2, 130.0, 86.3, 22.1};
  double cms_npart_err[nCms_highpT] = {0.,0.,0.,0.,0.,0.};
  double cms_npart_sys[nCms_highpT] = {6,6,6,6,6,6};
  double cms_raa[nCms_highpT]     = {0.24, 0.26, 0.31, 0.50, 0.70, 0.62};
  double cms_raa_err[nCms_highpT] = {0.03, 0.03, 0.04, 0.07, 0.11, 0.11};
  double cms_raa_sys[nCms_highpT] = {0.02, 0.02, 0.02, 0.05, 0.08, 0.10};
  grLhcRaaVsCent[1] = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_err, cms_raa_err);
  grLhcRaaVsCentSys[1] = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_sys, cms_raa_sys);
  */

  // CMS high pT
  // EPJC 77 (2017) 252
  // prompt Jpsi
  // global uncertainty = 7.7%
  // centrality bins: 0-5%, 5-10%, 10-15% ... 45-50%, 50-60%, 60-100%
  const int nCms_highpT = 12;
  double cms_npart[nCms_highpT] = {381.4, 329.8, 283.2, 240.9, 203.6, 171.3, 142.7, 117.6, 96.1, 76.7, 53.8, 14.3};
  double cms_npart_err[nCms_highpT] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double cms_npart_sys[nCms_highpT] = {6,6,6,6,6,6, 6,6,6,6,6,6};
  double cms_raa[nCms_highpT]     = {0.28, 0.30, 0.31, 0.36, 0.42, 0.46, 0.49, 0.47, 0.55, 0.56, 0.72, 0.72};
  double cms_raa_err[nCms_highpT] = {0.010, 0.011, 0.013, 0.014, 0.014, 0.014, 0.029, 0.029, 0.029, 0.043, 0.043, 0.043};
  double cms_raa_sys[nCms_highpT] = {0.023, 0.024, 0.033, 0.029, 0.043, 0.057, 0.057, 0.057, 0.072, 0.072, 0.115, 0.129};
  grLhcRaaVsCent[1] = new TGraphErrors(nCms_highpT, cms_npart, cms_raa, cms_npart_err, cms_raa_err);
  grLhcRaaVsCentSys[1] = new TGraphErrors(nCms_highpT, cms_npart, cms_raa, cms_npart_sys, cms_raa_sys);

  //==============================================
  // theory
  //==============================================
  TFile *fmodel = TFile::Open("Rootfiles/Paper/models/AuAu200.models.root","read");
  // tsinghua group
  TGraphErrors *gRaaVsCentModel[3][2];
  for(int i=0; i<nPtBins; i++)
    {
      for(int j=0; j<3; j++)
	{
	  if(i==1 && j==2) continue;
	  gRaaVsCentModel[j][i] = (TGraphErrors*)fmodel->Get(Form("gRaaVsNpart_Pt%1.0f_%s",ptBins_low[i],model_name[j]));
	}
    }
  
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


      gRaaVsCentModel[0][i]->Draw("samesE4");
      gRaaVsCentModel[0][i]->SetFillStyle(3444);
      gRaaVsCentModel[0][i]->SetFillColor(kGray+2);
      gRaaVsCentModel[1][i]->SetLineStyle(1);
      gRaaVsCentModel[1][i]->SetLineWidth(2);
      gRaaVsCentModel[1][i]->SetLineColor(kViolet);
      gRaaVsCentModel[1][i]->Draw("samesL");
      if(i==0)
	{
	  gRaaVsCentModel[2][i]->SetLineStyle(7);
	  gRaaVsCentModel[2][i]->SetLineColor(kGreen+2);
	  gRaaVsCentModel[2][i]->SetLineWidth(2);
	  gRaaVsCentModel[2][i]->Draw("samesL");
	}

      // LHC
      grLhcRaaVsCent[i]->SetMarkerStyle(21);
      grLhcRaaVsCent[i]->SetMarkerSize(1.5);
      grLhcRaaVsCent[i]->SetLineColor(marker_color[1]);
      grLhcRaaVsCent[i]->SetMarkerColor(marker_color[1]);
      grLhcRaaVsCentSys[i]->SetMarkerColor(marker_color[1]);
      grLhcRaaVsCentSys[i]->SetLineColor(marker_color[1]);
      grLhcRaaVsCentSys[i]->SetFillStyle(0);
      TBox *bLhc;
      if(i==0) bLhc = new TBox(360,1-0.13,365,1+0.13);
      else     bLhc = new TBox(360,1-0.077,365,1+0.077);
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


      TLegend *leg41 = new TLegend(0.3,0.67,0.5,0.8);
      leg41->SetBorderSize(0);
      leg41->SetFillColor(0);
      leg41->SetTextFont(62);
      leg41->SetTextSize(0.035);
      leg41->SetHeader(Form("p_{T} > %1.0f GeV/c",ptBins_low[i]));
      leg41->AddEntry(gRaaVsCentModel[1][i],"TM I: Tsinghua","L");
      leg41->AddEntry(gRaaVsCentModel[0][i],"TM II: TAMU","F");
      leg41->Draw();

      if(i==0)
	{
	  TLegend *leg42 = new TLegend(0.55,0.71,0.8,0.76);
	  leg42->SetBorderSize(0);
	  leg42->SetFillColor(0);
	  leg42->SetTextFont(62);
	  leg42->SetTextSize(0.035);
	  leg42->AddEntry(gRaaVsCentModel[2][i],"SHM","L");
	  leg42->Draw();
	}

      TLegend *leg_ncoll_2 = new TLegend(0.15,0.17,0.35,0.22);
      leg_ncoll_2->SetBorderSize(0);
      leg_ncoll_2->SetFillColor(0);
      leg_ncoll_2->SetTextFont(62);
      leg_ncoll_2->SetTextSize(0.035);;
      leg_ncoll_2->AddEntry(gPubNpartErr,"N_{coll} uncertainty","f");
      leg_ncoll_2->Draw();

      TLegend *leg3 = new TLegend(0.16,0.84,0.4,0.96);
      leg3->SetBorderSize(0);
      leg3->SetFillColor(0);
      leg3->SetTextFont(62);
      leg3->SetTextSize(0.04);
      if(i==0)
	{
	  leg3->AddEntry(grStar,"Inclusive: Au+Au @ 200 GeV, |y| < 0.5, p_{T} > 0.15 GeV/c","P");
	  leg3->AddEntry(grLhcRaaVsCent[i],"Inclusive: Pb+Pb @ 2.76 TeV, |y| < 0.8, p_{T} > 0 GeV/c","P");
	}
      else
	{
	  leg3->AddEntry(grStar,"Inclusive: Au+Au @ 200 GeV, |y| < 0.5, p_{T} > 5 GeV/c","P");
	  leg3->AddEntry(grLhcRaaVsCent[i],"Prompt: Pb+Pb @ 2.76 TeV, |y| < 2.4, p_{T} > 6.5 GeV/c","P");
	}
      leg3->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiRaaVsNpart_Pt%s_RHICvsLHC.pdf",run_type.Data(),run_cfg_name.Data(),pt_Name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiRaaVsNpart_Pt%s_RHICvsLHC.png",run_type.Data(),run_cfg_name.Data(),pt_Name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiRaaVsNpart_Pt%s_RHICvsLHC.eps",run_type.Data(),run_cfg_name.Data(),pt_Name[i]));
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

  // calculate chi2/NDF bewteen data and theory
  double chi2[nPtBins][3][3];
  double ndf[nPtBins][3][3];
  double pvalue[nPtBins][3][3];
  for(int i=0; i<nPtBins; i++)
    {
      for(int m=0; m<3; m++)
	{
	  if(i==1 && m==2) continue;
	  printf("+++ pt%s %s +++ \n",pt_Name[i],model_name[m]);
	  for(int j=0; j<3; j++)
	    {
	      chi2[i][m][j] = 0;
	      ndf[i][m][j] = 0;
	      pvalue[i][m][j] = 0;
	      for(int bin=1; bin<=raaVsNpart[i]->GetN(); bin++)
		{
		  raaVsNpart[i]->GetPoint(bin-1, x, y);
		  double pt = x;
		  double raa_data = y;
		  double err_data = 0;
		  if(i==0) err_data = sqrt( pow(raaVsNpart[i]->GetErrorY(bin-1), 2) + 
					    pow(raaVsNpartSys[i]->GetErrorY(bin-1), 2) +
					    pow(ncollErr[bin-1]/ncoll[bin-1],2));
		  if(i==1) err_data = sqrt( pow(raaVsNpart[i]->GetErrorY(bin-1), 2) + 
					    pow(raaVsNpartSys[i]->GetErrorY(bin-1), 2) +
					    pow(ncollErr2[bin-1]/ncoll2[bin-1],2));
		  if(j==1) 
		    {
		      raa_data = y * globalSys[i]->GetY2();
		      err_data = err_data/y*raa_data;
		    }
		  if(j==2) 
		    {
		      raa_data = y * globalSys[i]->GetY1();
		      err_data = err_data/y*raa_data;
		    }
	
		  double *model_x = gRaaVsCentModel[m][i]->GetX();
		  double *model_y = gRaaVsCentModel[m][i]->GetY();
		  double *model_ey = gRaaVsCentModel[m][i]->GetEY();
		  double raa_model = -1, err_model = -1;
		  if(m<2 && pt<model_x[gRaaVsCentModel[m][i]->GetN()-1]) continue;
		  if(m==2 && pt>model_x[gRaaVsCentModel[m][i]->GetN()-1]) continue;
		  ndf[i][m][j]++;
		  for(int ipoint=0; ipoint<gRaaVsCentModel[m][i]->GetN(); ipoint++)
		    {
		      if(m<2 && pt >= model_x[ipoint+1] && pt < model_x[ipoint] )
			{
			  raa_model = model_y[ipoint+1];
			  err_model = model_ey[ipoint+1];
			  break;
			}
		      if(m==2 && pt < model_x[ipoint+1] && pt >= model_x[ipoint] )
			{
			  raa_model = model_y[ipoint];
			  err_model = model_ey[ipoint];
			  break;
			}
		    }
		  if(m==2 && bin==1)
		    {
		      raa_model = 1;
		      err_model = 0;
		    }
		  double err_tot = sqrt(pow(err_data,2) + pow(err_model,2)); 
		  chi2[i][m][j] += pow((raa_data-raa_model)/err_tot,2);
		  //printf("pt = %4.2f, data = %4.2f #pm %4.2f, model = %4.2f #pm %4.2f\n",pt,raa_data,err_data,raa_model,err_model);
		}
	      pvalue[i][m][j] = TMath::Prob(chi2[i][m][j],ndf[i][m][j]);
	    }
	  for(int j=0; j<3; j++)
	    {
	      if(pvalue[i][m][j]>0.1)
		printf("%4.2f/%d & %4.2f & ",chi2[i][m][j],ndf[i][m][j],pvalue[i][m][j]);
	      else if(pvalue[i][m][j]>0.01)
		printf("%4.2f/%d & %4.3f & ",chi2[i][m][j],ndf[i][m][j],pvalue[i][m][j]);
	      else
		printf("%4.2f/%d & %4.2e &",chi2[i][m][j],ndf[i][m][j],pvalue[i][m][j]);
	    }
	  printf("\n");
	}
    }
  cout << "here" << endl;

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
    }
}

//================================================
void xsec(const int mode = 1, const bool savePlot = 0, const bool saveHisto = 0)
{
  // mode 0: get histogram from scratch
  // mode 1: get histogram from saved 

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
  else fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type.Data()),"read");
  TGraphAsymmErrors *hJpsipp = (TGraphAsymmErrors*)fout->Get("hpp200JpsiVsPtFinalSys");

  TGraphAsymmErrors *hJpsiXsec[nCentBins];
  TGraphAsymmErrors *hJpsiXsecSys[nCentBins];
  TGraphAsymmErrors *hJpsiRaa[nCentBins];
  TGraphAsymmErrors *hJpsiRaaSys[nCentBins];
  TGraphAsymmErrors *hJpsiRaaSys2[nCentBins];

  double x, y;
  double x1, y1;
  const double x_err = 0.25;
  if(mode==0)
    {
      TFile *fSys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type.Data()),"read");
      TH1F *hAuAuJpsiSys[nCentBins];
      for(int k=0; k<nCentBins; k++)
	{
	  hAuAuJpsiSys[k]= (TH1F*)fSys->Get(Form("JpsiSysVsPt_All_cent%s",cent_Title[k]));
	}

      TFile *fdata = TFile::Open(Form("Rootfiles/%s.JpsiXsec.pt%1.1f.pt%1.1f.root",run_type.Data(),pt1_cut,pt2_cut),"read");
      TH1F *hJpsiInvYield[nCentBins];

      TH1F *hJpsiPtPos = (TH1F*)fdata->Get("hJpsiPtPos_cent0080");
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
    }
  else
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiXsec[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
	  hJpsiXsecSys[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));

	  hJpsiRaa[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
	  hJpsiRaaSys[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
	  hJpsiRaaSys2[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_Raa_cent%s_sys_pp",cent_Title[k]));
	}
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
      if(k==1)
	{
	  hJpsiXsec[k]->GetPoint(0, x, y);
	  printf("[i] Jpsi yield = %4.2e +/- %4.2e +/- %4.2e\n",y,hJpsiXsec[k]->GetErrorYhigh(0),hJpsiXsecSys[k]->GetErrorYhigh(0));
	}
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
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.pdf",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.png",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.eps",run_type.Data(),run_cfg_name.Data()));
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

  /*
    // CMS: Eur. Phys. J. C 05 (2012) 063
    // Inclusive Jpsi, |y| < 2.4, 0-100%
    // global uncertainty = 8.3%
    const int ncms = 2;
    double cms_pt[ncms] = {8.11, 13.22};
    double cms_pt_err_low[ncms] = {1.61, 3.22};
    double cms_pt_err_high[ncms] = {1.89, 16.78};
    double cms_pt_sys_low[ncms] = {0.25, 0.25};
    double cms_pt_sys_high[ncms] = {0.25, 0.25};
    double cms_raa[ncms] = {0.32, 0.31};
    double cms_raa_err_low[ncms] = {0.03, 0.04};
    double cms_raa_err_high[ncms] = {0.03, 0.04};
    double cms_raa_sys_low[ncms] = {0.02, 0.01};
    double cms_raa_sys_high[ncms] = {0.02, 0.01};

   */
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
      gPhenixRaaVsPt[k]->SetMarkerStyle(kOpenCross);
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

      double gerr = sqrt(pow(ncollErr[k+1]/ncoll[k+1],2)+pow(ppInelasticErr/ppInelastic,2));
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.pdf",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.png",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.eps",run_type.Data(),run_cfg_name.Data()));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_RaaVsPt_CompPub.pdf"));
    }


  //==============================================
  // final plot

  // model calculation
  TFile *fmodel = TFile::Open("Rootfiles/Paper/models/AuAu200.models.root", "read");
  TGraphErrors *gRaaVsPtModel[2][nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      gRaaVsPtModel[0][k] = (TGraphErrors*)fmodel->Get(Form("gRaaVsPt_cent%s_TAMU",cent_Title[k]));
      gRaaVsPtModel[0][k]->SetFillStyle(3644);
      gRaaVsPtModel[0][k]->SetFillColor(kGray+2);
      gRaaVsPtModel[1][k] = (TGraphErrors*)fmodel->Get(Form("gRaaVsPt_cent%s_Tsinghua",cent_Title[k]));
      gRaaVsPtModel[1][k]->SetLineWidth(2);
      gRaaVsPtModel[1][k]->SetLineColor(kViolet);
      gRaaVsPtModel[1][k]->SetLineStyle(1);
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
	  gCmsRaaVsPt->SetMarkerStyle(27);
	  gCmsRaaVsPt->SetMarkerSize(2.5);
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

      else
	{
	  if(k<4)
	    {
	      gRaaLowPtSys[k-1]->Draw("sameE5");
	      gRaaHighPtSys[k-1]->Draw("sameE5");
	      gRaaLowPt[k-1]->Draw("sames PEZ");
	      gRaaHighPt[k-1]->Draw("sames PEZ");
	    }
	  gPhenixRaaVsPtSys[k-1]->Draw("sameE5");
	  gPhenixRaaVsPt[k-1]->Draw("sames PEZ");
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

      if(k<4) gRaaVsPtModel[0][k]->Draw("samesE4");
      else    gRaaVsPtModel[0][k]->Draw("samesL");
      gRaaVsPtModel[1][k]->Draw("samesL");

      hJpsiRaaSys2[k]->Draw("sameE5");
      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      if(k==0 || k==3) centLabel[k] = GetPaveText(0.27,0.32,0.9,0.95);
      else centLabel[k] = GetPaveText(0.19,0.24,0.9,0.95);
      centLabel[k]->SetTextFont(63);
      centLabel[k]->SetTextSize(22);
      centLabel[k]->AddText(Form("(%s) %s%%",alphabet[k],cent_Name[k]));
      centLabel[k]->Draw();

      // Global systematics
      double gerr = sqrt(pow(ncollErr[k]/ncoll[k],2)+pow(ppInelasticErr/ppInelastic,2));
      double x_pos = 14.25;
      if(k==2 || k==3) x_pos = 10.5;
      if(k==4) x_pos = 7.75;
      box_ncoll[k] = new TBox(x_pos,1-gerr,x_pos+0.25*x_max[k]/x_max[0],1+gerr);
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

      else
	{
	  TBox *box_phenix = new TBox(x_pos-0.25*x_max[k]/x_max[0],1-phenix_gsys[k-1],x_pos,1+phenix_gsys[k-1]);
	  box_phenix->SetFillStyle(1001);
	  box_phenix->SetFillColor(gPhenixRaaVsPt[0]->GetMarkerColor());
	  box_phenix->Draw();
	}
    }
  c->cd();
  pads[0]->cd();
  TLegend *leg = new TLegend(0.15,0.77,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(16);
  leg->SetHeader("Inclusive");
  leg->AddEntry(gAliceRaaVsPt,"ALICE: Pb+Pb @ 2.76 TeV, 0-40%, |y| < 0.8","P");
  leg->Draw();

  TLegend *leg = new TLegend(0.15,0.65,0.45,0.77);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(16);
  leg->SetHeader("Prompt");
  leg->AddEntry(gCmsRaaVsPt,"CMS: Pb+Pb @ 2.76 TeV, 0-100%, |y| < 2.4","P");
  leg->Draw();

  pads[1]->cd();
  TLegend *leg = new TLegend(0.05,0.65,0.4,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(16);
  leg->SetHeader("Inclusive");
  leg->AddEntry(gPhenixRaaVsPt[0],"PHENIX: |y| < 0.35 (PRL 98 (2007) 232301)","P");
  leg->AddEntry(gRaaLowPt[0],"STAR: |y| < 1 (PRC 90 (2014) 024906)","P");
  leg->AddEntry(gRaaHighPt[0],"STAR: |y| < 1 (PLB 733 (2013) 55)","P");
  leg->Draw();

  pads[4]->cd();
  TLegend *leg = new TLegend(0.5,0.9,0.8,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(18);
  leg->AddEntry(gPhenixRaaVsPt[3],"PHENIX: 60-92%","P");
  leg->Draw();

  c->cd();
  TPad *pad = GetSinglePad(Form("pad_cent6"), 0.68, 0.98, 0.01, 0.49);
  pad->Draw();
  pad->cd();
  TLegend *leg30 = new TLegend(0.25,0.45,0.7,0.87);
  leg30->SetBorderSize(0);
  leg30->SetFillColor(0);
  leg30->SetTextFont(63);
  leg30->SetTextSize(20);
  leg30->SetHeader("Au+Au @ 200 GeV");
  leg30->AddEntry(hJpsiRaa[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg30->AddEntry(hJpsiRaaSys[0],"Au+Au uncertainty","F");
  leg30->AddEntry(hJpsiRaaSys2[0],"p+p uncertainty","F");
  leg30->Draw();

  TLegend *leg32 = new TLegend(0.25,0.25,0.7,0.45);
  leg32->SetBorderSize(0);
  leg32->SetFillColor(0);
  leg32->SetTextFont(63);
  leg32->SetTextSize(20);
  leg32->AddEntry(gRaaVsPtModel[1][0],"TM I: Tsinghua","L");
  leg32->AddEntry(gRaaVsPtModel[0][0],"TM II: TAMU","F");
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.pdf",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.png",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.eps",run_type.Data(),run_cfg_name.Data()));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_RaaVsPt.pdf"));
    }
  if(gSavePaper)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Paper/%s/Figure_JpsiRaaVsPt.pdf",gPaperVersion));
    }

  // calculate chi2/NDF bewteen data and theory
  const int NDFs[nCentBins] = {9, 9, 8, 8, 6};
  double chi2[nCentBins][2][3];
  double ndf[nCentBins][2][3];
  double pvalue[nCentBins][2][3];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<2; i++)
	{
	  printf("+++ centrality %s%% - %s +++ \n",cent_Name[k],model_name[i]);
	  for(int j=0; j<3; j++)
	    {
	      chi2[k][i][j] = 0;
	      ndf[k][i][j] = 0;
	      for(int bin=1; bin<=NDFs[k]; bin++)
		{
		  hJpsiRaa[k]->GetPoint(bin-1, x, y);
		  double pt = x;
		  if(i==1)
		    {
		      if(k==0 || k==1 || k==4) { if(pt>6) continue; }
		      else { if(pt>8) continue; }
		    }
		  ndf[k][i][j]++;
		  double raa_data = y;
		  double err_data = sqrt( pow(hJpsiRaa[k]->GetErrorYhigh(bin-1), 2) + 
					  pow(hJpsiRaaSys[k]->GetErrorYhigh(bin-1), 2) + 
					  pow((hJpsiRaaSys2[k]->GetErrorYhigh(bin-1)+hJpsiRaaSys2[k]->GetErrorYlow(bin-1))/2, 2) );
		  if(j==1) 
		    {
		      raa_data = y * box_ncoll[k]->GetY2();
		      err_data = err_data/y*raa_data;
		    }
		  if(j==2) 
		    {
		      raa_data = y * box_ncoll[k]->GetY1();
		      err_data = err_data/y*raa_data;
		    }
	
		  double *model_x = gRaaVsPtModel[i][k]->GetX();
		  double *model_y = gRaaVsPtModel[i][k]->GetY();
		  double *model_ey = gRaaVsPtModel[i][k]->GetEY();
		  double raa_model, err_model;
		  for(int ipoint=0; ipoint<gRaaVsPtModel[i][k]->GetN(); ipoint++)
		    {
		      if(pt< model_x[ipoint+1] && pt >= model_x[ipoint] )
			{
			  raa_model = model_y[ipoint];
			  err_model = model_ey[ipoint];
			  break;
			}
		    }
		  double err_tot = sqrt(pow(err_data,2) + pow(err_model,2)); 
		  chi2[k][i][j] += pow((raa_data-raa_model)/err_tot,2);
		  //printf("pt = %4.2f, data = %4.2f #pm %4.2f, model = %4.2f #pm %4.2f\n",pt,raa_data,err_data,raa_model,err_model);
		}
	      pvalue[k][i][j] = TMath::Prob(chi2[k][i][j],ndf[k][i][j]);
	    }
	  
	  for(int j=0; j<3; j++)
	    {
	      if(pvalue[k][i][j]>0.1)
		printf("%4.2f/%d & %2.2f & ", chi2[k][i][j],ndf[k][i][j], pvalue[k][i][j]);
	      else if(pvalue[k][i][j]>0.01)
		printf("%4.2f/%d & %2.3f & ", chi2[k][i][j],ndf[k][i][j], pvalue[k][i][j]);
	      else
		printf("%4.2f/%d & %2.2e &", chi2[k][i][j],ndf[k][i][j], pvalue[k][i][j]);
	    }
	  printf("\n");
	}
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
void xsec_v2(const int mode = 1, const bool savePlot = 0, const bool saveHisto = 0)
{
  // mode 0: get histogram from scratch
  // mode 1: get histogram from saved 

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
  else fout = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type.Data()),"read");
  TGraphAsymmErrors *hJpsipp = (TGraphAsymmErrors*)fout->Get("hpp200JpsiVsPtFinalSys");

  TGraphAsymmErrors *hJpsiXsec[nCentBins];
  TGraphAsymmErrors *hJpsiXsecSys[nCentBins];
  TGraphAsymmErrors *hJpsiRaa[nCentBins];
  TGraphAsymmErrors *hJpsiRaaSys[nCentBins];
  TGraphAsymmErrors *hJpsiRaaSys2[nCentBins];

  double x, y;
  double x1, y1;
  const double x_err = 0.25;
  if(mode==0)
    {
      TFile *fSys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type.Data()),"read");
      TH1F *hAuAuJpsiSys[nCentBins];
      for(int k=0; k<nCentBins; k++)
	{
	  hAuAuJpsiSys[k]= (TH1F*)fSys->Get(Form("JpsiSysVsPt_All_cent%s",cent_Title[k]));
	}

      TFile *fdata = TFile::Open(Form("Rootfiles/%s.JpsiXsec.pt%1.1f.pt%1.1f.root",run_type.Data(),pt1_cut,pt2_cut),"read");
      TH1F *hJpsiInvYield[nCentBins];

      TH1F *hJpsiPtPos = (TH1F*)fdata->Get("hJpsiPtPos_cent0080");
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
    }
  else
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiXsec[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
	  hJpsiXsecSys[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));

	  hJpsiRaa[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
	  hJpsiRaaSys[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
	  hJpsiRaaSys2[k] = (TGraphAsymmErrors*)fout->Get(Form("Graph_Jpsi_Raa_cent%s_sys_pp",cent_Title[k]));
	}
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
      if(k==1)
	{
	  hJpsiXsec[k]->GetPoint(0, x, y);
	  printf("[i] Jpsi yield = %4.2e +/- %4.2e +/- %4.2e\n",y,hJpsiXsec[k]->GetErrorYhigh(0),hJpsiXsecSys[k]->GetErrorYhigh(0));
	}
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
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.pdf",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.png",run_type.Data(),run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/AuAu_JpsiInvYield.eps",run_type.Data(),run_cfg_name.Data()));
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

  /*
    // CMS: Eur. Phys. J. C 05 (2012) 063
    // Inclusive Jpsi, |y| < 2.4, 0-100%
    // global uncertainty = 8.3%
    const int ncms = 2;
    double cms_pt[ncms] = {8.11, 13.22};
    double cms_pt_err_low[ncms] = {1.61, 3.22};
    double cms_pt_err_high[ncms] = {1.89, 16.78};
    double cms_pt_sys_low[ncms] = {0.25, 0.25};
    double cms_pt_sys_high[ncms] = {0.25, 0.25};
    double cms_raa[ncms] = {0.32, 0.31};
    double cms_raa_err_low[ncms] = {0.03, 0.04};
    double cms_raa_err_high[ncms] = {0.03, 0.04};
    double cms_raa_sys_low[ncms] = {0.02, 0.01};
    double cms_raa_sys_high[ncms] = {0.02, 0.01};

   */
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
      gPhenixRaaVsPt[k]->SetMarkerStyle(kOpenCross);
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

      double gerr = sqrt(pow(ncollErr[k+1]/ncoll[k+1],2)+pow(ppInelasticErr/ppInelastic,2));
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.pdf",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.png",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_vs_pub.eps",run_type.Data(),run_cfg_name.Data()));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_RaaVsPt_CompPub.pdf"));
    }


  //==============================================
  // final plot

  // model calculation
  TFile *fmodel = TFile::Open("Rootfiles/Paper/models/AuAu200.models.root", "read");
  TGraphErrors *gRaaVsPtModel[2][nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      gRaaVsPtModel[0][k] = (TGraphErrors*)fmodel->Get(Form("gRaaVsPt_cent%s_TAMU",cent_Title[k]));
      gRaaVsPtModel[0][k]->SetFillStyle(1001);
      gRaaVsPtModel[0][k]->SetFillColor(29);
      gRaaVsPtModel[0][k]->SetLineColor(29);
      gRaaVsPtModel[1][k] = (TGraphErrors*)fmodel->Get(Form("gRaaVsPt_cent%s_Tsinghua",cent_Title[k]));
      gRaaVsPtModel[1][k]->SetLineWidth(2);
      gRaaVsPtModel[1][k]->SetLineColor(kViolet);
      gRaaVsPtModel[1][k]->SetLineStyle(1);
    }

  TCanvas *c = new TCanvas("Raa_Jpsi_dimuon","Raa_Jpsi_dimuon",1200,700);
  TBox *box_ncoll[nCentBins];
  TPaveText *centLabel[nCentBins];
  TPad *pads[nCentBins];
  const char* alphabet[nCentBins] = {"a", "b", "c", "d", "e"};
  for(int k=0; k<nCentBins; k++)
    {
      c->cd();
      if(k==0) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.345, 0.68, 0.51, 0.99);
      if(k==1) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.68, 0.98, 0.51, 0.99);
      if(k==2) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.05, 0.38, 0.01, 0.49);
      if(k==3) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.38, 0.68, 0.01, 0.49);
      if(k==4) pads[k] = GetSinglePad(Form("pad_cent%d",k), 0.68, 0.98, 0.01, 0.49);
      pads[k]->Draw();
      pads[k]->cd();
      SetPadMargin(gPad, 0.15, 0.13, 0.005, 0.005);

      TH1F *hRaa = new TH1F(Form("Raa_Jpsi_cent%d",k),";p_{T} [GeV/c];",150,0,15);
      hRaa->GetYaxis()->SetRangeUser(0.01,1.75);
      hRaa->GetYaxis()->CenterTitle();
      ScaleHistoTitle(hRaa,22,2.1,20,24,1.9,20,63);
      if(k!=0 && k!=2) 
	{
	  hRaa->GetYaxis()->SetLabelSize(0);
	  gPad->SetLeftMargin(0.02);
	}
      if(k==2 || k==3)  hRaa->GetXaxis()->SetRangeUser(0, 11);
      if(k==4)  hRaa->GetXaxis()->SetRangeUser(0, 8.1);
      hRaa->DrawCopy();

      if(k<4) gRaaVsPtModel[0][k]->Draw("samesE4");
      else    gRaaVsPtModel[0][k]->Draw("samesL");
      gRaaVsPtModel[1][k]->Draw("samesL");


      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();

      if(k==0)
	{
	  gCmsRaaVsPt->SetMarkerStyle(27);
	  gCmsRaaVsPt->SetMarkerSize(2.5);
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

      else
	{
	  if(k<4)
	    {
	      gRaaLowPtSys[k-1]->Draw("sameE5");
	      gRaaHighPtSys[k-1]->Draw("sameE5");
	      gRaaLowPt[k-1]->Draw("sames PEZ");
	      gRaaHighPt[k-1]->Draw("sames PEZ");
	    }
	  gPhenixRaaVsPtSys[k-1]->Draw("sameE5");
	  gPhenixRaaVsPt[k-1]->Draw("sames PEZ");
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

      hJpsiRaaSys2[k]->Draw("sameE5");
      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      if(k==0 || k==2) centLabel[k] = GetPaveText(0.27,0.32,0.9,0.95);
      else centLabel[k] = GetPaveText(0.19,0.24,0.9,0.95);
      centLabel[k]->SetTextFont(63);
      centLabel[k]->SetTextSize(22);
      centLabel[k]->AddText(Form("(%s) %s%%",alphabet[k],cent_Name[k]));
      centLabel[k]->Draw();

      // Global systematics
      double gerr = sqrt(pow(ncollErr[k]/ncoll[k],2)+pow(ppInelasticErr/ppInelastic,2));
      double x_pos = 14.25;
      if(k==2 || k==3) x_pos = 10.5;
      if(k==4) x_pos = 7.75;
      box_ncoll[k] = new TBox(x_pos,1-gerr,x_pos+0.25*x_max[k]/x_max[0],1+gerr);
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

      else
	{
	  TBox *box_phenix = new TBox(x_pos-0.25*x_max[k]/x_max[0],1-phenix_gsys[k-1],x_pos,1+phenix_gsys[k-1]);
	  box_phenix->SetFillStyle(1001);
	  box_phenix->SetFillColor(gPhenixRaaVsPt[0]->GetMarkerColor());
	  box_phenix->Draw();
	}
    }
  c->cd();
  TPad *pad = GetSinglePad(Form("pad_cent6"), 0.00, 0.34, 0.50, 0.99);
  pad->Draw();
  pad->cd();
  TLegend *leg30 = new TLegend(0.02,0.45,0.5,0.95);
  leg30->SetBorderSize(0);
  leg30->SetFillColor(0);
  leg30->SetTextFont(63);
  leg30->SetTextSize(17);
  leg30->SetHeader("Au+Au @ 200 GeV, Inclusive J/#psi");
  leg30->AddEntry(hJpsiRaa[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|< 0.5","P");
  leg30->AddEntry(hJpsiRaaSys[0],"Au+Au uncertainty","F");
  leg30->AddEntry(hJpsiRaaSys2[0],"p+p uncertainty","F");
  leg30->AddEntry(gPhenixRaaVsPt[0],"PHENIX: |y|<0.35 (PRL98 (2007) 232301)","P");
  leg30->AddEntry(gRaaLowPt[0],"STAR: |y|<1 (PRC 90 (2014) 024906)","P");
  leg30->AddEntry(gRaaHighPt[0],"STAR: |y|<1 (PLB 733 (2013) 55)","P");
  leg30->Draw();

  TLegend *leg = new TLegend(0.02,0.2,0.55,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(17);
  leg->SetHeader("Pb+Pb @ 2.76 TeV");
  leg->AddEntry(gAliceRaaVsPt,"ALICE: Inclusive J/#psi, 0-40%, |y|<0.8","P");
  leg->AddEntry(gCmsRaaVsPt,"CMS: Prompt J/#psi, 0-100%, |y|<2.4","P");
  leg->Draw();

  pads[0]->cd();
  TLegend *leg32 = new TLegend(0.5,0.78,0.9,0.95);
  leg32->SetBorderSize(0);
  leg32->SetFillColor(0);
  leg32->SetTextFont(63);
  leg32->SetTextSize(16);
  leg32->SetHeader("Au+Au @ 200 GeV, |y| < 0.5");
  leg32->AddEntry(gRaaVsPtModel[1][0],"TM I: Tsinghua","L");
  leg32->AddEntry(gRaaVsPtModel[0][0],"TM II: TAMU","F");
  leg32->Draw();


  TLegend *leg = new TLegend(0.05,0.65,0.4,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(16);
  leg->SetHeader("Inclusive");
  leg->AddEntry(gPhenixRaaVsPt[0],"PHENIX: |y| < 0.35 (PRL 98 (2007) 232301)","P");
  leg->AddEntry(gRaaLowPt[0],"STAR: |y| < 1 (PRC 90 (2014) 024906)","P");
  leg->AddEntry(gRaaHighPt[0],"STAR: |y| < 1 (PLB 733 (2013) 55)","P");
  //leg->Draw();

  pads[4]->cd();
  TLegend *leg = new TLegend(0.5,0.9,0.8,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(18);
  leg->AddEntry(gPhenixRaaVsPt[3],"PHENIX: 60-92%","P");
  leg->Draw();

  c->cd();
  
  TPad *pad = GetSinglePad(Form("pad_new"), 0.30, 0.34, 0.51, 0.99);
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.pdf",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.png",run_type.Data(),run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaaPt_dimuon.eps",run_type.Data(),run_cfg_name.Data()));
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

  // di-electron distribution
  TFile *fe = TFile::Open("Rootfiles/lowmidhigh_cent_0_8_cut1_Run14BHT2.root","read");
  TH1F *hSignalElec;
  for(int i=0; i<5; i++)
    {
      TH1F *htmp = (TH1F*)fe->Get(Form("hsmmix_sig_%d",3+i));
      if(i==0) hSignalElec = (TH1F*)htmp->Clone("hSignalElec");
      else     hSignalElec->Add(htmp);
    }
  draw1D(hSignalElec);

  // di-muon distribution
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
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type.Data()),"read");
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
  double g_bin_width[nPtBins] = {0.03, 0.04};
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

      // increase the error bars of the data points in psi2S mass range
      for(int bin=hhSignal->FindBin(3.6); bin<=hhSignal->FindBin(3.75); bin++)
	{
	  //hhSignal->SetBinError(bin, hhSignal->GetBinError(bin)*10);
	}
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

      if(i==1)
	{
	  hSignalElec->Scale(hhSignal->GetBinContent(hhSignal->FindBin(3.1))/hSignalElec->GetBinContent(hSignalElec->FindBin(3.1)));
	  hSignalElec->SetMarkerStyle(25);
	  hSignalElec->SetMarkerColor(kGreen+2);
	  hSignalElec->SetLineColor(hSignalElec->GetMarkerColor());
	  hSignalElec->SetMarkerSize(hhSignal->GetMarkerSize());
	  hSignalElec->Draw("sames");
	}

      TF1 *funcJpsi = new TF1(Form("Fit_Jpsi_pt%s",pt_Name[i]),"gausn",2.45,3.8);
      for(int ipar=0; ipar<3; ipar++)
	{
	  funcJpsi->SetParameter(ipar, funcSignal->GetParameter(ipar));
	}

      TLegend *leg = new TLegend(0.15,0.7,0.3,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(63);
      leg->SetTextSize(18);
      if(i==0)
	{
	  leg->AddEntry(hplot,"Unlike-sign (UL) pairs (#times 0.15)","P");
	  leg->AddEntry(hhSeLS,"Like-sign (LS) pairs (#times 0.15)","PL");
	  leg->AddEntry(hhMixBkg,"Mixed-event (ME) (#times 0.15)","L");
	}
      else
	{
	  leg->AddEntry(hplot,"Unlike-sign pairs","P");
	  leg->AddEntry(hhSeLS,"Like-sign pairs","PL");
	  leg->AddEntry(hhMixBkg,"Mixed-event","L");
	}
      leg->Draw();

      TLegend *leg1 = new TLegend(0.13,0.37,0.28,0.47);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextFont(63);
      leg1->SetTextSize(18);
      leg1->AddEntry(hhSignal,"Signal = UL - ME", "P");
      leg1->AddEntry(funcSignal, "Fit to signal","L");
      leg1->Draw();

      if(i==1)
	{
	  TLegend *leg2 = new TLegend(0.13,0.32,0.28,0.37);
	  leg2->SetBorderSize(0);
	  leg2->SetFillColor(0);
	  leg2->SetTextFont(63);
	  leg2->SetTextSize(18);
	  leg2->AddEntry(hSignalElec,"J/#psi#rightarrowe^{+}+e^{-}", "P");
	  leg2->Draw();
	}

      double nJpsiAll = funcSignal->GetParameter(0)/hhSignal->GetBinWidth(1);
      double nJpsiAll_err = funcSignal->GetParError(0)/hhSignal->GetBinWidth(1);
      double min_mass = 3.0, max_mass = 3.2;
      double nAll = hSeUL[i]->Integral(hSeUL[i]->FindFixBin(min_mass),hSeUL[i]->FindFixBin(max_mass));
      double nJpsiMass = funcJpsi->Integral(min_mass, max_mass)/funcJpsi->Integral(2.5, 3.7)*nJpsiAll;


      TPaveText *t1 = GetPaveText(0.60,0.77,0.65,0.92,0.035,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("J/#psi#rightarrow#mu^{+}+#mu^{-}, #it{L} ~ 14.2 nb^{-1}"));
      if(i==0) t1->AddText(Form("|y| < 0.5, p_{T} > 0.15 GeV/c"));
      else     t1->AddText(Form("|y| < 0.5, p_{T} > %1.0f GeV/c",ptBins_low[i]));
      t1->AddText(Form("N_{J/#psi} = %2.0f #pm %2.0f",nJpsiAll, nJpsiAll_err));
      t1->AddText(Form("S/B = 1:%2.1f",(nAll-nJpsiMass)/nJpsiMass));
      t1->AddText(Form("#chi^{2}/NDF = %2.1f/%d\n",funcSignal->GetChisquare(),funcSignal->GetNDF()));
      t1->Draw();

      TPaveText *t1 = GetPaveText(0.14,0.2,0.84,0.94,0.038,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("Au+Au @ 200 GeV, %s%%",cent_Name_pt[icent]));
      //t1->AddText(Form("%s%%"));
      t1->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiSig_pt%s_cent%s.pdf",run_type.Data(),run_cfg_name.Data(),pt_Name[i],cent_Title_pt[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiSig_pt%s_cent%s.png",run_type.Data(),run_cfg_name.Data(),pt_Name[i],cent_Title_pt[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiSig_pt%s_cent%s.eps",run_type.Data(),run_cfg_name.Data(),pt_Name[i],cent_Title_pt[icent]));
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
void publication(const bool saveHisto = 0)
{
  // 2011 data
  // https://drupal.star.bnl.gov/STAR/files/starpublications/250/data.html
  const int npoints = 5;
  double pt[npoints] = {0.5, 1.5, 2.5, 3.5, 4.5};
  double pt_err[npoints] = {0.5, 0.5, 0.5, 0.5, 0.5};
  double yield[4][npoints] = { {1.57e-5, 7.44e-6, 2.69e-6, 7.20e-7, 2.37e-7},
			       {2.96e-5, 1.55e-5, 5.18e-6, 1.27e-6, 5.23e-7},
			       {1.36e-5, 5.72e-6, 2.30e-6, 6.76e-7, 1.75e-7},
			       {5.06e-6, 2.03e-6, 7.97e-7, 2.62e-7, 8.00e-8}};
  double yield_err[4][npoints] = { {2.67e-6, 9.72e-7, 3.81e-7, 1.44e-7, 5.45e-8},
				   {8.96e-6, 3.29e-6, 1.29e-6, 4.80e-7, 1.84e-7},
				   {2.31e-6, 8.44e-7, 3.43e-7, 1.40e-7, 5.62e-8},
				   {6.75e-7, 2.48e-7, 1.05e-7, 4.36e-8, 1.91e-8}};
  double yield_sys[4][npoints] = { {8, 12, 10, 13, 11},
				   {10, 12, 11, 13, 9},
				   {8, 13, 10, 12, 10},
				   {8, 12, 10, 12, 10}};

  TGraphErrors *gJpsiYield[4];
  TGraphErrors *gJpsiYieldSys[4];
  const char *cent[4] = {"0060","0020","2040","4060"};
  for(int i=0; i<4; i++)
    {
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  yield_sys[i][ipoint] = yield_sys[i][ipoint]*0.01*yield[i][ipoint];
	}
      gJpsiYield[i] = new TGraphErrors(npoints, pt, yield[i], pt_err, yield_err[i]);
      gJpsiYield[i]->SetName(Form("Run11_JpsInvYield_cent%s",cent[i]));
      gJpsiYieldSys[i] = new TGraphErrors(npoints, pt, yield[i], pt_err, yield_sys[i]);
      gJpsiYieldSys[i]->SetName(Form("Run11_JpsInvYieldSys_cent%s",cent[i]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root", "update");
      for(int i=0; i<4; i++)
	{
	  gJpsiYield[i]->Write("",TObject::kOverwrite);
	  gJpsiYieldSys[i]->Write("",TObject::kOverwrite);
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

  const char *name[3] = {"STAR_2009","STAR_2012","PHENIX"};

  // check the y-distribution
  // PHYSICAL REVIEW C 93, 024919 (2016)
  double midRapFrac[3] = {0, 0, 0};
  double midRapFrac2[3] = {0, 0, 0};
  double midRapFracStat[3] = {0, 0, 0};
  double midRapFracSys[3] = {0, 0, 0};
  TCanvas *c = new TCanvas("Jpsi_yDis","Jpsi_yDis",800,600);
  TH1F *hplot = new TH1F("hplot",";y/y_{max};a.u.",-10,-1,1);
  hplot->GetYaxis()->SetRangeUser(0, 1.5);
  hplot->Draw();
  const double ymax = TMath::Log(200./3.097);
  TF1 *funcy1 = new TF1("func_y1","gausn",-0.6,0.6);
  funcy1->SetLineColor(2);
  funcy1->SetLineStyle(2);
  funcy1->SetParameters(9.99561e-01, 0, 3.67483e-01);
  funcy1->Draw("sames");
  TF1 *funcy2 = new TF1("func_y2","[1]*1.0/(1.0-x**2)*TMath::Exp(-[0]*pow(TMath::Log((1+x)/(1-x)),2))",-0.6,0.6);
  funcy2->SetParameters(9.69153e-01,1.04028);
  funcy2->Draw("sames");
  const double sigmaError = 1.09661e-02;
  TF1 *funTmp1 = (TF1*)funcy1->Clone(Form("%s_tmp1",funcy1->GetName()));
  funTmp1->SetParameter(2, funcy1->GetParameter(2)+sigmaError);
  TF1 *funTmp2 = (TF1*)funcy1->Clone(Form("%s_tmp2",funcy1->GetName()));
  funTmp2->SetParameter(2, funcy1->GetParameter(2)-sigmaError);
  for(int i=0; i<3; i++)
    {   
      // systematic uncertainty due to different fit functions
      if(i<2)
	{
	  midRapFrac[i] = funcy1->Integral(-0.5/ymax,0.5/ymax)/funcy1->Integral(-1/ymax,1/ymax) * 2;
	  midRapFrac2[i] = funcy2->Integral(-0.5/ymax,0.5/ymax)/funcy2->Integral(-1/ymax,1/ymax) * 2;
	  printf("[i] fraction = %4.4f, %4.4f\n",midRapFrac[i]/2, midRapFrac2[i]/2);
	}
      else
	{
	  midRapFrac[i] = 0.7 / (funcy1->Integral(-0.35/ymax,0.35/ymax)/funcy1->Integral(-0.5/ymax,0.5/ymax)) ;
	  midRapFrac2[i] = 0.7 / (funcy2->Integral(-0.35/ymax,0.35/ymax)/funcy2->Integral(-0.5/ymax,0.5/ymax));
	  printf("[i] fraction = %4.4f, %4.4f\n",0.7/midRapFrac[i], 0.7/midRapFrac2[i]);
	}
      midRapFracSys[i] = fabs(midRapFrac2[i]/midRapFrac[i]-1);

      // systematic uncertainty due to fit errors
      double fracTmp1 = funTmp1->Integral(-0.5/ymax,0.5/ymax)/funTmp1->Integral(-1/ymax,1/ymax) * 2;
      double fracTmp2 = funTmp2->Integral(-0.5/ymax,0.5/ymax)/funTmp2->Integral(-1/ymax,1/ymax) * 2;
      if(i==2)
	{
	  fracTmp1 = 0.7 / (funTmp1->Integral(-0.35/ymax,0.35/ymax)/funTmp1->Integral(-0.5/ymax,0.5/ymax));
	  fracTmp2 = 0.7 / (funTmp2->Integral(-0.35/ymax,0.35/ymax)/funTmp2->Integral(-0.5/ymax,0.5/ymax));
	}
      midRapFracStat[i] = fabs(fracTmp1/midRapFrac[i]-1) > fabs(fracTmp2/midRapFrac[i]-1) ? fabs(fracTmp1/midRapFrac[i]-1) : fabs(fracTmp2/midRapFrac[i]-1);
      printf("[i] %s: scale = %4.2f +/- %4.2f%%, %4.2f%%\n",name[i],midRapFrac[i],midRapFracSys[i]*100,midRapFracStat[i]*100);

      // get the total error
      midRapFracSys[i] = sqrt(midRapFracSys[i]*midRapFracSys[i] + midRapFracStat[i]*midRapFracStat[i]);
      printf("[i] total error = %4.3f%%\n",midRapFracSys[i]*100);
    }
  leg = new TLegend(0.2,0.65,0.4,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(funcy1,Form("#frac{%2.4f}{2#pi#times%2.4f}e^{-0.5#times(#frac{x}{%2.4f})^{2}}",funcy1->GetParameter(0),funcy1->GetParameter(2),funcy1->GetParameter(2)),"L");
  leg->AddEntry(funcy2,Form("#frac{%2.4f}{1-(y/y_{max})^{2}}e^{-%2.4f[ln(#frac{1+y/y_{max}}{1-y/y_{max}})]^{2}}",funcy2->GetParameter(1),funcy2->GetParameter(0)),"L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/JpsiYDis.pdf",run_type.Data(),run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Ydis.pdf"));
    }

  // check the pt dependence in different rapidity bins
  const char* modelName[3] = {"ICEM", "NRQCD-JXW", "NRQCD-ZHF"};
  TGraphErrors *hppXsecRatio[3][2]; // two models; two ratios |y|<0.5/|y|<1; |y|<0.35/|y|<0.5
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<2; j++)
	{
	  if(i==0) hppXsecRatio[i][j] = new TGraphErrors(60);
	  else
	    {
	      if(j==0) hppXsecRatio[i][j] = new TGraphErrors(7);
	      if(j==1) hppXsecRatio[i][j] = new TGraphErrors(8);
	    }
	  hppXsecRatio[i][j]->SetName(Form("hppXsecRatio_%s_Rap%d",modelName[i],j));
	}
    }
  /// ICEM from Ramona
  char tmp_str[256];
  double tmp_d;
  ifstream file;
  file.open("Rootfiles/Paper/pp/ptratio.agr");
  double pt_tmp, ratio_tmp;
  for(int i=0; i<346; i++)
    {
      file.getline(tmp_str, 256);
    }
  for(int i=0; i<60; i++)
    {
      file >> pt_tmp >> ratio_tmp;
      hppXsecRatio[0][0]->SetPoint(i, pt_tmp, ratio_tmp);
    }
  for(int i=0; i<4; i++)
    {
      file.getline(tmp_str, 256);
    }
  for(int i=0; i<60; i++)
    {
      file >> pt_tmp >> ratio_tmp;
      hppXsecRatio[0][1]->SetPoint(i, pt_tmp, ratio_tmp);
    }
  double x1_tmp, y1_tmp, x2_tmp, y2_tmp;
  for(int i=0; i<60; i++)
    {
      hppXsecRatio[0][0]->GetPoint(i, x1_tmp, y1_tmp);
      hppXsecRatio[0][1]->GetPoint(i, x2_tmp, y2_tmp);
      hppXsecRatio[0][0]->SetPoint(i, x1_tmp, y2_tmp/y1_tmp);
      hppXsecRatio[0][0]->SetPointError(i, 0.25, 1e-10);
      hppXsecRatio[0][1]->SetPoint(i, x1_tmp, 1./y2_tmp);
      hppXsecRatio[0][1]->SetPointError(i, 0.25, 1e-10);
    }
  file.close();
  
  // NRQCD
  double pt_model[2][3][10], yield_model[2][3][10], error_model[2][3][10];
  const char* modelNameTmp[2] = {"JXW", "ZHF"};
  const char* rapidityNameTmp[3] = {"035", "05", "10"};
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<3; j++)
	{
	  if(i==0) file.open(Form("Rootfiles/Paper/pp/UpdatedJXW/Jpsi-Y%s-dsig-dpt-%s.dat",rapidityNameTmp[j],modelNameTmp[i]));
	  else file.open(Form("Rootfiles/Paper/pp/dsig_dpt/Jpsi-Y%s-dsig-dpt-%s.dat",rapidityNameTmp[j],modelNameTmp[i]));
	  for(int l=0; l<5; l++)
	    {
	      file.getline(tmp_str, 256);
	    }
	  for(int l=0; l<9; l++)
	    {
	      if(j<2 && l==8) continue;
	      file >> pt_model[i][j][l] >>  yield_model[i][j][l] >>  error_model[i][j][l]; 
	    }
	  file.close();
	}

      for(int ipoint=0; ipoint<7; ipoint++)
	{
	  double xerror = 0.5;
	  if(ipoint>=3) xerror = 1;
	  double yerror =  yield_model[i][1][ipoint+1]/yield_model[i][2][ipoint] * fabs(error_model[i][1][ipoint+1]/yield_model[i][1][ipoint+1] - error_model[i][2][ipoint]/yield_model[i][2][ipoint]);

	  hppXsecRatio[i+1][0]->SetPoint(ipoint, pt_model[i][1][ipoint+1], yield_model[i][1][ipoint+1]/yield_model[i][2][ipoint]);
	  hppXsecRatio[i+1][0]->SetPointError(ipoint, xerror, yerror);
	}

      for(int ipoint=0; ipoint<8; ipoint++)
	{
	  double xerror = 0.5;
	  if(ipoint>=4) xerror = 1;
	  double yerror =  yield_model[i][0][ipoint]/yield_model[i][1][ipoint] * fabs(error_model[i][0][ipoint]/yield_model[i][0][ipoint] - error_model[i][1][ipoint]/yield_model[i][1][ipoint]);

	  hppXsecRatio[i+1][1]->SetPoint(ipoint, pt_model[i][0][ipoint], yield_model[i][0][ipoint]/yield_model[i][1][ipoint]);
	  hppXsecRatio[i+1][1]->SetPointError(ipoint, xerror, yerror);
	}
    }

  //|y|<0.5/|y|<1; |y|<0.35/|y|<0.5
  TCanvas *c = new TCanvas("Jpsi_ptDisInY","Jpsi_ptDisInY",800,600);
  TH1F *hplotInY = new TH1F("hplotInY",";p_{T} [GeV/c];J/#Psi yield ratio",150,0,15);
  hplotInY->GetYaxis()->SetRangeUser(0, 0.9);
  hplotInY->DrawCopy();
  TLine *linePar[3][2];
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<2; j++)
	{
	  if(i==0)
	    {
	      double fraction = midRapFrac[0]/2;
	      if(j==1) fraction = 0.7/midRapFrac[2];
	      linePar[i][j] = GetLine(0, fraction, 15, fraction, j+1);
	      linePar[i][j]->Draw();
	    }
	  hppXsecRatio[i][j]->SetMarkerStyle(24-i*2);
	  hppXsecRatio[i][j]->SetMarkerSize(1.5);
	  hppXsecRatio[i][j]->SetMarkerColor(j+1);
	  hppXsecRatio[i][j]->SetLineColor(j+1);
	  hppXsecRatio[i][j]->Draw("samesP");
	}
    }
  leg = new TLegend(0.12,0.2,0.3,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hppXsecRatio[2][0],"|y|<0.5/|y|<1","P");
  leg->AddEntry(hppXsecRatio[2][1],"|y|<0.35/|y|<0.5","P");
  leg->Draw();
  leg = new TLegend(0.35,0.15,0.6,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hppXsecRatio[0][0],"ICEM [PRD94.114029 (2016)]","P");
  leg->AddEntry(hppXsecRatio[1][0],"NRQCD [PRL110.042002(2013)]","P");
  leg->AddEntry(hppXsecRatio[2][0],"NRQCD [PRL114.092006(2015)]","P");
  leg->AddEntry(linePar[0][0],"World data [PRC93.024919(2016)]","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/JpsiPtDisInY.pdf",run_type.Data(),run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_PtDisInY.pdf"));
    }


  // fit the ICEM results to be used for scaling data
  TF1 *funcICEM[2];
  TCanvas *c = new TCanvas("fit_ICEM_ratio","fit_ICEM_ratio",1100,450);
  c->Divide(2,1);
  for(int j=0; j<2; j++)
    {
      TGraphErrors* graphTmp = (TGraphErrors*)hppXsecRatio[0][j]->Clone(Form("%s_clone",hppXsecRatio[0][j]->GetName()));
      for(int ipoint=0; ipoint<graphTmp->GetN(); ipoint++)
	{
	  graphTmp->SetPointError(ipoint, 0.25, 1e-3);
	}
      funcICEM[j] = new TF1(Form("funcICEM_ratio%d",j),"pol4",0,15);
      if(j==0) funcICEM[j]->SetParameters(0.5, -1e-2, 2e-3, -1.5e-4, 4e-6);
      if(j==1) funcICEM[j]->SetParameters(0.7, -4e-3, 1e-3, -9e-5, 2e-6);
      graphTmp->Fit(funcICEM[j], "R0Q");
      c->cd(j+1);
      if(j==0) hplotInY->GetYaxis()->SetRangeUser(0.3,0.7);
      if(j==1) hplotInY->GetYaxis()->SetRangeUser(0.55,0.8);
      hplotInY->DrawCopy();
      graphTmp->Draw("samesPE");
      funcICEM[j]->Draw("sames");
      linePar[0][j]->Draw();
      TPaveText *t1;
      if(j==0) t1 = GetTitleText("ICEM: |y|<0.5/|y|<1", 0.05);
      else t1 = GetTitleText("ICEM: |y|<0.35/|y|<0.5", 0.05);
      t1->Draw();
      leg = new TLegend(0.15,0.15,0.4,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(graphTmp,"ICEM [PRD94.114029 (2016)]","P");
      leg->AddEntry(funcICEM[j],"Fit","L");
      leg->AddEntry(linePar[0][0],"World data [PRC93.024919(2016)]","L");
      leg->Draw();
    }
  return;

  // pp data
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
  gHighPtPP->SetName("STAR_2009_xsec");
  gHighPtPP->SetMarkerStyle(20);
  gHighPtPP->SetMarkerSize(1.5);

  TGraphAsymmErrors *gHighPtPPSys = new TGraphAsymmErrors(npp, xpp, ypp, sxlpp, sxlpp, sylpp, syhpp);
  gHighPtPPSys->SetName("STAR_2009_xsec_sys");
  gHighPtPPSys->SetMarkerSize(0);
  gHighPtPPSys->SetFillStyle(0);

  // STAR Run12 results
  TFile *frun12 = TFile::Open("Rootfiles/jpsi_xsec_pp200_run12.root","read");
  TGraphAsymmErrors *gRun12Sys = (TGraphAsymmErrors*)frun12->Get("gJpsiXsecCombSys"); 
  gRun12Sys->SetName("STAR_2012_xsec_sys");
  gRun12Sys->SetMarkerSize(0);
  gRun12Sys->SetFillStyle(0);
  gRun12Sys->SetLineColor(4);
  gRun12Sys->SetLineWidth(1);

  TGraphAsymmErrors *gRun12 = (TGraphAsymmErrors*)frun12->Get("gJpsiXsecCombAsy");
  gRun12->SetName("STAR_2012_xsec");
  gRun12->SetMarkerStyle(24);
  gRun12->SetMarkerSize(1.5);
  gRun12->SetMarkerColor(4);
  gRun12->SetLineColor(4);

  // PHENIX measurements
  TFile *fjpsi = TFile::Open(Form("Rootfiles/2016sQM/jpsi_xsec_pp200_run12.root",run_cfg_name.Data()),"read");
  TGraphAsymmErrors *gPhenixSys = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix_Systematics");
  gPhenixSys->SetName("PHENIX_xsec_sys");
  gPhenixSys->SetMarkerSize(0);
  gPhenixSys->SetFillStyle(0);
  gPhenixSys->SetLineColor(6);
  gPhenixSys->SetLineWidth(1);

  TGraphAsymmErrors *gPhenix = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix");
  gPhenix->SetName("PHENIX_xsec");
  gPhenix->SetMarkerStyle(24);
  gPhenix->SetMarkerSize(1.5);
  gPhenix->SetMarkerColor(6);
  gPhenix->SetLineColor(6);

  TH1F *hpp = new TH1F("pp200_Jpsi",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hpp->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hpp,"Invariant J/psi cross section",kTRUE);
  c->SetName("pp200_Jpsi_original");
  //gHighPtPPSys->Draw("sameE5");
  //gHighPtPP->Draw("sames PE");
  gRun12Sys->Draw("sameE5");
  gRun12->Draw("samePE");  
  gPhenixSys->Draw("sameE5");
  gPhenix->Draw("sames PE");

  leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  //leg->AddEntry(gHighPtPP,"STAR 2009 |y|<1","PL");
  leg->AddEntry(gRun12,"STAR 2012 |y|<1","PL");
  leg->AddEntry(gPhenix,"PHENIX |y|<0.35","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/ppRef_Compare.pdf",run_type.Data(),run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Comp.pdf"));
    }

  // scale the measured x-section according to Jpsi y-distribution
  TGraphAsymmErrors *hJpsiXsecScale[3];
  TGraphAsymmErrors *hJpsiXsecSysScale[3];
  TGraphAsymmErrors *graph = 0x0, *gSys = 0x0;
  TF1 *funcICEMscale = 0x0;
  double xtmp, ytmp;
  for(int i=0; i<3; i++)
    {
      if(i==0)
	{
	  graph = gHighPtPP;
	  gSys = gHighPtPPSys;
	  funcICEMscale = funcICEM[0];
	}
      else if(i==1)
	{
	  graph = gRun12;
	  gSys = gRun12Sys;
	  funcICEMscale = funcICEM[0];
	}
      else if(i==2)
	{
	  graph = gPhenix;
	  gSys = gPhenixSys;
	  funcICEMscale = funcICEM[1];
	}
      hJpsiXsecScale[i] = new TGraphAsymmErrors(*graph);
      hJpsiXsecScale[i]->SetName(Form("%s_scaled",graph->GetName()));
      for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
	{
	  graph->GetPoint(ipoint, xtmp, ytmp);
	  double scale = funcICEMscale->Eval(xtmp) * 2;
	  if(i==2) scale = 0.7/funcICEMscale->Eval(xtmp);
	  hJpsiXsecScale[i]->SetPoint(ipoint, xtmp, ytmp*scale);
	  hJpsiXsecScale[i]->SetPointError(ipoint, graph->GetErrorXlow(ipoint), graph->GetErrorXhigh(ipoint),
					   graph->GetErrorYlow(ipoint)*midRapFrac[i], graph->GetErrorYhigh(ipoint)*scale);
	}
      hJpsiXsecSysScale[i] = new TGraphAsymmErrors(gSys->GetN());
      hJpsiXsecSysScale[i]->SetName(Form("%s_scaled",gSys->GetName()));
      hJpsiXsecSysScale[i]->SetMarkerSize(0);
      hJpsiXsecSysScale[i]->SetFillStyle(0);
      hJpsiXsecSysScale[i]->SetLineColor(hJpsiXsecScale[i]->GetMarkerColor());
      for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
      	{
      	  gSys->GetPoint(ipoint, xtmp, ytmp);
	  double scale = funcICEMscale->Eval(xtmp) * 2;
	  if(i==2) scale = 0.7/funcICEMscale->Eval(xtmp);
      	  hJpsiXsecSysScale[i]->SetPoint(ipoint, xtmp, ytmp*scale);
      	  hJpsiXsecSysScale[i]->SetPointError(ipoint, gSys->GetErrorXlow(ipoint), gSys->GetErrorXhigh(ipoint),
					      sqrt( pow(gSys->GetErrorYlow(ipoint),2) + pow(midRapFracSys[i]*ytmp,2) )*scale, 
					      sqrt( pow(gSys->GetErrorYhigh(ipoint),2) + pow(midRapFracSys[i]*ytmp,2) )*scale);
      	}
    }
  TCanvas *c = draw1D(hpp,"Invariant J/psi cross section scaled to |y| < 0.5",kTRUE);
  hJpsiXsecSysScale[0]->Draw("sameE5");
  hJpsiXsecScale[0]->Draw("sames PE");
  hJpsiXsecSysScale[1]->Draw("sameE5");
  hJpsiXsecScale[1]->Draw("sames PE");
  hJpsiXsecSysScale[2]->Draw("sameE5");
  hJpsiXsecScale[2]->Draw("sames PE");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/ppRef_Compare_scaled.pdf",run_type.Data(),run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Comp_scaled.pdf"));
    }

  const int nResults = 3;
  TGraphAsymmErrors *hJpsiRebin[nResults];
  TGraphAsymmErrors *hJpsiRebinSys[nResults];
  const int markercolor[3] = {1, 4, 6};
  const int markerstyle[3] = {20, 24, 28};
  double fraction1[nResults];
  double fraction2[nResults];

  TF1 *funcJpsiXsec[nResults];
  double pT, y1, xh1, xl1, yh1, yl1, syh1, syl1;
  for(int i=0; i<nResults; i++)
    {
      if(i==0) continue;
      cout << name[i] << endl;
      hJpsiRebin[i] = new TGraphAsymmErrors(nbins);
      hJpsiRebin[i]->SetName(Form("Jpsi_xsec_%s",name[i]));
      hJpsiRebinSys[i] = new TGraphAsymmErrors(nbins);
      hJpsiRebinSys[i]->SetName(Form("Jpsi_xsec_%s_Sys",name[i]));
      graph = hJpsiXsecScale[i];
      gSys = hJpsiXsecSysScale[i];

      TCanvas *c = new TCanvas(Form("fit_xsec_%d",i),Form("fit_xsec_%s",name[i]),800,600);
      hpp->DrawCopy("");
      gPad->SetLogy();
      graph->GetXaxis()->SetRangeUser(0,15);
      graph->GetYaxis()->SetRangeUser(1e-6,10);
      funcJpsiXsec[i] = new TF1(Form("Func_Jpsi_xsec_%s",name[i]),"[0]*(([1]-1)*([1]-2))/(2*3.14*[1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]))**(-1*[1])",0,14);
      if(i==2) funcJpsiXsec[i]->SetRange(0,6);
      funcJpsiXsec[i]->SetParameters(1, 3, 1, 4);
      graph->Fit(funcJpsiXsec[i], "R0Q");
      graph->Draw("samesPE");
      funcJpsiXsec[i]->Draw("sames");
      TF1 *funcTmp = new TF1(Form("%s_yield",funcJpsiXsec[i]->GetName()),Form("%s*x",funcJpsiXsec[i]->GetExpFormula().Data()));
      funcTmp->SetParameters(funcJpsiXsec[i]->GetParameters());
      if(i==2) fraction1[i] = funcTmp->Integral(0.15,0.25)/funcTmp->Integral(0,0.25);
      else     fraction1[i] = funcTmp->Integral(0.15,0.5)/funcTmp->Integral(0,0.5);
      fraction2[i] = funcTmp->Integral(10,14)/funcTmp->Integral(10,15);
      TPaveText *t1 = GetTitleText(name[i],0.05);
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Fit_%s.pdf",run_type.Data(),run_cfg_name.Data(),name[i]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Fit%s.pdf",name[i]));
	}
      printf("[i] Fraction of [0.15,0.5]/[0,0.5] = %4.4f, [10,14]/[10,15] = %4.4f\n",fraction1[i],fraction2[i]);

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
	      if(i==1 && ipoint==15) continue; // reject unpublished data point
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
	      if(i==0)
		{
		  // STAR 2009
		  if(pT<4) dpT = 0.5;
		  else if(pT<8) dpT = 1;
		  else dpT = 2;
		}
	      else if(i==1)
		{
		  // STAR 2012
		  if(pT<4) dpT = 0.5;
		  else if(pT<8) dpT = 1;
		  else dpT = 2;
		}
	      else if(i==2)
		{
		  // Phenix
		  if(pT<5) dpT = 0.25;
		  else     dpT = 1;
		}

	      double yield = y1 * dpT * pT;
	      if(bin==1 && ipoint==0) 
		{
		  yield = yield * fraction1[i];
		}
	      y += yield;
	      yh += pow(yh1/y1 * yield, 2);
	      yl += pow(yl1/y1 * yield, 2);
	      syh += syh1/y1 * yield;
	      syl += syl1/y1 * yield;
	    }
	  y  = y / (xbins[bin]-xbins[bin-1]) / x;
	  yh = sqrt(yh) / (xbins[bin]-xbins[bin-1]) / x;
	  yl = sqrt(yl) / (xbins[bin]-xbins[bin-1]) / x;
	  syh = syh / (xbins[bin]-xbins[bin-1]) / x;
	  syl = syl / (xbins[bin]-xbins[bin-1]) / x;
	  if(y==0) continue;
	  if(i==2 && x>8) continue;
	  if(xbins[bin]>14) 
	    {
	      y/ = fraction2[i];
	      yl/ = fraction2[i];
	      yh/ = fraction2[i];
	      syl/ = fraction2[i];
	      syh/ = fraction2[i];
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

  // PHENIX fit results
  TF1 *funcPhenix = new TF1("PHENIX_paper","[0]*x*1./((1+(x/[1])**2)**[2])",0,8);
  funcPhenix->SetParameters(28.7, 3.41, 4.6);
  printf("[i] PHENIX paper: fraction of [0.15, 0.25]/[0, 0.25] = %4.3f\n",funcPhenix->Integral(0.15,0.25)/funcPhenix->Integral(0,0.25));

  TH1F *hppRebin = new TH1F("pp200_Jpsi_Rebin",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hppRebin->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hppRebin,"",kTRUE);
  for(int i=0; i<nResults; i++)
    {
      if(i==0) continue;
      hJpsiRebinSys[i]->Draw("sameE5");
      hJpsiRebin[i]->Draw("samePE");
    }
  leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  leg->AddEntry(gRun12,"STAR 2012 |y|<1","PL");
  leg->AddEntry(gPhenix,"PHENIX |y|<0.35","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/paper/ppRef_Compare_rebin.pdf",run_type.Data(),run_cfg_name.Data()));

  // take the average
  TH1F *hPPJpsiFinal = new TH1F("hpp200JpsiVsPtFinal","",nbins,xbins);
  TGraphAsymmErrors *hPPJpsiFinalSys = new TGraphAsymmErrors(nbins);
  hPPJpsiFinalSys->SetName("hpp200JpsiVsPtFinalSys");

  double xsec_0 = 0, xsec_0_errh = 0, xsec_0_errl = 0, xsec_0_sysh = 0, xsec_0_sysl = 0;
  double xsec_5 = 0, xsec_5_errh = 0, xsec_5_errl = 0, xsec_5_sysh = 0, xsec_5_sysl = 0;
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
	  if(i==0) continue;
	  hJpsiRebin[i]->GetPoint(ipoint, pT, y1);
	  if(y1<=0) continue;
	  yh1 = hJpsiRebin[i]->GetErrorYhigh(ipoint);
	  yl1 = hJpsiRebin[i]->GetErrorYlow(ipoint);
	  syh1 = hJpsiRebinSys[i]->GetErrorYhigh(ipoint);
	  syl1 = hJpsiRebinSys[i]->GetErrorYlow(ipoint);
	  double weight = 1./(pow((yh1+yl1)/2,2) + pow((syh1+syl1)/2,2));
	  total_weight += weight;
	  y += y1 * weight;
	  yh += pow(yh1*weight,2);
	  yl += pow(yl1*weight,2);
	  syh += syh1*weight;
	  syl += syl1*weight;
	}
      y /= total_weight;
      yh = sqrt(yh)/total_weight;
      yl = sqrt(yl)/total_weight;
      syh = syh/total_weight;
      syl = syl/total_weight;
      double tot_err_l = sqrt(yl*yl+syl*syl);
      double tot_err_h = sqrt(yh*yh+syh*syh);

      hPPJpsiFinal->SetBinContent(ipoint+1, y);
      hPPJpsiFinal->SetBinError(ipoint+1, (yh+yl)/2);
      hPPJpsiFinalSys->SetPoint(ipoint, x, y);
      hPPJpsiFinalSys->SetPointError(ipoint, xl, xh, tot_err_l, tot_err_h);

      double yield = y * 2 * pi * x * (xbins[ipoint+1]-xbins[ipoint]);
      double yield_errh = yh/y * yield;
      double yield_errl = yl/y * yield;
      double yield_sysh = syh/y * yield;
      double yield_sysl = syl/y * yield;
      if(x>0) 
	{
	  xsec_0 += yield;
	  xsec_0_errh += pow(yield_errh, 2);
	  xsec_0_errl += pow(yield_errl, 2);	  
	  xsec_0_sysh += yield_sysh;
	  xsec_0_sysl += yield_sysl;
	}
      if(x>5) 
	{
	  xsec_5 += yield;
	  xsec_5_errh += pow(yield_errh, 2);
	  xsec_5_errl += pow(yield_errh, 2);
	  xsec_5_sysh += yield_sysl;
	  xsec_5_sysl += yield_sysl; // to take out the small asymmetric in phenix error
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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Final.pdf",run_type.Data(),run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Comb.pdf"));
    }

  // compare the error bar
  TCanvas *c = new TCanvas("ppRatioToFit", "ppRatioToFit", 800, 600);
  TH1F *hRatio = new TH1F("hRatio",";p_{T} [GeV/c];Data/Combine",100,0,15);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.8);
  hRatio->Draw();
  TF1 *funcTmp = new TF1(Form("%s_yield2",funcJpsiXsec[1]->GetName()),Form("%s*x",funcJpsiXsec[1]->GetExpFormula().Data()));
  funcTmp->SetParameters(funcJpsiXsec[1]->GetParameters());
  TGraphAsymmErrors *gRatioToFit[3];
  TGraphAsymmErrors *gRatioToFitSys[3];
  double x,y;
  leg = new TLegend(0.6,0.7,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(int i=0; i<3; i++)
    {
      if(i==0)
	{
	  gRatioToFit[i] = new TGraphAsymmErrors(hPPJpsiFinal);
	  gRatioToFitSys[i] = new TGraphAsymmErrors(hPPJpsiFinal->GetNbinsX());
	  for(int ipoint=0; ipoint<gRatioToFitSys[i]->GetN(); ipoint++)
	    {
	      gRatioToFitSys[i]->SetPointError(ipoint, 0.2, 0.2, 
					       sqrt( pow(hPPJpsiFinalSys->GetErrorYhigh(ipoint), 2) - pow(hPPJpsiFinal->GetBinError(ipoint+1), 2) ),
					       sqrt( pow(hPPJpsiFinalSys->GetErrorYhigh(ipoint), 2) - pow(hPPJpsiFinal->GetBinError(ipoint+1), 2) ) );
	    }
	  gRatioToFit[i]->SetMarkerColor(kGreen+2);
	  gRatioToFit[i]->SetLineColor(kGreen+2);
	  gRatioToFit[i]->SetMarkerStyle(21);
	  gRatioToFit[i]->SetMarkerSize(1.5);
	  gRatioToFitSys[i]->SetFillStyle(0);
	  gRatioToFitSys[i]->SetMarkerColor(kGreen+2);
	  gRatioToFitSys[i]->SetLineColor(kGreen+2);
	}
      else
	{
	  gRatioToFit[i] = new TGraphAsymmErrors(*hJpsiRebin[i]);
	  gRatioToFitSys[i] = new TGraphAsymmErrors(*hJpsiRebinSys[i]);
	}

      for(int ipoint=0; ipoint<gRatioToFitSys[0]->GetN(); ipoint++)
	{
	  double scale = hPPJpsiFinal->GetBinContent(ipoint+1);
	  gRatioToFit[i]->GetPoint(ipoint, x, y);
	  x = x + i*0.2;
	  gRatioToFit[i]->SetPoint(ipoint, x, y/scale);
	  gRatioToFit[i]->SetPointError(ipoint, gRatioToFit[i]->GetErrorXlow(ipoint), gRatioToFit[i]->GetErrorXhigh(ipoint),
					gRatioToFit[i]->GetErrorYlow(ipoint)/scale, gRatioToFit[i]->GetErrorYhigh(ipoint)/scale);
	  gRatioToFitSys[i]->SetPoint(ipoint, x, y/scale);
	  gRatioToFitSys[i]->SetPointError(ipoint, 0.1, 0.1,
					   gRatioToFitSys[i]->GetErrorYlow(ipoint)/scale, gRatioToFitSys[i]->GetErrorYhigh(ipoint)/scale);
	  
	}
      gRatioToFitSys[i]->Draw("samesE5");
      gRatioToFit[i]->Draw("samesPEZ");
      if(i==0) leg->AddEntry(gRatioToFit[i], "Combined", "P");
      if(i==1) leg->AddEntry(gRatioToFit[i], "STAR", "P");
      if(i==2) leg->AddEntry(gRatioToFit[i], "PHENIX", "P");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Final_Ratio.pdf",run_type.Data(),run_cfg_name.Data()));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiXsecPP_Comb_Ratio.pdf"));
    }

  xsec_0_errh = sqrt(xsec_0_errh); xsec_0_errl = sqrt(xsec_0_errl);
  xsec_5_errh = sqrt(xsec_5_errh); xsec_5_errl = sqrt(xsec_5_errl);

  double xsec_0_toth = sqrt(xsec_0_errh*xsec_0_errh + xsec_0_sysh*xsec_0_sysh);
  double xsec_0_totl = sqrt(xsec_0_errl*xsec_0_errl + xsec_0_sysl*xsec_0_sysl);
  double xsec_5_toth = sqrt(xsec_5_errh*xsec_5_errh + xsec_5_sysh*xsec_5_sysh);
  double xsec_5_totl = sqrt(xsec_5_errl*xsec_5_errl + xsec_5_sysl*xsec_5_sysl);
  printf("[i] pT > 0, xsec = %4.2f +/- %4.3f +%4.3f - %4.3f\n",xsec_0,xsec_0_errh,xsec_0_sysh,xsec_0_sysl);
  printf("[i] pT > 5, xsec = %4.2e +/- %4.3e +%4.3e - %4.3e\n",xsec_5,xsec_5_errh,xsec_5_sysh,xsec_5_sysl);
  printf("[i] pT > 0, xsec = %4.2f +/- %4.3f%% +%4.3f%% - %4.3f%%\n",xsec_0,xsec_0_errh/xsec_0*100,xsec_0_sysh/xsec_0*100,xsec_0_sysl/xsec_0*100);
  printf("[i] pT > 5, xsec = %4.2e +/- %4.3f%% +%4.3f%% - %4.3f%%\n",xsec_5,xsec_5_errh/xsec_5*100,xsec_5_sysh/xsec_5*100,xsec_5_sysl/xsec_5*100);
  TH1F *hJpsiIntXsec = new TH1F("hpp200JpsiVsCentFinal",";p_{T} (GeV/c);Bd#sigma/dy [nb]",2,0,2);
  hJpsiIntXsec->SetBinContent(1, xsec_0);
  hJpsiIntXsec->SetBinError(1, (xsec_0_errh+xsec_0_errl)/2);
  hJpsiIntXsec->SetBinContent(2, xsec_5);
  hJpsiIntXsec->SetBinError(2, (xsec_5_errh+xsec_5_errl)/2);

  TGraphAsymmErrors *hJpsiIntXsecSys = new TGraphAsymmErrors(2);
  hJpsiIntXsecSys->SetName("hpp200JpsiVsCentFinalSys");
  hJpsiIntXsecSys->SetPoint(0, 0, xsec_0);
  hJpsiIntXsecSys->SetPointError(0, 0, 0, xsec_0_totl, xsec_0_toth);
  hJpsiIntXsecSys->SetPoint(1, 5, xsec_5);
  hJpsiIntXsecSys->SetPointError(1, 0, 0, xsec_5_totl, xsec_5_toth);
  
  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","update");
      fout->cd();
      for(int i=0; i<3; i++)
	{
	  if(i==0)
	    {
	      graph = gHighPtPP;
	      gSys = gHighPtPPSys;
	    }
	  else if(i==1)
	    {
	      graph = gRun12;
	      gSys = gRun12Sys;
	    }
	  else if(i==2)
	    {
	      graph = gPhenix;
	      gSys = gPhenixSys;
	    }
	  graph->Write("",TObject::kOverwrite);
	  gSys->Write("",TObject::kOverwrite);
	}
      hJpsiIntXsec->Write("",TObject::kOverwrite);
      hJpsiIntXsecSys->Write("",TObject::kOverwrite);
      hPPJpsiFinal->Write("",TObject::kOverwrite);
      hPPJpsiFinalSys->Write("",TObject::kOverwrite);
    }
}


//================================================
void ppRef2(const int savePlot = 0, const int saveHisto = 0)
{
  gStyle->SetOptStat(0);
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


  const int nHistos = 2;
  const char* data[nHistos] = {"Run15_dimuon","Run12_dielectron"};
  TFile *file[nHistos];
  file[0] = TFile::Open("Rootfiles/Run15pp200_pt1.5_pt1.3JpsiXsec.root", "read");
  file[1] = TFile::Open("Rootfiles/jpsi_xsec_pp200_run12_20180128.root", "read"); 
  TH1F *hJpsiXsec[nHistos];
  TGraphErrors *gJpsiXsecSys[nHistos];
  double x, y;
  const int nbins1 = 15;
  const double xbins1[nbins1+1] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0};
  for(int i=0; i<nHistos; i++)
    {
      if(i==0)
	{
	  hJpsiXsec[i] = (TH1F*)file[i]->Get("Runpp200_hxsec");
	  TH1D *hJpsiXsecSys1 = (TH1D*)file[i]->Get("Runpp200_hxsectotsys");
	  gJpsiXsecSys[i] = new TGraphErrors(hJpsiXsecSys1);
	}
      else
	{
	  hJpsiXsec[i] = new TH1F("hJpsiXsecComb_run12", "", nbins1, xbins1);
	  gJpsiXsecSys[i] = (TGraphErrors*)file[i]->Get("gJpsiXsecCombSys_run12");
	  TGraphErrors *gJpsiXsec = (TGraphErrors*)file[i]->Get("gJpsiXsecComb_run12");
	  for(int ipoint=0; ipoint<nbins1; ipoint++)
	    {
	      gJpsiXsec->GetPoint(ipoint, x, y);
	      double yield = y * x / hJpsiXsec[i]->GetBinCenter(ipoint+1);
	      double err =  gJpsiXsec->GetErrorY(ipoint)/y * yield;
	      hJpsiXsec[i]->SetBinContent(ipoint+1, y);
	      hJpsiXsec[i]->SetBinError(ipoint+1, err);
	      gJpsiXsecSys[i]->SetPoint(ipoint, hJpsiXsec[i]->GetBinCenter(ipoint+1), hJpsiXsec[i]->GetBinContent(ipoint+1));
	      gJpsiXsecSys[i]->SetPointError(ipoint, hJpsiXsec[i]->GetBinWidth(ipoint+1)/2, gJpsiXsecSys[i]->GetErrorY(ipoint)/y * yield);
	    }
	}
      int clr = color[i+1];
      hJpsiXsec[i]->SetMarkerStyle(20+i);
      hJpsiXsec[i]->SetMarkerSize(1.2);
      hJpsiXsec[i]->SetMarkerColor(clr);
      hJpsiXsec[i]->SetLineColor(clr);
      gJpsiXsecSys[i]->SetLineColor(clr);
      gJpsiXsecSys[i]->SetLineWidth(1);
    }
  TCanvas *c = draw1D(hJpsiXsec[1],"Differential cross-section for inclusive J/psi;p_{T} (GeV/c);B#times#frac{d^{2}#sigma}{2#pip_{T}dp_{T}dy}",true);
  hJpsiXsec[0]->Draw("sames");
  TLegend *leg = new TLegend(0.5,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  for(int i=0; i<nHistos; i++)
    {
      gJpsiXsecSys[i]->SetFillStyle(0);
      gJpsiXsecSys[i]->Draw("sameE5");
      leg->AddEntry(hJpsiXsec[i], data[i], "PL");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Run12vs15.pdf",run_type.Data(),run_cfg_name.Data()));

  // rebin
  TF1 *funcJpsiXsec[nHistos];
  TH1F *hppRef[nHistos];
  TGraphErrors *gppRefSys[nHistos];
  TGraphAsymmErrors *gppRefSysTot[nHistos];
  double fraction1[nHistos];
  double fraction2[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("fit_xsec_%d",data[i]),Form("fit_xsec_%d",data[i]),800,600);
      gPad->SetLogy();
      funcJpsiXsec[i] = new TF1(Form("Func_Jpsi_xsec_%s",data[i]),"[0]*(([1]-1)*([1]-2))/(2*3.14*[1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]))**(-1*[1])",0,14);
      funcJpsiXsec[i]->SetParameters(1, 3, 1, 3.1);
      hJpsiXsec[i]->Fit(funcJpsiXsec[i], "IR0");
      hJpsiXsec[i]->SetTitle(Form("%s: inclusive J/psi cross-section;p_{T} (GeV/c);B#times#frac{d^{2}#sigma}{2#pip_{T}dp_{T}dy}",data[i]));
      hJpsiXsec[i]->Draw("");
      funcJpsiXsec[i]->Draw("sames");
      TF1 *funcTmp = new TF1(Form("%s_yield",funcJpsiXsec[i]->GetName()),Form("%s*x",funcJpsiXsec[i]->GetExpFormula().Data()));
      funcTmp->SetParameters(funcJpsiXsec[i]->GetParameters());
      fraction1[i] = funcTmp->Integral(0.15,1)/funcTmp->Integral(0,1);
      fraction2[i] = funcTmp->Integral(10,14)/funcTmp->Integral(10,15);
      printf("[i] low pt fraction = %4.2f%%, high pt fraction = %4.2f%%\n",fraction1[i]*100,fraction2[i]*100);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Fit_%s.pdf",run_type.Data(),run_cfg_name.Data(),data[i]));
     
      hppRef[i] = new TH1F(Form("hpp200JpsiVsPtFinal_%s",data[i]),"",nbins,xbins);
      gppRefSys[i] = new TGraphErrors(nbins);
      gppRefSys[i]->SetName(Form("hpp200JpsiVsPtFinalSys_%s",data[i]));
      gppRefSysTot[i] = new TGraphAsymmErrors(nbins);
      gppRefSysTot[i]->SetName(Form("hpp200JpsiVsPtFinalSysTot_%s",data[i]));
      for(int ibin=1; ibin<=nbins; ibin++)
	{
	  double yield = 0, yerr = 0, ysys = 0;
	  for(int jbin=1; jbin<=hJpsiXsec[i]->GetNbinsX(); jbin++)
	    {
	      double pt1 = hJpsiXsec[i]->GetBinCenter(jbin);
	      if(pt1<xbins[ibin-1] || pt1>xbins[ibin]) continue;
	      double dpt1 = hJpsiXsec[i]->GetBinWidth(jbin);
	      double y1 = hJpsiXsec[i]->GetBinContent(jbin);
	      double yerr1 = hJpsiXsec[i]->GetBinError(jbin);
	      double ysys1 = gJpsiXsecSys[i]->GetErrorY(jbin-1);

	      yield += y1 * pt1 * dpt1;
	      yerr = pow(yerr1*pt1*dpt1, 2);
	      ysys = ysys1*pt1*dpt1;
	    }
	  ysys = ysys / yield;
	  double pt = hppRef[i]->GetBinCenter(ibin);
	  double dpt = hppRef[i]->GetBinWidth(ibin);
	  yield = yield / pt / dpt;
	  yerr = sqrt(yerr) / pt / dpt;
	  ysys = ysys * yield;
	  //if(ibin==1) cout << pt << "  " << dpt << endl;
	  if(ibin==1)
	    {
	      yield *= fraction1[i];
	      yerr  *= fraction1[i];
	      ysys  *= fraction1[i];
	    }
	  if(ibin==nbins)
	    {
	      yield /= fraction1[i];
	      yerr  /= fraction1[i];
	      ysys  /= fraction1[i];
	    }
	  hppRef[i]->SetBinContent(ibin, yield);
	  hppRef[i]->SetBinError(ibin, yerr);
	  gppRefSys[i]->SetPoint(ibin-1, pt, yield);
	  gppRefSys[i]->SetPointError(ibin-1, dpt/2, ysys);
	  gppRefSysTot[i]->SetPoint(ibin-1, pt, yield);
	  gppRefSysTot[i]->SetPointError(ibin-1, dpt/2, dpt/2, sqrt(ysys*ysys+yerr*yerr), sqrt(ysys*ysys+yerr*yerr));
	}
      int clr = color[i+1];
      hppRef[i]->SetMarkerStyle(20+i);
      hppRef[i]->SetMarkerSize(1.2);
      hppRef[i]->SetMarkerColor(clr);
      hppRef[i]->SetLineColor(clr);
      gppRefSys[i]->SetLineColor(clr);
      gppRefSys[i]->SetLineWidth(1);
    }
  TCanvas *c = draw1D(hppRef[1],"Differential cross-section for inclusive J/psi;p_{T} (GeV/c);B#times#frac{d^{2}#sigma}{2#pip_{T}dp_{T}dy}",true);
  hppRef[0]->Draw("sames");
  for(int i=0; i<nHistos; i++)
    {
      gppRefSys[i]->SetFillStyle(0);
      gppRefSys[i]->Draw("sameE5");
    }
  leg->Draw();

  TFile *fout = 0x0;
  if(saveHisto)  fout = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","update");
  else           fout = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","read");
  TH1F *hppRefPub = (TH1F*)fout->Get("hpp200JpsiVsPtFinal");
  TGraphAsymmErrors *gppRefPubSys = (TGraphAsymmErrors*)fout->Get("hpp200JpsiVsPtFinalSys");
  hppRefPub->SetMarkerStyle(29);
  hppRefPub->SetMarkerSize(2);
  hppRefPub->Draw("sames");
  gppRefPubSys->SetLineWidth(1);
  gppRefPubSys->SetFillStyle(0);
  gppRefPubSys->Draw("sameE5");
  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.65);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(hppRefPub, "Published STAR+PHENIX", "PL");
  leg2->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/ppRef_Run12vs15vsPub.pdf",run_type.Data(),run_cfg_name.Data()));

  if(saveHisto)
    {
      fout->cd();
      for(int i=0; i<nHistos; i++)
	{
	  hppRef[i]->Write("",TObject::kOverwrite);
	  gppRefSys[i]->Write("",TObject::kOverwrite);
	  gppRefSysTot[i]->Write("",TObject::kOverwrite);
	}
    }

  /*
    = (TH1F*)file1->Get("Runpp200_hxsec");
    TGraphAsymmErrors *gJpsiXsecSys1 = (TGraphAsymmErrors*)file1->Get("Runpp200_xsectotsys");
    hJpsiXsec1->SetMarkerStyle(20);
    TCanvas *c = draw1D(hJpsiXsec1,"Run15_pp200: J/psi invariant yield;p_{T} (GeV/c);B#times#frac{d^{2}#sigma}{2#pip_{T}dp_{T}dy}",true);
    gJpsiXsecSys1->SetMarkerColor(1);
    gJpsiXsecSys1->SetLineColor(1);
    gJpsiXsecSys1->SetFillStyle(0);
    gJpsiXsecSys1->Draw("sameE5");
    TH1F *hppRef1    = new TH1F("hpp200JpsiVsPtFinal_Run15","",nbins,xbins);
    TGraphErrors *hppRefSys1 = new TGraphErrors(nbins);
    hppRefSys1->SetName("hpp200JpsiVsPtFinalSys_Run15");
    for(int ibin=1; ibin<=nbins-1; ibin++)
    {
    hppRef1->SetBinContent(ibin, hJpsiXsec1->GetBinContent(ibin)*fraction1);
    hppRef1->SetBinError(ibin, hJpsiXsec1->GetBinError(ibin)*fraction1);
    hppRefSys1->SetPoint(ibin-1, hppRef1->GetBinCenter(ibin), hppRef1->GetBinContent(ibin));
    hppRefSys1->SetPointError(ibin-1, 0.2, gJpsiXsecSys1->GetErrorYhigh(ibin-1)*fraction1);
    }
    hppRef1->SetMarkerStyle(24);
    hppRef1->SetMarkerColor(2);
    hppRef1->SetLineColor(2);
    hppRef1->Draw("sames");
  */
  
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

