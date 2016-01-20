const char *run_type = "Run14_AuAu200";
TString run_cfg_name = "2015QM";
const double luminosity = 3.8;
const double ppInelastic = 42.; // mb
const double ppInelasticErr = 3.; // mb
const double ncoll[nCentBins] = {393., 785., 300., 95.}; 
const double ncollErr[nCentBins] = {27., 29., 31., 21.};
const double npart[nCentBins] = {161, 280, 142, 62};

//================================================
void plot_AuAu200()
{  
  rawSignal();
  //xsec();
  //upsilon();
  //publication();
  //v2();
  //efficiency();
}

//================================================
void efficiency(const bool savePlot = 0)
{
  gStyle->SetOptStat(0);
  TFile *fmuon = TFile::Open("Rootfiles/Run14.AuAu200.JpsiEff.pt1.0.pt1.0.root","read");
  TH1F *hJpsiPt[2], *hJpsiPtRebin[2];
  hJpsiPt[0] = (TH1F*)fmuon->Get("MCinput_Jpsi_pT_0060_WeightPt");
  hJpsiPt[1] = (TH1F*)fmuon->Get("MTDreco_Jpsi_pT_0060_WeightPt");
  const int nbins = 10;
  double xbins[nbins+1] = {0,1,2,3,4,5,6,7,8,9,10};
  for(int i=0; i<2; i++)
    {
      hJpsiPtRebin[i] = (TH1F*)hJpsiPt[i]->Rebin(nbins,Form("%s_rebin",hJpsiPt[i]->GetName()),xbins);
    }
  TH1F *hJpsiEffMuon = (TH1F*)hJpsiPtRebin[1]->Clone("hJpsiEffMuon");
  hJpsiEffMuon->Divide(hJpsiPtRebin[0]);
  hJpsiEffMuon->Scale(0.8/0.5);

  TFile *ftrig = TFile::Open("Rootfiles/Run14.AuAu200.JpsiTrigEff.pt1.0.pt1.0.root","read");
  TH1F *htrig = (TH1F*)ftrig->Get("JpsiTrigEff_cent0060_rebin");
  TH1F *hresp = (TH1F*)ftrig->Get("JpsiRespEff_cent0060_rebin");
  hJpsiEffMuon->Multiply(htrig);
  hJpsiEffMuon->Multiply(hresp);

  TFile *felec = TFile::Open("Rootfiles/2015QM/HT15forMrr.root","read");
  TH1F *hJpsiEffElec = (TH1F*)felec->Get("jpsiRcEff");
  hJpsiEffElec->SetMarkerStyle(20);
  hJpsiEffElec->SetMarkerColor(4);
  hJpsiEffElec->SetLineColor(4);
  hJpsiEffElec->SetMarkerSize(2);
  

  TCanvas *c1 = new TCanvas("Jpsi_eff","Jpsi_eff",800,650);
  gPad->SetLogy();
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hJpsiEffMuon,0.05,1,0.04,0.05,1.1,0.04,62);
  hJpsiEffMuon->GetYaxis()->SetRangeUser(1e-3,1);
  hJpsiEffMuon->GetYaxis()->SetNdivisions(505);
  hJpsiEffMuon->SetTitle(";p_{T} (GeV/c);J/#psi efficiency");
  hJpsiEffMuon->SetMarkerStyle(21);
  hJpsiEffMuon->SetMarkerColor(2);
  hJpsiEffMuon->SetLineColor(2);
  hJpsiEffMuon->SetMarkerSize(2);
  hJpsiEffMuon->Draw();
  hJpsiEffElec->Draw("sames");

  TLegend *leg = new TLegend(0.2,0.7,0.4,0.92);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(26);
  leg->SetHeader("AuAu @ 200 GeV");
  leg->AddEntry(hJpsiEffElec,"Electron: BEMC trigger with E_{T} > 3.4 GeV","PL");
  leg->AddEntry(hJpsiEffMuon,"Muon: MTD trigger","PL");
  leg->Draw();

  TPaveText *star = GetPaveText(0.7,0.75,0.2,0.25,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiEff.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiEff.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiEff.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiEff.eps",run_type,run_cfg_name.Data()));
    }
}

//================================================
void v2(const bool savePlot = 0)
{
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("v2_combined","v2_combined",800,500);
  TH1F *hv2Draw = new TH1F("hv2Draw",";p_{T} (GeV/c);v_{2}               ",10,0,10);
  hv2Draw->GetYaxis()->SetRangeUser(-0.2,0.25);
  ScaleHistoTitle(hv2Draw,0.06,1,0.05,0.06,1,0.05,62);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hv2Draw->DrawCopy();
  
  // different models

  // prompt jpsi Yan et al (prl 97:232301)
  float pt_yan_init[11] = {0.07, 0.49, 1.00, 1.50, 2.01, 2.50, 3.01, 3.50, 4.00, 4.51, 4.84};
  float v2_yan_init[11] = {0.0000, 0.0004, 0.0017, 0.0030, 0.0047, 0.0061, 0.0074, 0.0084, 0.0091, 0.0097, 0.0101};
  TGraph* line5 = new TGraph(11, pt_yan_init, v2_yan_init);
  line5->SetLineColor(kBlue+3);
  line5->SetLineStyle(1);
  line5->SetLineWidth(3.);
  line5->Draw("lsame");

  // Coalecence at freezeout Greco et al (plb 595:202)
  float pt_greco[24] = {0.26, 0.51, 0.74, 0.99, 1.23, 1.45, 1.65, 1.85, 2.00, 2.18, 2.34, 2.52,2.68, 2.88, 3.06, 3.27, 3.49, 3.72, 3.93, 4.16, 4.40, 4.65, 4.88, 5.00};
  float v2_greco[24] = {0.0018, 0.0018, 0.0027, 0.0044, 0.0071, 0.0169, 0.0267, 0.0391, 0.0516, 0.0658, 0.0782, 0.0924, 0.1058, 0.1173, 0.1271, 0.1378, 0.1458, 0.1547, 0.1609, 0.1671, 0.1698, 0.1716, 0.1733, 0.1733};
  TGraph* line1 = new TGraph(24, pt_greco, v2_greco);
  line1->SetLineColor(kGreen+3);
  line1->SetLineStyle(3);
  line1->SetLineWidth(3.);
  line1->Draw("lsame");

  // coalescence + initial mix Zaho et al (0806.1239[nucl-th])
  float pt_zhao[11]= {0.00, 0.51, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00};
  float v2_zhao[11]= {0.0000, 0.0035, 0.0081, 0.0151, 0.0249, 0.0350, 0.0417, 0.0434, 0.0427, 0.0396, 0.0361};
  float er_zhao[11];
  TGraph* line4 = new TGraph(11, pt_zhao, v2_zhao);
  line4->SetLineColor(kBlue-4);
  line4->SetLineStyle(5);
  line4->SetLineWidth(3.);
  line4->Draw("lsame");

  double pt_transModelMix[51] = {0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8,  2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, 6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8, 9., 9.2, 9.4, 9.6, 9.8, 10.};
  double v2_transModelMix[51] = {0.0000, 0.0312, 0.1044, 0.2510, 0.4609, 0.6470, 0.8453, 1.0247, 1.1745, 1.2634, 1.3348, 1.3585, 1.3664, 1.3536, 1.3316, 1.2945, 1.2647, 1.2231, 1.1903, 1.1542, 1.1315, 1.1048, 1.0800, 1.0566, 1.0318, 1.0147, 0.9990, 0.9771, 0.9559, 0.9435, 0.9232, 0.9013, 0.8913, 0.8724, 0.8541, 0.8448, 0.8256, 0.8157, 0.7974, 0.7884, 0.7692, 0.7594, 0.7397, 0.7300, 0.7195, 0.7099, 0.6904, 0.6809, 0.6707, 0.6612, 0.6505};
  for(int i=0; i<51; i++)
    v2_transModelMix[i] *= 0.01;
  TGraph* line8 = new TGraph(51, pt_transModelMix, v2_transModelMix);
  line8->SetLineColor(35);
  line8->SetLineStyle(10);
  line8->SetLineWidth(3.);
  line8->Draw("lsame");

  ifstream ifs;
  string tempstr;
  float pt_Ulrich_visHydro_T120_00_100[80], v2_Ulrich_visHydro_T120_00_100[80];
  ifs.open("v2/JpsiSpv2020T120Tau6-00100.dat");
  for(int i=0; i<80; i++)
    {
      ifs>>pt_Ulrich_visHydro_T120_00_100[i]>>tempstr>>v2_Ulrich_visHydro_T120_00_100[i]>>tempstr;
    }
  ifs.close();
  TGraph* line9 = new TGraph(40, pt_Ulrich_visHydro_T120_00_100, v2_Ulrich_visHydro_T120_00_100);
  line9->SetLineColor(kPink+9);
  line9->SetLineStyle(9);
  line9->SetLineWidth(3.);
  line9->Draw("lsame");

  // data point
  double v2_pt_cent[5][4];
  double v2StaError_pt_cent[5][4];
  double v2SysErrorLow_pt_cent[5][4];
  double v2SysErrorHigh_pt_cent[5][4];
  double pt[4];
  double pts[6][4];
  double ptSmear[6] = {-0.2, 0., 0.2, -0.1, 0.1, 0.3};
  for(int i=0; i<6; i++)
    for(int j=0; j<4; j++)
        pts[i][j] = pt[j]+ptSmear[i];
  double ptError[4] = {0.014, 0.022, 0.094, 0.1};
  double nonFlow[5][4];

  char fileName[256];
  char aa[256];

  ifstream ifs;
  ifs.open("Rootfiles/Published/Jpsi_v2_200/v2Results_4ptBins.dat");
  for(int i=0; i<5; i++)
    {
      ifs.getline(aa, 256);
      ifs.getline(aa, 256);
      ifs.getline(aa, 256);
      for(int j=0; j<4; j++)
	{
	  ifs>>pt[j]>>tempstr>>ptError[j]>>v2_pt_cent[i][j]>>tempstr>>v2StaError_pt_cent[i][j]>>tempstr>>v2SysErrorHigh_pt_cent[i][j]>>tempstr>>v2SysErrorLow_pt_cent[i][j]>>tempstr>>nonFlow[i][j];
	  //cout<<pt[j]<<"  "<<ptError[j]<<"   "<<v2_pt_cent[i][j]<<"   "<<v2StaError_pt_cent[i][j]<<"  "<<v2SysErrorHigh_pt_cent[i][j]<<"   "<<v2SysErrorLow_pt_cent[i][j]<<"   "<<nonFlow[i][j]<<endl;
	}
      ifs.getline(aa, 256);
    }
  ifs.close();

 
  
  v2_pt_cent[4][0] = 0.05; v2_pt_cent[4][1] = 0.01;  v2_pt_cent[4][2] = 0.05;  v2_pt_cent[4][3] = 0.08;  
  v2StaError_pt_cent[4][0] = 0.03; v2StaError_pt_cent[4][1] = 0.04; v2StaError_pt_cent[4][2] = 0.04; v2StaError_pt_cent[4][3] = 0.05;
  v2SysErrorLow_pt_cent[4][0] = 0.02;   v2SysErrorLow_pt_cent[4][1] = 0.01;  v2SysErrorLow_pt_cent[4][2] = 0.02;  v2SysErrorLow_pt_cent[4][3] = 0.01; 
  

  
  v2_pt_cent[4][0] = 0.08; v2_pt_cent[4][1] = 0.015; v2_pt_cent[4][2] = 0.06; v2_pt_cent[4][3] = 0.10459; 
  v2StaError_pt_cent[4][0] = 0.0212132; v2StaError_pt_cent[4][1] = 0.0282843; v2StaError_pt_cent[4][2] = 0.0282843; v2StaError_pt_cent[4][3] = 0.0384111;
  v2SysErrorLow_pt_cent[4][0] = sqrt(v2SysErrorLow_pt_cent[4][0]*v2SysErrorLow_pt_cent[4][0]+0.04*0.04);
  v2SysErrorLow_pt_cent[4][1] = sqrt(v2SysErrorLow_pt_cent[4][1]*v2SysErrorLow_pt_cent[4][1]+0.02*0.02);
  v2SysErrorLow_pt_cent[4][2] = sqrt(v2SysErrorLow_pt_cent[4][2]*v2SysErrorLow_pt_cent[4][2]+0.02*0.02);
  v2SysErrorLow_pt_cent[4][3] = sqrt(v2SysErrorLow_pt_cent[4][3]*v2SysErrorLow_pt_cent[4][3]+0.07*0.07);
  

  for(int i=0; i<4; i++)
    {
      v2SysErrorHigh_pt_cent[4][i] = v2SysErrorLow_pt_cent[4][i];
    }
  
 
  printf("%10s %10s %10s %10s %10s %10s %10s\n","pt","v2","+stat.","-stat.","+sys.","-sys.","non-flow");
  for(int j=0; j<4; j++)
    {
      printf("%10.2f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",pt[j],v2_pt_cent[4][j],v2StaError_pt_cent[4][j],-1*v2StaError_pt_cent[4][j],v2SysErrorLow_pt_cent[4][j],-1*v2SysErrorLow_pt_cent[4][j], -1*nonFlow[3][j]);
    }

  gr1 = new TGraphErrors(4, pt, v2_pt_cent[4], ptError, v2StaError_pt_cent[4]);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(1);
  gr1->SetMarkerSize(1.5);
  gr1->SetLineColor(1);
  gr1->SetLineWidth(2.);

  TBox *box[4];
  for(int i=0; i<4; i++)
    {
      box[i] = new TBox(pt[i]-0.1, v2_pt_cent[4][i], pt[i]+0.1, v2_pt_cent[4][i]-nonFlow[3][i]);
      box[i]->SetLineColor(kBlack);
      box[i]->SetFillColor(kBlack);
      box[i]->SetLineWidth(2.);
      box[i]->SetFillStyle(0);
      box[i]->Draw("fsame");
    }
  TBox *boxLabel = new TBox(1, -0.05, 1.5, -0.03);
  boxLabel->SetLineColor(kBlack);
  boxLabel->SetFillColor(kBlack);
  boxLabel->SetLineWidth(2.);
  boxLabel->SetFillStyle(0);
  boxLabel->Draw("fsame");
  TLatex* lat3 = new TLatex(1.8, -0.045, "maximum non-flow");
  lat3->SetTextSize(0.045);
  lat3->Draw("same");
  TGraph* bracketUp[4];
  TGraph* bracketDown[4];
  
  for(int i=0; i<4; i++)
    {
      float bracketUpX[4] = {pt[i]-0.15, pt[i]-0.1, pt[i]+0.1, pt[i]+0.15};
      float v2Up = v2_pt_cent[4][i]+v2SysErrorHigh_pt_cent[4][i];
      float bracketUpY[4] = {v2Up-0.004, v2Up, v2Up, v2Up-0.004};
      bracketUp[i] = new TGraph(4, bracketUpX, bracketUpY);
      bracketUp[i]->SetLineColor(1);
      bracketUp[i]->SetLineWidth(2.);
      bracketUp[i]->Draw("lsame");
      float bracketDownX[4] = {pt[i]-0.15, pt[i]-0.1, pt[i]+0.1, pt[i]+0.15};
      float v2Down = v2_pt_cent[4][i]-v2SysErrorLow_pt_cent[4][i];
      float bracketDownY[4] = {v2Down+0.004, v2Down, v2Down, v2Down+0.004};
      bracketDown[i] = new TGraph(4, bracketDownX, bracketDownY);
      bracketDown[i]->SetLineColor(1);
      bracketDown[i]->SetLineWidth(2.);
      bracketDown[i]->Draw("lsame");
    }

  gr1->Draw("Psame");
  TLegend* leg1 = new TLegend(0.16, 0.14, 0.84, 0.38);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);

  leg1->AddEntry(line5, "initially produced", "l");
  leg1->AddEntry(line1, "coalescence from thermalized c#bar{c}", "l");
  leg1->AddEntry(line4, "initial + coalescence", "l");
  leg1->AddEntry(line8, "initial + coalescence", "l");
  leg1->AddEntry(line9, "hydrodynamic", "l");
  leg1->SetTextSize(0.045);
  leg1->Draw("same");


  TPaveText *t1 = GetPaveText(0.2,0.4,0.9,0.95,0.05,62);
  t1->AddText("Au+Au 200 GeV 0-80 %");
  t1->Draw();

  TLegend* leg1 = new TLegend(0.5, 0.9, 0.7, 0.95);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.05);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr1,"STAR J/#psi#rightarrowe^{+}e^{-} Run10+11","P");
  leg1->Draw();

  TPaveText *star = GetPaveText(0.75,0.8,0.4,0.5,0.05);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_Run10+11.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_Run10+11.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_Run10+11.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_Run10+11.eps",run_type,run_cfg_name.Data()));
    }


  // add MTD result
  double pt_mtd[3] = {1.03, 3.2, 6.49};
  double ptError_mtd[3] = {0,0,0};
  double v2_mtd[3] = {-0.122102, 0.00276895, 0.127219};
  double v2StatErr_mtd[3] = {0.156197, 0.134758, 0.10168};
  double v2SysErr_mtd[3] = {0.0557383, 0.0115725, 0.114443};
  double ptSysErr_mtd[3] = {0.2,0.2,0.2};

  TCanvas *c2 = new TCanvas("v2_MTD","v2_MTD",800,500);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hv2Draw->GetYaxis()->SetRangeUser(-0.4,0.4);
  hv2Draw->Draw();
  for(int i=0; i<4; i++)
    {
      box[i]->Draw("fsame");
    }
  TBox *boxLabel = new TBox(3.5, -0.33, 3.8, -0.3);
  boxLabel->SetLineColor(kBlack);
  boxLabel->SetFillColor(kBlack);
  boxLabel->SetLineWidth(2.);
  boxLabel->SetFillStyle(0);
  boxLabel->Draw("fsame");
  TLatex* lat3 = new TLatex(3.9, -0.325, "maximum non-flow");
  lat3->SetTextSize(0.045);
  lat3->Draw("same");
  for(int i=0; i<4; i++)
    {
      bracketUp[i]->Draw("lsame");
      bracketDown[i]->Draw("lsame");
    }
  gr1->Draw("Psame");

  gr3 = new TGraphErrors(3, pt_mtd, v2_mtd, ptSysErr_mtd, v2SysErr_mtd);
  gr3->SetMarkerStyle(29);
  gr3->SetMarkerColor(4);
  gr3->SetMarkerSize(2);
  gr3->SetLineColor(4);
  gr3->SetLineWidth(2.);
  gr3->SetFillStyle(0);
  gr3->Draw("samesE5");

  gr2 = new TGraphErrors(3, pt_mtd, v2_mtd, ptError_mtd, v2StatErr_mtd);
  gr2->SetMarkerStyle(29);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerSize(2.5);
  gr2->SetLineColor(4);
  gr2->SetLineWidth(2.);
  gr2->Draw("samesPEZ");

  TLegend* leg1 = new TLegend(0.4, 0.25, 0.6, 0.4);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.05);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr2,"J/#psi#rightarrow#mu^{+}#mu^{-} |y| < 0.5, Run14","P");
  leg1->AddEntry(gr1,"J/#psi#rightarrowe^{+}e^{-} |y| < 1, Run10+11","P");
  leg1->Draw();

  TPaveText *t1 = GetPaveText(0.2,0.4,0.9,0.95,0.05,62);
  t1->AddText("Au+Au 200 GeV 0-80 %");
  t1->Draw();

  TPaveText *star = GetPaveText(0.75,0.8,0.9,0.95,0.05);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_uu_Run14.eps",run_type,run_cfg_name.Data()));
    }


  TCanvas *c1 = new TCanvas("v2_combined_vs_LHC","v2_combined_vs_LHC",800,500);
  TH1F *hv2Draw = new TH1F("hv2Draw2",";p_{T} (GeV/c);v_{2}               ",10,0,15);
  hv2Draw->GetYaxis()->SetRangeUser(-0.1,0.4);
  ScaleHistoTitle(hv2Draw,0.06,1,0.05,0.06,1,0.05,62);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hv2Draw->DrawCopy();

  // coalescence + initial mix Zhao et al (0806.1239[nucl-th])
  TGraph* line4 = new TGraph(11, pt_zhao, v2_zhao);
  line4->SetLineColor(2);
  line4->SetLineStyle(1);
  line4->SetLineWidth(2);
  line4->Draw("lsame");

  // Liu
  TGraph* line8 = new TGraph(51, pt_transModelMix, v2_transModelMix);
  line8->SetLineColor(2);
  line8->SetLineStyle(2);
  line8->SetLineWidth(2);
  line8->Draw("lsame");

  TGraphErrors *gr1 = new TGraphErrors(4, pt, v2_pt_cent[4], ptError, v2StaError_pt_cent[4]);
  gr1->SetMarkerStyle(29);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerSize(3);
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);

  TBox *box[4];
  for(int i=0; i<4; i++)
    {
      box[i] = new TBox(pt[i]-0.1, v2_pt_cent[4][i], pt[i]+0.1, v2_pt_cent[4][i]-nonFlow[3][i]);
      box[i]->SetLineColor(2);
      box[i]->SetFillColor(2);
      box[i]->SetLineWidth(2.);
      box[i]->SetFillStyle(0);
      box[i]->Draw("fsame");
    }
  TBox *boxLabel = new TBox(1, -0.06, 1.5, -0.04);
  boxLabel->SetLineColor(2);
  boxLabel->SetFillColor(2);
  boxLabel->SetLineWidth(2.);
  boxLabel->SetFillStyle(0);
  boxLabel->Draw("fsame");
  TLatex* lat3 = new TLatex(1.8, -0.055, "STAR maximum non-flow");
  lat3->SetTextSize(0.045);
  lat3->Draw("same");
  TGraph* bracketUp[4];
  TGraph* bracketDown[4];
  
  for(int i=0; i<4; i++)
    {
      float bracketUpX[4] = {pt[i]-0.15, pt[i]-0.1, pt[i]+0.1, pt[i]+0.15};
      float v2Up = v2_pt_cent[4][i]+v2SysErrorHigh_pt_cent[4][i];
      float bracketUpY[4] = {v2Up-0.004, v2Up, v2Up, v2Up-0.004};
      bracketUp[i] = new TGraph(4, bracketUpX, bracketUpY);
      bracketUp[i]->SetLineColor(2);
      bracketUp[i]->SetLineWidth(2.);
      bracketUp[i]->Draw("lsame");
      float bracketDownX[4] = {pt[i]-0.15, pt[i]-0.1, pt[i]+0.1, pt[i]+0.15};
      float v2Down = v2_pt_cent[4][i]-v2SysErrorLow_pt_cent[4][i];
      float bracketDownY[4] = {v2Down+0.004, v2Down, v2Down, v2Down+0.004};
      bracketDown[i] = new TGraph(4, bracketDownX, bracketDownY);
      bracketDown[i]->SetLineColor(2);
      bracketDown[i]->SetLineWidth(2.);
      bracketDown[i]->Draw("lsame");
    }

  gr1->Draw("Psame");
  TLegend *leg = new TLegend(0.16,0.8,0.4,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);;
  leg->AddEntry(gr1,"STAR, 200 GeV Au+Au, |y| < 1","P");
  leg->AddEntry(line8,"#it{Liu et al.}, RHIC","L");
  leg->AddEntry(line4,"#it{Zhao et al.}, RHIC","L");
  leg->Draw();

  TPaveText *star = GetPaveText(0.8,0.85,0.9,0.95,0.05);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee.eps",run_type,run_cfg_name.Data()));
    }

  double PtCMS24[]={7.17,8.86,13.18};
  double v2CMS24[]={0.045,0.081,0.030};
  double dv2CMS24[]={0.0227,0.0221,0.0207};
  double dv21CMS24[]={0.0100,0.00936,0.00668};
  TGraphErrors *grCMS24 = new TGraphErrors(3,PtCMS24, v2CMS24, 0, dv2CMS24);
  grCMS24->SetMarkerStyle(20);
  grCMS24->SetMarkerSize(2);
  grCMS24->SetMarkerColor(4);
  grCMS24->SetLineColor(grCMS24->GetMarkerColor());
  grCMS24->SetLineWidth(1);
  double PtCMS16[]={4.5,8.7};
  double v2CMS16[]={0.051,0.065};
  double dv2CMS16[]={0.021,0.023};
  double dv21CMS16[]={0.023,0.0087};
  TGraphErrors *grCMS16 = new TGraphErrors(2,PtCMS16, v2CMS16, 0, dv2CMS16);
  grCMS16->SetMarkerStyle(24);
  grCMS16->SetMarkerSize(2);
  grCMS16->SetMarkerColor(4);
  grCMS16->SetLineColor(grCMS16->GetMarkerColor());
  grCMS16->SetLineWidth(1);

  TGraph *grZhuang = new TGraph("Rootfiles/Published/Jpsi_Raa_200/model/Jpsi_v2_Zhuang_LHC_B_Thermalized.dat");
  grZhuang->SetLineColor(grCMS16->GetLineColor());
  grZhuang->SetLineWidth(2);
  grZhuang->SetLineStyle(2);
  grZhuang->Draw("c");

  TGraph *grZhuang2 = new TGraph("Rootfiles/Published/Jpsi_Raa_200/model/Jpsi_v2_Zhuang_LHC_B_NotThermalized.dat");
  grZhuang2->SetLineColor(grZhuang->GetLineColor());
  grZhuang2->SetLineWidth(2);
  grZhuang2->SetLineStyle(3);
  grZhuang2->Draw("c");


  grCMS24->Draw("psame");
  drawBoxes(grCMS24->GetN(), PtCMS24, 0.2, v2CMS24, dv21CMS24,2,grCMS24->GetLineColor());
  grCMS16->Draw("psame");
  drawBoxes(grCMS16->GetN(), PtCMS16, 0.2, v2CMS16, dv21CMS16,2,grCMS16->GetLineColor());

  TLegend *leg = new TLegend(0.16,0.6,0.4,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);;
  leg->AddEntry(grCMS24,"CMS Preliminary, 2.76 TeV Pb+Pb, |y|<2.4","P");
  leg->AddEntry(grCMS16,"CMS Preliminary, 2.76 TeV Pb+Pb, 1.6<|y|<2.4","P");
  leg->AddEntry(grZhuang,"#it{Liu et al.}, LHC, B thermalized","L");
  leg->AddEntry(grZhuang2,"#it{Liu et al.}, LHC, B not thermalized","L");
  leg->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_vs_CMS.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_vs_CMS.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_vs_CMS.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/JpsiV2_ee_vs_CMS.eps",run_type,run_cfg_name.Data()));
    }
}

//================================================
void xsec(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  const int marker_style[nCentBins] = {kFullCircle, kFullStar, kFullSquare, kFullCross};
  const double marker_size[nCentBins] = {1.5,2,1.5,2};
  const int marker_color[nCentBins] = {1,2,4,6};
  const int pub_marker_style[nCentBins] = {kOpenCircle, kOpenStar, kOpenSquare, kOpenCross};
  const int pub_marker_color_low[nCentBins] = {kGray+2, kRed+2, kBlue+2, kMagenta+2};
  const int pub_marker_color_high[nCentBins] = {1,2,4,6};
  const double scale_factor[nCentBins] = {10,1,0.2,0.1};
  const int lowpt_color = 9;
  const int highpt_color = 12;

  // published data
  TFile *fpub = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","read");
  TGraphAsymmErrors *gAuAuLowPt[nCentBins];
  TGraphAsymmErrors *gAuAuLowPtSys[nCentBins];
  TGraphAsymmErrors *gAuAuHighPt[nCentBins];
  TGraphAsymmErrors *gAuAuHighPtSys[nCentBins];
  TH1F *hTBW[nCentBins];
  TGraphAsymmErrors *gRaaLowPt[nCentBins];
  TGraphAsymmErrors *gRaaLowPtSys[nCentBins];
  TGraphAsymmErrors *gRaaHighPt[nCentBins];
  TGraphAsymmErrors *gRaaHighPtSys[nCentBins];
  TH1F *hJpsipp = (TH1F*)fpub->Get("Jpsi_xsec_pp200_combined");
  double x,y;
  for(int k=0; k<nCentBins; k++)
    {
      gAuAuLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[k]));
      gAuAuLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[k]));
      gAuAuHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[k]));
      gAuAuHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[k]));
      hTBW[k] = (TH1F*)fpub->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[k]));
      scaleGraph(gAuAuLowPt[k], scale_factor[k]);
      scaleGraph(gAuAuLowPtSys[k], scale_factor[k]);
      scaleGraph(gAuAuHighPt[k], scale_factor[k]);
      scaleGraph(gAuAuHighPtSys[k], scale_factor[k]);
      hTBW[k]->Scale(scale_factor[k]);

      gRaaLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",cent_Title[k]));
      gRaaLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",cent_Title[k]));
      gRaaHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",cent_Title[k]));
      gRaaHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",cent_Title[k]));
    }

  // data
  TFile *fSys = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hSys = (TH1F*)fSys->Get(Form("Sys_all_%s",cent_Title[0]));
  TFile *fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiInvYield[nCentBins];
  TGraphErrors *hJpsiXsec[nCentBins];
  TGraphErrors *hJpsiXsecSys[nCentBins];
  TGraphErrors *hJpsiRaa[nCentBins];
  TGraphErrors *hJpsiRaaSys[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)fdata->Get(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiXsec[k] = new TGraphErrors(hJpsiInvYield[k]);
      hJpsiXsec[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiXsecSys[k] = new TGraphErrors(hJpsiInvYield[k]);
      hJpsiXsecSys[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));
      for(int i=0; i<hJpsiXsec[k]->GetN(); i++)
	{
	  hJpsiXsec[k]->GetPoint(i,x,y);
	  hJpsiXsec[k]->SetPoint(i,x,y*scale_factor[k]);
	  hJpsiXsec[k]->SetPointError(i,0,hJpsiXsec[k]->GetErrorY(i)*scale_factor[k]);

	  hJpsiXsecSys[k]->SetPoint(i,x,y*scale_factor[k]);
	  hJpsiXsecSys[k]->SetPointError(i,0.2,hSys->GetBinContent(i+1)*y*scale_factor[k]);
	}

      hJpsiRaa[k] = new TGraphErrors();
      hJpsiRaa[k]->SetName(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
      hJpsiRaaSys[k] = new TGraphErrors();
      hJpsiRaaSys[k]->SetName(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
      for(int i=0; i<hJpsiXsec[k]->GetN(); i++)
	{
	  hJpsiXsec[k]->GetPoint(i,x,y);
	  double AuAu_val = y/scale_factor[k];
	  double AuAu_err = hJpsiXsec[k]->GetErrorY(i)/scale_factor[k];
	  double pp_val = hJpsipp->GetBinContent(i+1);
	  double pp_err = hJpsipp->GetBinError(i+1);
	  double prefix = ppInelastic/ncoll[k] * 1e6;
	  double val = prefix * AuAu_val / pp_val;
	  double err = val * sqrt(AuAu_err*AuAu_err/AuAu_val/AuAu_val + pp_err*pp_err/pp_val/pp_val);
	  hJpsiRaa[k]->SetPoint(i,x,val);
	  hJpsiRaa[k]->SetPointError(i,0,err);
	  printf("%s: raa -> %2.2f +/- %2.2f%%, AuAu -> %2.2e +/- %2.2f%%, pp -> %2.2e +/-%2.2f%%\n",pt_Name[i+1],val,err/val*100,AuAu_val,AuAu_err/AuAu_val*100,pp_val,pp_err/pp_val*100);

	  hJpsiRaaSys[k]->SetPoint(i,x,val);
	  hJpsiRaaSys[k]->SetPointError(i,0.2,val*hSys->GetBinContent(i+1));
	}
    }



  TCanvas *c1 = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",800,700);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);B_{ll}d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{-2}]",10,0,10);
  hAuAu->GetYaxis()->SetRangeUser(1e-11,5e-3);
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
  leg->AddEntry(hJpsiXsec[0],"0-60%#times10","P");
  leg->AddEntry(hJpsiXsec[1],"0-20%","P");
  leg->AddEntry(hJpsiXsec[2],"20-40%/5","P");
  leg->AddEntry(hJpsiXsec[3],"40-60%/10","P");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.15,0.3,0.85,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  t1->Draw();
  TPaveText *star = GetPaveText(0.25,0.5,0.25,0.3,0.04);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c2 = new TCanvas("AuAu200_Jpsi_vs_pub","AuAu200_Jpsi_vs_pub",800,700);
  gPad->SetLogy();
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hAuAu->Draw();
  TGraphErrors *gJpsi[nCentBins], *gJpsiSys[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      gAuAuLowPt[k]->SetMarkerStyle(pub_marker_style[k]);
      gAuAuLowPt[k]->SetMarkerSize(marker_size[k]);
      gAuAuLowPt[k]->SetMarkerColor(pub_marker_color_low[k]);
      gAuAuLowPt[k]->SetLineColor(pub_marker_color_low[k]);
      gAuAuLowPt[k]->Draw("sames PE");
      gAuAuLowPtSys[k]->SetMarkerColor(pub_marker_color_low[k]);
      gAuAuLowPtSys[k]->SetLineColor(pub_marker_color_low[k]);
      gAuAuLowPtSys[k]->Draw("sameE5");
      gAuAuHighPt[k]->SetMarkerStyle(pub_marker_style[k]);
      gAuAuHighPt[k]->SetMarkerSize(marker_size[k]);
      gAuAuHighPt[k]->SetMarkerColor(pub_marker_color_high[k]);
      gAuAuHighPt[k]->SetLineColor(pub_marker_color_high[k]);
      gAuAuHighPt[k]->Draw("sames PE");
      gAuAuHighPtSys[k]->SetMarkerColor(pub_marker_color_high[k]);
      gAuAuHighPtSys[k]->SetLineColor(pub_marker_color_high[k]);
      gAuAuHighPtSys[k]->Draw("sameE5");
      hTBW[k]->Draw("sames");

      // gJpsi[k] = new TGraphErrors(*hJpsiXsec[k]);
      // gJpsiSys[k] = new TGraphErrors(*hJpsiXsecSys[k]);
      // gJpsi[k]->SetMarkerStyle(marker_style[k]);
      // gJpsi[k]->SetMarkerColor(2);
      // gJpsi[k]->SetLineColor(2);
      // gJpsi[k]->SetMarkerSize(marker_size[k]+0.5);
      // gJpsiSys[k]->SetFillStyle(0);
      // gJpsiSys[k]->SetLineColor(gJpsi[k]->GetLineColor());
      hJpsiXsecSys[k]->Draw("samesE5");
      hJpsiXsec[k]->Draw("samesP");
    }
  TLegend *leg = new TLegend(0.7,0.65,0.9,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hJpsiXsec[0],"0-60%#times10","P");
  leg->AddEntry(hJpsiXsec[1],"0-20%","P");
  leg->AddEntry(hJpsiXsec[2],"20-40%/5","P");
  leg->AddEntry(hJpsiXsec[3],"40-60%/10","P");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.15,0.3,0.9,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV"));
  t1->Draw();

  /*
  TLegend *leg = new TLegend(0.14,0.29,0.3,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gAuAuHighPt[0],"star","");
  leg->Draw();
  TLegend *leg = new TLegend(0.17,0.17,0.3,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gAuAuLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y| < 1","P");
  leg->AddEntry(gJpsi[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg->AddEntry(hTBW[0],"TBW fit (#beta=0)","L");
  leg->Draw();
  */

  TPaveText *t1 = GetPaveText(0.29,0.39,0.3,0.35,24,63);
  t1->AddText("Open: J/#psi#rightarrowe^{+}e^{-}, |y| < 1");
  t1->Draw();

  TPaveText *t1 = GetPaveText(0.31,0.41,0.25,0.3,24,63);
  t1->AddText("Filled: J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  t1->Draw();

  TLegend *leg = new TLegend(0.16,0.2,0.3,0.25);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(24);
  leg->SetTextFont(63);
  leg->AddEntry(hTBW[0],"TBW fit (#beta=0)","L");
  leg->Draw();
  

  TPaveText *star = GetPaveText(0.4,0.6,0.8,0.85,0.04);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiInvYield_vs_pub.eps",run_type,run_cfg_name.Data()));
    }

  // RAA plot

  TCanvas *c3 = new TCanvas("Raa_Jpsi_vs_pub","Raa_Jpsi_vs_pub",1100,700);
  //  c->Divide(2,2,0,0);
  TPad *pad = GetSinglePad("jet_trigger",0.02,0.98,0.02,0.98);
  pad->Divide(2,2,0,0);
  pad->SetTickx();
  pad->SetTicky();
  pad->Draw();
  //c3->Divide(2,2);

  TH1F *hRaa = new TH1F("Raa_Jpsi",";p_{T} (GeV/c);R_{AA}",10,-0.1,9.8);
  hRaa->GetYaxis()->SetRangeUser(0.05,1.95);
  hRaa->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hRaa,24,1.9,20,24,1.9,20,63);
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      if(k%2==1) gPad->SetRightMargin(0.003);
      else gPad->SetLeftMargin(0.13);
      if(k>1) gPad->SetBottomMargin(0.15);
      hRaa->DrawCopy();
      TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
      line->Draw();
      gRaaLowPt[k]->SetMarkerStyle(kOpenStar);
      gRaaLowPt[k]->SetMarkerSize(2);
      gRaaLowPt[k]->SetMarkerColor(lowpt_color);
      gRaaLowPt[k]->SetLineColor(lowpt_color);
      gRaaLowPtSys[k]->SetMarkerColor(lowpt_color);
      gRaaLowPtSys[k]->SetLineColor(lowpt_color);

      gRaaHighPt[k]->SetMarkerStyle(kOpenStar);
      gRaaHighPt[k]->SetMarkerSize(2);
      gRaaHighPt[k]->SetMarkerColor(highpt_color);
      gRaaHighPt[k]->SetLineColor(highpt_color);
      gRaaHighPtSys[k]->SetMarkerColor(highpt_color);
      gRaaHighPtSys[k]->SetLineColor(highpt_color);

      hJpsiRaa[k]->SetMarkerStyle(29);
      hJpsiRaa[k]->SetMarkerColor(2);
      hJpsiRaa[k]->SetLineColor(2);
      hJpsiRaa[k]->SetMarkerSize(2.5);
      hJpsiRaaSys[k]->SetFillStyle(0);
      hJpsiRaaSys[k]->SetLineColor(hJpsiRaa[k]->GetLineColor());

      gRaaLowPtSys[k]->Draw("sameE5");
      gRaaHighPtSys[k]->Draw("sameE5");
      gRaaLowPt[k]->Draw("sames PEZ");
      gRaaHighPt[k]->Draw("sames PEZ");
      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");

      if(k%2==0) t1 = GetPaveText(0.22,0.25,0.9,0.95);
      else       t1 = GetPaveText(0.1,0.12,0.9,0.95);
      t1->SetTextFont(63);
      t1->SetTextSize(24);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();

      // Global systematics
      double gerr = sqrt(ppInelasticErr*ppInelasticErr/ppInelastic/ppInelastic+ncollErr[k]*ncollErr[k]/ncoll[k]/ncoll[k]);
      cout << gerr << endl;
      TBox *box = new TBox(0.1,1-gerr,0.2,1+gerr);
      box->SetFillStyle(3001);
      box->SetFillColor(kBlack);
      box->Draw();
    }
  pad->cd(1);

  /*
  TLegend *leg = new TLegend(0.36,0.65,0.55,0.97);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);;
  leg->AddEntry(gRaaHighPt[0],"J/#psi#rightarrowe^{+}e^{-}","P");
  leg->Draw();
  TLegend *leg = new TLegend(0.4,0.65,0.6,0.97);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);
  leg->SetHeader("Au+Au @ 200 GeV");
  leg->AddEntry(gRaaLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y| < 1","P");
  leg->AddEntry(hJpsiRaa[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg->Draw();
  */

  TPaveText *t1 = GetPaveText(0.596,0.696,0.9,0.95,20,63);
  t1->AddText("STAR Au+Au @ 200 GeV");
  t1->Draw();

  TPaveText *t1 = GetPaveText(0.58,0.68,0.78,0.88,20,63);
  t1->AddText("Open: J/#psi#rightarrowe^{+}e^{-}, |y| < 1");
  t1->Draw();

  TPaveText *t1 = GetPaveText(0.6,0.7,0.68,0.78,20,63);
  t1->AddText("Filled: J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  t1->SetTextColor(2);
  t1->Draw();


  pad->cd(2);
  TPaveText *star = GetPaveText(0.6,0.8,0.85,0.9,24,23);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_pub.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_pub.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_pub.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_pub.eps",run_type,run_cfg_name.Data()));
    }

  ifstream ifin;
  char tmp[256];
  double dtmp;

  double phenix_pt[4][5];
  double phenix_pt_err[4][5];
  double phenix_pt_sys[4][5];
  double phenix_raa[4][5];
  double phenix_raa_stat[4][5];
  double phenix_raa_sys[4][5];
  const char *phenix_name[4] = {"060","020","2040","4060"};
  TGraphErrors *gPhenixRaa[4];
  TGraphErrors *gPhenixRaaSys[4];
  for(int k=0; k<4; k++)
    {
      ifin.open(Form("model/raa_phenix_%s.dat",phenix_name[k]));
      for(int j=0; j<5; j++)
	{
	  ifin >> dtmp >> dtmp >> dtmp >> dtmp >> phenix_pt[k][j] >> dtmp >> phenix_raa[k][j] >> phenix_raa_stat[k][j] >> phenix_raa_sys[k][j] >> dtmp;
	  phenix_pt_err[k][j] = 0;
	  phenix_pt_sys[k][j] = 0.2;
	}
      ifin.close();
      gPhenixRaa[k] = new TGraphErrors(5,phenix_pt[k], phenix_raa[k],phenix_pt_err[k],phenix_raa_stat[k]);
      gPhenixRaaSys[k] = new TGraphErrors(5,phenix_pt[k], phenix_raa[k],phenix_pt_sys[k],phenix_raa_sys[k]);
    }
  
  for(int k=0; k<nCentBins; k++)
    {
      pad->cd(k+1);
      gPhenixRaa[k]->SetMarkerStyle(20);
      gPhenixRaa[k]->SetMarkerSize(2);
      gPhenixRaa[k]->SetMarkerColor(kGreen+2);
      gPhenixRaa[k]->SetLineColor(kGreen+2);
      gPhenixRaaSys[k]->SetMarkerColor(kGreen+2);
      gPhenixRaaSys[k]->SetLineColor(kGreen+2);
      gPhenixRaaSys[k]->SetFillStyle(0);

      gPhenixRaaSys[k]->Draw("sameE5");
      gPhenixRaa[k]->Draw("sames PEZ");

      hJpsiRaaSys[k]->Draw("sameE5");
      hJpsiRaa[k]->Draw("samesPEZ");
    }
  pad->cd(2);
  TLegend *leg = new TLegend(0.05,0.6,0.3,0.7);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);
  leg->AddEntry(gPhenixRaa[0],"PHENIX: J/#psi#rightarrowe^{+}e^{-}, |y| < 0.35","P");
  leg->Draw();
  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_PHENIX.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_PHENIX.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_PHENIX.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_PHENIX.eps",run_type,run_cfg_name.Data()));
    }

  // Raa vs Npart
  TCanvas *c3 = new TCanvas("Jpsi_raa_vs_npart","Jpsi_raa_vs_npart",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  TH1F *hRaaVsNpart = new TH1F("hRaaVsNpart",";N_{part};R_{AA}                   ",100,-1,370);
  ScaleHistoTitle(hRaaVsNpart,28,1,24,28,1,24,63);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,2);
  hRaaVsNpart->DrawCopy();
  double pub_npart[5] = {325, 235.6, 167.7, 115.9, 62.44};
  double pub_npart_err[5] = {4, 8.8, 10.6, 11.1, 10.0};
  double pub_ncoll[5] = {962.2, 608.6, 377.3, 223.8, 91.33};
  double pub_ncoll_err[5] = {27.6, 31.1, 33.3, 30.5, 20.0};
  double pub_npart_xerr[5];
  double pub_npart_xsys[5];
  double pub_npart_relerr[5];
  double pub_npart_y[5];
  for(int i=0; i<5; i++)
    {
      pub_npart_xerr[i] = 0;
      pub_npart_xsys[i] = 6;
      pub_npart_y[i] = 1;
      pub_npart_relerr[i] = pub_ncoll_err[i]/pub_ncoll[i];
    }
  TGraphErrors *gPubNpartErr = new TGraphErrors(5, pub_npart, pub_npart_y, pub_npart_err, pub_npart_relerr);
  gPubNpartErr->SetFillColor(kGreen-10);
  gPubNpartErr->SetLineColor(kGreen-10);
  gPubNpartErr->Draw("e3sames");
  TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  // tsinghua group
  ifin.open("model/data_and_figures_from_Yunpeng/data_pt_5-10_GeV.txt");
  ifin.getline(tmp,256);
  ifin.getline(tmp,256);
  double raa_5_tsu[11];
  double npart_tsu[11];
  for(int i=0; i<11; i++)
    {
      ifin >> npart_tsu[i] >> raa_5_tsu[i];
    }
  ifin.close();
  TGraph *gTsuRaaVsNpart = new TGraph(11,npart_tsu,raa_5_tsu);
  gTsuRaaVsNpart->SetMarkerStyle(21);
  gTsuRaaVsNpart->SetMarkerColor(1);
  gTsuRaaVsNpart->SetLineColor(1);
  gTsuRaaVsNpart->SetLineWidth(2);
  gTsuRaaVsNpart->Draw("plsames");
  
  // tamu group
  ifin.open("model/data_from_Xingbo/raa_centra_rhic_pt4.5to10_tozebo_110513.dat");
  double raa_5_tamu[15];
  double npart_tamu[15];
  for(int i=0; i<15; i++)
    {
      ifin >> npart_tamu[i] >> raa_5_tamu[i];
      cout << npart_tamu[i] << "  " <<   raa_5_tamu[i] << endl;
    }
  ifin.close();
  TGraph *gTamuRaaVsNpart = new TGraph(15,npart_tamu,raa_5_tamu);
  gTamuRaaVsNpart->SetMarkerStyle(20);
  gTamuRaaVsNpart->SetMarkerColor(1);
  gTamuRaaVsNpart->SetLineColor(1);
  gTamuRaaVsNpart->SetLineWidth(2);
  gTamuRaaVsNpart->SetLineStyle(2);
  gTamuRaaVsNpart->Draw("plsames");
  
  double pub_raa_vs_npart[5] = {0.64, 0.68, 0.72, 0.95, 1.08};
  double pub_raa_vs_npart_err[5] = {0.14, 0.14, 0.14, 0.16, 0.14};
  double pub_raa_vs_npart_sys[5] = {0.03, 0.03, 0.06, 0.06, 0.08};
  TGraphErrors *pubRaaVsNpart = new TGraphErrors(5, pub_npart, pub_raa_vs_npart, pub_npart_xerr, pub_raa_vs_npart_err);
  pubRaaVsNpart->SetMarkerStyle(20);
  pubRaaVsNpart->SetMarkerColor(4);
  pubRaaVsNpart->SetMarkerSize(2);
  pubRaaVsNpart->SetLineColor(4);
  pubRaaVsNpart->Draw("samesPEZ");
  TGraphErrors *pubRaaVsNpartSys = new TGraphErrors(5, pub_npart, pub_raa_vs_npart, pub_npart_xsys, pub_raa_vs_npart_sys);
  pubRaaVsNpartSys->SetMarkerStyle(20);
  pubRaaVsNpartSys->SetMarkerColor(4);
  pubRaaVsNpartSys->SetMarkerSize(0);
  pubRaaVsNpartSys->SetLineColor(4);
  pubRaaVsNpartSys->SetFillStyle(0);
  pubRaaVsNpartSys->Draw("samesE5");
  TBox *boxLabel = new TBox(360,1-0.14,365,1+0.14);
  boxLabel->SetLineColor(4);
  boxLabel->SetFillColor(4);
  boxLabel->SetLineWidth(2.);
  boxLabel->SetFillStyle(3001);
  boxLabel->Draw("fsame");

  const int start_bin = hJpsiInvYield[0]->FindBin(4.5);
  const int end_bin = hJpsiInvYield[0]->FindBin(9.5);
  double pp_yield = 0;
  double pp_err = 0;
  double raa_5GeV[3];
  double raa_5GeV_stat[3];
  double raa_5GeV_sys[3];
  double npart_raa_5GeV[3];
  double npart_raa_5GeV_err[3];
  double npart_raa_5GeV_sys[3];
  for(int k=1; k<nCentBins; k++)
    {
      double auau_yield = 0;
      double auau_err = 0;
      double auau_sys = 0;
      double auau_weight = 0;
      pp_yield = 0;
      pp_err = 0;
      for(int bin=start_bin; bin<=end_bin; bin++)
	{
	  double pt = hJpsiInvYield[k]->GetBinCenter(bin);
	  double dpt = hJpsiInvYield[k]->GetBinWidth(bin);
	  auau_yield += hJpsiInvYield[k]->GetBinContent(bin) * pt * dpt;
	  auau_err += TMath::Power(hJpsiInvYield[k]->GetBinError(bin) * pt * dpt,2);
	  //cout << hJpsiInvYield[k]->GetBinError(bin)/hJpsiInvYield[k]->GetBinContent(bin) << endl;
	  double weight = 1/TMath::Power(hJpsiInvYield[k]->GetBinError(bin)/hJpsiInvYield[k]->GetBinContent(bin),2);
	  auau_weight += weight;
	  auau_sys += weight * hSys->GetBinContent(bin);
	  pp_yield += hJpsipp->GetBinContent(bin) * pt * dpt;
	  pp_err += TMath::Power(hJpsipp->GetBinError(bin) * pt *dpt, 2);
	  //cout << "bin = " << bin << " with raa = " << ppInelastic/ncoll[k] * 1e6 * hJpsiInvYield[k]->GetBinContent(bin)/hJpsipp->GetBinContent(bin) << endl; 
	}
      auau_err = sqrt(auau_err);
      pp_err = sqrt(pp_err);
      auau_sys /= auau_weight;
      auau_yield /= 7.;
      auau_yield /= 6.;
      pp_yield /= 7.;
      pp_yield /= 6.;
      auau_err /= 7.;
      auau_err /= 6.;
      pp_err /= 7.;
      pp_err /= 6.;

      double raa = ppInelastic/ncoll[k] * 1e6 * auau_yield/pp_yield;
      double raa_err = auau_err / auau_yield * raa;
      double raa_sys = auau_sys * raa;
      printf("Cent %s: raa = %1.2f +/- %1.2f +/- %1.2f\n",cent_Name[k],raa,raa_err,raa_sys);

      raa_5GeV[k-1] = raa;
      raa_5GeV_stat[k-1] = raa_err;
      raa_5GeV_sys[k-1] = raa_sys;
      npart_raa_5GeV[k-1] = npart[k];
      npart_raa_5GeV_err[k-1] = 0;
      npart_raa_5GeV_sys[k-1] = 6;
    }

  TGraphErrors *raaVsNpart = new TGraphErrors(3, npart_raa_5GeV, raa_5GeV, npart_raa_5GeV_err, raa_5GeV_stat);
  TGraphErrors *raaVsNpartSys = new TGraphErrors(3, npart_raa_5GeV, raa_5GeV, npart_raa_5GeV_sys, raa_5GeV_sys);

  raaVsNpart->SetMarkerStyle(29);
  raaVsNpart->SetMarkerColor(2);
  raaVsNpart->SetMarkerSize(3);
  raaVsNpart->SetLineColor(2);
  raaVsNpart->Draw("samesPEZ");
  raaVsNpartSys->SetMarkerStyle(20);
  raaVsNpartSys->SetMarkerColor(2);
  raaVsNpartSys->SetMarkerSize(0);
  raaVsNpartSys->SetLineColor(2);
  raaVsNpartSys->SetFillStyle(0);
  raaVsNpartSys->Draw("samesE5");
  double gSys = TMath::Sqrt(ppInelasticErr/ppInelastic*ppInelasticErr/ppInelastic + pp_err*pp_err/pp_yield/pp_yield);
  TBox *boxLabel = new TBox(350,1-gSys,355,1+gSys);
  boxLabel->SetLineColor(kRed);
  boxLabel->SetFillColor(kRed);
  boxLabel->SetLineWidth(2.);
  boxLabel->SetFillStyle(3001);
  boxLabel->Draw("fsame");

  TPaveText *t1 = GetPaveText(0.15,0.3,0.9,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV"));
  t1->Draw();
  TLegend *leg = new TLegend(0.16,0.7,0.4,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);;
  leg->AddEntry(raaVsNpart,"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5, p_{T} > 4 GeV/c","P");
  leg->AddEntry(pubRaaVsNpart,"STAR J/#psi#rightarrowe^{+}e^{-}, |y| < 1, p_{T} > 5 GeV/c","P");
  leg->AddEntry(gTsuRaaVsNpart,"#it{Liu et al.}, p_{T} > 5 GeV/c","PL");
  leg->AddEntry(gTamuRaaVsNpart,"#it{Zhao et al.}, p_{T} > 5 GeV/c","PL");
  leg->Draw();
  TLegend *leg = new TLegend(0.16,0.2,0.4,0.25);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);;
  leg->AddEntry(gPubNpartErr,"N_{coll} uncertainty","f");
  leg->Draw();
  TPaveText *star = GetPaveText(0.7,0.8,0.2,0.25,26,23);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c3 = new TCanvas("Jpsi_raa_vs_npart_LHC","Jpsi_raa_vs_npart_LHC",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  hRaaVsNpart->GetYaxis()->SetRangeUser(0,2);
  hRaaVsNpart->DrawCopy();
  TLine *line = GetLine(hRaaVsNpart->GetXaxis()->GetXmin(),1,hRaaVsNpart->GetXaxis()->GetXmax(),1,1);
  line->Draw();
  TGraphErrors *gPubNpartErr = new TGraphErrors(5, pub_npart, pub_npart_y, pub_npart_err, pub_npart_relerr);
  gPubNpartErr->SetFillColor(kGreen-10);
  gPubNpartErr->SetLineColor(kGreen-10);
  gPubNpartErr->Draw("e3sames");
  TGraph *gTsuRaaVsNpart = new TGraph(11,npart_tsu,raa_5_tsu);
  gTsuRaaVsNpart->SetLineColor(2);
  gTsuRaaVsNpart->SetLineWidth(2);
  gTsuRaaVsNpart->SetLineStyle(2);
  gTsuRaaVsNpart->Draw("lsames");
  TGraph *gTamuRaaVsNpart = new TGraph(15,npart_tamu,raa_5_tamu);
  gTamuRaaVsNpart->SetLineColor(2);
  gTamuRaaVsNpart->SetLineWidth(2);
  gTamuRaaVsNpart->SetLineStyle(1);
  gTamuRaaVsNpart->Draw("lsames");

  TGraphErrors *raaVsNpart = new TGraphErrors(3, npart_raa_5GeV, raa_5GeV, npart_raa_5GeV_err, raa_5GeV_stat);
  TGraphErrors *raaVsNpartSys = new TGraphErrors(3, npart_raa_5GeV, raa_5GeV, npart_raa_5GeV_sys, raa_5GeV_sys);
  raaVsNpart->SetMarkerStyle(29);
  raaVsNpart->SetMarkerColor(2);
  raaVsNpart->SetMarkerSize(3);
  raaVsNpart->SetLineColor(2);
  raaVsNpart->Draw("samesPEZ");
  raaVsNpartSys->SetMarkerStyle(20);
  raaVsNpartSys->SetMarkerColor(2);
  raaVsNpartSys->SetMarkerSize(0);
  raaVsNpartSys->SetLineColor(2);
  raaVsNpartSys->SetFillStyle(0);
  raaVsNpartSys->Draw("samesE5");

  TBox *boxLabel = new TBox(350,1-gSys,355,1+gSys);
  boxLabel->SetLineColor(kRed);
  boxLabel->SetFillColor(kRed);
  boxLabel->SetLineWidth(2.);
  boxLabel->SetFillStyle(3001);
  boxLabel->Draw("fsame");

  TLegend *leg = new TLegend(0.16,0.8,0.4,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);;
  leg->AddEntry(raaVsNpart,"STAR, 200 GeV Au+Au, |y| < 0.5, p_{T} > 4 GeV/c","P");
  leg->AddEntry(gTsuRaaVsNpart,"#it{Liu et al.}, RHIC, p_{T} > 5 GeV/c","PL");
  leg->AddEntry(gTamuRaaVsNpart,"#it{Zhao et al.}, RHIC, p_{T} > 5 GeV/c","PL");
  leg->Draw();
  TLegend *leg = new TLegend(0.16,0.15,0.4,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("STAR");
  leg->AddEntry(gPubNpartErr,"N_{coll} uncertainty","f");
  leg->AddEntry(boxLabel,"p+p uncertainty","f");
  leg->Draw();
  TPaveText *star = GetPaveText(0.7,0.8,0.6,0.7,26,23);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();
  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_muon.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_muon.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_muon.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_muon.eps",run_type,run_cfg_name.Data()));
    }


  double cms_npart[6] = {355.4, 261.4, 187.2, 130.0, 86.3, 22.1};
  double cms_npart_err[6] = {0.,0.,0.,0.,0.,0.};
  double cms_npart_sys[6] = {6,6,6,6,6,6};
  double cms_raa[6] = {0.24, 0.26, 0.31, 0.5, 0.7, 0.62};
  double cms_raa_err[6] = {0.03, 0.03, 0.04, 0.07, 0.11, 0.11};
  double cms_raa_sys[6] = {0.02, 0.02, 0.03, 0.06, 0.09, 0.11};
  TGraphErrors *cmsRaaVsNpart = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_err, cms_raa_err);
  TGraphErrors *cmsRaaVsNpartSys = new TGraphErrors(6, cms_npart, cms_raa, cms_npart_sys, cms_raa_sys);
  cmsRaaVsNpart->SetMarkerStyle(20);
  cmsRaaVsNpart->SetMarkerColor(4);
  cmsRaaVsNpart->SetMarkerSize(2);
  cmsRaaVsNpart->SetLineColor(cmsRaaVsNpart->GetMarkerColor());
  cmsRaaVsNpart->Draw("samesPEZ");
  cmsRaaVsNpartSys->SetMarkerStyle(20);
  cmsRaaVsNpartSys->SetMarkerColor(cmsRaaVsNpart->GetMarkerColor());
  cmsRaaVsNpartSys->SetMarkerSize(0);
  cmsRaaVsNpartSys->SetLineColor(cmsRaaVsNpart->GetMarkerColor());
  cmsRaaVsNpartSys->SetFillStyle(0);
  cmsRaaVsNpartSys->Draw("samesE5");
  ifin.open("model/Jpsi_Rapp_LHC_mid_highPt.dat");
  double npart_tamu_lhc[15];
  double raa_tamu_lhc[15];
  for(int i=0; i<15; i++)
    {
      ifin >> npart_tamu_lhc[i]  >> raa_tamu_lhc[i];
    }
  ifin.close();
  TGraph *gTamuRaaVsNpartLHC = new TGraph(15,npart_tamu_lhc, raa_tamu_lhc);
  gTamuRaaVsNpartLHC->SetLineColor(cmsRaaVsNpart->GetMarkerColor());
  gTamuRaaVsNpartLHC->SetLineWidth(2);
  gTamuRaaVsNpartLHC->SetLineStyle(1);
  gTamuRaaVsNpartLHC->Draw("lsames");


  TLegend *leg = new TLegend(0.16,0.69,0.4,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);;
  leg->AddEntry(cmsRaaVsNpart,"CMS, 2.76 TeV Pb+Pb, |y| < 2.4, 6.5 < p_{T} < 30 GeV/c","P");
  leg->AddEntry(gTamuRaaVsNpartLHC,"#it{Zhao et al.}, LHC","PL");
  leg->Draw();


  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_vs_CMS.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_vs_CMS.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_vs_CMS.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_npart_vs_CMS.eps",run_type,run_cfg_name.Data()));
    }


  // Tsinghua group
  ifin.open("model/data_and_figures_from_Yunpeng/data_Raa_Pt.dat");
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
	   >> v2_init[0][i] >> v2_gen[0][i] >> v2_all[0][i]
	   >> v2_init[1][i] >> v2_gen[1][i] >> v2_all[1][i]
	   >> v2_init[2][i] >> v2_gen[2][i] >> v2_all[2][i]
	   >> v2_init[3][i] >> v2_gen[3][i] >> v2_all[3][i];
    }
  TGraph *grTHU = new TGraph(11,pt_thu,v2_all[0]);
  grTHU->SetName("THU_v2_all_0020");
  ifin.close();

  // TAMU group
  ifin.open("model/data_from_Xingbo/raa_pt_rhic_020_tozebo_110507.dat");
  double pt_tamu[51];
  double v2_tamu[51];
  for(int i=0; i<51; i++)
    {
      ifin >> pt_tamu[i] >> v2_tamu[i];
      //cout << pt_tamu[i] << endl;
    }
  TGraph *grTAMU = new TGraph(51,pt_tamu,v2_tamu);
  grTAMU->SetName("TAMU_v2_all_0020");
  ifin.close();

  TCanvas *c3 = new TCanvas("Raa_Jpsi_vs_pub_0-10","Raa_Jpsi_vs_pub_0-10",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  ScaleHistoTitle(hRaa,28,1,24,28,1,24,63);
  hRaa->SetMaximum(1.7);
  hRaa->DrawCopy();
  TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  gRaaLowPtSys[1]->Draw("sameE5");
  gRaaHighPtSys[1]->Draw("sameE5");
  gRaaLowPt[1]->Draw("sames PEZ");
  gRaaHighPt[1]->Draw("sames PEZ");
  hJpsiRaaSys[1]->Draw("sameE5");
  hJpsiRaa[1]->Draw("samesPEZ");

  double gerr = sqrt(ppInelasticErr*ppInelasticErr/ppInelastic/ppInelastic+ncollErr[1]*ncollErr[1]/ncoll[1]/ncoll[1]);
  TBox *box = new TBox(0.1,1-gerr,0.2,1+gerr);
  box->SetFillStyle(3001);
  box->SetFillColor(kBlack);
  box->Draw();

  TPaveText *t1 = GetPaveText(0.15,0.3,0.9,0.95,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV 0-20%%"));
  t1->Draw();
  TLegend *leg = new TLegend(0.16,0.825,0.4,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);;
  leg->AddEntry(gRaaHighPt[1],"J/#psi#rightarrowe^{+}e^{-}","P");
  leg->Draw();
  TLegend *leg = new TLegend(0.2,0.75,0.46,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);
  leg->AddEntry(gRaaLowPt[1],"STAR J/#psi#rightarrowe^{+}e^{-}, |y| < 1","P");
  leg->AddEntry(hJpsiRaa[1],"STAR J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5","P");
  leg->Draw();

  TPaveText *star = GetPaveText(0.6,0.8,0.2,0.25,24,23);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_0020.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_0020.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_0020.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_0020.eps",run_type,run_cfg_name.Data()));
    }

  grTHU->SetLineColor(kGreen+1);
  grTHU->SetLineWidth(2);
  grTHU->SetLineStyle(2);
  grTHU->Draw("samesCL");

  grTAMU->SetLineColor(kGreen+1);
  grTAMU->SetLineWidth(2);
  grTAMU->SetLineStyle(1);
  grTAMU->Draw("samesCL");

  TLegend *leg = new TLegend(0.65,0.76,0.8,0.97);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(20);
  leg->SetHeader("Model calculation");
  leg->AddEntry(grTAMU,"#it{Zhao et al.}","L");
  leg->AddEntry(grTHU,"#it{Liu et al.}","L");
  leg->Draw();

  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_model_0020.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_model_0020.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_model_0020.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_JpsiRaa_vs_model_0020.eps",run_type,run_cfg_name.Data()));
    }


  if(0)
    {
      // comparison
      TCanvas *c = new TCanvas("AuAu200_Jpsi_ratio","AuAu200_Jpsi_ratio",1100,700);
      c->Divide(2,2);
      TH1F *hRatio = new TH1F("hRatio",";p_{T} (GeV/c);Ratio/TBW",10,0,10);
      hRatio->GetYaxis()->SetRangeUser(0,2);
      ScaleHistoTitle(hRatio,0.06,1,0.05,0.06,1,0.05,62);

      double x,y;
      for(int k=0; k<nCentBins; k++)
	{
	  TGraphErrors *gJpsiYieldRatio = new TGraphErrors(hJpsiInvYield[k]);
	  TGraphErrors *gJpsiYieldSysRatio = new TGraphErrors(hJpsiInvYield[k]);
	  for(int i=0; i<gJpsiYieldRatio->GetN(); i++)
	    {
	      gJpsiYieldRatio->GetPoint(i,x,y);
	      double scale = hTBW[k]->GetBinContent(hTBW[k]->FindBin(x)) / scale_factor[k];

	      gJpsiYieldRatio->SetPoint(i,x,y/scale);
	      gJpsiYieldRatio->SetPointError(i,0,gJpsiYieldRatio->GetErrorY(i)/scale);
	      
	      gJpsiYieldSysRatio->SetPoint(i,x,y/scale);
	      gJpsiYieldSysRatio->SetPointError(i,0.2,hSys->GetBinContent(i+1)*y/scale);
	    }

	  TGraphErrors *gAuAuLowPtRatio = new TGraphErrors(gAuAuLowPt[k]->GetN());
	  gAuAuLowPtRatio->SetName(Form("gAuAuLowPtRatio_Cent%d",k));
	  TGraphErrors *gAuAuLowPtSysRatio = new TGraphErrors(gAuAuLowPtSys[k]->GetN());
	  gAuAuLowPtSysRatio->SetName(Form("gAuAuLowPtSysRatio_Cent%d",k));

	  TGraphErrors *gAuAuHighPtRatio = new TGraphErrors(gAuAuHighPt[k]->GetN());
	  gAuAuHighPtRatio->SetName(Form("gAuAuHighPtRatio_Cent%d",k));
	  TGraphErrors *gAuAuHighPtSysRatio = new TGraphErrors(gAuAuHighPtSys[k]->GetN());
	  gAuAuHighPtSysRatio->SetName(Form("gAuAuHighPtSysRatio_Cent%d",k));

	  for(int i=0; i<gAuAuLowPt[k]->GetN(); i++)
	    {
	      gAuAuLowPt[k]->GetPoint(i,x,y);
	      double scale = hTBW[k]->GetBinContent(hTBW[k]->FindBin(x));

	      gAuAuLowPtRatio->SetPoint(i,x,y/scale);
	      gAuAuLowPtRatio->SetPointError(i,0,gAuAuLowPt[k]->GetErrorY(i)/scale);
	      
	      gAuAuLowPtSysRatio->SetPoint(i,x,y/scale);
	      gAuAuLowPtSysRatio->SetPointError(i,0.2,gAuAuLowPtSys[k]->GetErrorY(i)/scale);
	    }

	  for(int i=0; i<gAuAuHighPt[k]->GetN(); i++)
	    {
	      gAuAuHighPt[k]->GetPoint(i,x,y);
	      double scale = hTBW[k]->GetBinContent(hTBW[k]->FindBin(x));
	      gAuAuHighPtRatio->SetPoint(i,x,y/scale);
	      gAuAuHighPtRatio->SetPointError(i,0,gAuAuHighPt[k]->GetErrorY(i)/scale);
	      
	      gAuAuHighPtSysRatio->SetPoint(i,x,y/scale);
	      gAuAuHighPtSysRatio->SetPointError(i,0.2,gAuAuHighPtSys[k]->GetErrorY(i)/scale);
	    }

	  c->cd(k+1);
	  gPad->SetGridy();
	  SetPadMargin(gPad,0.15,0.15,0.05,0.02);
	  hRatio->Draw();

	  gAuAuLowPtSysRatio->SetFillStyle(0);
	  gAuAuLowPtSysRatio->SetFillColor(4);
	  gAuAuLowPtSysRatio->SetLineColor(4);
	  gAuAuLowPtSysRatio->Draw("samesE5");
	  gAuAuLowPtRatio->SetMarkerStyle(25);
	  gAuAuLowPtRatio->SetMarkerColor(4);
	  gAuAuLowPtRatio->SetLineColor(4);
	  gAuAuLowPtRatio->SetMarkerSize(1.5);
	  gAuAuLowPtRatio->Draw("samesPEZ");

	  gAuAuHighPtSysRatio->SetFillStyle(0);
	  gAuAuHighPtSysRatio->Draw("samesE5");
	  gAuAuHighPtRatio->SetMarkerStyle(25);
	  gAuAuHighPtRatio->SetMarkerSize(1.5);
	  gAuAuHighPtRatio->Draw("samesPEZ");

	  gJpsiYieldSysRatio->SetFillStyle(0);
	  gJpsiYieldSysRatio->SetFillColor(2);
	  gJpsiYieldSysRatio->Draw("samesE5");
	  gJpsiYieldRatio->Draw("samesPEZ");


	  TPaveText *t1 = GetPaveText(0.2,0.3,0.85,0.9,0.06,62);
	  t1->AddText(Form("%s%%",cent_Name[k]));
	  t1->Draw();
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/Run14_JpsiYield_vs_pub.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/Run14_JpsiYield_vs_pub.png",run_type));
	}
    }
}

//================================================
void rawSignal(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(kFALSE);
  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut));
  TH1F *hJpsiSignal[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiSignal[i] = (TH1F*)fin->Get(Form("Jpsi_Signal_cent0020_pt%s_save",pt_Name[i]));
    }
  TH1F *hJpsiCount = (TH1F*)fin->Get("Jpsi_BinCountYield_cent0020");
  TH1F *hSigToBkg = (TH1F*)fin->Get("Jpsi_SigToBkg_cent0020");

  /*
  // raw signal
  TGaxis::SetMaxDigits(3); 
  TH1F *hSeUL = (TH1F*)fin->Get("hJpsiInfoRaw_di_mu_InvMass_jpsi_PtBin0_CentBin0");
  TH1F *hSeLS = (TH1F*)fin->Get("hBkgLSPosRaw_di_mu_InvMass_jpsi_PtBin0_CentBin0");
  TH1F *hMixBkg = (TH1F*)fin->Get("mix_bkg_pt0-10_cent0-60");
  TCanvas *c = new TCanvas("Jpsi_All","Jpsi_All",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hSeUL,0.05,1,0.04,0.05,1.1,0.04,62);
  hSeUL->SetMaximum(1.3*hSeUL->GetMaximum());
  hSeUL->GetYaxis()->SetNdivisions(505);
  hSeUL->SetYTitle("Counts");
  hSeUL->Draw();
  hSeLS->SetMarkerStyle(25);
  hSeLS->Draw("sames");
  hMixBkg->Draw("samesHIST");
  hJpsiSignal[0]->Draw("sames");
  TPaveText *t1 = GetPaveText(0.12,0.5,0.72,0.92,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5, p_{T} > 1 GeV/c");
  t1->AddText("p_{T}^{#mu} > 1.5 GeV/c");
  t1->Draw();

  TPaveText *star = GetPaveText(0.3,0.4,0.2,0.25,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  TLegend *leg = new TLegend(0.65,0.65,0.85,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(24);
  leg->AddEntry(hSeUL,"Unlike-sign pairs","P");
  leg->AddEntry(hSeLS,"Like-sign pairs","P");
  leg->AddEntry(hMixBkg,"Mixed-event","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt0-10_cent0060.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt0-10_cent0060.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt0-10_cent0060.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_pt0-10_cent0060.eps",run_type,run_cfg_name.Data()));
    }

  // all signal
  TCanvas *c1 = new TCanvas("Fit_Jpsi_All","Fit_Jpsi_All",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hJpsiSignal[0],0.05,1,0.04,0.05,1.1,0.04,62);
  hJpsiSignal[0]->SetMaximum(1.5*hJpsiSignal[0]->GetMaximum());
  hJpsiSignal[0]->SetYTitle("Counts");
  hJpsiSignal[0]->SetMarkerSize(1.5);
  hJpsiSignal[0]->Draw();
  TF1 *funcSignal = new TF1("Fit_JpsiSignal_AllPt","gausn(0)+pol3(3)",2.3,4);
  funcSignal->SetParameter(1,3.09);
  funcSignal->SetParameter(2,0.06);
  TFitResultPtr ptr = hJpsiSignal[0]->Fit(funcSignal,"IR0");
  funcSignal->SetLineColor(2);
  funcSignal->Draw("same");

  TPaveText *t1 = GetPaveText(0.55,0.75,0.5,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText("Au+Au @ 200 GeV 0-60%");
  t1->AddText(Form("#it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}");
  t1->AddText("|y|<0.5, p_{T} > 1 GeV/c");
  t1->AddText(Form("S/B = 1:%2.1f",1./hSigToBkg->GetBinContent(0)));
  t1->AddText(Form("N_{J/#psi} = %2.0f",hJpsiCount->GetBinContent(0)));
  t1->AddText(Form("Significance = %2.1f#sigma",hJpsiCount->GetBinContent(0)/hJpsiCount->GetBinError(0)));
  t1->Draw();

  TPaveText *star = GetPaveText(0.22,0.32,0.85,0.9,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_pt0-10_cent0060.eps",run_type,run_cfg_name.Data()));
    }
  */
  TCanvas *c2 = new TCanvas("Fit_jpsi_ptbins","Fit_jpsi_ptbins",1100,650);
  c2->Divide(3,2);
  for(int i=1; i<nPtBins; i++)
    {
      c2->cd(i);
      SetPadMargin(gPad,0.13,0.13,0.05,0.05);
      ScaleHistoTitle(hJpsiSignal[i],0.06,1,0.05,0.06,1.1,0.05,62);
      hJpsiSignal[i]->SetMaximum(1.8*hJpsiSignal[i]->GetMaximum());
      hJpsiSignal[i]->SetYTitle("Counts");
      hJpsiSignal[i]->SetMarkerSize(1.5);
      hJpsiSignal[i]->Draw();
      TF1 *funcSignal = new TF1(Form("Fit_JpsiSignal_pt%s",pt_Name[i]),"gausn(0)+pol1(3)",2.3,4);
      funcSignal->SetParameter(1,3.09);
      funcSignal->SetParameter(2,0.06);
      hJpsiSignal[i]->Fit(funcSignal,"IR0");
      funcSignal->SetLineColor(2);
      funcSignal->Draw("same");
      if(i==1)
	{
	  TPaveText *t1 = GetPaveText(0.17,0.37,0.65,0.92,0.06,62);
	  t1->SetTextAlign(11);
	  t1->AddText("Au+Au @ 200 GeV");
	  t1->AddText("0-20%");
	  t1->AddText(Form("#it{L} ~ %1.1f nb^{-1}",luminosity));
	  t1->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}");
	  t1->Draw();
	}

      t1 = GetPaveText(0.53,0.75,0.5,0.83,0.06,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]));
      t1->AddText(Form("S/B = 1:%2.1f",1./hSigToBkg->GetBinContent(i)));
      t1->AddText(Form("N_{J/#psi} = %2.0f",hJpsiCount->GetBinContent(i)));
      t1->AddText(Form("Significance: %2.1f",hJpsiCount->GetBinContent(i)/hJpsiCount->GetBinError(i)));
      t1->Draw();
    }
  /*
  c2->cd(2);
  TPaveText *star = GetPaveText(0.25,0.5,0.85,0.93,0.06);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  */
  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0020.pdf",run_type,run_cfg_name.Data()));
      /*
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0060.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0060.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/FitJpsi_ptbins_cent0060.eps",run_type,run_cfg_name.Data()));
      */
    }

}

//================================================
void upsilon(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("Rootfiles/UpsFinal.root","read");
  TH1F *hUL = (TH1F*)f->Get("SE_UL");
  TH1F *hLS = (TH1F*)f->Get("SE_LS");
  TH1F *hSub = (TH1F*)f->Get("SE_Ups");
  TF1 *func[4];
  func[0] = (TF1*)f->Get("fun");
  func[1] = (TF1*)f->Get("r1s");
  func[2] = (TF1*)f->Get("r2s");
  func[3] = (TF1*)f->Get("r3s");

  TCanvas *c1 = new TCanvas("Fit_Jpsi_All","Fit_Jpsi_All",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hUL,0.05,1,0.04,0.05,1.1,0.04,62);
  hUL->SetMaximum(120);
  hUL->SetTitle("");
  hUL->SetLineColor(2);
  hUL->Draw();
  hLS->Draw("sames HIST");

  int lowbin = hUL->FindBin(9);
  int highbin = hUL->FindBin(12);
  printf("[i] %d counts in unlike-sign\n",hUL->Integral(lowbin,highbin));
  printf("[i] %d counts in like-sign\n",hLS->Integral(lowbin,highbin));
  TPaveText *t1 = GetPaveText(0.15,0.25,0.85,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->Draw();
  TLegend *leg = new TLegend(0.5,0.55,0.75,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("#Upsilon#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hUL,"Unlike-sign pairs","P");
  leg->AddEntry(hLS,"Like-sign pairs","L");
  leg->Draw();
  TPaveText *star = GetPaveText(0.15,0.4,0.15,0.2,0.04,20);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c2 = new TCanvas("Fit_Upsilon","Fit_Upsilon",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  //gPad->SetGrid(1,1);
  ScaleHistoTitle(hSub,0.05,1,0.04,0.05,1.1,0.04,62);
  hSub->SetTitle("");
  hSub->SetMarkerColor(2);
  hSub->SetLineColor(1);
  hSub->GetYaxis()->SetRangeUser(-20,40);
  hSub->Draw();
  for(int i=0; i<4; i++)
    {
      func[i]->SetLineStyle(i+1);
      func[i]->Draw("sames");
    }
  TF1 *funcBkg = new TF1("funcBkg","[0]*expo(1)",8,12);
  for(int i=0; i<3; i++)
    {
      funcBkg->SetParameter(i,func[0]->GetParameter(9+i));
    }
  funcBkg->SetLineColor(kGreen+1);
  funcBkg->SetLineStyle(5);
  funcBkg->Draw("sames");
  printf("Fit: chi2/ndf = %1.1f/%1.0f\n",func[0]->GetChisquare(), func[0]->GetNDF());

  TPaveText *t1 = GetPaveText(0.15,0.25,0.85,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->Draw();
  TLegend *leg = new TLegend(0.65,0.63,0.85,0.92);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("#Upsilon#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hSub,"UL-LS","P");
  leg->AddEntry(func[0],"Combined fit","L");
  leg->AddEntry(func[1],"#Upsilon(1S)","L");
  leg->AddEntry(func[2],"#Upsilon(2S)","L");
  leg->AddEntry(func[3],"#Upsilon(3S)","L");
  leg->AddEntry(funcBkg,"Drell-Yan+b#bar{b}","L");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.65,0.85,0.55,0.60,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("#chi^{2}/NDF = %1.1f/%d\n",func[0]->GetChisquare(), func[0]->GetNDF()));
  t1->Draw();
  TPaveText *star = GetPaveText(0.7,0.85,0.15,0.2,0.04,20);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.eps",run_type,run_cfg_name.Data()));
    }
 
}

//================================================
void publication(const bool savePlot = 1, const bool saveHisto = 1)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TFile *fout = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","update");

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TFile *f = TFile::Open("Rootfiles/Publication.Jpsi.pp200GeV.root","read");
  gStyle->SetOptStat(0);

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

  TH1F *hLowPtPP = (TH1F*)f->Get("hInvariantYield_pp");
  hLowPtPP->SetMarkerStyle(24);
  hLowPtPP->SetMarkerSize(1.5);
  hLowPtPP->Draw("sames");

  TH1F *hJpsiPP = new TH1F("Jpsi_xsec_pp200_combined",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",nbins,xbins);
  for(int ibin=1; ibin<=3; ibin++)
    {
      hJpsiPP->SetBinContent(ibin,hLowPtPP->GetBinContent(ibin+1));
      hJpsiPP->SetBinError(ibin,hLowPtPP->GetBinError(ibin+1));
    }
  double val = (ypp[4]*xpp[4]+ypp[5]*xpp[5])*1./(2.*5.);
  double err = sqrt(eylpp[4]*eylpp[4]*xpp[4]*xpp[4]+eylpp[5]*eylpp[5]*xpp[5]*xpp[5])*1./(2.*5.);
  double syl = (sylpp[4]*xpp[4]+sylpp[5]*xpp[5])*1./(2.*5.);
  double syh = (syhpp[4]*xpp[4]+syhpp[5]*xpp[5])*1./(2.*5.);
  double sys = (syl+syh)/2.;
  hJpsiPP->SetBinContent(4,val);
  hJpsiPP->SetBinError(4,sqrt(err*err+sys*sys));
  val = (ypp[6]*xpp[6]*1.+ypp[7]*xpp[7]*1.)/(2.*7.);
  err = sqrt(eylpp[6]*eylpp[6]*xpp[6]*xpp[6]+eylpp[7]*eylpp[7]*xpp[7]*xpp[7])/(2.*7.);
  syl = (sylpp[6]*xpp[6]*1.+sylpp[7]*xpp[7]*1.)/(2.*7.);
  syh = (syhpp[6]*xpp[6]*1.+syhpp[7]*xpp[7]*1.)/(2.*7.);
  sys = (syl+syh)/2.;
  hJpsiPP->SetBinContent(5,val);
  hJpsiPP->SetBinError(5,sqrt(err*err+sys*sys));
  val = ypp[8];
  err = eylpp[8];
  sys = (sylpp[8]+syhpp[8])/2;
  hJpsiPP->SetBinContent(6,val);
  hJpsiPP->SetBinError(6,sqrt(err*err+sys*sys));


  hJpsiPP->SetMarkerStyle(21);
  hJpsiPP->SetMarkerSize(1.5);
  hJpsiPP->SetMarkerColor(2);
  hJpsiPP->SetLineColor(2);
  hJpsiPP->Draw("sames");
  leg = new TLegend(0.5,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  leg->AddEntry(gHighPtPP,"STAR high p_{T}","PL");
  leg->AddEntry(hJpsiPP,"PHENIX + STAR","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_xsec_pp.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_xsec_pp.png",run_type));
    }

  // AuAu
  const int nAuAuLowPt = 5;
  double xAuAuLowPt[nCentBins][nAuAuLowPt];
  double exlAuAuLowPt[nCentBins][nAuAuLowPt];
  double sxlAuAuLowPt[nCentBins][nAuAuLowPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nAuAuLowPt; j++)
	{
	  xAuAuLowPt[i][j] = 0.5 + j; 
	  exlAuAuLowPt[i][j] = 0; 
	  sxlAuAuLowPt[i][j] = 0.2; 
	}
    }
  double yAuAuLowPt[nCentBins][nAuAuLowPt];
  double eylAuAuLowPt[nCentBins][nAuAuLowPt];
  double sylAuAuLowPt[nCentBins][nAuAuLowPt];
  double syhAuAuLowPt[nCentBins][nAuAuLowPt];

  double xRaaLowPt[nCentBins][nAuAuLowPt];
  double exlRaaLowPt[nCentBins][nAuAuLowPt];
  double sxlRaaLowPt[nCentBins][nAuAuLowPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nAuAuLowPt; j++)
	{
	  xRaaLowPt[i][j] = 0.5 + j; 
	  exlRaaLowPt[i][j] = 0; 
	  sxlRaaLowPt[i][j] = 0.2; 
	}
    }
  double yRaaLowPt[nCentBins][nAuAuLowPt];
  double eylRaaLowPt[nCentBins][nAuAuLowPt];
  double sylRaaLowPt[nCentBins][nAuAuLowPt];
  double syhRaaLowPt[nCentBins][nAuAuLowPt];

  ifstream data_in;
  data_in.open("/Users/admin/Work/STAR/analysis/Rootfiles/Publication.Jpsi.AuAu200.LowPt.txt");
  char data[10][256];
  double meanx, meany, staty, sysyhigh, sysylow, raa, statraa, sysraahigh, sysraalow, sysglobal;
  int counter = 0;
  while(!data_in.eof())
    {
      for(int i=0; i<10; i++)
	data_in >> data[i];

      int cent = counter/nAuAuLowPt;
      int ptbin = counter%nAuAuLowPt;
      yAuAuLowPt[cent][ptbin] = atof(data[1]);
      eylAuAuLowPt[cent][ptbin] = atof(data[2]);
      syhAuAuLowPt[cent][ptbin] = atof(data[3]);
      sylAuAuLowPt[cent][ptbin] = fabs(atof(data[4]));

      yRaaLowPt[cent][ptbin] = atof(data[5]);
      eylRaaLowPt[cent][ptbin] = atof(data[6]);
      syhRaaLowPt[cent][ptbin] = atof(data[7]);
      sylRaaLowPt[cent][ptbin] = fabs(atof(data[8]));
      //cout << sylAuAuLowPt[cent][ptbin] << endl;
      //cout << counter << "  " << cent << "  " << ptbin << endl;
      counter++;
      if(counter==nAuAuLowPt*nCentBins) break; 
    }
  data_in.close();

  const int nAuAuHighPt = 6;
  double xAuAuHighPt[nCentBins][nAuAuHighPt];
  double exlAuAuHighPt[nCentBins][nAuAuHighPt];
  double sxlAuAuHighPt[nCentBins][nAuAuHighPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nAuAuHighPt; j++)
	{
	  exlAuAuHighPt[i][j] = 0; 
	  sxlAuAuHighPt[i][j] = 0.2; 
	}
    }
  double yAuAuHighPt[nCentBins][nAuAuHighPt];
  double eylAuAuHighPt[nCentBins][nAuAuHighPt];
  double sylAuAuHighPt[nCentBins][nAuAuHighPt];
  double syhAuAuHighPt[nCentBins][nAuAuHighPt];

  data_in.open("/Users/admin/Work/STAR/analysis/Rootfiles/Publication.Jpsi.AuAu200.HighPt.txt");
  int counter = 0;
  while(!data_in.eof())
    {
      for(int i=0; i<6; i++)
	data_in >> data[i];

      int cent = counter/nAuAuHighPt;
      int ptbin = counter%nAuAuHighPt;
      xAuAuHighPt[cent][ptbin] = atof(data[1]);
      yAuAuHighPt[cent][ptbin] = atof(data[2]);
      eylAuAuHighPt[cent][ptbin] = atof(data[3]);
      sylAuAuHighPt[cent][ptbin] = atof(data[4]);
      syhAuAuHighPt[cent][ptbin] = fabs(atof(data[5]));

      counter++;
      if(counter==nAuAuHighPt*nCentBins) break;
    }
  data_in.close();


  const int nRaaHighPt = 5;
  double xRaaHighPt[nCentBins][nRaaHighPt];
  double exlRaaHighPt[nCentBins][nRaaHighPt];
  double sxlRaaHighPt[nCentBins][nRaaHighPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nRaaHighPt; j++)
	{
	  exlRaaHighPt[i][j] = 0; 
	  sxlRaaHighPt[i][j] = 0.2; 
	}
    }
  double yRaaHighPt[nCentBins][nRaaHighPt];
  double eylRaaHighPt[nCentBins][nRaaHighPt];
  double sylRaaHighPt[nCentBins][nRaaHighPt];
  double syhRaaHighPt[nCentBins][nRaaHighPt];
  data_in.open("/Users/admin/Work/STAR/analysis/Rootfiles/Publication.Jpsi.AuAu200.Raa.HighPt.txt");
  int counter = 0;
  while(!data_in.eof())
    {
      for(int i=0; i<6; i++)
	data_in >> data[i];

      int cent = counter/nRaaHighPt;
      int ptbin = counter%nRaaHighPt;
      xRaaHighPt[cent][ptbin] = atof(data[1]);
      yRaaHighPt[cent][ptbin] = atof(data[2]);
      eylRaaHighPt[cent][ptbin] = atof(data[3]);
      sylRaaHighPt[cent][ptbin] = atof(data[4]);
      syhRaaHighPt[cent][ptbin] = fabs(atof(data[5]));
      cout << xRaaHighPt[cent][ptbin] << endl;
      counter++;
      if(counter==nRaaHighPt*nCentBins) break;
    }
  data_in.close();
  
  // Jpsi Raa
  TCanvas *c = new TCanvas("AuAu200_Raa","AuAu200_Raa",1100,700);
  c->Divide(2,2);
  TH1F *hRaa = new TH1F("AuAu200_Raa",";p_{T} (GeV/c);R_{AA}",10,0,10);
  hRaa->GetYaxis()->SetRangeUser(0,1.8);
  ScaleHistoTitle(hRaa,0.06,1,0.05,0.06,1,0.05,62);

  TGraphAsymmErrors *gRaaLowPt[nCentBins];
  TGraphAsymmErrors *gRaaLowPtSys[nCentBins];
  TGraphAsymmErrors *gRaaHighPt[nCentBins];
  TGraphAsymmErrors *gRaaHighPtSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
      hRaa->Draw();
      
      gRaaLowPt[i] = new TGraphAsymmErrors(nAuAuLowPt, xRaaLowPt[i], yRaaLowPt[i], exlRaaLowPt[i], exlRaaLowPt[i], eylRaaLowPt[i], eylRaaLowPt[i]);
      gRaaLowPt[i]->SetName(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",cent_Title[i]));
      gRaaLowPt[i]->SetMarkerStyle(20);
      gRaaLowPt[i]->SetMarkerSize(1.5);
      gRaaLowPt[i]->SetMarkerColor(2);
      gRaaLowPt[i]->SetLineColor(2);
      gRaaLowPt[i]->Draw("sames PE");

      gRaaLowPtSys[i] = new TGraphAsymmErrors(nAuAuLowPt, xRaaLowPt[i], yRaaLowPt[i], sxlRaaLowPt[i], sxlRaaLowPt[i], sylRaaLowPt[i], syhRaaLowPt[i]);
      gRaaLowPtSys[i]->SetName(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",cent_Title[i]));
      gRaaLowPtSys[i]->SetMarkerStyle(20);
      gRaaLowPtSys[i]->SetMarkerSize(0);
      gRaaLowPtSys[i]->SetFillStyle(0);
      gRaaLowPtSys[i]->SetMarkerColor(2);
      gRaaLowPtSys[i]->SetLineColor(2);
      gRaaLowPtSys[i]->Draw("sameE5");

      gRaaHighPt[i] = new TGraphAsymmErrors(nRaaHighPt, xRaaHighPt[i], yRaaHighPt[i], exlRaaHighPt[i], exlRaaHighPt[i], eylRaaHighPt[i], eylRaaHighPt[i]);
      gRaaHighPt[i]->SetName(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",cent_Title[i]));
      gRaaHighPt[i]->SetMarkerStyle(21);
      gRaaHighPt[i]->SetMarkerSize(1.5);
      gRaaHighPt[i]->SetMarkerColor(4);
      gRaaHighPt[i]->SetLineColor(4);
      gRaaHighPt[i]->Draw("sames PE");

      gRaaHighPtSys[i] = new TGraphAsymmErrors(nRaaHighPt, xRaaHighPt[i], yRaaHighPt[i], sxlRaaHighPt[i], sxlRaaHighPt[i], sylRaaHighPt[i], syhRaaHighPt[i]);
      gRaaHighPtSys[i]->SetName(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",cent_Title[i]));
      gRaaHighPtSys[i]->SetMarkerStyle(21);
      gRaaHighPtSys[i]->SetMarkerSize(0);
      gRaaHighPtSys[i]->SetFillStyle(0);
      gRaaHighPtSys[i]->SetMarkerColor(4);
      gRaaHighPtSys[i]->SetLineColor(4);
      gRaaHighPtSys[i]->Draw("sameE5");

      TPaveText *t1 = GetPaveText(0.56,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[i]));
      t1->Draw();

      TLine *line = GetLine(0,1,10,1,1);
      line->Draw();
    }
  c->cd(1);
  leg = new TLegend(0.2,0.6,0.4,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Au+Au @ 200 GeV");
  leg->AddEntry(gRaaLowPt[0],"STAR: Run10","P");
  leg->AddEntry(gRaaHighPt[0],"STAR: Run10","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Raa.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Raa.png",run_type));
    }

  // combined data point for fit
  TCanvas *c = new TCanvas("AuAu200_Jpsi_fit","AuAu200_Jpsi_fit",1100,700);
  c->Divide(2,2);

  const int nbinsAuAu = 9;
  const double xbinsAuAu[nbinsAuAu+1] = {0,1,2,3,4,5,6,7,8,10};
  TH1F *hJpsiAuAu[nCentBins];
  TF1 *funcJpsiAuAu[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hJpsiAuAu[i] = new TH1F(Form("Jpsi_InvYield_AuAu200_Combined_cent%s",cent_Title[i]),Form("Invariant yield of J/#psi (%s%%);p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",cent_Name[i]), nbinsAuAu, xbinsAuAu);
      for(int bin=1; bin<=nbinsAuAu; bin++)
	{
	  if(bin<=4)
	    {
	      hJpsiAuAu[i]->SetBinContent(bin, yAuAuLowPt[i][bin-1] * xAuAuLowPt[i][bin-1]);
	      hJpsiAuAu[i]->SetBinError(bin, eylAuAuLowPt[i][bin-1] * xAuAuLowPt[i][bin-1]);
	    }
	  else
	    {
	      hJpsiAuAu[i]->SetBinContent(bin,  yAuAuHighPt[i][bin-4] * xAuAuHighPt[i][bin-4]);
	      hJpsiAuAu[i]->SetBinError(bin, eylAuAuHighPt[i][bin-4] * xAuAuHighPt[i][bin-4]);
	    }
	}
      
      hJpsiAuAu[i]->SetMarkerStyle(21);
      hJpsiAuAu[i]->GetYaxis()->SetRangeUser(1e-10,1e-4);
      ScaleHistoTitle(hJpsiAuAu[i],0.06,1,0.05,0.06,1,0.05,62);
      c->cd(i+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);

      funcJpsiAuAu[i] = new TF1(Form("Fit_Jpsi_InvYield_AuAu200_Combined_cent%s",cent_Title[i]),"[0]*exp([1]*x+[2]*x*x+[3]*x*x*x)",0.1,10);
      if(i==3)funcJpsiAuAu[i]->SetParameter(0,1e-5);
      hJpsiAuAu[i]->Fit(funcJpsiAuAu[i],"IR0");
      hJpsiAuAu[i]->Draw();
      funcJpsiAuAu[i]->SetLineColor(2);
      funcJpsiAuAu[i]->Draw("sames");
    }


  // Jpsi yield
  TH1F *hTBW[4];
  TCanvas *c = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",1100,700);
  c->Divide(2,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",10,0,10);
  hAuAu->GetYaxis()->SetRangeUser(1e-10,1e-4);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);

  TGraphAsymmErrors *gAuAuLowPt[nCentBins];
  TGraphAsymmErrors *gAuAuLowPtSys[nCentBins];
  TGraphAsymmErrors *gAuAuHighPt[nCentBins];
  TGraphAsymmErrors *gAuAuHighPtSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      c->cd(i+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
      hAuAu->Draw();
      
      gAuAuLowPt[i] = new TGraphAsymmErrors(nAuAuLowPt, xAuAuLowPt[i], yAuAuLowPt[i], exlAuAuLowPt[i], exlAuAuLowPt[i], eylAuAuLowPt[i], eylAuAuLowPt[i]);
      gAuAuLowPt[i]->SetName(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[i]));
      gAuAuLowPt[i]->SetMarkerStyle(20);
      gAuAuLowPt[i]->SetMarkerSize(1.5);
      gAuAuLowPt[i]->SetMarkerColor(2);
      gAuAuLowPt[i]->SetLineColor(2);
      gAuAuLowPt[i]->Draw("sames PE");

      gAuAuLowPtSys[i] = new TGraphAsymmErrors(nAuAuLowPt, xAuAuLowPt[i], yAuAuLowPt[i], sxlAuAuLowPt[i], sxlAuAuLowPt[i], sylAuAuLowPt[i], syhAuAuLowPt[i]);
      gAuAuLowPtSys[i]->SetName(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[i]));
      gAuAuLowPtSys[i]->SetMarkerStyle(20);
      gAuAuLowPtSys[i]->SetMarkerSize(0);
      gAuAuLowPtSys[i]->SetFillStyle(0);
      gAuAuLowPtSys[i]->SetMarkerColor(2);
      gAuAuLowPtSys[i]->SetLineColor(2);
      gAuAuLowPtSys[i]->Draw("sameE5");

      gAuAuHighPt[i] = new TGraphAsymmErrors(nAuAuHighPt, xAuAuHighPt[i], yAuAuHighPt[i], exlAuAuHighPt[i], exlAuAuHighPt[i], eylAuAuHighPt[i], eylAuAuHighPt[i]);
      gAuAuHighPt[i]->SetName(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[i]));
      gAuAuHighPt[i]->SetMarkerStyle(21);
      gAuAuHighPt[i]->SetMarkerSize(1.5);
      gAuAuHighPt[i]->SetMarkerColor(4);
      gAuAuHighPt[i]->SetLineColor(4);
      gAuAuHighPt[i]->Draw("sames PE");

      gAuAuHighPtSys[i] = new TGraphAsymmErrors(nAuAuHighPt, xAuAuHighPt[i], yAuAuHighPt[i], sxlAuAuHighPt[i], sxlAuAuHighPt[i], sylAuAuHighPt[i], syhAuAuHighPt[i]);
      gAuAuHighPtSys[i]->SetName(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[i]));
      gAuAuHighPtSys[i]->SetMarkerStyle(21);
      gAuAuHighPtSys[i]->SetMarkerSize(0);
      gAuAuHighPtSys[i]->SetFillStyle(0);
      gAuAuHighPtSys[i]->SetMarkerColor(4);
      gAuAuHighPtSys[i]->SetLineColor(4);
      gAuAuHighPtSys[i]->Draw("sameE5");

      TPaveText *t1 = GetPaveText(0.56,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[i]));
      t1->Draw();

      hTBW[i] = (TH1F*)fout->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[i]));
      hTBW[i]->SetLineColor(1);
      hTBW[i]->SetLineStyle(2);
      hTBW[i]->SetLineWidth(2);
      hTBW[i]->Draw("sames");
    }
  c->cd(1);
  leg = new TLegend(0.2,0.2,0.4,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Au+Au @ 200 GeV");
  leg->AddEntry(gAuAuLowPt[0],"STAR: Run10","P");
  leg->AddEntry(gAuAuHighPt[0],"STAR: Run10","P");
  leg->AddEntry(hTBW[0],"TBW (#beta=0)","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Yield.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Yield.png",run_type));
    }

  TH1F *hTBW2[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hTBW2[i] = (TH1F*)hTBW[i]->Clone(Form("TBW_Jpsi_Yield_cent%s",cent_Title[i]));
      for(int bin=1; bin<=hTBW2[i]->GetNbinsX(); bin++)
	{
	  double pt = hTBW2[i]->GetBinCenter(bin);
	  hTBW2[i]->SetBinContent(bin,hTBW2[i]->GetBinContent(bin)*pt);
	  cout << hTBW2[i]->GetBinError(bin) << endl;
	  //hTBW2[i]->SetBinError(bin,hTBW2[i]->GetBinError(bin)*pt);
	}
    }

  if(saveHisto)
    {
      fout->cd();

      gHighPtPP->Write("",TObject::kOverwrite);
      gHighPtPPSys->Write("",TObject::kOverwrite);
      hLowPtPP->Write("",TObject::kOverwrite);
      hJpsiPP->Write("",TObject::kOverwrite);

      for(int i=0; i<nCentBins; i++)
	{
	  gAuAuLowPt[i]->Write("",TObject::kOverwrite);
	  gAuAuLowPtSys[i]->Write("",TObject::kOverwrite);
	  gAuAuHighPt[i]->Write("",TObject::kOverwrite);
	  gAuAuHighPtSys[i]->Write("",TObject::kOverwrite);
	  hTBW[i]->Write("",TObject::kOverwrite);
	  hTBW2[i]->Write("",TObject::kOverwrite);
	  hJpsiAuAu[i]->Write("",TObject::kOverwrite);
	  funcJpsiAuAu[i]->Write("",TObject::kOverwrite);
	}

      for(int i=0; i<nCentBins; i++)
	{
	  gRaaLowPt[i]->Write("",TObject::kOverwrite);
	  gRaaLowPtSys[i]->Write("",TObject::kOverwrite);
	  gRaaHighPt[i]->Write("",TObject::kOverwrite);
	  gRaaHighPtSys[i]->Write("",TObject::kOverwrite);
	}
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

