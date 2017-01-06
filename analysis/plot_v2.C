//================================================
void plot_v2()
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
  ifs.open("Rootfiles/Published/Jpsi_v2_200/JpsiSpv2020T120Tau6-00100.dat");
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
}


//-----------------------------------------
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size = 0.04, const Int_t font = 42)
{
  TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}
