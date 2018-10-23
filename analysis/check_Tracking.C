//================================================
void check_Tracking()
{
  gStyle->SetOptStat(0);

  //trackQA();
  //tightCutsEmb();
  trackYield();
  //gRefMult();
}

//================================================
void trackYield(const int savePlot = 0)
{
  const char* trigTitle = "P15ic";
  TFile *fin = TFile::Open(Form("output/Run14_AuAu200.MB.%s.Study.TpcTracking.root",trigTitle), "read");
  char *trigName = "mb";

  TFile *femb = NULL;
  if(iMbEmb==0) femb = TFile::Open(Form("output/%s.Embed.Jpsi.root",run_type.Data()),"read");
  if(iMbEmb==1) femb = TFile::Open(Form("output/%s.Embed_MB.Jpsi.root",run_type.Data()),"read");

  // vertex and centrality
  const int gVzCut = 30;
  THnSparseF *hnCent = (THnSparse*)fin->Get(Form("mhnCent_%s",trigName));
  hnCent->GetAxis(1)->SetRange(1,16);
  const int nHisto = 4;
  const char* histoName[nHisto] = {"gRefMultCorr", "cent", "ZDC", "TpcVz"};
  TH1F *histo[nHisto][gNTrgSetup];
  TH1F *hCent[gNTrgSetup];
  TH2F *hCentVsZdc[gNTrgSetup];
  for(int k=2; k<gNTrgSetup; k++)
    {
      hnCent->GetAxis(4)->SetRange(k, k);
      
      hnCent->GetAxis(3)->SetRangeUser(-1*gVzCut, gVzCut);
      hCent[k] = (TH1F*)hnCent->Projection(1);
      hCent[k]->SetName(Form("hCent%s",gTrgSetupTitle[k]));
      hCentVsZdc[k] = (TH2F*)hnCent->Projection(1,2);
      hCentVsZdc[k]->SetName(Form("hCentVsZdc%s",gTrgSetupTitle[k]));
      hnCent->GetAxis(3)->SetRange(0,-1);

      for(int i=0; i<nHisto; i++)
	{
	  histo[i][k] = (TH1F*)hnCent->Projection(i);
	  histo[i][k]->Sumw2();
	  histo[i][k]->SetName(Form("%s%s_0080",histoName[i], gTrgSetupTitle[k]));
	  if(i!=2) histo[i][k]->Scale(1./histo[i][k]->Integral());
	}
    }
  hnCent->GetAxis(4)->SetRange(0, -1);

  TCanvas *c = new TCanvas("cCent", "cCent", 1100, 700);
  c->Divide(2,2);
  TLegend *leg =  new TLegend(0.25, 0.25, 0.5, 0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int i=0; i<nHisto; i++)
    {
      c->cd(i+1);
      for(int k=2; k<gNTrgSetup; k++)
	{
	  histo[i][k]->SetLineColor(gColor[k-1]);
	  histo[i][k]->SetTitle("");
	  if(i==0)
	    {
	      gPad->SetLogy();
	      histo[i][k]->GetXaxis()->SetRangeUser(0, 700);
	      histo[i][k]->GetYaxis()->SetRangeUser(1e-7, 1e-1);
	      histo[i][k]->SetXTitle("gRefMultCorr");
	      leg->AddEntry(histo[i][k], Form("%s",gLegNameTrg[k].Data()), "PL");
	    }
	  else if(i==1)
	    {
	      histo[i][k]->GetYaxis()->SetRangeUser(0, 0.08);
	    }
	  else if(i==2)
	    {
	      histo[i][k]->GetYaxis()->SetRangeUser(0, 7e5);
	    }
	  if(k==0) histo[i][k]->Draw("HIST");
	  else     histo[i][k]->Draw("samesHIST");
	}
      TPaveText *t1 = GetTitleText(Form("%s_%s (%s, 0-80%%)",run_type.Data(),trigName,trigTitle), 0.06);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_CentVsLumi.pdf",run_type.Data(),trigName));

  const int nCentBins = 4; 
  int centBins_low[nCentBins]  = {16,12,8,4};
  int centBins_high[nCentBins] = {16,12,8,4};
  TString cent_Name[nCentBins];
  TString cent_Title[nCentBins];
  for(int j=0; j<nCentBins; j++)
    {
      cent_Name[j] = Form("%d%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
      cent_Title[j] = Form("%d-%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
    }
  
  // compare track quality
  THnSparseF *hnTrk = (THnSparseF*)fin->Get(Form("mhnTrkPt_%s",trigName));
  hnTrk->GetAxis(8)->SetRangeUser(-1*gVzCut, gVzCut);

  THnSparseF *hnTrkEtaPhi = (THnSparseF*)fin->Get(Form("mhnTrkEtaPhi_%s",trigName));
  hnTrkEtaPhi->GetAxis(8)->SetRangeUser(-1*gVzCut, gVzCut);
  //hnTrkEtaPhi->GetAxis(5)->SetRange(2,2);
  //hnTrkEtaPhi->GetAxis(8)->SetRangeUser(-30, 30);
  hnTrkEtaPhi->GetAxis(7)->SetRange(1,1); // DCA < 1 cm

  
  const int nVar = 4;
  const char* varName[nVar] = {"NHitsFit", "DCA", "eta", "phi"};
  const char* varTitle[nVar] = {"NHitsFit", "DCA", "#eta", "#varphi"};
  TH1F *hTrkVar[nVar][nCentBins][gNTrgSetup];
  TH2F *hTrkVarVsZdc[nVar][nCentBins][gNTrgSetup];
  TProfile *hTrkVarPro[nVar][nCentBins][gNTrgSetup];
  THnSparseF* hntmp = 0x0;
  for(int v=0; v<nVar; v++)
    {
      if(v<2) hntmp = hnTrk;
      else    hntmp = hnTrkEtaPhi;
      //hntmp->GetAxis(6-v/2)->SetRange(2,2); // match to MTD
      //if(v<2) hntmp->GetAxis(5)->SetRange(2,2); // muon candidates
      for(int k=2; k<gNTrgSetup; k++)
	{
	  hntmp->GetAxis(4)->SetRange(k, k);

	  TH3F *h3tmp = (TH3F*)hntmp->Projection(3, 7-v/2, v+1-v/2*2);
	  h3tmp->SetName(Form("h3tmp_%d%d",v,k));
	  for(int j=0; j<nCentBins; j++)
	    {
	      hTrkVar[v][j][k] = (TH1F*)h3tmp->ProjectionZ(Form("hTrk%s_%s%s",varName[v],cent_Name[j].Data(),gTrgSetupTitle[k]), centBins_low[j], centBins_high[j], 0, -1);

	      h3tmp->GetXaxis()->SetRange(centBins_low[j], centBins_high[j]);
	      hTrkVarVsZdc[v][j][k] = (TH2F*)h3tmp->Project3D("zy");
	      hTrkVarVsZdc[v][j][k]->SetName(Form("hTrk%sVsZdc_%s%s",varName[v],cent_Name[j].Data(),gTrgSetupTitle[k]));
	      hTrkVarPro[v][j][k] = (TProfile*)hTrkVarVsZdc[v][j][k]->ProfileX(Form("hTrk%sPro_%s%s",varName[v],cent_Name[j].Data(),gTrgSetupTitle[k]));
	    }
	  hntmp->GetAxis(4)->SetRange(0, -1);
	}
    }
  TCanvas *cVar[nVar];
  for(int v=0; v<nVar; v++)
    {
      cVar[v] = new TCanvas(Form("c%s",varName[v]),Form("c%s",varName[v]), 1100, 700);
      cVar[v]->Divide(4,2);
      TLegend *leg = NULL;
      if(v==0) leg = new TLegend(0.2, 0.6, 0.4, 0.85);
      else if(v==1) leg = new TLegend(0.4, 0.6, 0.6, 0.85);
      else     leg = new TLegend(0.3, 0.2, 0.5, 0.45);
      leg->SetHeader(Form("|vz| < %d cm, p_{T} > 1 GeV/c",gVzCut));
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      for(int j=0; j<nCentBins; j++)
	{
	  cVar[v]->cd(j+1);
	  if(v==1) gPad->SetLogy();
	  for(int k=2; k<gNTrgSetup; k++)
	    {
	      hTrkVar[v][j][k]->Sumw2();
	      hTrkVar[v][j][k]->SetMarkerStyle(18+k);
	      hTrkVar[v][j][k]->SetMarkerColor(gColor[k-2]);
	      hTrkVar[v][j][k]->SetLineColor(gColor[k-2]);
	      hTrkVar[v][j][k]->Scale(1./hCent[k]->Integral(centBins_low[j], centBins_high[j]));
	      hTrkVar[v][j][k]->SetTitle("");
	      hTrkVar[v][j][k]->SetMaximum(1.2*hTrkVar[v][j][k]->GetMaximum());
	      if(v!=1) hTrkVar[v][j][k]->SetMinimum(0);
	      if(v==1) hTrkVar[v][j][k]->GetXaxis()->SetRangeUser(0,3);
	      if(k==2) hTrkVar[v][j][k]->Draw();
	      else     hTrkVar[v][j][k]->Draw("sames");
	      if(j==0) leg->AddEntry(hTrkVar[v][j][k], gLegNameTrg[k].Data(), "PL");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s_%s (%s, %s%%)",run_type.Data(),trigName, trigTitle,cent_Title[j].Data()), 0.045);
	  t1->Draw();

	  cVar[v]->cd(j+5);
	  for(int k=3; k<gNTrgSetup; k++)
	    {
	      TH1F *hRatio = (TH1F*)hTrkVar[v][j][k]->Clone(Form("%s_ratio",hTrkVar[v][j][k]->GetName()));
	      hRatio->Divide(hTrkVar[v][j][2]);
	      if(v<=1) hRatio->GetYaxis()->SetRangeUser(0,2);
	      else     hRatio->GetYaxis()->SetRangeUser(0.7,1.1);
	      if(k==3) hRatio->Draw();
	      else     hRatio->Draw("sames");
	    }
	}
      cVar[v]->cd(1);
      leg->Draw();
      if(savePlot) cVar[v]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_Comp%s.pdf",run_type.Data(),trigTitle,varName[v]));
    }
  TCanvas *cVarVsZdc[nVar];
  for(int v=0; v<nVar; v++)
    {
      cVarVsZdc[v] = new TCanvas(Form("c%sVsZdc",varName[v]),Form("c%sVzZdc",varName[v]), 1100, 700);
      cVarVsZdc[v]->Divide(2,2);
      TLegend *leg = NULL;
      if(v==0 || v==3) leg = new TLegend(0.5, 0.6, 0.7, 0.85);
      else leg = new TLegend(0.5, 0.2, 0.7, 0.45);
      leg->SetHeader(Form("|vz| < %d cm, p_{T} > 1 GeV/c",gVzCut));
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      for(int j=0; j<nCentBins; j++)
	{
	  cVarVsZdc[v]->cd(j+1);
	  for(int k=2; k<gNTrgSetup; k++)
	    {
	      hTrkVarPro[v][j][k]->Sumw2();
	      hTrkVarPro[v][j][k]->SetMarkerStyle(18+k);
	      hTrkVarPro[v][j][k]->SetMarkerColor(gColor[k-2]);
	      hTrkVarPro[v][j][k]->SetLineColor(gColor[k-2]);
	      hTrkVarPro[v][j][k]->SetTitle(Form(";ZDC [kHz];%s",varTitle[v]));
	      hTrkVarPro[v][j][k]->SetMaximum(1.2*hTrkVarPro[v][j][k]->GetMaximum());
	      if(v==0) 
		{
		  if(j==0) hTrkVarPro[v][j][k]->GetYaxis()->SetRangeUser(28, 32);
		  else     hTrkVarPro[v][j][k]->GetYaxis()->SetRangeUser(30, 40);
		}
	      if(v==1) hTrkVarPro[v][j][k]->GetYaxis()->SetRangeUser(0.1, 0.7);
	      if(v==2) hTrkVarPro[v][j][k]->GetYaxis()->SetRangeUser(0, 0.03);
	      if(v==3) hTrkVarPro[v][j][k]->GetYaxis()->SetRangeUser(3.0, 3.2);
	      if(k==2) hTrkVarPro[v][j][k]->Draw();
	      else     hTrkVarPro[v][j][k]->Draw("sames");
	      if(j==0) leg->AddEntry(hTrkVarPro[v][j][k], gLegNameTrg[k].Data(), "PL");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s_%s (%s, %s%%)",run_type.Data(),trigName, trigTitle,cent_Title[j].Data()), 0.055);
	  t1->Draw();
	}
      cVarVsZdc[v]->cd(1);
      leg->Draw();
      if(savePlot) cVarVsZdc[v]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_Comp%sVsZdc.pdf",run_type.Data(),trigTitle,varName[v]));
    }

  // check track phi and eta vs. ZDC rate
  const int nCentBins3 = 4; 
  int centBins_low3[nCentBins3]  = {15,11,7,3};
  int centBins_high3[nCentBins3] = {16,12,8,4};
  TString cent_Name3[nCentBins3];
  TString cent_Title3[nCentBins3];
  for(int j=0; j<nCentBins3; j++)
    {
      cent_Name3[j] = Form("%d%d",(16-centBins_high3[j])*5,(17-centBins_low3[j])*5);
      cent_Title3[j] = Form("%d-%d",(16-centBins_high3[j])*5,(17-centBins_low3[j])*5);
    }
  TH3F *hTrkPtVsEtaPhiVsZdc[2][nCentBins3];
  TH2F *hTrkEtaVsPhi[nCentBins3][4];
  hnTrkEtaPhi->GetAxis(0)->SetRangeUser(1.5,5);
  for(int j=0; j<nCentBins3; j++)
    {
      hnTrkEtaPhi->GetAxis(3)->SetRange(centBins_low3[j], centBins_high3[j]);
      TH3F* h3tmp = (TH3F*)hnTrkEtaPhi->Projection(2,1,6);
      h3tmp->SetName(Form("%s_%d",h3tmp->GetName(),j));
      h3tmp->Rebin3D(1,1,2);
      for(int z=0; z<4; z++)
	{
	  h3tmp->GetZaxis()->SetRange(z+2,z+2);
	  hTrkEtaVsPhi[j][z] = (TH2F*)h3tmp->Project3D("yx");
	  hTrkEtaVsPhi[j][z]->SetName(Form("hTrkEtaVsPhi_%s_Zdc%d",cent_Name3[j].Data(),z));
	}
      for(int v=0; v<2; v++)
	{
	  if(v==0) hnTrkEtaPhi->GetAxis(2)->SetRangeUser(0,5);
	  hTrkPtVsEtaPhiVsZdc[v][j] = (TH3F*)hnTrkEtaPhi->Projection(0,v+1,6);
	  hTrkPtVsEtaPhiVsZdc[v][j]->SetName(Form("hTrkPtVs%sVsZdc_%s",varName[v+2],cent_Name3[j].Data()));
	  hTrkPtVsEtaPhiVsZdc[v][j]->Sumw2();
	  hnTrkEtaPhi->GetAxis(2)->SetRange(0,-1);
	}
      hnTrkEtaPhi->GetAxis(3)->SetRange(0,-1);
    }
  hnTrkEtaPhi->GetAxis(7)->SetRange(0,-1);

  TCanvas *c = new TCanvas("cTrkPhiVsEta","cTrkPhiVsEta",1100,800);
  c->Divide(4,4);
  for(int j=0; j<nCentBins3; j++)
    {
      for(int z=0; z<4; z++)
	{
	  cTrkPhiVsEta->cd(z*nCentBins3+j+1);
	  double scale = hTrkEtaVsPhi[j][z]->Integral(0,-1,18,18)/hTrkEtaVsPhi[j][z]->GetNbinsX(); // 0.6 < eta < 0.8
	  //hTrkEtaVsPhi[j][z]->Scale(1./scale);
	  //hTrkEtaVsPhi[j][z]->GetZaxis()->SetRangeUser(0,1.05);
	  hTrkEtaVsPhi[j][z]->SetTitle("");
	  hTrkEtaVsPhi[j][z]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("%s%%, %d < ZDC < %d kHz",cent_Title3[j].Data(),20+z*20,40+z*20), 0.07);
	  t1->Draw();
	}
    }
  return;

  // compare with MC
  THnSparseF *hnEmbTrkEtaPhi = (THnSparseF*)femb->Get("hTrkEtaPhi_MCreco_mb");
  TH2F *hEmbTrkEtaVsPhi[nCentBins3][gNTrgSetup];
  for(int j=0; j<nCentBins3; j++)
    {
      hnEmbTrkEtaPhi->GetAxis(3)->SetRange(centBins_low3[j], centBins_high3[j]);
      TH3F *h3tmp = (TH3F*)hnEmbTrkEtaPhi->Projection(1,2,4);
      h3tmp->SetName(Form("h3tmp_%d",j));
      for(int k=2; k<gNTrgSetup; k++)
	{
	  if(k>0) h3tmp->GetZaxis()->SetRange(k,k);
	  hEmbTrkEtaVsPhi[j][k] = (TH2F*)h3tmp->Project3D("xy");
	  hEmbTrkEtaVsPhi[j][k]->SetName(Form("hEmbTrkEtaVsPhi_%s%s",cent_Name3[j].Data(),gTrgSetupTitle[k]));
	}
    }
  TCanvas *c = new TCanvas("cEmbTrkPhiVsEta","cEmbTrkPhiVsEta",1100,800);
  c->Divide(4,3);
  for(int j=0; j<nCentBins3; j++)
    {
      for(int k=2; k<gNTrgSetup; k++)
	{
	  c->cd((k-2)*nCentBins3+j+1);
	  hEmbTrkEtaVsPhi[j][k]->Rebin2D(10,1);
	  double scale = hEmbTrkEtaVsPhi[j][k]->Integral(0,-1,15,15)/hEmbTrkEtaVsPhi[j][k]->GetNbinsX(); // 0.6 < eta < 0.8
	  hEmbTrkEtaVsPhi[j][k]->Scale(1./scale);
	  hEmbTrkEtaVsPhi[j][k]->GetZaxis()->SetRangeUser(0,1.05);
	  hEmbTrkEtaVsPhi[j][k]->SetTitle("");
	  hEmbTrkEtaVsPhi[j][k]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("%s%%, %s",cent_Title3[j].Data(),gLegNameTrg[k].Data()), 0.07);
	  t1->Draw();
	}
    }

  const int nPtBins = 3;
  const double xPtBins[nPtBins+1] = {1, 2, 3, 5};
  const int gNZdcBins = 11;
  TH1F *hTrkEtaPhi[2][nCentBins3][nPtBins][gNZdcBins];
  TCanvas *cTrkEtaVsPtVsZdc[2];
  for(int v=0; v<2; v++)
    {
      cTrkEtaVsPtVsZdc[v] = new TCanvas(Form("cTrk%sVsPtVsZdc",varName[v+2]),Form("cTrk%sVsPtVsZdc",varName[v+2]),1100,800);
      cTrkEtaVsPtVsZdc[v]->Divide(4,3);
      for(int j=0; j<nCentBins3; j++)
	{
	  for(int p=0; p<nPtBins; p++)
	    {
	      int low_pt_bin  = hTrkPtVsEtaPhiVsZdc[v][j]->GetXaxis()->FindBin(xPtBins[p]+1e-4);
	      int high_pt_bin = hTrkPtVsEtaPhiVsZdc[v][j]->GetXaxis()->FindBin(xPtBins[p+1]-1e-4);
	      cTrkEtaVsPtVsZdc[v]->cd(p*nCentBins+j+1);
	      for(int z=0; z<gNZdcBins; z++)
		{
		  hTrkEtaPhi[v][j][p][z] = (TH1F*)hTrkPtVsEtaPhiVsZdc[v][j]->ProjectionY(Form("hTrk%s_%d%d%d",varName[v+2],j,p,z),low_pt_bin,high_pt_bin,z+1,z+1);
		  hTrkEtaPhi[v][j][p][z]->SetMarkerStyle(20+z);
		  hTrkEtaPhi[v][j][p][z]->SetMarkerColor(gColor[z%6]);
		  hTrkEtaPhi[v][j][p][z]->SetLineColor(gColor[z%6]);
		  double scale = hTrkEtaPhi[v][j][p][z]->Integral();
		  if(scale>0) hTrkEtaPhi[v][j][p][z]->Scale(1./scale);
		  if(z==0) hTrkEtaPhi[v][j][p][z]->Draw();
		  else     hTrkEtaPhi[v][j][p][z]->Draw("sames");
		}
	    }
	}
    }
  
  // compare track yield
  TH2F *hTrkPtVsCent[gNTrgSetup];
  TH1F *hTrkYieldVsZdc[nCentBins][gNTrgSetup];
  hnTrkEtaPhi->GetAxis(7)->SetRange(3,3); // DCA < 1 cm
  hnTrkEtaPhi->GetAxis(2)->SetRangeUser(0,5); // 0 < phi < 5
  hnTrkEtaPhi->GetAxis(1)->SetRangeUser(0.4,0.8); // 0 < phi < 5
  hnTrkEtaPhi->GetAxis(0)->SetRangeUser(1.5,5); // pt>2
  hnTrkEtaPhi->GetAxis(5)->SetRange(1,2);
  hnTrkEtaPhi->GetAxis(8)->SetRangeUser(-1*gVzCut, gVzCut);
  for(int k=2; k<gNTrgSetup; k++)
    {
      hnTrkEtaPhi->GetAxis(4)->SetRange(k, k);
      hTrkPtVsCent[k] = (TH2F*)hnTrkEtaPhi->Projection(0,3);
      hTrkPtVsCent[k]->Sumw2();
      hTrkPtVsCent[k]->SetName(Form("hTrkPtVsCent%s",gTrgSetupTitle[k]));
      for(int j=0; j<nCentBins; j++)
	{
	  hnTrkEtaPhi->GetAxis(3)->SetRange(centBins_low[j], centBins_high[j]);
	  hTrkYieldVsZdc[j][k] = (TH1F*)hnTrkEtaPhi->Projection(6);
	  hTrkYieldVsZdc[j][k]->Sumw2();
	  hTrkYieldVsZdc[j][k]->SetName(Form("hTrkYieldVsZdc_%s%s",cent_Name[j].Data(),gTrgSetupTitle[k]));
	  for(int bin=1; bin<=hTrkYieldVsZdc[j][k]->GetNbinsX(); bin++)
	    {
	      double nevent = hCentVsZdc[k]->Integral((bin-1)*2+1,(bin-1)*2+2,centBins_low[j], centBins_high[j]);
	      if(nevent<=0) continue;
	      hTrkYieldVsZdc[j][k]->SetBinContent(bin, hTrkYieldVsZdc[j][k]->GetBinContent(bin)/nevent);
	      hTrkYieldVsZdc[j][k]->SetBinError(bin, hTrkYieldVsZdc[j][k]->GetBinError(bin)/nevent);
	    }
	  hnTrkEtaPhi->GetAxis(3)->SetRange(0, -1);
	}
      hnTrkEtaPhi->GetAxis(4)->SetRange(0, -1);
    }
 
  // raw distribution
  TCanvas *c = new TCanvas(Form("cTrkRawYield_%s",trigName), Form("cTrkRawYield_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.5, 0.6, 0.7, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader(Form("|vz| < %d cm, p_{T} > 1 GeV/c",gVzCut));
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      for(int k=2; k<gNTrgSetup; k++)
	{
	  hTrkYieldVsZdc[j][k]->SetMarkerStyle(18+k);
	  hTrkYieldVsZdc[j][k]->SetMarkerColor(gColor[k-2]);
	  hTrkYieldVsZdc[j][k]->SetLineColor(gColor[k-2]);
	  hTrkYieldVsZdc[j][k]->SetMarkerSize(1.2);
	  hTrkYieldVsZdc[j][k]->SetTitle(";ZDC [kHz];dN/N_{evt}");
	  hTrkYieldVsZdc[j][k]->SetMaximum(hTrkYieldVsZdc[j][4]->GetBinContent(2)*1.2);
	  hTrkYieldVsZdc[j][k]->SetMinimum(hTrkYieldVsZdc[j][4]->GetBinContent(11)/1.2);
	  if(k==0) hTrkYieldVsZdc[j][k]->Draw();
	  else     hTrkYieldVsZdc[j][k]->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkYieldVsZdc[j][k], gLegNameTrg[k].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s_%s (%s, %s%%)",run_type.Data(),trigName,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkRawVsZdc.pdf",run_type.Data(),trigTitle));

  // tracking efficiency from embedding
  TH1F *hTrkEffVsZdc[nCentBins][gNTrgSetup];
  THnSparseF *hnEmbTrk[2];
  hnEmbTrk[0] = (THnSparseF*)femb->Get(Form("mhMcTrkPtEff_Tpc_%s",trigName));
  hnEmbTrk[1] = (THnSparseF*)femb->Get(Form("mhMcTrkPtEff_MC_%s",trigName));
  TH3F *hEmbTrkPtVsCentVsZdc[2][gNTrgSetup];
  for(int i=0; i<2; i++)
    {
      hnEmbTrk[i]->GetAxis(1)->SetRangeUser(-0.8, 0.8); // phi cut
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hnEmbTrk[i]->GetAxis(4)->SetRange(k, k);
	  hEmbTrkPtVsCentVsZdc[i][k] = (TH3F*)hnEmbTrk[i]->Projection(0,2,3);
	  hEmbTrkPtVsCentVsZdc[i][k]->SetName(Form("hEmbTrkPtVsCentVsZdc%s_%d",gTrgSetupTitle[k],i));
	  hEmbTrkPtVsCentVsZdc[i][k]->Sumw2();
	  hEmbTrkPtVsCentVsZdc[i][k]->Rebin3D(10,1,1);

	  for(int j=0; j<nCentBins; j++)
	    {
	      TH1F *h1tmp = (TH1F*)hEmbTrkPtVsCentVsZdc[i][k]->ProjectionZ(Form("hEmbTrkVsZdc_%s%s_%d",cent_Name[j].Data(),gTrgSetupTitle[k],i),2,10,centBins_low[j],centBins_high[j]);
	      if(i==0) hTrkEffVsZdc[j][k] = (TH1F*)h1tmp->Clone(Form("hTrkEffVsZdc_%s%s",cent_Name[j].Data(),gTrgSetupTitle[k]));
	      else     hTrkEffVsZdc[j][k]->Divide(h1tmp);
	    }
	  hnEmbTrk[i]->GetAxis(4)->SetRange(0, -1);
	}
    }
  //return;

  TCanvas *c = new TCanvas(Form("cTrkEff_%s",trigName), Form("cTrkEff_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.2, 0.65, 0.8, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("p_{T} > 1 GeV/c");
  leg->SetNColumns(3);
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  hTrkEffVsZdc[j][k]->SetMarkerStyle(20+k);
	  hTrkEffVsZdc[j][k]->SetMarkerColor(gColor[k]);
	  hTrkEffVsZdc[j][k]->SetLineColor(gColor[k]);
	  hTrkEffVsZdc[j][k]->SetMarkerSize(1.2);
	  hTrkEffVsZdc[j][k]->GetXaxis()->SetRangeUser(0,120);
	  hTrkEffVsZdc[j][k]->GetYaxis()->SetRangeUser(0.5,1);
	  hTrkEffVsZdc[j][k]->SetTitle(";ZDC [kHz];Efficiency");
	  if(k==0) hTrkEffVsZdc[j][k]->Draw();
	  else     hTrkEffVsZdc[j][k]->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkEffVsZdc[j][k], gLegNameTrg[k].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s embedding, %s (%s%%)",run_type.Data(),trigName,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkEffVsZdc.pdf",run_type.Data(),trigTitle));

  // corrected yield
  TCanvas *c = new TCanvas(Form("cTrkYieldCorr_%s",trigName), Form("cTrkYieldCorr_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.5, 0.6, 0.7, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader(Form("|vz| < %d cm, p_{T} > 1 GeV/c",gVzCut));
  TH1F *hTrkPtCorr[nCentBins][gNTrgSetup];
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      for(int k=2; k<gNTrgSetup; k++)
	{
	  hTrkPtCorr[j][k] = (TH1F*)hTrkYieldVsZdc[j][k]->Clone(Form("hTrkPtCorr_%s%s",cent_Name[j].Data(),gTrgSetupTitle[k]));
	  for(int bin=1; bin<=hTrkPtCorr[j][k]->GetNbinsX(); bin++)
	    {
	      double eff = hTrkEffVsZdc[j][0]->GetBinContent(bin);
	      if(eff<=0 || hTrkPtCorr[j][k]->GetBinContent(bin)<=0) continue;
	      if(j==1) 
		{
		  printf("[i] %s, bin %d has %4.2f%% eff\n",hTrkEffVsZdc[j][k]->GetName(),bin, eff*100);
		}
	      double err_eff = hTrkEffVsZdc[j][0]->GetBinError(bin)/eff;
	      double err_data = hTrkPtCorr[j][k]->GetBinError(bin)/hTrkPtCorr[j][k]->GetBinContent(bin);
	      hTrkPtCorr[j][k]->SetBinContent(bin, hTrkPtCorr[j][k]->GetBinContent(bin)/eff);
	      hTrkPtCorr[j][k]->SetBinError(bin, sqrt(err_eff*err_eff+err_data*err_data)*hTrkPtCorr[j][k]->GetBinContent(bin));
	    }
	}
      for(int k=2; k<gNTrgSetup; k++)
	{
	  hTrkPtCorr[j][k]->SetMaximum(hTrkPtCorr[j][4]->GetBinContent(2)*1.2);
	  hTrkPtCorr[j][k]->SetMinimum(hTrkPtCorr[j][4]->GetBinContent(11)/1.2);
	  if(k==0) hTrkPtCorr[j][k]->Draw();
	  else     hTrkPtCorr[j][k]->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkPtCorr[j][k], gLegNameTrg[k].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s_%s (%s, %s%%)",run_type.Data(),trigName,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkYieldCorrVsZdc.pdf",run_type.Data(),trigTitle));

  // compare pt distribution
  const int nCentBins2 = 4;
  // const int centBinsLow2[nCentBins2]  = {16,12,8,4};
  // const int centBinsHigh2[nCentBins2] = {16,12,8,4};
  const int centBinsLow2[nCentBins2]  = {13,9,5,1};
  const int centBinsHigh2[nCentBins2] = {16,12,8,4};
  TString centName2[nCentBins2];
  TString centTitle2[nCentBins2];
  for(int j=0; j<nCentBins2; j++)
    {
      centName2[j] = Form("%d%d",(16-centBinsHigh2[j])*5,(17-centBinsLow2[j])*5);
      centTitle2[j] = Form("%d-%d",(16-centBinsHigh2[j])*5,(17-centBinsLow2[j])*5);
    }
  const double ncoll[16] = {10.54, 15.98, 23.96, 35.64, 51.62, 72.79, 101.51, 138.71, 186.69, 246.95, 320.78, 411.86, 524.31, 663.03, 838.41, 1048.11}; 
  TH1F *hNEvents[gNTrgSetup];
  TH1F *hEmbTrkPt[nCentBins2][gNTrgSetup][2];
  TH1F *hEmbTrkEff[nCentBins2][gNTrgSetup];
  TF1  *funcEmbTrkEff[nCentBins2][gNTrgSetup];
  TH1F *hDataTrkPt[nCentBins2][gNTrgSetup];
  for(int k=2; k<gNTrgSetup; k++)
    {
      // efficiency from embedding
      TH2F *hWeight = (TH2F*)hCentVsZdc[k]->Clone(Form("%s_clone",hCentVsZdc[k]->GetName()));
      hWeight->Rebin2D(2,1);
      hWeight->Scale(1./hWeight->Integral());
      for(int j=0; j<nCentBins2; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hEmbTrkPt[j][k][i] = (TH1F*)hEmbTrkPtVsCentVsZdc[i][k]->ProjectionX(Form("hEmbTrkPt_%s%s%d",centName2[j].Data(),gTrgSetupTitle[k],i));
	      hEmbTrkPt[j][k][i]->Reset("AC");
	    }
	  double totalScale = 0;
	  const int nbinsZ = hWeight->GetNbinsX();
	  cout << hEmbTrkPt[j][k][0]->GetName() << endl;
	  for(int icent=centBinsLow2[j]; icent<=centBinsHigh2[j]; icent++)
	    {
	      for(int zbin=1; zbin<=nbinsZ; zbin++)
		{
		  double weight = hWeight->GetBinContent(zbin, icent) * ncoll[icent-1];
		  totalScale += weight;
		  for(int i=0; i<2; i++)
		    {
		      TH1F *htmp = (TH1F*)hEmbTrkPtVsCentVsZdc[i][0]->ProjectionX(Form("htmp_%d%d%d%d",k,icent,zbin,i),icent,icent,zbin,zbin);
		      hEmbTrkPt[j][k][i]->Add(htmp, weight);
		    }
		}
	    }
	  for(int i=0; i<2; i++)
	    {
	      hEmbTrkPt[j][k][i]->Scale(1./totalScale);
	    }
	  hEmbTrkEff[j][k] = (TH1F*)hEmbTrkPt[j][k][0]->Clone(Form("hEmbTrkEff_%s%s",centName2[j].Data(),gTrgSetupTitle[k]));
	  hEmbTrkEff[j][k]->Divide(hEmbTrkPt[j][k][1]);
	  funcEmbTrkEff[j][k] = new TF1(Form("funcEmbTrkEff_%s%s",centName2[j].Data(),gTrgSetupTitle[k]), "pol0", 2, 10);
	  hEmbTrkEff[j][k]->Fit(funcEmbTrkEff[j][k], "0RQ");
	  if(j==1)
	    {
	      c = draw1D(hEmbTrkEff[j][k]);
	      funcEmbTrkEff[j][k]->SetLineColor(2);
	      funcEmbTrkEff[j][k]->Draw("sames");
	    }
	}

      // # of events
      hNEvents[k] = (TH1F*)hCentVsZdc[k]->ProjectionY(Form("hNEvents%s",gTrgSetupTitle[k]));

      // track pt distribution in data
      for(int j=0; j<nCentBins2; j++)
	{
	  hDataTrkPt[j][k] = (TH1F*)hTrkPtVsCent[k]->ProjectionY(Form("hDataTrkPt_%s%s",centName2[j].Data(),gTrgSetupTitle[k]),centBinsLow2[j],centBinsHigh2[j]);
	  hDataTrkPt[j][k]->Scale(1./hNEvents[k]->Integral(centBinsLow2[j],centBinsHigh2[j]));
	  if(j==1)
	    {
	      printf("[i] %s has %4.2f entries with %4.2f%% eff\n",hDataTrkPt[j][k]->GetName(),hDataTrkPt[j][k]->GetEntries(),funcEmbTrkEff[j][k]->GetParameter(0)*100);
	    }
	  for(int bin=1; bin<=hDataTrkPt[j][k]->GetNbinsX(); bin++)
	    {
	      double eff = funcEmbTrkEff[j][k]->GetParameter(0);
	      double val = hDataTrkPt[j][k]->GetBinContent(bin);
	      hDataTrkPt[j][k]->SetBinContent(bin, val/eff);
	      hDataTrkPt[j][k]->SetBinError(bin, hDataTrkPt[j][k]->GetBinError(bin)/eff);
	    }
	}
    }

  // corrected pt distribution
  TCanvas *c = new TCanvas(Form("cTrkPtCorr_%s",trigName), Form("cTrkPtCorr_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TH1F *hTrkPtCorr[nCentBins][gNTrgSetup];
  for(int j=0; j<nCentBins2; j++)
    {
      TLegend *leg = new TLegend(0.15, 0.15, 0.35, 0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader(Form("|vz| < %d cm",gVzCut));
      c->cd(j+1);
      for(int k=3; k<gNTrgSetup; k++)
	{
	  TH1F *hRatio = (TH1F*)hDataTrkPt[j][k]->Clone(Form("%s_ratio",hDataTrkPt[j][k]->GetName()));
	  hRatio->Divide(hDataTrkPt[j][2]);
	  hRatio->SetMarkerStyle(18+k);
	  hRatio->SetMarkerColor(gColor[k-2]);
	  hRatio->SetLineColor(gColor[k-2]);
	  hRatio->GetYaxis()->SetRangeUser(0.7, 1.2);
	  hRatio->GetXaxis()->SetRangeUser(1,8);
	  hRatio->SetTitle("");
	  if(k==3) hRatio->Draw();
	  else     hRatio->Draw("sames");

	   TF1 *func = new TF1(Form("%s_fit",hDataTrkPt[j][k]->GetName()), "pol0", 1, 8);
	   hRatio->Fit(func, "IR0Q");
	   func->SetLineColor(gColor[k-2]);
	   func->SetLineStyle(2);
	   func->Draw("sames");
	   leg->AddEntry(func, Form("%s: %2.3f #pm %2.3f",gLegNameTrg[k].Data(),func->GetParameter(0),func->GetParError(0)), "L");
	}
      TPaveText *t1 = GetTitleText(Form("%s_%s (%s, %s%%)",run_type.Data(),trigName,trigTitle,centTitle2[j].Data()), 0.055);
      t1->Draw();
      leg->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkCorrPtVsZdc.pdf",run_type.Data(),trigTitle));
}


//================================================
void gRefMult(const int mode = 0, const int savePlot = 0)
{
  /// compare the gRefMult distribution between 
  /// MB and dimuon triggers
  const char* dataName[2] = {"VPD-NoVtx", "Dimuon"};
  const char* dataTitle[2] = {"mb", "dm"};

  const int nVzBin = 6;
  const double lowVzBin[nVzBin] = {-100,-100, -60, -20, 20, 60};
  const double upVzBin[nVzBin]  = {100, -60, -20, 20, 60, 100};
  const int nZdcBin = 6;
  const double lowZdcBin[nZdcBin] = {10, 10, 30, 50, 70, 90};
  const double upZdcBin[nZdcBin]  = {110, 30, 50, 70, 90, 110};

  TH1F *hgRefMult[2][gNTrgSetup-1][2];  
  TH1F *hgRefMultDep[2][gNTrgSetup-1][nVzBin][nZdcBin];
  if(mode==0)
    {
      TFile *fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyTracking.root", "read");
      for(int i=0; i<2; i++)
	{
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hgRefMult[i][k][j] = (TH1F*)fin->Get(Form("hgRefMult%s_DM%d_%s",gTrgSetupTitle[k+1],j,dataName[i]));
		}
	      for(int v=0; v<nVzBin; v++)
		{
		  for(int z=0; z<nZdcBin; z++)
		    {
		      hgRefMultDep[i][k][v][z] = (TH1F*)fin->Get(Form("hgRefMult_Vz%d_Zdc%d%s_%s",v,z,gTrgSetupTitle[k+1],dataName[i]));
		    }
		}
	    }
	}

      // check the vz and luminosity dependence in MB
      const char* depName[2] = {"Vz", "Zdc"};
      const char* depUnit[2] = {"cm", "kHz"};
      TCanvas *cRefMult[2][2];
      TH1F *htmp = 0x0;
      TH1F *hplot = new TH1F("hplot",";gRefMultCorr;",750,0,750);
      for(int i=0; i<2; i++)
	{
	  for(int d=0; d<2; d++)
	    {
	      cRefMult[i][d] = new TCanvas(Form("cRefMult_%s_%s",depName[d],dataName[i]), Form("cRefMult_%s_%s",depName[d],dataName[i]), 1100, 700);
	      cRefMult[i][d]->Divide(2,2);
	      TLegend *leg =  new TLegend(0.3, 0.15, 0.5, 0.4);
	      leg->SetBorderSize(0);
	      leg->SetFillColor(0);
	      leg->SetTextSize(0.045);
	      for(int k=0; k<gNTrgSetup-1; k++)
		{
		  cRefMult[i][d]->cd(k+1);
		  gPad->SetLogy();
		  hplot->GetYaxis()->SetRangeUser(1e-6, 5e-2);
		  hplot->DrawCopy();
		  int nHisto = nVzBin-1;
		  if(d==1) nHisto = nZdcBin-1;
		  for(int ih=0; ih<nHisto; ih++)
		    {
		      if(d==0) htmp = (TH1F*)fin->Get(Form("hgRefMult_Vz%d_Zdc%d%s_%s",ih+1,0,gTrgSetupTitle[k+1],dataName[i]));
		      if(d==1) htmp = (TH1F*)fin->Get(Form("hgRefMult_Vz%d_Zdc%d%s_%s",0,ih+1,gTrgSetupTitle[k+1],dataName[i]));
		      htmp = (TH1F*)htmp->Clone(Form("%s_clone",htmp->GetName()));
		      htmp->Rebin(5);
		      if(i==1)
			{
			  double hvalue = 550.46;
			  if(k==3) hvalue = 558.41;
			  htmp->Scale(5e-3/htmp->GetBinContent(htmp->FindBin(hvalue)));
			}
		      else
			{
			  if(htmp->Integral()>0) htmp->Scale(1./htmp->Integral());
			}
		      htmp->SetMarkerStyle(20+ih);
		      htmp->SetMarkerColor(gColor[ih]);
		      htmp->SetLineColor(gColor[ih]);
		      htmp->DrawCopy("sames");
		      if(k==0)
			{
			  if(d==0) leg->AddEntry(htmp, Form("%1.0f < v_{z} < %1.0f cm",lowVzBin[ih+1],upVzBin[ih+1]));
			  if(d==1) leg->AddEntry(htmp, Form("%1.0f < ZDC < %1.0f kHz",lowZdcBin[ih+1],upZdcBin[ih+1]));
			}
		    }
		  TPaveText *t1 = GetTitleText(Form("%s: %s%s",dataName[i],run_type.Data(),gTrgSetupTitle[k+1]), 0.055);
		  t1->Draw();
		}
	      cRefMult[i][d]->cd(1);
	      leg->Draw();
	      if(savePlot) cRefMult[i][d]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_RefMultIn%s_%s.pdf",run_type.Data(),depName[d],dataTitle[i]));
	    } 
	}
      //return;

      // compare MB and MTD
      TH1F *hgRefMultRebin[2][gNTrgSetup-1][2];
      TH1F *hgRefMultScale[gNTrgSetup-1];
      TCanvas *cRefMultMtd = new TCanvas("cRefMultMtd", "cRefMultMtd", 1100, 700);
      cRefMultMtd->Divide(2,2);
      TString legName[4] = {"VPD-NoVtx", "VPD-NoVtx + dimuon", "Dimuon", "Dimuon scaled"};
      TLegend *leg =  new TLegend(0.3, 0.15, 0.5, 0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader("|v_{z}| < 100 cm");
      const double hvalue1 = 530.206, hvalue2 = 550.46, hvalue3 = 558.41;
      TF1 *func = 0x0;
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  cRefMultMtd->cd(k+1);
	  gPad->SetLogy();
	  hplot->DrawCopy();
	  double hvalue = 0;
	  for(int i=0; i<2; i++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hgRefMultRebin[i][k][j] = (TH1F*)hgRefMult[i][k][j]->Clone(Form("%s_clone",hgRefMult[i][k][j]->GetName()));
		  hgRefMultRebin[i][k][j]->Rebin(5);
		  hgRefMultRebin[i][k][j]->Scale(0.2);
		  if(hgRefMultRebin[i][k][j]->Integral()>0) hgRefMultRebin[i][k][j]->Scale(1./hgRefMultRebin[i][k][j]->Integral());
		  if(i==0 && j==0)
		    {
		      func = new TF1(Form("fit_%s",hgRefMultRebin[i][k][j]->GetName()),"[0]*TMath::Erf(-1*[1]*(x-[2]))+[0]", 450,700);
		      func->SetParameters(0.001, 10, 565);
		      func->SetLineColor(8);
		      func->SetLineStyle(2);
		      hgRefMultRebin[i][k][j]->Fit(func,"IR0");
		    }
		  if(i==1 && j==1) continue;
		  if(i!=0 || j!=0)
		    {
		      int bin = hgRefMultRebin[i][k][j]->FindBin(hvalue1);
		      if(k==gNTrgSetup-2) bin = hgRefMultRebin[i][k][j]->FindBin(hvalue3);
		      hgRefMultRebin[i][k][j]->Scale(hgRefMultRebin[0][k][0]->GetBinContent(bin)/hgRefMultRebin[i][k][j]->GetBinContent(bin));
		    }
		  hgRefMultRebin[i][k][j]->SetMarkerStyle(20+(i*2+j)*2);
		  hgRefMultRebin[i][k][j]->SetMarkerColor(gColor[i*2+j]);
		  hgRefMultRebin[i][k][j]->SetLineColor(gColor[i*2+j]);
		  hgRefMultRebin[i][k][j]->DrawCopy("sames");
		  if(k==0)
		    {
		      leg->AddEntry(hgRefMultRebin[i][k][j], legName[i*2+j].Data(), "P");
		    }
		  if(i==0 && j==0) func->Draw("sames");
		  if(i==1 && k<gNTrgSetup-2)
		    {
		      hgRefMultScale[k] = (TH1F*)hgRefMult[i][k][j]->Clone(Form("%s_scale",hgRefMult[i][k][j]->GetName()));
		      hgRefMultScale[k]->Reset();
		      int nbins = hgRefMult[i][k][j]->GetNbinsX();
		      for(int ibin=1; ibin<=nbins; ibin++)
			{
			  double counts = hgRefMult[i][k][j]->GetBinContent(ibin)/100;
			  double low = hgRefMult[i][k][j]->GetXaxis()->GetBinLowEdge(ibin);
			  double up  = hgRefMult[i][k][j]->GetXaxis()->GetBinUpEdge(ibin);
			  for(int iexpr=0; iexpr<counts; iexpr++)
			    {
			      hgRefMultScale[k]->Fill(gRandom->Uniform(low,up)*hvalue1/hvalue2);
			    }
			}
	
		      hgRefMultScale[k]->Rebin(5);
		      int bin = hgRefMultRebin[0][k][0]->FindBin(hvalue1);
		      hgRefMultScale[k]->Scale(hgRefMultRebin[0][k][0]->GetBinContent(bin)/hgRefMultScale[k]->GetBinContent(bin));
		      hgRefMultScale[k]->SetMarkerStyle(34);
		      hgRefMultScale[k]->SetMarkerColor(6);
		      hgRefMultScale[k]->SetLineColor(6);
		      hgRefMultScale[k]->DrawCopy("sames");
		      if(k==0)
			{
			  leg->AddEntry(hgRefMultScale[k], legName[3].Data(), "P");
			}
		    }
		}
	  TPaveText *t1 = GetTitleText(Form("%s%s",run_type.Data(),gTrgSetupTitle[k+1]), 0.055);
	  t1->Draw();
	    }
	}
      cRefMultMtd->cd(1);
      leg->Draw();
      if(savePlot) cRefMultMtd->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_RefMultComp_DmVsMb.pdf",run_type.Data()));
    }

  if(mode==1)
    {
      TFile *fin[2];
      fin[0] = TFile::Open("output/Run14_AuAu200.MB.Study.TpcTracking.root", "read");
      fin[1] = TFile::Open("output/Run14_AuAu200.Study.TpcTracking.root", "read");
      THnSparseF *hnRefMult[2];
      hnRefMult[0] = (THnSparseF*)fin[0]->Get("mhRefMultStudy_mb");
      hnRefMult[1] = (THnSparseF*)fin[1]->Get("mhRefMultStudy_di_mu");
      for(int i=0; i<2; i++)
	{
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hnRefMult[i]->GetAxis(4)->SetRange(k+1, k+1);

	      for(int j=0; j<2; j++)
		{
		  if(j==1) hnRefMult[i]->GetAxis(5)->SetRange(2,2);
		  hgRefMult[i][k][j] = (TH1F*)hnRefMult[i]->Projection(3);
		  hgRefMult[i][k][j]->SetName(Form("hgRefMult%s_DM%d_%s",gTrgSetupTitle[k+1],j,dataName[i]));
		  hgRefMult[i][k][j]->Sumw2();
		  hgRefMult[i][k][j]->SetTitle("");
		  hnRefMult[i]->GetAxis(5)->SetRange(0,-1);
		}

	      for(int v=0; v<nVzBin; v++)
		{
		  hnRefMult[i]->GetAxis(0)->SetRangeUser(lowVzBin[v]+1e-4, upVzBin[v]-1e-4);
		  for(int z=0; z<nZdcBin; z++)
		    {
		      hnRefMult[i]->GetAxis(1)->SetRangeUser(lowZdcBin[z]+1e-4, upZdcBin[z]-1e-4);
		      hgRefMultDep[i][k][v][z] = (TH1F*)hnRefMult[i]->Projection(3);
		      hgRefMultDep[i][k][v][z]->SetName(Form("hgRefMult_Vz%d_Zdc%d%s_%s",v,z,gTrgSetupTitle[k+1],dataName[i]));
		      hgRefMultDep[i][k][v][z]->Sumw2();
		      hgRefMultDep[i][k][v][z]->SetTitle("");
		      hnRefMult[i]->GetAxis(1)->SetRange(0,-1);
		    }
		  hnRefMult[i]->GetAxis(0)->SetRange(0,-1);
		}
	      hnRefMult[i]->GetAxis(4)->SetRange(0,-1);
	    }
	}
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyTracking.root", "update");
      for(int i=0; i<2; i++)
	{
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hgRefMult[i][k][j]->Write("",TObject::kOverwrite);
		}
	      for(int v=0; v<nVzBin; v++)
		{
		  for(int z=0; z<nZdcBin; z++)
		    {
		      hgRefMultDep[i][k][v][z]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
      fout->Close();
    }
  
  
}


//================================================
void trackQA(const int savePlot = 0)
{
  const int nCentBins       = 8; 
  TString cent_Name[nCentBins];
  TString cent_Title[nCentBins];
  const int centBins_low[nCentBins]  = {1,5,7,9,11,13,15,16};
  const int centBins_high[nCentBins] = {4,6,8,10,12,14,15,16}; 
  for(int j=0; j<nCentBins; j++)
    {
      cent_Name[j]  = Form("%d-%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
      cent_Title[j] = Form("%d%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
    }
  const int nType = 2;
  const char* typeName[nType] = {"NHitsFit", "Dca"};
  const double min_pt = 1;
  const double max_pt = 2;


  //----------------------------------------------
  const int nData = 4;
  const char* dataTitle[4] = {"Data_dimu", "Data_mb", "Embed_dimu", "Embed_mb"};
  const char* histoName[4] = {"di_mu", "mb", "di_mu", "mb"};
  const char* saveName = histoName[0];
  TString legNameData[4] = {"Data dimuon: #mu cand.", "Data mb: #mu cand.", "Embed: #mu", "MB_embed: #mu"};
  TFile *fin[4];
  fin[0] = TFile::Open("output/Run14_AuAu200.Study.TpcTracking.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.MB.Study.prod_mid.P15ic.root", "read");
  fin[2] = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");
  fin[3] = TFile::Open("output/Run14_AuAu200.Embed_MB.Jpsi.root", "read");

  /*
  const int nData = 3;
  const char* dataTitle[3] = {"Data_dimu", "Data_mb", "Embed_mb"};
  const char* histoName[3] = {"di_mu", "mb", "mb"};
  const char* saveName = histoName[0];
  TString legNameData[3] = {"Data dimuon: #mu cand.", "Data mb: #mu cand.", "Embed: #mu"};
  TFile *fin[3];
  fin[0] = TFile::Open("output/Run14_AuAu200.Study.TpcTracking.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.MB.Study.TpcTracking.root", "read");
  fin[2] = TFile::Open("output/Run14_AuAu200.Embed_MB.Jpsi.root", "read");
   */

  /*
  const int nData = 2;
  const char* dataTitle[3] = {"Data", "Embed_MC", "Embed_data"};
  const char* histoName = "mb";
  const char* legTitle = "vpd30";
  TString legNameData[3] = {"Data: hadron", "Embed: MC #pi", "Embed: data #pi"};
  TFile *fin[3];
  fin[0] = TFile::Open("output/Run14_AuAu200.MB.P17id.Study.TpcTracking.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
  fin[2] = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
  */
  

  THnSparseF *hnQa[nData][nType];
  TH1F *hTrkDis[nData][nType][nCentBins][gNTrgSetup-1];
  THnSparseF* hnTmp;
  for(int d=0; d<nData; d++)
    {
      for(int i=0; i<nType; i++)
	{
	  if(d<2)      hnQa[d][i] = (THnSparseF*)fin[d]->Get(Form("mhnTrkPt_%s",histoName[d]));
	  else         hnQa[d][i] = (THnSparseF*)fin[d]->Get(Form("hTrk%s_MCreco_%s",typeName[i],histoName[d]));
	  hnQa[d][i]->SetName(Form("%s_%s",dataTitle[d], typeName[i]));
	  if(d<2)
	    {
	      hnQa[d][i]->GetAxis(5)->SetRange(2,2);
	    }

	  hnQa[d][i]->GetAxis(0)->SetRangeUser(min_pt, max_pt);
	  for(int j=0; j<nCentBins; j++)
	    {
	      hnQa[d][i]->GetAxis(3)->SetRange(centBins_low[j], centBins_high[j]);
	      for(int k=0; k<gNTrgSetup-1; k++)
		{
		  hnQa[d][i]->GetAxis(4)->SetRange(k+1, k+1);
		  if(d<2) hTrkDis[d][i][j][k] = (TH1F*)hnQa[d][i]->Projection(i+1);
		  else    hTrkDis[d][i][j][k] = (TH1F*)hnQa[d][i]->Projection(1);
		  hTrkDis[d][i][j][k]->SetName(Form("%s_Trk%s_cent%s%s",dataTitle[d], typeName[i],cent_Title[j].Data(), gTrgSetupTitle[k+1]));
		  if(d>1) hTrkDis[d][i][j][k]->Sumw2();
		  if(hTrkDis[d][i][j][k]->Integral()>0)
		    hTrkDis[d][i][j][k]->Scale(1./hTrkDis[d][i][j][k]->Integral());
		}
	      hnQa[d][i]->GetAxis(4)->SetRange(0, -1);
	    }
	  hnQa[d][i]->GetAxis(3)->SetRange(0, -1);
	  hnQa[d][i]->GetAxis(0)->SetRange(0, -1);
	  if(d<2)
	    {
	      hnQa[d][i]->GetAxis(5)->SetRange(0,-1);
	    }
	}
    }

  /// study the muon candidate vs. charged hadrons
  TH1F *hTrkDisInMb[nType][nCentBins][gNTrgSetup-1][2];
  const char* partName[2] = {"Hadron", "Muon"};
  const int iData = 1;
  for(int i=0; i<nType; i++)
    {
      hnQa[iData][i]->GetAxis(0)->SetRangeUser(min_pt, max_pt);
      for(int j=0; j<nCentBins; j++)
	{
	  hnQa[iData][i]->GetAxis(3)->SetRange(centBins_low[j], centBins_high[j]);
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hnQa[iData][i]->GetAxis(4)->SetRange(k+1, k+1);
	      for(int m=0; m<2; m++)
		{
		  hnQa[iData][i]->GetAxis(5)->SetRange(m+1, m+1);
		  hTrkDisInMb[i][j][k][m] = (TH1F*)hnQa[iData][i]->Projection(i+1);
		  hTrkDisInMb[i][j][k][m]->SetName(Form("%s_%s%s_cent%s%s",dataTitle[iData],partName[m],typeName[i],cent_Title[j].Data(), gTrgSetupTitle[k+1]));
		  if(hTrkDisInMb[i][j][k][m]->Integral()>0)  
		    hTrkDisInMb[i][j][k][m]->Scale(1./hTrkDisInMb[i][j][k][m]->Integral());
		}
	      hnQa[iData][i]->GetAxis(5)->SetRange(0,-1);
	    }
	  hnQa[iData][i]->GetAxis(4)->SetRange(0,-1);
	}
      hnQa[iData][i]->GetAxis(3)->SetRange(0,-1);
      hnQa[iData][i]->GetAxis(0)->SetRange(0,-1);
    }

  TCanvas *cHadronVsMuon[nType][nCentBins/2];
  TH1F *hplot = new TH1F("hplot", "", 100, 0, 50);
  for(int i=0; i<nType; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  if(j%2==0) continue;
	  cHadronVsMuon[i][j/2] = new TCanvas(Form("cHadronVsMuon_%s_%s",typeName[i], cent_Name[j].Data()), Form("cHadronVsMuon_%s_%s",typeName[i], cent_Name[j].Data()), 1100, 700);
	  cHadronVsMuon[i][j/2]->Divide(2, 2);
	  TLegend *leg =  0x0;
	  if(i==0) leg = new TLegend(0.15, 0.6, 0.3, 0.85);
	  if(i==1) leg = new TLegend(0.5, 0.6, 0.7, 0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      cHadronVsMuon[i][j/2]->cd(k+1);
	      if(i==0)
		{
		  hplot->GetYaxis()->SetRangeUser(0,0.1);
		}
	      if(i==1)
		{
		  gPad->SetLogy();
		  hplot->GetXaxis()->SetRangeUser(0, 3.5);
		  hplot->GetYaxis()->SetRangeUser(1e-4,1);
		}
	      hplot->DrawCopy();
	      for(int m=0; m<2; m++)
		{
		  hTrkDisInMb[i][j][k][m]->SetMarkerStyle(20+2*m);
		  hTrkDisInMb[i][j][k][m]->SetMarkerColor(gColor[m]);
		  hTrkDisInMb[i][j][k][m]->SetLineColor(gColor[m]);
		  hTrkDisInMb[i][j][k][m]->SetMarkerSize(1.2);
		  hTrkDisInMb[i][j][k][m]->SetTitle("");
		  hTrkDisInMb[i][j][k][m]->DrawCopy("samesP");
		  if(k==0) leg->AddEntry(hTrkDisInMb[i][j][k][m], Form("%s",partName[m]), "P");
		}
	      TPaveText *t1 = 0x0;
	      if(iData==0) t1 = GetTitleText(Form("Dimuon: %s%% (%s%s)",cent_Name[j].Data(),run_type.Data(),gTrgSetupTitle[k+1]), 0.05);
	      if(iData==1) t1 = GetTitleText(Form("VPD-noVtx: %s%% (%s%s)",cent_Name[j].Data(),run_type.Data(),gTrgSetupTitle[k+1]), 0.05);
	      t1->Draw();
	    }
	  cHadronVsMuon[i][j/2]->cd(1);
	  leg->Draw();
	  if(savePlot) cHadronVsMuon[i][j/2]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sHadronVsMuon_%s_%s.pdf",run_type.Data(),typeName[i],cent_Title[j].Data(),dataTitle[iData]));
	}
    }

  // different data sets
  TCanvas *cData[nType][nCentBins];
  TH1F *hCutAll[nType][nData][gNTrgSetup-1];
  TH1F *hCutAcc[nType][nData][gNTrgSetup-1];
  for(int i=0; i<nType; i++)
    {
      for(int d=0; d<nData; d++)
	{
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hCutAll[i][d][k] = new TH1F(Form("h%s_All_%s%s",typeName[i], dataTitle[d], gTrgSetupTitle[k+1]), "", nCentBins, 0, nCentBins);
	      hCutAcc[i][d][k] = new TH1F(Form("h%s_Acc_%s%s",typeName[i], dataTitle[d], gTrgSetupTitle[k+1]), "", nCentBins, 0, nCentBins);
	    }
	}
      for(int j=0; j<nCentBins; j++)
	{
	  cData[i][j] = new TCanvas(Form("cData_%s_%s",typeName[i], cent_Name[j].Data()), Form("cData_%s_%s",typeName[i], cent_Name[j].Data()), 1100, 700);
	  cData[i][j]->Divide(2, 2);
	  TLegend *leg =  0x0;
	  if(i==0) leg = new TLegend(0.15, 0.6, 0.3, 0.85);
	  if(i==1) leg = new TLegend(0.5, 0.6, 0.7, 0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      cData[i][j]->cd(k+1);
	      if(i==1) gPad->SetLogy();
	      double all, all_err, acc, acc_err;
	      for(int d=0; d<nData; d++)
		{
		  TH1F *hplot = (TH1F*)hTrkDis[d][i][j][k]->Clone(Form("%s_clone",hTrkDis[d][i][j][k]->GetName()));
		  hplot->SetMarkerStyle(20+2*d);
		  hplot->SetMarkerColor(gColor[d]);
		  hplot->SetLineColor(gColor[d]);
		  hplot->SetMarkerSize(1.2);
		  hplot->SetTitle("");
		  if(i==0) 
		    {
		      hplot->GetYaxis()->SetRangeUser(0,0.12);
		      all = hTrkDis[d][i][j][k]->IntegralAndError(hTrkDis[d][i][j][k]->FindFixBin(15),
								  hTrkDis[d][i][j][k]->GetNbinsX(),
								  all_err);
		      acc = hTrkDis[d][i][j][k]->IntegralAndError(hTrkDis[d][i][j][k]->FindFixBin(20),
								  hTrkDis[d][i][j][k]->GetNbinsX(),
								  acc_err);
		    }
		  if(i==1) 
		    {
		      hplot->GetXaxis()->SetRangeUser(0,3.5);
		      hplot->GetYaxis()->SetRangeUser(5e-4, 1);
		      all = hTrkDis[d][i][j][k]->IntegralAndError(1, hTrkDis[d][i][j][k]->FindFixBin(3), all_err);
		      acc = hTrkDis[d][i][j][k]->IntegralAndError(1, hTrkDis[d][i][j][k]->FindFixBin(1), acc_err);
		    }
		  hCutAll[i][d][k]->SetBinContent(j+1, all);
		  hCutAll[i][d][k]->SetBinError(j+1, all_err);
		  hCutAcc[i][d][k]->SetBinContent(j+1, acc);
		  hCutAcc[i][d][k]->SetBinError(j+1, acc_err);

		  if(d==0) hplot->DrawCopy("P");
		  else     hplot->DrawCopy("samesP");
		  if(k==0) leg->AddEntry(hplot, legNameData[d].Data(), "P");
		}
	      TPaveText *t1 = GetTitleText(Form("%s%s: %s%%",run_type.Data(),gTrgSetupTitle[k+1],cent_Name[j].Data()), 0.055);
	      t1->Draw();
	    }
	  cData[i][j]->cd(1);
	  leg->Draw();
	  if(savePlot) cData[i][j]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sVsData_%s_%s.pdf",run_type.Data(),typeName[i],cent_Title[j].Data(),saveName));
	}
    }

  TCanvas *cEffVsCent[nData];
  TCanvas *cEffVsData[nData];
  TGraphAsymmErrors *gCutEff[nType][nCentBins-1][nData];
  for(int i=0; i<nType; i++)
    {
      cEffVsCent[i] = new TCanvas(Form("cEffVsCent_%s",typeName[i]), Form("cEffVsCent_%s",typeName[i]), 1100, 700);
      cEffVsCent[i]->Divide(2,2);
      TLegend *leg = new TLegend(0.15, 0.15, 0.3, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  cEffVsCent[i]->cd(k+1);
	  for(int d=0; d<nData; d++)
	    {
	      gCutEff[i][d][k] = new TGraphAsymmErrors(hCutAcc[i][d][k], hCutAll[i][d][k],"cl=0.683 b(1,1) mode");
	      if(d==0)
		{
		  TH1F *hplot = (TH1F*)hCutAcc[i][d][k]->Clone(Form("plot_%s",hCutAcc[i][d][k]->GetName()));
		  for(int bin=1; bin<=hplot->GetNbinsX(); bin++)
		    {
		      hplot->GetXaxis()->SetBinLabel(bin, Form("%s%%",cent_Name[bin-1].Data()));
		    }
		  hplot->GetXaxis()->SetLabelSize(0.06);
		  hplot->Reset();
		  if(i==0) hplot->GetYaxis()->SetRangeUser(0.85, 1.0);
		  if(i==1) hplot->GetYaxis()->SetRangeUser(0.7, 1.0);
		  hplot->Draw();
		}
	      gCutEff[i][d][k]->SetMarkerStyle(20+2*d);
	      gCutEff[i][d][k]->SetMarkerColor(gColor[d]);
	      gCutEff[i][d][k]->SetLineColor(gColor[d]);
	      gCutEff[i][d][k]->SetMarkerSize(1.2);
	      gCutEff[i][d][k]->SetTitle("");
	      gCutEff[i][d][k]->Draw("samesPE");
	      if(j==0) leg->AddEntry(gCutEff[i][d][k], legNameData[d].Data(), "P");
	    }
	  if(i==0)
	    {
	      TPaveText *t1 = GetTitleText(Form("eff = (NHits>=20)/(NHits>=10) (%s)",gTrgSetupTitle[k+1]),0.055);
	      t1->Draw();
	    }
	  if(i==1)
	    {
	      TPaveText *t1 = GetTitleText(Form("eff = (DCA<=3)/(DCA<=1) (%s)",gTrgSetupTitle[k+1]),0.055);
	      t1->Draw();
	    }
	}
      cEffVsCent[i]->cd(3);
      leg->Draw();
      if(savePlot) cEffVsCent[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sEffComp_%s.pdf",run_type.Data(),typeName[i],saveName));
    }
  //return;

  // in luminosity bins
  TCanvas *cCent[nData][nType];
  for(int d=0; d<2; d++)
    {
      for(int i=0; i<nType; i++)
	{
	  cCent[d][i] = new TCanvas(Form("cCent_%s_%s",dataTitle[d], typeName[i]), Form("cCent_%s_%s",dataTitle[d], typeName[i]), 1100, 700);
	  cCent[d][i]->Divide(2, 2);
	  TLegend *leg =  0x0;
	  if(i==0) leg = new TLegend(0.2, 0.45, 0.4, 0.8);
	  if(i==1) leg = new TLegend(0.6, 0.45, 0.8, 0.8);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
	  for(int j=0; j<nCentBins-1; j++)
	    {
	      if(j%2==1) continue;
	      cCent[d][i]->cd(j/2+1);
	      if(i==1) gPad->SetLogy();
	      for(int k=0; k<gNTrgSetup-1; k++)
		{
		  TH1F *hplot = (TH1F*)hTrkDis[d][i][j][k]->Clone(Form("%s_clone",hTrkDis[d][i][j][k]->GetName()));
		  hplot->SetMarkerStyle(20+k);
		  hplot->SetMarkerColor(gColor[k]);
		  hplot->SetLineColor(gColor[k]);
		  hplot->SetMarkerSize(1.2);
		  hplot->SetTitle("");
		  if(i==0) 
		    {
		      hplot->GetYaxis()->SetRangeUser(0,0.1);
		    }
		  if(i==1) 
		    {
		      hplot->GetXaxis()->SetRangeUser(0,3.5);
		      if(d==1) hplot->GetYaxis()->SetRangeUser(1e-4, 1);
		      if(d==0) hplot->GetYaxis()->SetRangeUser(1e-3, 1);
		    }
		  if(k==0) hplot->DrawCopy("P");
		  else     hplot->DrawCopy("samesP");
		  if(j==0) leg->AddEntry(hplot, Form("%s",gTrgSetupTitle[k+1]), "P");
		}
	      TPaveText *t1 = GetTitleText(Form("%s_%s (%d-%d%%, %s)",run_type.Data(),dataTitle[d],(16-centBins_high[j+1])*5,(17-centBins_low[j+1])*5,legTitle), 0.055);
	      t1->Draw();
	    }
	  cCent[d][i]->cd(1);
	  leg->Draw();
	  if(savePlot) cCent[d][i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sVsCent_%s.pdf",run_type.Data(),typeName[i],dataTitle[d]));
	}
    }
}

//================================================
void tightCutsEmb(const int savePlot = 0)
{
  const char* trigTitle = "dimuon";
  TFile *fin = 0x0;
  char *trigName = "";
  char *libVer = "";

  if(trigTitle=="vpd30")
    {
      fin = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
      trigName = "mb";
      libVer = "P17id";
    }
  else if(trigTitle=="dimuon")
    {
      fin = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");
      trigName = "di_mu";
      libVer = "P15ie";
    }

  const int nType = 3;
  const char* typeName[nType] = {"MC", "Tpc", "ElecCheck"};
  THnSparseF *hnTpc[nType];
  TH1F *hTrkPtVsCent[nType][gNTrgSetup-1];
  for(int i=0; i<nType; i++)
    {
      hnTpc[i] = (THnSparseF*)fin->Get(Form("mhMcTrkPtEff_%s_%s",typeName[i],trigName));

      hnTpc[i]->GetAxis(0)->SetRangeUser(2, 20);
      hnTpc[i]->GetAxis(1)->SetRangeUser(-0.5, 0.5);
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hnTpc[i]->GetAxis(4)->SetRange(k+1, k+1);
	  hTrkPtVsCent[i][k] = (TH1F*)hnTpc[i]->Projection(2);
	  hTrkPtVsCent[i][k]->SetName(Form("hTrkPtVsCent_%s%s",typeName[i],gTrgSetupTitle[k+1]));
	  hTrkPtVsCent[i][k]->Sumw2();
	  hTrkPtVsCent[i][k]->Rebin(2);
	}
      hnTpc[i]->GetAxis(0)->SetRange(0,-1);
      hnTpc[i]->GetAxis(1)->SetRange(0,-1);
      hnTpc[i]->GetAxis(4)->SetRange(0,-1);
    }

  TList *list = new TList;
  TString legName_trg[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      legName_trg[k] = Form("AuAu_200%s",gTrgSetupTitle[k+1]);
    }
  TH1F *hTrkEffVsCent[2][gNTrgSetup-1];
  const char *typeTitle[2] = {"standard cuts", "tighter cuts"};
  for(int i=0; i<2; i++)
    {    
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hTrkEffVsCent[i][k] = (TH1F*)hTrkPtVsCent[i+1][k]->Clone(Form("hTrkEffVsCent_%s%s",typeName[i+1],gTrgSetupTitle[k+1]));
	  hTrkEffVsCent[i][k]->Divide(hTrkPtVsCent[0][k]);
	  hTrkEffVsCent[i][k]->SetMarkerStyle(20+k);
	  hTrkEffVsCent[i][k]->SetMarkerColor(gColor[k]);
	  hTrkEffVsCent[i][k]->SetLineColor(gColor[k]);
	  for(int bin=1; bin<=hTrkEffVsCent[i][k]->GetNbinsX(); bin++)
	    {
	      hTrkEffVsCent[i][k]->GetXaxis()->SetBinLabel(bin,Form("%d-%d%%",80-bin*10,90-bin*10));
	    }
	  list->Add(hTrkEffVsCent[i][k]);
	}
      c = drawHistos(list,Form("TpcEff_vs_Cent_%s",typeName[i+1]),Form("%s: TPC tracking efficiency above 2 GeV/c (%s);;Efficiency",trigTitle,typeTitle[i]),false,2.0,3.8,true,0,1,kFALSE,kTRUE,legName_trg,true,"|#eta| < 0.5",0.5,0.7,0.2,0.45,kTRUE,0.04,0.04); 
      c->SetGrid(0,1);
      list->Clear();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TpcEffVsCent_Cut%d.pdf",run_type.Data(),trigName,i));
    }
}
