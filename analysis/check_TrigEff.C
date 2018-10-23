const int year = YEAR;
const int nCentBins = nCentBins_pt; 
const char** cent_Name = cent_Name_pt;
const char** cent_Title = cent_Title_pt;
const int* centBins_low = centBins_low_pt;
const int* centBins_high = centBins_high_pt;
const TString legName[5] = {"All","prod","prod_low","prod_mid","prod_high"};
//================================================
void check_TrigEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //makeTrigEff();
  //anaTrigEff();
  //genMuonTrigEffLS();
  //compareTrigEffLS();

  respAndMthEff();
  //mtdTrigTime();
  //dmMthEff();
}

//================================================
void dmMthEff(const int savePlot = 0)
{
  TFile* fstudy = TFile::Open("output/Run14_AuAu200.Study.root", "read");


  // Get # of events
  double nDmEvt[gNTrgSetup][24];
  THnSparseF *hnCent = (THnSparseF*)fstudy->Get("mhMtdQaCentAll_di_mu");
  hnCent->GetAxis(2)->SetRange(2,2); // select good vertex
  hnCent->GetAxis(1)->SetRange(1,16); // 0-80%
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hnCent->GetAxis(0)->SetRange(k,k);
      TH1F* h1tmp = (TH1F*)hnCent->Projection(3);
      h1tmp->SetName(Form("hnCent_%d",k));
      for(int bin=1; bin<=h1tmp->GetNbinsX(); bin++)
	{
	  nDmEvt[k][bin-1] = h1tmp->GetBinContent(bin);
	}
      hnCent->GetAxis(0)->SetRange(0,-1);
    }
  hnCent->GetAxis(1)->SetRange(0,-1);
  hnCent->GetAxis(2)->SetRange(0,-1);


  // Get # of matched tracks
  //mhnMtdMthEff_mb;1	Track p_{T} vs isMth vs. TrigSetup vs ZDC vs v_{z} vs #varphi vs. isTrigMth
  THnSparseF *hnTrk = (THnSparseF*)fstudy->Get("mhnMtdMthEff_di_mu");
  hnTrk->GetAxis(0)->SetRangeUser(2,10); // pt
  hnTrk->GetAxis(1)->SetRange(2,2); // matched tracks
  TH1F *hNMthTrkVsZdc[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hnTrk->GetAxis(2)->SetRange(k,k);
      hNMthTrkVsZdc[k] = (TH1F*)hnTrk->Projection(3);
      hNMthTrkVsZdc[k]->SetName(Form("hNMthTrkVsZdc%s",gTrgSetupTitle[k]));
      hNMthTrkVsZdc[k]->Sumw2();
      for(int bin=1; bin<=hNMthTrkVsZdc[k]->GetNbinsX(); bin++)
	{
	  if(nDmEvt[k][bin-1]<=0) continue;
	  hNMthTrkVsZdc[k]->SetBinContent(bin, hNMthTrkVsZdc[k]->GetBinContent(bin)/nDmEvt[k][bin-1]);
	  hNMthTrkVsZdc[k]->SetBinError(bin, hNMthTrkVsZdc[k]->GetBinError(bin)/nDmEvt[k][bin-1]);
	}
      hnTrk->GetAxis(5)->SetRange(0,-1);
    }

  TCanvas *cMthTrk = new TCanvas("cMthTrk", "cMthTrk", 800, 600);
  leg = new TLegend(0.5,0.6,0.7,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(int k=2; k<gNTrgSetup; k++)
    {
      hNMthTrkVsZdc[k]->SetTitle("");
      hNMthTrkVsZdc[k]->SetMarkerStyle(20+k);
      hNMthTrkVsZdc[k]->SetMarkerColor(gColor[k]);
      hNMthTrkVsZdc[k]->GetYaxis()->SetRangeUser(0.13, 0.26);
      if(k==0) hNMthTrkVsZdc[k]->Draw("PE");
      else     hNMthTrkVsZdc[k]->Draw("samesPE");
      leg->AddEntry(hNMthTrkVsZdc[k], gTrgSetupTitle[k], "PL");
    }
  TPaveText *t1 = GetTitleText("di_muon: # of matched tracks");
  t1->Draw();
  leg->Draw();
  if(savePlot)
    cMthTrk->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_DM_MthTrkVsZdc.pdf",run_type.Data()));
}

//================================================
void respAndMthEff(const int savePlot = 1)
{
  const int iMB = 1;
  const char* trig_name[2] = {"di_mu", "mb"};
  TFile *fstudy = NULL;
  if(iMB==1) fstudy = TFile::Open("output/Run14_AuAu200.MB.P15ic.Study.root", "read");
  if(iMB==0) fstudy = TFile::Open("output/Run14_AuAu200.Study.root", "read");

  const int nHistos = 4;
  const char* histoName[nHistos] = {"TrigElec", "Resp", "Mth", "TrigMth"};
  TH1F *hTrkPtAll[nHistos][gNTrgSetup];
  TH1F *hTrkPtAcc[nHistos][gNTrgSetup];
  TH1F *hTrkPtEff[nHistos][gNTrgSetup];

  const int nPtBins = 8;
  const double xPtBins[nPtBins+1] = {0, 1, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0};

  // trigger electronics efficiency
  THnSparseF *hnTrigElecEff = (THnSparseF*)fstudy->Get(Form("mhMtdTrigElecEff_%s",trig_name[iMB]));
  hnTrigElecEff->GetAxis(5)->SetRange(1,16); // 0-80% centrality
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hnTrigElecEff->GetAxis(6)->SetRange(k,k);
      TH1F *h1tmp = (TH1F*)hnTrigElecEff->Projection(0);
      h1tmp->SetName(Form("hTrkPtAll_%s%s_tmp",histoName[0],gTrgSetupTitle[k]));
      h1tmp->Sumw2();
      hTrkPtAll[0][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAll_%s%s",histoName[0],gTrgSetupTitle[k]), xPtBins);

      hnTrigElecEff->GetAxis(1)->SetRange(2,2);
      hnTrigElecEff->GetAxis(2)->SetRange(2,2);
      h1tmp = (TH1F*)hnTrigElecEff->Projection(0);
      h1tmp->SetName(Form("hTrkPtAcc_%s%s_tmp",histoName[0],gTrgSetupTitle[k]));
      h1tmp->Sumw2();
      hTrkPtAcc[0][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAcc_%s%s",histoName[0],gTrgSetupTitle[k]), xPtBins);

      hTrkPtEff[0][k] = (TH1F*)hTrkPtAcc[0][k]->Clone(Form("hTrkPtEff_%s%s",histoName[0],gTrgSetupTitle[k]));
      hTrkPtEff[0][k]->Divide(hTrkPtAll[0][k]);

      hnTrigElecEff->GetAxis(1)->SetRange(0,-1);
      hnTrigElecEff->GetAxis(2)->SetRange(0,-1);
      hnTrigElecEff->GetAxis(6)->SetRange(0,-1);
    }

  // response & matching efficiency
  double nMbEvt[gNTrgSetup];
  TH2F *hMbEvt = (TH2F*)fstudy->Get(Form("mhMtdQaCent_%s",trig_name[iMB]));
  nMbEvt[0] = hMbEvt->GetEntries();
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMbEvt->ProjectionY("",k,k);
      nMbEvt[k] = htmp->Integral(1,16);
    }

  const int gVzCut = 30;
  const int gNHitsFitCut = 15;
  THnSparseF *hnMthEff = (THnSparseF*)fstudy->Get(Form("mhnMtdMthEff_%s",trig_name[iMB]));
  hnMthEff->GetAxis(4)->SetRangeUser(-1*gVzCut, gVzCut);
  hnMthEff->GetAxis(7)->SetRangeUser(gNHitsFitCut, 45);
  THnSparseF *hnRespEff = (THnSparseF*)fstudy->Get(Form("mhnMtdRespEff_%s",trig_name[iMB]));
  hnRespEff->GetAxis(4)->SetRangeUser(-1*gVzCut, gVzCut);
  hnRespEff->GetAxis(7)->SetRangeUser(gNHitsFitCut, 45);
  THnSparseF *hntmp = NULL;
  const char* histoName[nHistos] = {"TrigElec", "Resp", "Mth", "TrigMth"};
  for(int i=0; i<3; i++)
    {
      if(i==0) hntmp = hnRespEff;
      else     hntmp = hnMthEff;
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hntmp->GetAxis(2)->SetRange(k,k);
	  TH1F *h1tmp = (TH1F*)hntmp->Projection(0);
	  h1tmp->SetName(Form("hTrkPtAll_%s%s_tmp",histoName[i+1],gTrgSetupTitle[k]));
	  h1tmp->Sumw2();
	  hTrkPtAll[i+1][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAll_%s%s",histoName[i+1],gTrgSetupTitle[k]), xPtBins);
	  printf("[i] %s has %1.0f entries\n",hTrkPtAll[i+1][k]->GetName(),hTrkPtAll[i+1][k]->GetEntries());

	  if(i<2) hntmp->GetAxis(1)->SetRange(2,2);
	  else hntmp->GetAxis(6)->SetRange(2,2);
	  h1tmp = (TH1F*)hntmp->Projection(0);
	  h1tmp->SetName(Form("hTrkPtAcc_%s%s_tmp",histoName[i+1],gTrgSetupTitle[k]));
	  h1tmp->Sumw2();
	  hTrkPtAcc[i+1][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAcc_%s%s",histoName[i+1],gTrgSetupTitle[k]), xPtBins);

	  hTrkPtEff[i+1][k] = (TH1F*)hTrkPtAcc[i+1][k]->Clone(Form("hTrkPtEff_%s%s",histoName[i+1],gTrgSetupTitle[k]));
	  hTrkPtEff[i+1][k]->Divide(hTrkPtAll[i+1][k]);

	  hntmp->GetAxis(1)->SetRange(0,-1);
	  hntmp->GetAxis(2)->SetRange(0,-1);
	  if(i==2) hntmp->GetAxis(6)->SetRange(0,-1);
	}
    }

  // show the efficiency
  const char* histoTitle[nHistos] = {"Trigger electronics efficiency (MB, muon candidates)", 
				     "MTD response efficiency (MB, charged tracks)",
				     "MTD matching efficiency (MB, tof electronics)",
				     "MTD matching efficiency (MB, trigger electronics)"};
  TCanvas *cEff[nHistos];
  TF1 *funcEff[nHistos][gNTrgSetup];
  TF1 *funcEffRatio[nHistos][gNTrgSetup];
  for(int i=0; i<nHistos; i++)
    {
      cEff[i] = new TCanvas(Form("cEff_%s",histoName[i]), Form("cEff_%s",histoName[i]), 1100, 500);
      cEff[i]->Divide(2,1);

      cEff[i]->cd(1);
      TLegend *leg = new TLegend(0.15,0.65,0.35,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("|v_{z}| < %d",gVzCut));
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k==1) continue;
	  hTrkPtEff[i][k]->SetTitle(";p_{T} (GeV/c);");
	  hTrkPtEff[i][k]->SetMarkerStyle(20+k);
	  hTrkPtEff[i][k]->SetMarkerColor(gColor[k]);
	  hTrkPtEff[i][k]->SetLineColor(gColor[k]);
	  if(i==0) hTrkPtEff[i][k]->GetYaxis()->SetRangeUser(0.8,1.2);
	  if(i==1) hTrkPtEff[i][k]->GetYaxis()->SetRangeUser(0,0.15);
	  if(i==2) hTrkPtEff[i][k]->GetYaxis()->SetRangeUser(0,0.1);
	  if(i==3) hTrkPtEff[i][k]->GetYaxis()->SetRangeUser(0,0.2);

	  if(k==0) hTrkPtEff[i][k]->Draw("PE");
	  else     hTrkPtEff[i][k]->Draw("samesPE");
	  if(k>0) leg->AddEntry(hTrkPtEff[i][k], gTrgSetupTitle[k], "PL");
	}
      TPaveText *t1 = GetTitleText(Form("%s",histoTitle[i]));
      t1->Draw();
      leg->Draw();

      cEff[i]->cd(2);
      leg = new TLegend(0.15,0.15,0.8,0.25);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetNColumns(2);
      for(int k=1; k<gNTrgSetup; k++)
	{
	  if(k==1) continue;
	  TH1F *hRatio = (TH1F*)hTrkPtEff[i][k]->Clone(Form("%s_ratio",hTrkPtEff[i][k]->GetName()));
	  hRatio->Divide(hTrkPtEff[i][0]);
	  funcEffRatio[i][k] = new TF1(Form("func%sEffRatio%s",histoName[i],gTrgSetupTitle[k]), "pol0", 1, 10);
	  hRatio->Fit(funcEffRatio[i][k],"0RQ");
	  hRatio->SetTitle(";p_{T} (GeV/c);Ratio to inclusive");
	  hRatio->GetYaxis()->SetRangeUser(0.8,1.2);
	  if(k==1) hRatio->Draw("PE");
	  else     hRatio->Draw("samesPE");
	  leg->AddEntry(hRatio, Form("%4.3f +/- %4.3f",funcEffRatio[i][k]->GetParameter(0),funcEffRatio[i][k]->GetParError(0)), "PL");
	  funcEffRatio[i][k]->SetLineColor(gColor[k]);
	  funcEffRatio[i][k]->SetLineStyle(2);
	  funcEffRatio[i][k]->Draw("sames");
	}
      t1 = GetTitleText(Form("Ratio to inclusive sample"));
      t1->Draw();
      leg->Draw();
      if(savePlot)
	cEff[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_CompMtd%sEff.pdf",run_type.Data(),histoName[i]));
    }
  
  // response efficiency vs. ZDC rate, vz, phi
  const int nDepen = 4;
  const char* depenName[4] = {"Zdc", "Vz", "Phi", "Cent"};
  const char* depenTitle[4] = {"Zdc", "v_{z}", "#varphi", "Cent"};
  const char* depenUnit[4] = {" (kHz)", " (cm)", "", ""};
  TH1F * hTrkPtVsZdcProj[3][gNTrgSetup][4];
  TH1F * hTrkPtVsZdcMth[3][gNTrgSetup][4];
  TH1F * hTrkPtVsZdcEff[3][gNTrgSetup][4];
  for(int i=0; i<3; i++)
    {
      if(i==0) hntmp = hnRespEff;
      else     hntmp = hnMthEff;
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hntmp->GetAxis(2)->SetRange(k,k);
	  for(int j=0; j<nDepen; j++)
	    {
	      //if(j>0) hntmp->GetAxis(3)->SetRangeUser(40,60);
	      //if(j!=1) hntmp->GetAxis(4)->SetRangeUser(-1*gVzCut, gVzCut);
	      if(j<3) hTrkPtVsZdcProj[i][k][j] = (TH1F*)hntmp->Projection(3+j);
	      else    hTrkPtVsZdcProj[i][k][j] = (TH1F*)hntmp->Projection(8);
	      hTrkPtVsZdcProj[i][k][j]->SetName(Form("hTrkPtVs%sProj_%s%s",depenName[j],histoName[i+1],gTrgSetupTitle[k]));

	      if(i<2) hntmp->GetAxis(1)->SetRange(2,2);
	      else    hntmp->GetAxis(6)->SetRange(2,2);
	      if(j<3) hTrkPtVsZdcMth[i][k][j] = (TH1F*)hntmp->Projection(3+j);
	      else    hTrkPtVsZdcMth[i][k][j] = (TH1F*)hntmp->Projection(8);
	      hTrkPtVsZdcMth[i][k][j]->SetName(Form("hTrkPtVs%sMth_%s%s",depenName[j],histoName[i+1],gTrgSetupTitle[k]));
	      if(i<2) hntmp->GetAxis(1)->SetRange(0,-1);
	      else    hntmp->GetAxis(6)->SetRange(0,-1);

	      hTrkPtVsZdcEff[i][k][j] = (TH1F*)hTrkPtVsZdcMth[i][k][j]->Clone(Form("hTrkPtVs%sEff_%s%s",depenName[j],histoName[i+1],gTrgSetupTitle[k]));
	      hTrkPtVsZdcEff[i][k][j]->Sumw2();
	      hTrkPtVsZdcEff[i][k][j]->Divide(hTrkPtVsZdcProj[i][k][j]);
	      hTrkPtVsZdcEff[i][k][j]->SetMarkerStyle(20+k);
	      hTrkPtVsZdcEff[i][k][j]->SetMarkerColor(gColor[k]);
	      hTrkPtVsZdcEff[i][k][j]->SetLineColor(gColor[k]);
	      
	      hntmp->GetAxis(3)->SetRange(0,-1);
	      hntmp->GetAxis(4)->SetRange(0,-1);
	    }
	  hntmp->GetAxis(2)->SetRange(0,-1);
	}
    }

  TCanvas *cEffDep[3];
  for(int i=0; i<3; i++)
    {
      cEffDep[i] = new TCanvas(Form("c%sEffDep",histoName[i+1]), Form("c%sEffDep",histoName[i+1]), 1100, 700);
      cEffDep[i]->Divide(2,2);
      for(int j=0; j<nDepen; j++)
	{
	  cEffDep[i]->cd(j+1);
	  leg = new TLegend(0.15,0.15,0.35,0.4);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  if(j==0) leg->SetHeader(Form("|v_{z}| < %d",gVzCut));
	  //if(j==1) leg->SetHeader("40 < Zdc < 60 kHz");
	  //if(j==2) leg->SetHeader(Form("|v_{z}| < %d, 40 < Zdc < 60 kHz",gVzCut));
	  for(int k=2; k<gNTrgSetup; k++)
	    {
	      hTrkPtVsZdcEff[i][k][j]->SetTitle(Form(";%s%s;Eff.",depenTitle[j],depenUnit[j]));
	      if(i==0) hTrkPtVsZdcEff[i][k][j]->GetYaxis()->SetRangeUser(0,0.035);
	      if(i==1) hTrkPtVsZdcEff[i][k][j]->GetYaxis()->SetRangeUser(0,0.03);
	      if(i==2) hTrkPtVsZdcEff[i][k][j]->GetYaxis()->SetRangeUser(0.02,0.09);

	      if(k==2) hTrkPtVsZdcEff[i][k][j]->Draw();
	      else     hTrkPtVsZdcEff[i][k][j]->Draw("sames");
	      leg->AddEntry(hTrkPtVsZdcEff[i][k][j], gTrgSetupTitle[k], "PL");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s",histoTitle[i+1]),0.05);
	  t1->Draw();
	  leg->Draw();
	}
      if(savePlot)
	cEffDep[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_Mtd%sEffDep.pdf",run_type.Data(),histoName[i+1]));
    }
  return;

  if(1)
    {
      TFile *fstudyP15ic = TFile::Open("output/Run14_AuAu200.MB.Study.prod_mid.P15ic.root", "read");
      THnSparseF *hnMthEffP15ic = (THnSparseF*)fstudyP15ic->Get(Form("mhnMtdMthEff_%s",trig_name[1]));
      THnSparseF *hnRespEffP15ic = (THnSparseF*)fstudyP15ic->Get(Form("mhnMtdRespEff_%s",trig_name[1]));
      hnMthEffP15ic->GetAxis(4)->SetRangeUser(-1*gVzCut, gVzCut);
      hnRespEffP15ic->GetAxis(4)->SetRangeUser(-1*gVzCut, gVzCut);
      hnMthEffP15ic->GetAxis(7)->SetRangeUser(gNHitsFitCut, 45);
      hnRespEffP15ic->GetAxis(7)->SetRangeUser(gNHitsFitCut, 45);
      TH1F *hTrkPtAllP15ic[nHistos][2];
      TH1F *hTrkPtAccP15ic[nHistos][2];
      TH1F *hTrkPtEffP15ic[nHistos][2];
      for(int i=0; i<3; i++)
	{
	  if(i==0) hntmp = hnRespEffP15ic;
	  else     hntmp = hnMthEffP15ic;
	  for(int k=0; k<2; k++)
	    {
	      hntmp->GetAxis(2)->SetRange(k+2, k+2);
	      TH1F *h1tmp = (TH1F*)hntmp->Projection(0);
	      h1tmp->SetName(Form("hTrkPtAll_%s%s_P15ic_tmp",histoName[i+1],gTrgSetupTitle[k+2]));
	      h1tmp->Sumw2();
	      hTrkPtAllP15ic[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAll_%s%s_P15ic",histoName[i+1],gTrgSetupTitle[k+2]), xPtBins);
	      printf("[i] %s has %1.0f entries\n",hTrkPtAllP15ic[i][k]->GetName(),hTrkPtAllP15ic[i][k]->GetEntries());

	      if(i<2) hntmp->GetAxis(1)->SetRange(2,2);
	      else hntmp->GetAxis(6)->SetRange(2,2);
	      h1tmp = (TH1F*)hntmp->Projection(0);
	      h1tmp->SetName(Form("hTrkPtAcc_%s%s_P15ic_tmp",histoName[i+1],gTrgSetupTitle[k+2]));
	      h1tmp->Sumw2();
	      hTrkPtAccP15ic[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAcc_%s%s_P15ic",histoName[i+1],gTrgSetupTitle[k+2]), xPtBins);

	      hTrkPtEffP15ic[i][k] = (TH1F*)hTrkPtAccP15ic[i][k]->Clone(Form("hTrkPtEff_%s%s_P15ic",histoName[i+1],gTrgSetupTitle[k+2]));
	      hTrkPtEffP15ic[i][k]->Divide(hTrkPtAllP15ic[i][k]);
	  
	      hntmp->GetAxis(1)->SetRange(0,-1);
	      hntmp->GetAxis(2)->SetRange(0,-1);
	      if(i==2) hntmp->GetAxis(6)->SetRange(0,-1);
	    }

	  cEff[i+1]->cd(1);
	  TLegend *leg = new TLegend(0.15,0.55,0.35,0.65);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  for(int k=0; k<2; k++)
	    {
	      hTrkPtEffP15ic[i][k]->SetMarkerSize(1.5);
	      hTrkPtEffP15ic[i][k]->SetMarkerStyle(29+k);
	      hTrkPtEffP15ic[i][k]->SetMarkerColor(2+k*7);
	      hTrkPtEffP15ic[i][k]->SetLineColor(2+k*7);
	      hTrkPtEffP15ic[i][k]->Draw("samesPE");
	      leg->AddEntry(hTrkPtEffP15ic[i][k], Form("%s (P15ic)",gTrgSetupTitle[k+2]), "PL");
	    }
	  leg->Draw();

	  cEff[i+1]->cd(2);
	  leg = new TLegend(0.475,0.11,0.8,0.20);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  for(int k=0; k<2; k++)
	    {
	      TH1F *hRatio = (TH1F*)hTrkPtEffP15ic[i][k]->Clone(Form("%s_ratio",hTrkPtEffP15ic[i][k]->GetName()));
	      hRatio->Divide(hTrkPtEff[i+1][0]);
	      TF1 *funcRatio = new TF1(Form("fit_%s",hRatio->GetName()), "pol0", 1, 10);
	      hRatio->Fit(funcRatio,"0RQ");
	      hRatio->Draw("samesPE");
	      funcRatio->SetLineColor(hRatio->GetMarkerColor());
	      funcRatio->SetLineStyle(2);
	      funcRatio->Draw("sames");
	      leg->AddEntry(hRatio, Form("%4.3f +/- %4.3f",funcRatio->GetParameter(0),funcRatio->GetParError(0)), "PL");
	      leg->Draw();
	    }
	  if(savePlot)
	    cEff[i+1]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_CompMtd%sEff_P15ic.pdf",run_type.Data(),histoName[i+1]));
	}
      hnMthEffP15ic->GetAxis(4)->SetRange(0,-1);
      hnRespEffP15ic->GetAxis(4)->SetRange(0,-1);

      TH1F * hTrkPtVsZdcProjP15ic[3][2][nDepen];
      TH1F * hTrkPtVsZdcMthP15ic[3][2][nDepen];
      TH1F * hTrkPtVsZdcEffP15ic[3][2][nDepen];
      const char* depenTitle[3] = {"Zdc", "v_{z}", "#varphi"};
      for(int i=0; i<3; i++)
	{
	  if(i==0) hntmp = hnRespEffP15ic;
	  else     hntmp = hnMthEffP15ic;
	  for(int k=0; k<2; k++)
	    {
	      hntmp->GetAxis(2)->SetRange(k+2,k+2); // prod_mid
	      for(int j=0; j<nDepen; j++)
		{
		  if(j>0) hntmp->GetAxis(3)->SetRangeUser(40,60);
		  if(j!=1) hntmp->GetAxis(4)->SetRangeUser(-1*gVzCut, gVzCut);
	      
		  hTrkPtVsZdcProjP15ic[i][k][j] = (TH1F*)hntmp->Projection(3+j);
		  hTrkPtVsZdcProjP15ic[i][k][j]->SetName(Form("hTrkPtVs%sProj_%s%s_P15ic",depenName[j],histoName[i+1],gTrgSetupTitle[k+2]));

		  if(i<2) hntmp->GetAxis(1)->SetRange(2,2);
		  else    hntmp->GetAxis(6)->SetRange(2,2);
		  hTrkPtVsZdcMthP15ic[i][k][j] = (TH1F*)hntmp->Projection(3+j);
		  hTrkPtVsZdcMthP15ic[i][k][j]->SetName(Form("hTrkPtVs%sMth_%s%s_P15ic",depenName[j],histoName[i+1],gTrgSetupTitle[k+2]));
		  if(i<2) hntmp->GetAxis(1)->SetRange(0,-1);
		  else    hntmp->GetAxis(6)->SetRange(0,-1);

		  hntmp->GetAxis(3)->SetRange(0,-1);
		  hntmp->GetAxis(4)->SetRange(0,-1);

		  hTrkPtVsZdcEffP15ic[i][k][j] = (TH1F*)hTrkPtVsZdcMthP15ic[i][k][j]->Clone(Form("hTrkPtVs%sEff_%s%s_P15ic",depenName[j],histoName[i+1],gTrgSetupTitle[k+2]));
		  hTrkPtVsZdcEffP15ic[i][k][j]->Sumw2();
		  hTrkPtVsZdcEffP15ic[i][k][j]->Divide(hTrkPtVsZdcProjP15ic[i][k][j]);
		  hTrkPtVsZdcEffP15ic[i][k][j]->SetMarkerStyle(29+k);
		  hTrkPtVsZdcEffP15ic[i][k][j]->SetMarkerSize(1.5);
		  hTrkPtVsZdcEffP15ic[i][k][j]->SetMarkerColor(2+7*k);
		  hTrkPtVsZdcEffP15ic[i][k][j]->SetLineColor(2+7*k);
		  
		  cEffDep[i]->cd(j+1);
		  hTrkPtVsZdcEffP15ic[i][k][j]->Draw("sames");
		}
	      hntmp->GetAxis(2)->SetRange(0,-1);
	    }
	  if(savePlot)
	    cEffDep[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_Mtd%sEffDep_P15ic.pdf",run_type.Data(),histoName[i+1]));
	}
 
    }
}

//================================================
void compareTrigEffLS(const int savePlot = 1)
{
  TFile *fdata = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root", "read");
  THnSparseF *hnInvMass = (THnSparseF*)fdata->Get("mhInvMassEff_di_mu");
  TH1F *hPairPtAll[gNTrgSetup-1];
  TH1F *hPairPtAcc[gNTrgSetup-1];
  TH1F *hPairPtEff[gNTrgSetup-1];

  hnInvMass->GetAxis(2)->SetRange(1,16);
  hnInvMass->GetAxis(4)->SetRange(1,2);
  hnInvMass->GetAxis(0)->SetRangeUser(2.6, 3.6);
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      hnInvMass->GetAxis(3)->SetRange(k+1, k+1);
      hPairPtAll[k] = (TH1F*)hnInvMass->Projection(1);
      hPairPtAll[k]->SetName(Form("hPairPtAll%s",gTrgSetupTitle[k+1]));
      //hPairPtAll[k]->Rebin(2);

      hnInvMass->GetAxis(5)->SetRange(2, 2);
      hPairPtAcc[k] = (TH1F*)hnInvMass->Projection(1);
      hPairPtAcc[k]->SetName(Form("hPairPtAcc%s",gTrgSetupTitle[k+1]));
      //hPairPtAcc[k]->Rebin(2);
      hnInvMass->GetAxis(5)->SetRange(0, -1);

      hnInvMass->GetAxis(3)->SetRange(0, -1);
      printf("[i] LS pair: %1.0f/%1.0f\n",hPairPtAll[k]->Integral(0,-1),hPairPtAcc[k]->Integral(0,-1));
    }

  TList *list = new TList;
  TString legName1[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      hPairPtEff[k] = (TH1F*)hPairPtAcc[k]->Clone(Form("hPairPtEff%s",gTrgSetupTitle[k+1]));
      hPairPtEff[k]->Divide(hPairPtAll[k]);
      list->Add(hPairPtEff[k]);
      legName1[k] = Form("%s",gTrgSetupTitle[k+1]);
    }
  TCanvas *c = drawHistos(list,"TrigEffVsPt_LS","",false,0,10,false,0.8,1.05,kFALSE,true,legName1,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  
  // compare LS efficiency using different methods
  TCanvas *cTrigEff = new TCanvas("cTrigEffLS", "cTrigEffLS", 1100, 700);
  cTrigEff->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.2,0.35,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  TFile *fEff = TFile::Open("Rootfiles/Run14_AuAu200.StudyTrigEff.root", "read");
  TH1F *hLSeff[3][gNTrgSetup-1];
  TH1F *hLSeffWeight[2];
  const double weights[gNTrgSetup-1] = {0.0744, 0.1194, 0.2655, 0.5408};
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      cTrigEff->cd(k+1);
      for(int i=0; i<3; i++)
	{
	  if(i<2) hLSeff[i][k] = (TH1F*)fEff->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_M%d_TrigStudy",gTrgSetupTitle[k+1],i));
	  else    hLSeff[i][k] = (TH1F*)hPairPtEff[k]->Clone(Form("%s_clone",hPairPtEff[k]->GetName()));
	  hLSeff[i][k]->SetTitle(";p_{T} (GeV/c);Eff.");
	  hLSeff[i][k]->SetMarkerStyle(20+i);
	  hLSeff[i][k]->SetMarkerColor(gColor[i]);
	  hLSeff[i][k]->SetLineColor(gColor[i]);
	  hLSeff[i][k]->GetYaxis()->SetRangeUser(0,1);
	  if(i==0) hLSeff[i][k]->Draw();
	  else     hLSeff[i][k]->Draw("sames");
	  if(k==0)
	    {
	      if(i==0) leg->AddEntry(hLSeff[i][k], "Run15 p+p dimuon", "P");
	      if(i==1) leg->AddEntry(hLSeff[i][k], "Run14 Au+Au MB (single track)", "P");
	      if(i==2) leg->AddEntry(hLSeff[i][k], "Run14 Au+Au MB (LS+US, 2.6 < M < 3.6)", "P");
	    }
	  if(i<2)
	    {
	      if(k==0) 
		{
		  hLSeffWeight[i] = (TH1F*)hLSeff[i][k]->Clone(Form("Run14_AuAu200_JpsiTrigEffVsPt_M%d",i));
		  hLSeffWeight[i]->Reset("AC");
		}
	      hLSeffWeight[i]->Add(hLSeff[i][k], weights[k]);
	    }
	}
      TPaveText *t1 = GetTitleText(Form("MTD trigger efficiency for J/psi (%s)",gTrgSetupTitle[k+1]),0.055);
      t1->Draw();
    }
  cTrigEff->cd(1);
  leg->Draw();
  if(savePlot) cTrigEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_Comp_ppVsMBAuAu.pdf",run_type.Data()));

  hLSeffWeight[1]->Divide(hLSeffWeight[0]);
  c = draw1D(hLSeffWeight[1],"Ratio of MTD trigger efficiency for LS pairs: AuAuMB/pp;p_{T} (GeV/c);Ratio");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_Sys_ppVsMBAuAu.pdf",run_type.Data()));
}


//================================================
void genMuonTrigEffLS(const int mode = 2, const int savePlot = 1, const int saveHisto = 1)
{
  if(mode<2)
    {
      const int nPtBins = 6;
      const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 
    }
  else
    {
      const int nPtBins = 10;
      const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0}; 
    }
  const int nLumi = 4;
  const char *lumi_name[nLumi] = {"prod","prod_low","prod_mid","prod_high"};
  const double min_TacDiffCut[nLumi] = {785+1,785+1,788+1,788+1};
  const double max_TacDiffCut[nLumi] = {980, 980, 980, 980};
  TGraphAsymmErrors* gTacDiffEff[nLumi];
  double pt_arr[2][nPtBins], all_arr[2][nPtBins], all_err_arr[2][nPtBins], acc_arr[2][nPtBins], acc_err_arr[2][nPtBins];

  if(mode==0)
    {
      // use muon candidates from Run15 data
      TH2F *hppMuonTacDiffVsTrigUnit[nPtBins];
      TH1F *hppMuonTacDiff[nPtBins];
      TH2F *hppTacDiffVsTrigUnit[nPtBins];
      TH1F *hppTacDiff[nPtBins];
      TF1  *funcppTacDiff[nPtBins];
      TFile *fpp = TFile::Open("Rootfiles/Run15_pp200.MtdTrigEff.root", "read");
      TCanvas *c = new TCanvas("fit_TacDiff_LS", "fit_TacDiff_LS", 1100, 700);
      c->Divide(3,2);
      for(int i=0; i<nPtBins; i++)
	{
	  hppMuonTacDiffVsTrigUnit[i] = (TH2F*)fpp->Get(Form("Run15_pp200_hTacDiffVsTrigUnit_Muon_PtBin%d",i+1));
	  hppMuonTacDiff[i] = (TH1F*)hppMuonTacDiffVsTrigUnit[i]->ProjectionY(Form("Run15_pp200_hTacDiff_Muon_PtBin%d",i+1));
	  hppMuonTacDiff[i]->Rebin(2);
	  hppMuonTacDiff[i]->SetMarkerStyle(24);
	  hppMuonTacDiff[i]->SetMarkerColor(2);
	  hppMuonTacDiff[i]->SetLineColor(2);
	  hppMuonTacDiff[i]->Scale(1./hppMuonTacDiff[i]->Integral()*0.5);

	  hppTacDiffVsTrigUnit[i] = (TH2F*)fpp->Get(Form("Run15_pp200_hTacDiffVsTrigUnit_LS_PtBin%d",i+1));
	  hppTacDiff[i] = (TH1F*)hppTacDiffVsTrigUnit[i]->ProjectionY(Form("Run15_pp200_hTacDiff_LS_PtBin%d",i+1));
	  hppTacDiff[i]->Rebin(2);
	  hppTacDiff[i]->Scale(1./hppTacDiff[i]->Integral()*0.5);
	  if(i<nPtBins-1)
	  //if(i<0)
	    {
	      funcppTacDiff[i] = new TF1(Form("fit_%s",hppTacDiff[i]->GetName()), CrystalBall, 880, 960, 5);
	      funcppTacDiff[i]->SetParameters(1,925,8,1,0.2);
	    }
	  else
	    {
	      funcppTacDiff[i] = new TF1(Form("fit_%s",hppTacDiff[i]->GetName()), "gaus", 910, 960);
	    }

	  hppTacDiff[i]->Fit(funcppTacDiff[i], "IR0");
	  c->cd(i+1);
	  hppTacDiff[i]->SetMarkerStyle(20);
	  hppTacDiff[i]->GetXaxis()->SetRangeUser(860, 980);
	  hppTacDiff[i]->SetMaximum(1.4*hppTacDiff[i]->GetMaximum());
	  hppTacDiff[i]->Draw();
	  hppMuonTacDiff[i]->Draw("sames");
	  funcppTacDiff[i]->SetLineStyle(2);
	  funcppTacDiff[i]->SetLineColor(4);
	  funcppTacDiff[i]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("Run15_pp200: %1.1f < p_{T} < %1.1f GeV/c",xPtBins[i],xPtBins[i+1]),0.05);
	  t1->Draw();
	}

      // Get the mean and sigma for TacDiff distributions from Au+Au
      TFile *fin = 0x0;
      if(saveHisto) fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyTrigEff.root", "update");
      else          fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyTrigEff.root", "read");
      TH1F *hTacDiffPtMean[gNTrgSetup];
      TH1F *hTacDiffPtSigma[gNTrgSetup];  
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffPtMean[j]  = (TH1F*)fin->Get(Form("hTacDiffMeanVsPt%s",gTrgSetupName[j]));
	  hTacDiffPtSigma[j] = (TH1F*)fin->Get(Form("hTacDiffSigmaVsPt%s",gTrgSetupName[j]));
	}

      double min_TacDiffCut_new[nLumi][nPtBins];
      for(int l=0; l<nLumi; l++)
	{
	  printf("+++ %s +++\n",lumi_name[l]);
	  // calculate new cut values
	  for(int bin=1; bin<=nPtBins; bin++)
	    {
	      double pt = hTacDiffPtMean[0]->GetBinCenter(bin);
	      double sigma_auau = hTacDiffPtSigma[l+1]->GetBinContent(bin);
	      double mean_auau = hTacDiffPtMean[l+1]->GetBinContent(bin);
	      double sigma_pp = funcppTacDiff[bin-1]->GetParameter(2);
	      double mean_pp = funcppTacDiff[bin-1]->GetParameter(1);
	      min_TacDiffCut_new[l][bin-1] = mean_pp - (mean_auau - min_TacDiffCut[l]) * sigma_pp/sigma_auau;
	      if(bin==nPtBins) min_TacDiffCut_new[l][bin-1] = min_TacDiffCut_new[l][bin-2];
	      if(l==1 || l==3) printf("[i] %s: auau (%4.2f, %4.2f), pp (%4.2f, %4.2f), old_cut = %4.2f, new_cut = %4.2f\n",gTrgSetupTitle[l+1],mean_auau,sigma_auau,mean_pp,sigma_pp,min_TacDiffCut[l],min_TacDiffCut_new[l][bin-1]);
	    }

	  // bin counting
	  for(int bin=1; bin<=nPtBins; bin++)
	    {
	      TH1F *hppTac = hppTacDiff[bin-1];
	      int xbins = hppTac->GetNbinsX(); 
	      pt_arr[0][bin-1] = hppTac->GetBinCenter(bin);
	      all_arr[0][bin-1] = hppTac->IntegralAndError(1, xbins, all_err_arr[0][bin-1]);
	      int low_bin =    hppTac->FindFixBin(min_TacDiffCut_new[l][bin-1]);
	      acc_arr[0][bin-1] = hppTac->IntegralAndError(low_bin+1,
							   hppTac->FindFixBin(max_TacDiffCut[l]-0.1),
							   acc_err_arr[0][bin-1]);
	      double fraction = (hppTac->GetXaxis()->GetBinUpEdge(low_bin)-min_TacDiffCut_new[l][bin-1])/hppTac->GetBinWidth(low_bin);
	      acc_arr[0][bin-1] += hppTac->GetBinContent(low_bin)*fraction;
	      acc_err_arr[0][bin-1] = sqrt(pow(acc_err_arr[0][bin-1],2)+pow(hppTac->GetBinError(low_bin)*fraction,2));
	      if(l==1 || l==3) printf("[i] %s: all = %4.2f, acc = %4.2f, eff = %4.2f%%\n",gTrgSetupTitle[l+1],all_arr[0][bin-1],acc_arr[0][bin-1],acc_arr[0][bin-1]/all_arr[0][bin-1]*100);
	    }

	  // fitting
	  for(int bin=1; bin<=nPtBins; bin++)
	    {
	      TF1 *funcppTac    = funcppTacDiff[bin-1];
	      pt_arr[1][bin-1]  = pt_arr[0][bin-1];
	      all_arr[1][bin-1] = funcppTac->Integral(840,max_TacDiffCut[l]);
	      acc_arr[1][bin-1] = funcppTac->Integral(min_TacDiffCut_new[l][bin-1],max_TacDiffCut[l]);
	      all_err_arr[1][bin-1] = all_err_arr[0][bin-1]/all_arr[0][bin-1] * all_arr[1][bin-1];
	      acc_err_arr[1][bin-1] = acc_err_arr[0][bin-1]/acc_arr[0][bin-1] * acc_arr[1][bin-1];	  
	    }
      
	  gTacDiffEff[l] = GetEfficiencyCurve(nPtBins, pt_arr[1], all_arr[1], all_err_arr[1], acc_arr[1], acc_err_arr[1]);
	  gTacDiffEff[l]->SetName(Form("%s_gTacDiffEff_LS_%s_M%d",run_type.Data(),lumi_name[l],mode));
	}
    }
  if(mode==1)
    {
      // use muon candidates from Run14 MB data
       TFile *fin = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root", "read");
       TH3F *hTacDiffVsPtVsLumi = (TH3F*)fin->Get("mhMT101TacDiffMuon_di_mu");
       //TH3F *hTacDiffVsPtVsLumi = (TH3F*)fin->Get("mhMT101TacDiffMatched_di_mu");
       hTacDiffVsPtVsLumi->Sumw2();
       TH1F *hTacDiff[nLumi][nPtBins];
       for(int l=0; l<nLumi; l++)
	 {
	   printf("+++ %s +++\n",lumi_name[l]);
	  for(int bin=1; bin<=nPtBins; bin++)
	    {
	      int lowbin = hTacDiffVsPtVsLumi->GetZaxis()->FindFixBin(xPtBins[bin-1]+1e-4);
	      int upbin  = hTacDiffVsPtVsLumi->GetZaxis()->FindFixBin(xPtBins[bin]-1e-4);
	      hTacDiff[l][bin-1] = (TH1F*)hTacDiffVsPtVsLumi->ProjectionX(Form("hTacDiff_Pt%d_%s",bin,lumi_name[l]),l+1,l+1,lowbin,upbin);
	      
	      int xbins = hTacDiff[l][bin-1]->GetNbinsX(); 
	      pt_arr[0][bin-1] = (xPtBins[bin-1]+xPtBins[bin])/2.;
	      all_arr[0][bin-1] = hTacDiff[l][bin-1]->IntegralAndError(1, xbins, all_err_arr[0][bin-1]);
	      acc_arr[0][bin-1] = hTacDiff[l][bin-1]->IntegralAndError(hTacDiff[l][bin-1]->FindBin(min_TacDiffCut[l]+0.1),
								       hTacDiff[l][bin-1]->FindBin(max_TacDiffCut[l]-0.1),
								       acc_err_arr[0][bin-1]);
	    }
	  gTacDiffEff[l] = GetEfficiencyCurve(nPtBins, pt_arr[0], all_arr[0], all_err_arr[0], acc_arr[0], acc_err_arr[0]);
	  gTacDiffEff[l]->SetName(Form("%s_gTacDiffEff_LS_%s_M%d",run_type.Data(),lumi_name[l],mode));
	 }
       TCanvas *c = new TCanvas("MT101TacDiffMuon","MT101TacDiffMuon", 1100, 700);
       c->Divide(3,2);
       for(int bin=1; bin<=nPtBins; bin++)
	 {
	   c->cd(bin);
	   for(int l=0; l<nLumi; l++)
	     {
	       hTacDiff[l][bin-1]->Scale(1./hTacDiff[l][bin-1]->Integral());
	       hTacDiff[l][bin-1]->SetMarkerStyle(20+l);
	       hTacDiff[l][bin-1]->SetMarkerColor(gColor[l]);
	       hTacDiff[l][bin-1]->SetLineColor(gColor[l]);
	       hTacDiff[l][bin-1]->GetXaxis()->SetRangeUser(750, 840);
	       if(l==0) hTacDiff[l][bin-1]->Draw();
	       else     hTacDiff[l][bin-1]->Draw("sames");
	     }
	 } 
    }
  if(mode==2)
    {
      // use matched tracks from Run14 MB data
       TFile *fin = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root", "read");
       TH3F *hTacDiffVsPtVsLumi = (TH3F*)fin->Get("mhMT101TacDiffMatched_di_mu");
       hTacDiffVsPtVsLumi->Sumw2();
       TH1F *hTacDiff[nLumi][nPtBins];
       for(int l=0; l<nLumi; l++)
	 {
	   printf("+++ %s +++\n",lumi_name[l]);
	  for(int bin=1; bin<=nPtBins; bin++)
	    {
	      int lowbin = hTacDiffVsPtVsLumi->GetZaxis()->FindFixBin(xPtBins[bin-1]+1e-4);
	      int upbin  = hTacDiffVsPtVsLumi->GetZaxis()->FindFixBin(xPtBins[bin]-1e-4);
	      hTacDiff[l][bin-1] = (TH1F*)hTacDiffVsPtVsLumi->ProjectionX(Form("hTacDiff_Pt%d_%s",bin,lumi_name[l]),l+1,l+1,lowbin,upbin);
	      
	      int xbins = hTacDiff[l][bin-1]->GetNbinsX(); 
	      pt_arr[0][bin-1] = (xPtBins[bin-1]+xPtBins[bin])/2.;
	      all_arr[0][bin-1] = hTacDiff[l][bin-1]->IntegralAndError(1, xbins, all_err_arr[0][bin-1]);
	      acc_arr[0][bin-1] = hTacDiff[l][bin-1]->IntegralAndError(hTacDiff[l][bin-1]->FindBin(min_TacDiffCut[l]+0.1),
								       hTacDiff[l][bin-1]->FindBin(max_TacDiffCut[l]-0.1),
								       acc_err_arr[0][bin-1]);
	    }
	  gTacDiffEff[l] = GetEfficiencyCurve(nPtBins, pt_arr[0], all_arr[0], all_err_arr[0], acc_arr[0], acc_err_arr[0]);
	  gTacDiffEff[l]->SetName(Form("%s_gTacDiffEff_LS_%s_M%d",run_type.Data(),lumi_name[l],mode));
	 }
       TCanvas *c = new TCanvas("MT101TacDiffMuon","MT101TacDiffMuon", 1100, 700);
       c->Divide(5,2);
       for(int bin=1; bin<=nPtBins; bin++)
	 {
	   c->cd(bin);
	   for(int l=0; l<nLumi; l++)
	     {
	       hTacDiff[l][bin-1]->Scale(1./hTacDiff[l][bin-1]->Integral());
	       hTacDiff[l][bin-1]->SetMarkerStyle(20+l);
	       hTacDiff[l][bin-1]->SetMarkerColor(gColor[l]);
	       hTacDiff[l][bin-1]->SetLineColor(gColor[l]);
	       hTacDiff[l][bin-1]->GetXaxis()->SetRangeUser(650, 840);
	       if(l==0) hTacDiff[l][bin-1]->Draw();
	       else     hTacDiff[l][bin-1]->Draw("sames");
	     }
	 } 
    }

  TFile *fRun14 = TFile::Open("Rootfiles/Run14_AuAu200.MtdTrigEff.root","read");
  TGraphAsymmErrors *gRun14 = (TGraphAsymmErrors*)fRun14->Get("Run14_AuAu200_gTacDiffEffFinal_BinCount_prod_Run15_pp200");
  double x, y, x1, y1;
  for(int l=0; l<nLumi; l++)
    {
      int npoints = gTacDiffEff[l]->GetN();
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  gTacDiffEff[l]->GetPoint(ipoint,x,y);
	  if(mode<2)
	    {
	      gRun14->GetPoint(ipoint, x1, y1);
	      gTacDiffEff[l]->SetPoint(ipoint,x1,y);
	      gTacDiffEff[l]->SetPointEXhigh(ipoint,gRun14->GetErrorXhigh(ipoint));
	      gTacDiffEff[l]->SetPointEXlow(ipoint,gRun14->GetErrorXlow(ipoint));
	    }
	  else
	    {
	      gTacDiffEff[l]->SetPointEXhigh(ipoint,0.1);
	      gTacDiffEff[l]->SetPointEXlow(ipoint,0.1);
	    }
	}
    }

  TCanvas *cFitTrigEff = new TCanvas("cFitTrigEff", "cFitTrigEff", 1100, 700);
  cFitTrigEff->Divide(2,2);
  TF1 *funcTacEff[nLumi];
  TH1F *hplot = new TH1F("hplot",";p_{T}^{#mu} [GeV/c];Efficiency",100,1,20);
  if(mode<2) hplot->GetXaxis()->SetRangeUser(1,7);
  else hplot->GetXaxis()->SetRangeUser(1,10);
  for(int l=0; l<nLumi; l++)
    {
      funcTacEff[l] = new TF1(Form("%s_gTacDiffEff_LS_%s_func_M%d",run_type.Data(),lumi_name[l],mode),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
      if(mode==2) funcTacEff[l]->SetRange(1.3, 10);
      funcTacEff[l]->SetParameters(0.9, 3.5, 0.8);
      TGraphAsymmErrors *gFitTmp = (TGraphAsymmErrors*)gTacDiffEff[l]->Clone(Form("%s_clone",gTacDiffEff[l]->GetName()));
      gFitTmp->Fit(funcTacEff[l],"IR0");
      cFitTrigEff->cd(l+1);
      hplot->GetYaxis()->SetRangeUser(0.2,1.1);
      hplot->GetYaxis()->SetNdivisions(505);
      hplot->Draw();
      gTacDiffEff[l]->SetMarkerStyle(21);
      gTacDiffEff[l]->SetMarkerColor(1);
      gTacDiffEff[l]->SetLineColor(1);
      gTacDiffEff[l]->Draw("PEsame");
      funcTacEff[l]->SetLineColor(2);
      funcTacEff[l]->SetLineStyle(2);
      funcTacEff[l]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s: estimated MTD trigger efficiency for muon candidates",run_type.Data()));
      if(mode==2) t1 = GetTitleText(Form("%s: estimated MTD trigger efficiency for matched tracks",run_type.Data()));
      t1->Draw();
      if(l==0)
	{
	  TLegend *leg = new TLegend(0.45,0.15,0.65,0.35);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  leg->SetHeader(lumi_name[l]);
	  leg->AddEntry(gTacDiffEff[l],"Data: bin counting","P");
	  leg->AddEntry(funcTacEff[l],"Fit: p0-e^{-p1*(x-p2)}","L");
	  leg->Draw();
	}
    }
  if(savePlot)
    cFitTrigEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_FitTacDiffEff_M%d.pdf",run_type.Data(),mode));

  if(saveHisto)
    {
      TFile *fout = NULL;
      if(mode==0) fout = fin;
      else fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyTrigEff.root", "update");
      fout->cd();
      if(mode==0)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      funcppTacDiff[i]->Write("",TObject::kOverwrite);
	    }
	}
      for(int l=0; l<nLumi; l++)
	{
	  gTacDiffEff[l]->Write("",TObject::kOverwrite);
	  funcTacEff[l]->Write("",TObject::kOverwrite);
	} 
    }
  
}

//================================================
void anaTrigEff(const int savePlot = 1, const int saveHisto = 1)
{
  //==============================================
  // compare the TacDiff in MT101
  TFile *fStudy = TFile::Open("output/Run14_AuAu200.Study.MtdTrig.root", "read");
  TH2F *h2TacDiff = (TH2F*)fStudy->Get("mhMT101TacDiff_di_mu");
  TH1F *hTacDiffMT101[4];
  TF1 *funcTacDiffMT101[4];
  TCanvas *c = new TCanvas("cTacDiffMT101", "cTacDiffMT101", 800, 600);
  leg = new TLegend(0.15,0.6,0.4,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader(run_type.Data());
  for(int j=0; j<gNTrgSetup-1; j++)
    {
      hTacDiffMT101[j] = (TH1F*)h2TacDiff->ProjectionX(Form("hTacDiffMT101_%s",gTrgSetupName[j+1]), j+1, j+1);
      hTacDiffMT101[j]->Sumw2();
      hTacDiffMT101[j]->SetXTitle("TacSum_{MT101} - TacSum_{VPD}/8 + 1024");
      hTacDiffMT101[j]->SetMarkerStyle(20+j+1);
      hTacDiffMT101[j]->SetMarkerColor(gColor[j+1]);
      hTacDiffMT101[j]->SetLineColor(gColor[j+1]);
      hTacDiffMT101[j]->SetMarkerSize(1.2);
      hTacDiffMT101[j]->Scale(1./hTacDiffMT101[j]->GetBinContent(hTacDiffMT101[j]->FindBin(790)));
      hTacDiffMT101[j]->GetYaxis()->SetRangeUser(1e-2,1.5);
      hTacDiffMT101[j]->GetXaxis()->SetRangeUser(760, 820);
      if(j==0) hTacDiffMT101[j]->DrawCopy();
      else     hTacDiffMT101[j]->DrawCopy("sames");
      leg->AddEntry(hTacDiffMT101[j], legName[j+1].Data(), "P");
      funcTacDiffMT101[j] = new TF1(Form("func_%s",hTacDiffMT101[j]->GetName()), "gaus", 790, 805);
      if(j<2) funcTacDiffMT101[j]->SetRange(786, 805);
      if(j==3) funcTacDiffMT101[j]->SetRange(790, 805);
      funcTacDiffMT101[j]->SetParameters(1, 795, 8);
      hTacDiffMT101[j]->Fit(funcTacDiffMT101[j], "IR0Q");
      funcTacDiffMT101[j]->SetLineColor(gColor[j+1]);
      funcTacDiffMT101[j]->SetLineStyle(2);
      funcTacDiffMT101[j]->Draw("sames");
      cout << hTacDiffMT101[j]->Integral(787,789)/hTacDiffMT101[j]->Integral(787,1024) << endl;
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_TacDiffInLumi_All.pdf",run_type.Data()));


  //==============================================
  // matched tracks
  const char* type_name[3] = {"Muon","UL", "LS"};
  const char* type_title[3] = {"Muon = UL - LS","Unlike-sign", "Like-sign"};
  const int nTrigUnit = 28;
  const int nPtBins = 6;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 
  const double fit_min = 789;
  const double fit_max = 806;

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.StudyTrigEff.root",run_type.Data()), "update");
  else          fin = TFile::Open(Form("Rootfiles/%s.StudyTrigEff.root",run_type.Data()), "read");

  TCanvas *c = new TCanvas("cTacDiff", "cTacDiff", 1100, 700);
  c->Divide(2, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type.Data());
  TH1F *hTacDiffAll[3][gNTrgSetup];
  for(int i=0; i<3; i++)
    {
      c->cd(i+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffAll[i][j] = (TH1F*)fin->Get(Form("hTacDiffAll_%s%s",type_name[i],gTrgSetupName[j]));
	  hTacDiffAll[i][j]->SetXTitle("#DeltaTacSum");
	  scaleHisto(hTacDiffAll[i][j], 1, 1, true, false, false);
	  hTacDiffAll[i][j]->SetMarkerStyle(20+j);
	  hTacDiffAll[i][j]->SetMarkerColor(gColor[j]);
	  hTacDiffAll[i][j]->SetLineColor(gColor[j]);
	  hTacDiffAll[i][j]->SetMarkerSize(1.2);
	  hTacDiffAll[i][j]->Scale(1./hTacDiffAll[i][j]->GetBinContent(hTacDiffAll[i][j]->FindBin(795)));
	  hTacDiffAll[i][j]->GetYaxis()->SetRangeUser(1e-2,1.2);
	  hTacDiffAll[i][j]->GetXaxis()->SetRangeUser(780, 820);
	  if(j==0) hTacDiffAll[i][j]->DrawCopy();
	  else     hTacDiffAll[i][j]->DrawCopy("sames");
	  if(i==0) leg->AddEntry(hTacDiffAll[i][j], legName[j].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(type_title[i],0.05);
      t1->Draw();
    }
  c->cd(4);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_TacDiffInLumi_Muon.pdf",run_type.Data()));
  for(int i=0; i<3; i++)
    {
      c->cd(i+1);
      gPad->SetLogy();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_TacDiffInLumi_Muon_Log.pdf",run_type.Data()));

  //==============================================
  // matched tracks vs MT101 vs pp shifted muons
  TFile *flumi = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  TF1 *funcppLS = (TF1*)flumi->Get("fit_Run15_pp200_hTacDiff_LS_PtBin1");
  TFile *fMtdTrig = TFile::Open("Rootfiles/Run14_AuAu200.MtdTrigEff.root", "read");
  TH1F *hppShifted = (TH1F*)fMtdTrig->Get("Run15_pp200_hMuonTacDiffCombined_PtBin0");
  hppShifted->Scale(1./hppShifted->GetBinContent(hppShifted->FindBin(795)));
  TCanvas *c = new TCanvas("TacSum_LSvsMT101", "TacSum_LSvsMT101", 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.6,0.35,0.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  leg->SetHeader(run_type.Data());
  for(int j=0; j<gNTrgSetup-1; j++)
    {
      double tacSumCutMin = 786;
      double tacSumCutMax = 837;
      if(j>1) tacSumCutMin = 789;

      c->cd(j+1);
      hTacDiffMT101[j]->SetMarkerStyle(20);
      hTacDiffMT101[j]->SetMarkerColor(1);
      hTacDiffMT101[j]->SetLineColor(1);
      hTacDiffMT101[j]->SetTitle(";#DeltaTacSum;a.u.");
      hTacDiffMT101[j]->DrawCopy();
      funcTacDiffMT101[j]->SetLineColor(kGray);
      funcTacDiffMT101[j]->Draw("sames");
      double tacSumCutMinNew = funcppLS->GetParameter(1) - (funcTacDiffMT101[j]->GetParameter(1)-tacSumCutMin)/funcTacDiffMT101[j]->GetParameter(2)*funcppLS->GetParameter(2);
      double tacSumCutMaxNew = funcppLS->GetParameter(1) + (tacSumCutMax - funcTacDiffMT101[j]->GetParameter(1))/funcTacDiffMT101[j]->GetParameter(2)*funcppLS->GetParameter(2);
      double eff = funcppLS->Integral(tacSumCutMinNew,tacSumCutMaxNew)/funcppLS->Integral(900,1200);
      TPaveText *t1 = GetPaveText(0.65,0.85,0.75,0.85,0.05);
      t1->AddText(Form("#mu = %4.2f, #sigma = %4.2f",funcTacDiffMT101[j]->GetParameter(1),funcTacDiffMT101[j]->GetParameter(2)));
      t1->AddText(Form("eff = %4.3f",eff));
      t1->Draw();

      hppShifted->SetMarkerStyle(25);
      hppShifted->SetMarkerColor(4);
      hppShifted->SetLineColor(4);
      hppShifted->Draw("sames");

      hTacDiffAll[2][j+1]->SetMarkerStyle(24);
      hTacDiffAll[2][j+1]->SetMarkerColor(2);
      hTacDiffAll[2][j+1]->SetLineColor(2);
      hTacDiffAll[2][j+1]->DrawCopy("sames");
      TF1 *functmp = new TF1(Form("fit_%s",hTacDiffAll[2][j+1]->GetName()),"gaus",786,805);
      if(j>1) functmp->SetRange(790,805);
      functmp->SetParameters(1, 795, 8);
      hTacDiffAll[2][j+1]->Fit(functmp, "IR0Q");
      functmp->SetLineStyle(2);
      functmp->SetLineColor(2);
      functmp->Draw("sames");
      tacSumCutMinNew = funcppLS->GetParameter(1) - (functmp->GetParameter(1)-tacSumCutMin)/functmp->GetParameter(2)*funcppLS->GetParameter(2);
      tacSumCutMaxNew = funcppLS->GetParameter(1) + (tacSumCutMax - functmp->GetParameter(1))/functmp->GetParameter(2)*funcppLS->GetParameter(2);
      eff = funcppLS->Integral(tacSumCutMinNew,tacSumCutMaxNew)/funcppLS->Integral(900,1200);
      t1 = GetPaveText(0.65,0.85,0.62,0.72,0.05);
      t1->AddText(Form("#mu = %4.2f, #sigma = %4.2f",functmp->GetParameter(1),functmp->GetParameter(2)));
      t1->AddText(Form("eff = %4.3f",eff));
      t1->SetTextColor(2);
      t1->Draw();
      TPaveText *t1 = GetTitleText(Form("%s: #DeltaTacSum distribution",legName[j+1].Data()),0.055);
      t1->Draw();
      if(j==0)
	{
	  leg->AddEntry(hTacDiffMT101[j], "All hits", "P");
	  leg->AddEntry(hTacDiffAll[2][j+1], "Like-sign pairs", "P");
	  leg->AddEntry(hppShifted, "pp muons (shifted)", "P");
	}
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacDiff_AuAuvspp.pdf",run_type.Data()));

  //==============================================
  // fit TacSum distributions
  TH1F *hTacDiffPt[nPtBins][gNTrgSetup];
  TF1 *funcTacDiffPt[nPtBins][gNTrgSetup];
  TH1F *hTacDiffPtMean[gNTrgSetup];
  TH1F *hTacDiffPtSigma[gNTrgSetup];  
  for(int j=0; j<gNTrgSetup; j++)
    {
      hTacDiffPtMean[j]  = new TH1F(Form("hTacDiffMeanVsPt%s",gTrgSetupName[j]),"Mean of #DeltaTacSum;p_{T} (GeV/c);<#DeltaTacSum>", nPtBins, xPtBins);
      hTacDiffPtSigma[j] = new TH1F(Form("hTacDiffSigmaVsPt%s",gTrgSetupName[j]),"Width of #DeltaTacSum;p_{T} (GeV/c);#sigma(#DeltaTacSum)", nPtBins, xPtBins);      
    }
  for(int t=0; t<nPtBins; t++)
    {
      TCanvas *c = new TCanvas(Form("cFitTacDiff_%d",t), Form("cFitTacDiff_Pt%d",t), 1100, 700);
      c->Divide(3, 2);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffPt[t][j] = (TH1F*)fin->Get(Form("hTacDiff_Pt%d_%s%s",t+1,type_name[2],gTrgSetupName[j]));
	  if(t==nPtBins-1) hTacDiffPt[t][j]->Rebin(2);
	  scaleHisto(hTacDiffPt[t][j], 1, 1, true, false, false);
	  double fit_min_tmp = fit_min;
	  double fit_max_tmp = fit_max;
	  //if(j==1 || j==2) fit_min_tmp = 786;
	  if(t>=nPtBins-1) fit_max_tmp = 815;
	  funcTacDiffPt[t][j] = new TF1(Form("fit_%s",hTacDiffPt[t][j]->GetName()),"gaus",fit_min_tmp,fit_max_tmp);
	  funcTacDiffPt[t][j]->SetParameters(1, 795, 8);
	  hTacDiffPt[t][j]->Fit(funcTacDiffPt[t][j], "IR0Q");
	  c->cd(j+1);
	  hTacDiffPt[t][j]->GetXaxis()->SetRangeUser(780, 820);
	  hTacDiffPt[t][j]->SetMarkerStyle(20);
	  hTacDiffPt[t][j]->Draw();
	  funcTacDiffPt[t][j]->SetLineColor(2);
	  funcTacDiffPt[t][j]->SetLineStyle(2);
	  funcTacDiffPt[t][j]->Draw("sames");
	  TPaveText *t1 = GetPaveText(0.6, 0.8, 0.7, 0.85, 0.045);
	  t1->AddText(Form("#mu = %4.2f#pm%4.2f",funcTacDiffPt[t][j]->GetParameter(1),funcTacDiffPt[t][j]->GetParError(1)));
	  t1->AddText(Form("#sigma = %4.2f#pm%4.2f",funcTacDiffPt[t][j]->GetParameter(2),funcTacDiffPt[t][j]->GetParError(2)));
	  t1->Draw();
	  hTacDiffPtMean[j]->SetBinContent(t+1, funcTacDiffPt[t][j]->GetParameter(1));
	  hTacDiffPtMean[j]->SetBinError(t+1, funcTacDiffPt[t][j]->GetParError(1));
	  hTacDiffPtMean[j]->SetMarkerStyle(20+j);
	  hTacDiffPtMean[j]->SetMarkerColor(gColor[j]);
	  hTacDiffPtMean[j]->SetLineColor(gColor[j]);
	  hTacDiffPtMean[j]->SetMarkerSize(1.2);
	  hTacDiffPtSigma[j]->SetBinContent(t+1, funcTacDiffPt[t][j]->GetParameter(2));
	  hTacDiffPtSigma[j]->SetBinError(t+1, funcTacDiffPt[t][j]->GetParError(2));
	  hTacDiffPtSigma[j]->SetMarkerStyle(20+j);
	  hTacDiffPtSigma[j]->SetMarkerColor(gColor[j]);
	  hTacDiffPtSigma[j]->SetLineColor(gColor[j]);
	  hTacDiffPtSigma[j]->SetMarkerSize(1.2);
	  TPaveText *t1 = GetTitleText(Form("%s: %1.1f < p_{T} < %1.1f GeV/c",legName[j].Data(),xPtBins[t],xPtBins[t+1]),0.05);
	  t1->Draw();
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_FitTacDiffLS_pt%d.pdf",run_type.Data(),t));
    }  

  TCanvas *c = new TCanvas("cTacDiffVsPt", "cTacDiffVsPt", 1100, 500);
  c->Divide(2,1);
  leg = new TLegend(0.5,0.15,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type.Data());
  c->cd(1);
  SetPadMargin(gPad, 0.11, 0.13, 0.1, 0.1);
  for(int j=0; j<gNTrgSetup; j++)
    {
      ScaleHistoTitle(hTacDiffPtMean[j], 0.045, 1.1, 0.04, 0.045, 1.2, 0.04, 62);
      hTacDiffPtMean[j]->SetTitle("");
      hTacDiffPtMean[j]->GetYaxis()->SetRangeUser(792, 802);
      if(j==0) hTacDiffPtMean[j]->Draw();
      else     hTacDiffPtMean[j]->Draw("sames");
      leg->AddEntry(hTacDiffPtMean[j], legName[j].Data(), "P");
    }
  TPaveText *t1 = GetTitleText("Mean of #DeltaTacSum distribution",0.05);
  t1->Draw();
  leg->Draw();
  c->cd(2);
  SetPadMargin(gPad, 0.11, 0.13, 0.1, 0.1);
  for(int j=0; j<gNTrgSetup; j++)
    {
      ScaleHistoTitle(hTacDiffPtSigma[j], 0.045, 1.1, 0.04, 0.045, 1.2, 0.04, 62);
      hTacDiffPtSigma[j]->SetTitle("");
      hTacDiffPtSigma[j]->GetYaxis()->SetRangeUser(4, 9);
      if(j==0) hTacDiffPtSigma[j]->Draw();
      else     hTacDiffPtSigma[j]->Draw("sames");
    }
  TPaveText *t1 = GetTitleText("Width of #DeltaTacSum distribution",0.05);
  t1->Draw();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_TacDiffLS_FitVsPt.pdf",run_type.Data()));

  TH1F *hTacDiffTrigUnit[nPtBins][nTrigUnit][gNTrgSetup];
  TF1 *funcTacDiffTrigUnit[nPtBins][nTrigUnit][gNTrgSetup];
  TH1F *hTacDiffMean[nPtBins][gNTrgSetup];
  TH1F *hTacDiffSigma[nPtBins][gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {
      for(int t=0; t<nPtBins-1; t++)
	{
	  hTacDiffMean[t][j]  = new TH1F(Form("hTacDiffMean_Pt%d%s",t,gTrgSetupName[j]), "Mean of #DeltaTacSum;TrigUnit;<#DeltaTacSum>", nTrigUnit, 0, nTrigUnit);
	  hTacDiffSigma[t][j] = new TH1F(Form("hTacDiffSigma_Pt%d%s",t,gTrgSetupName[j]), "Width of #DeltaTacSum;TrigUnit;#sigma(#DeltaTacSum)", nTrigUnit, 0, nTrigUnit);	  
	  TCanvas *c = new TCanvas(Form("cFitTacDiff_%d%d",j,t), Form("cFitTacDiff_Pt%d%s",t,gTrgSetupName[j]), 1100, 700);
	  c->Divide(6, 5);
	  for(int k=0; k<nTrigUnit; k++)
	    {
	      hTacDiffTrigUnit[t][k][j] = (TH1F*)fin->Get(Form("hTacDiff_Pt%d_Unit%d_%s%s",t+1,k+1,type_name[2],gTrgSetupName[j]));
	      scaleHisto(hTacDiffTrigUnit[t][k][j], 1, 1, true, false, false);
	      double fit_min_tmp = 789;
	      double fit_max_tmp = 810;
	      if(j==1 || j==2) fit_min_tmp = 786;
	      if(k==1) fit_max_tmp = 799;
	      if(k==4 || k==5) fit_max_tmp = 803;
	      funcTacDiffTrigUnit[t][k][j] = new TF1(Form("fit_%s",hTacDiffTrigUnit[t][k][j]->GetName()),"gaus",fit_min_tmp,fit_max_tmp);
	      funcTacDiffTrigUnit[t][k][j]->SetParameters(1, 795, 8);
	      hTacDiffTrigUnit[t][k][j]->Fit(funcTacDiffTrigUnit[t][k][j], "IR0Q");
	      c->cd(k+1);
	      hTacDiffTrigUnit[t][k][j]->GetXaxis()->SetRangeUser(780, 820);
	      hTacDiffTrigUnit[t][k][j]->SetMarkerStyle(24);
	      hTacDiffTrigUnit[t][k][j]->SetMarkerSize(1.2);
	      hTacDiffTrigUnit[t][k][j]->Draw();
	      funcTacDiffTrigUnit[t][k][j]->Draw("sames");
	      hTacDiffMean[t][j]->SetBinContent(k+1, funcTacDiffTrigUnit[t][k][j]->GetParameter(1));
	      hTacDiffMean[t][j]->SetBinError(k+1, funcTacDiffTrigUnit[t][k][j]->GetParError(1));
	      hTacDiffMean[t][j]->SetMarkerStyle(20+j);
	      hTacDiffMean[t][j]->SetMarkerColor(gColor[j]);
	      hTacDiffMean[t][j]->SetLineColor(gColor[j]);
	      hTacDiffMean[t][j]->SetMarkerSize(1.2);
	      hTacDiffSigma[t][j]->SetBinContent(k+1, funcTacDiffTrigUnit[t][k][j]->GetParameter(2));
	      hTacDiffSigma[t][j]->SetBinError(k+1, funcTacDiffTrigUnit[t][k][j]->GetParError(2));
	      hTacDiffSigma[t][j]->SetMarkerStyle(20+j);
	      hTacDiffSigma[t][j]->SetMarkerColor(gColor[j]);
	      hTacDiffSigma[t][j]->SetLineColor(gColor[j]);
	      hTacDiffSigma[t][j]->SetMarkerSize(1.2);
	    }
	}
    }
  
  TCanvas *c = new TCanvas("cTacDiffMean", "cTacDiffMean", 1100, 700);
  c->Divide(3, 2);
  leg = new TLegend(0.2,0.4,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(run_type.Data());
  for(int t=0; t<nPtBins-1; t++)
    {
      c->cd(t+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffMean[t][j]->SetTitle("");
	  hTacDiffMean[t][j]->GetYaxis()->SetRangeUser(785, 810);
	  if(j==0) hTacDiffMean[t][j]->Draw();
	  else     hTacDiffMean[t][j]->Draw("sames");
	  if(t==0) leg->AddEntry(hTacDiffMean[t][j], legName[j].Data(), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xPtBins[t],xPtBins[t+1]),0.05);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_TacDiffLS_MeanVsTrigUnit.pdf",run_type.Data()));

  TCanvas *c = new TCanvas("cTacDiffSigma", "cTacDiffSigma", 1100, 700);
  c->Divide(3, 2);
  for(int t=0; t<nPtBins-1; t++)
    {
      c->cd(t+1);
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffSigma[t][j]->SetTitle("");
	  hTacDiffSigma[t][j]->GetYaxis()->SetRangeUser(0, 20);
	  if(j==0) hTacDiffSigma[t][j]->Draw();
	  else     hTacDiffSigma[t][j]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xPtBins[t],xPtBins[t+1]),0.05);
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuDm_TacDiffLS_WidthVsTrigUnit.pdf",run_type.Data()));

  if(saveHisto)
    {
      fin->cd();
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hTacDiffPtMean[j]->Write("",TObject::kOverwrite);
	  hTacDiffPtSigma[j]->Write("",TObject::kOverwrite);
	}
    }	  
}

//================================================
void makeTrigEff(const int saveHisto = 1)
{
  const int nBinsTacDiff1 = 35;
  const double xBinsTacDiff1[nBinsTacDiff1+1] = {760,765,770,775,780,782,784,786,787,788,789,790,791,792,793,795,797,799,801,803,805,807,809,811,813,815,817,819,821,823,825,827,829,833,837,841};

  const int nBinsTacDiff2 = 19;
  const double xBinsTacDiff2[nBinsTacDiff2+1] = {760,765,770,775,780,782,784,786,788,789,790,792,795,801,807,813,819,825,833,841};
  const char* type_name[3] = {"Muon","UL", "LS"};
  const int nTrigUnit = 28;
  const int nPtBins = 6;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 


  TFile *fin = TFile::Open("output/Run14_AuAu200.JpsiMuon.root", "read");
  THnSparseF *hn = (THnSparseF*)fin->Get("mhJpsiMuonTrigEff_di_mu");
  hn->GetAxis(0)->SetRangeUser(3.0, 3.2); // cut in Jpsi region
  
  TH1F *hTacDiffAll[3][gNTrgSetup];
  TH1F *hTacDiffPt[2][nPtBins][gNTrgSetup];
  TH1F *hTacDiffTrigUnit[2][nPtBins][nTrigUnit][gNTrgSetup];

  // get UL and LS
  for(int j=0; j<gNTrgSetup; j++)
    {
      if(j>0) hn->GetAxis(7)->SetRange(j, j);
      for(int i=0; i<2; i++)
	{
	  hn->GetAxis(5)->SetRange(i+1,i+1);
	  hn->GetAxis(2)->SetRangeUser(1.3, 10);
	  TH1F *h1tmp = (TH1F*)hn->Projection(1);
	  h1tmp->SetName(Form("hTacDiffAll_%d%d",i,j));
	  //hTacDiffAll[i+1][j] = (TH1F*)h1tmp->Rebin(nBinsTacDiff1, Form("hTacDiffAll_%s%s",type_name[i+1],gTrgSetupName[j]), xBinsTacDiff1);
	  hTacDiffAll[i+1][j] = (TH1F*)h1tmp->Clone(Form("hTacDiffAll_%s%s",type_name[i+1],gTrgSetupName[j]));
	  hTacDiffAll[i+1][j]->Sumw2();
	  for(int t=0; t<nPtBins; t++)
	    {
	      hn->GetAxis(2)->SetRangeUser(xPtBins[t], xPtBins[t+1]);
	      TH1F *h1tmp = (TH1F*)hn->Projection(1);
	      h1tmp->SetName(Form("hTacDiffPt%d%d%d",t,i,j));
	      hTacDiffPt[i][t][j] = (TH1F*)h1tmp->Rebin(nBinsTacDiff1, Form("hTacDiff_Pt%d_%s%s",t+1,type_name[i+1],gTrgSetupName[j]), xBinsTacDiff1);
	      TH2F *h2tmp = (TH2F*)hn->Projection(1, 3);
	      h2tmp->SetName(Form("hTacDiffVsTrigUnit_%d%d%d",j,i,t));
	      h2tmp->Sumw2();
	      for(int k=0; k<nTrigUnit; k++)
		{
		  TH1F *h1tmp = (TH1F*)h2tmp->ProjectionY(Form("hTacDiffTrigUnit_%d%d%d%d",i,j,t,k), k+2, k+2);
		  hTacDiffTrigUnit[i][t][k][j] = (TH1F*)h1tmp->Rebin(nBinsTacDiff1, Form("hTacDiff_Pt%d_Unit%d_%s%s",t+1,k+1,type_name[i+1],gTrgSetupName[j]), xBinsTacDiff1);
		  //hTacDiffTrigUnit[i][t][k][j] = (TH1F*)h2tmp->ProjectionY(Form("hTacDiff_Pt%d_Unit%d_%s%s",t+1,k+1,type_name[i+1],gTrgSetupName[j]), k+2, k+2);
		}
	    }
	  hn->GetAxis(2)->SetRange(0,-1);
	}
      hn->GetAxis(5)->SetRange(0,-1);
    }
  hn->GetAxis(7)->SetRange(0,-1);

  // get pure muon
  for(int j=0; j<gNTrgSetup; j++)
    {
      hTacDiffAll[0][j] = (TH1F*)hTacDiffAll[1][j]->Rebin(nBinsTacDiff2, Form("hTacDiffAll_%s%s",type_name[0],gTrgSetupName[j]), xBinsTacDiff2);
      TH1F *h1Rebin = (TH1F*)hTacDiffAll[2][j]->Rebin(nBinsTacDiff2, Form("%s_rebin",hTacDiffAll[2][j]->GetName()), xBinsTacDiff2);
      hTacDiffAll[0][j]->Add(h1Rebin, -1);
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.StudyTrigEff.root",run_type.Data()), "update");
      for(int j=0; j<gNTrgSetup; j++)
	{
	  for(int i=0; i<3; i++)
	    {
	      hTacDiffAll[i][j]->SetTitle("");
	      hTacDiffAll[i][j]->Write("",TObject::kOverwrite);
	      if(i==0) continue;
	      for(int t=0; t<nPtBins; t++)
		{
		  hTacDiffPt[i-1][t][j]->SetTitle("");
		  hTacDiffPt[i-1][t][j]->Write("",TObject::kOverwrite);
		  for(int k=0; k<nTrigUnit; k++)
		    {
		      hTacDiffTrigUnit[i-1][t][k][j]->SetTitle("");
		      hTacDiffTrigUnit[i-1][t][k][j]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
      fout->Close();
    }
  
}


//================================================
void mtdTrigTime(const int savePlot = 0)
{
  TFile *fstudy = TFile::Open("output/Run14_AuAu200.MB.Study.root", "read");
  THnSparseF *hnMtdTrigTime = (THnSparseF*)fstudy->Get("mhMtdTrigTime_mb");


  TH1F *hTrigTime[30][gNTrgSetup];
  for(int i=0; i<30; i++)
    {
      hnMtdTrigTime->GetAxis(0)->SetRange(i+1, i+1);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hnMtdTrigTime->GetAxis(3)->SetRange(k, k);
	  hTrigTime[i][k] = (TH1F*)hnMtdTrigTime->Projection(1);
	  hTrigTime[i][k]->SetName(Form("hTrigTime_BL%d%s",i+1,gTrgSetupTitle[k]));
	  hTrigTime[i][k]->Sumw2();
	  if(hTrigTime[i][k]->Integral()>0)
	    hTrigTime[i][k]->Scale(1./hTrigTime[i][k]->Integral(0,-1));
	  hnMtdTrigTime->GetAxis(3)->SetRange(0,-1);
	}
      hnMtdTrigTime->GetAxis(0)->SetRange(0,-1);
    }

  TCanvas *cBL[3];
  for(int i=0; i<3; i++)
    {
      cBL[i] = new TCanvas(Form("cBL_%d-%d",1+i*12,13+i*12), Form("cBL_%d-%d",1+i*12,13+i*12), 1100, 700);
      cBL[i]->Divide(4,3);
    }
  for(int i=0; i<30; i++)
    {
      cBL[i/12]->cd(i%12+1);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  hTrigTime[i][k]->SetTitle("");
	  hTrigTime[i][k]->SetMarkerStyle(20+k);
	  hTrigTime[i][k]->SetMarkerColor(gColor[k]);
	  hTrigTime[i][k]->SetLineColor(gColor[k]);
	  hTrigTime[i][k]->SetMaximum(1.2*hTrigTime[i][k]->GetMaximum());
	  if(i<15) hTrigTime[i][k]->GetXaxis()->SetRangeUser(2810, 2920);
	  else hTrigTime[i][k]->GetXaxis()->SetRangeUser(2870, 2970);
	  if(k==0) hTrigTime[i][k]->Draw("P");
	  else hTrigTime[i][k]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("BL = %d",i+1), 0.08);
      t1->Draw();

      if(i%12==0)
	{
	  leg = new TLegend(0.6,0.6,0.8,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  for(int k=1; k<gNTrgSetup; k++)
	    {
	      leg->AddEntry(hTrigTime[i][k], gTrgSetupTitle[k], "PL");
	    }
	  leg->Draw();
	}
    }

  if(savePlot)
    {
      for(int i=0; i<3; i++)
	{
	  cBL[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_MtdTrigTime%d.pdf",run_type.Data(),i));
	}
    }
	  
}
