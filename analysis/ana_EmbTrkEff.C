const char *run_config = "";
const int year = YEAR;
TString run_cfg_name;
const int nDet = 4;
const char *det_name[nDet] = {"","Tpc","Mtd","Final"};
const char *charge_name[3] = {"","_pos","_neg"};
const char *charge_title[3] = {" "," positive "," negative "};

TFile *f;

//================================================
void ana_EmbTrkEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
  TString fileName;
  TString outName;
  TString outPDF;

  //const char *particle = "Upsilon";
  const char *particle = "Jpsi";

  if(year==2013)
    {
      fileName = Form("Run13.pp500.jpsi.Embed.root");
    }
  else if(year==2014)
    {
      fileName = Form("Run14_AuAu200.Embed.%s.root",particle);
    }
  outName = Form("%s.EmbTrkEff.root",run_type);
  outPDF = Form("%s.EmbTrkEff.%s.pdf",run_type,particle);

  run_cfg_name = run_config;

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("# of events: %4.4e\n",hStat->GetBinContent(3));

  //efficiency(outName);
  //effVsCent();
  effVsEta();
  //TrkEff3D(outName, outPDF);
}

//================================================
void effVsCent(const int savePlot = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  TFile *fin  = TFile::Open(Form("Rootfiles/%s.EmbTrkEff.root",run_type),"read");
  TList *list = new TList;

  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }
  TString legName_trg[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      legName_trg[k] = Form("AuAu_200%s",gTrgSetupTitle[k+1]);
    }

  // tracking efficiency
  const int nEffType = 5;
  const char *trkEffType[nEffType] = {"MC","Tpc","MtdMth","MuonPid","MtdTrig"};
  TH1F *hMcTrkPt[nEffType][gNTrgSetup][nCentBins];
  TH1F *hMcTrkPtEff[nEffType][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsZdc[nEffType][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsCent[nEffType][gNTrgSetup];

  for(int i=0; i<nEffType; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hMcTrkPtVsCent[i][j] = (TH2F*)fin->Get(Form("McTrkPtVsCent_%s%s",trkEffType[i],gTrgSetupTitle[j]));
	  hMcTrkPtVsCent[i][j]->RebinX(2);
	  for(int k=0; k<nCentBins; k++)
	    {
	      hMcTrkPtVsZdc[i][j][k] = (TH2F*)fin->Get(Form("McTrkPtVsZdc_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      hMcTrkPt[i][j][k] = (TH1F*)fin->Get(Form("McTrkPt_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      if(k<4) hMcTrkPt[i][j][k]->Rebin(5);
	      else hMcTrkPt[i][j][k]->Rebin(20);
	      if(i>0) 
		{
		  hMcTrkPtEff[i][j][k] = DivideTH1ForEff(hMcTrkPt[i][j][k], hMcTrkPt[i-1][j][k], Form("hMcTrkPtEff_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
		}	
	    }
	}
    }

  TCanvas *c = new TCanvas("TpcEff_Lumi","TpcEff_Lumi",1100,700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.4,0.2,0.7,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int k=1; k<nCentBins; k++)
    {
      c->cd(k);
      SetPadMargin(gPad,0.13, 0.13, 0.03,0.1);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  TH1F *htmp = (TH1F*)hMcTrkPtEff[1][j][k]->Clone(Form("%s_clone",hMcTrkPtEff[1][j][k]->GetName()));
	  htmp->SetMarkerStyle(20+j-1);
	  htmp->SetLineColor(color[j-1]);
	  htmp->SetMarkerColor(color[j-1]);
	  htmp->GetYaxis()->SetRangeUser(0,1);
	  htmp->GetYaxis()->SetTitle("Efficiency");
	  ScaleHistoTitle(htmp,0.06,1,0.05,0.06,0.9,0.05,62);
	  if(j==1) htmp->Draw();
	  else     htmp->Draw("sames");
	  if(k==1)
	    {
	      leg->AddEntry(htmp,Form("AuAu_200%s",gTrgSetupTitle[j]),"PL");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency (%s%%)",cent_Name[k]),0.06);
      t1->Draw();
      TLine *line = GetLine(1.3,0,1.3,0.8,1);
      line->Draw();
    }
  c->cd(1);
  leg->SetHeader("|#eta_{MC}| < 0.5");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/McTpcPtEff_in_Lumi.pdf",run_type));

  // efficiency vs luminosity
  const double pt_cut = 2;
  TH1F *hTpcEffVsLumi[gNTrgSetup][nCentBins];
  for(int j=0; j<gNTrgSetup; j++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  int bin_cut = hMcTrkPtVsZdc[1][j][k]->GetYaxis()->FindFixBin(pt_cut);
	  hTpcEffVsLumi[j][k] = (TH1F*)hMcTrkPtVsZdc[1][j][k]->ProjectionX(Form("TpcEffVsLumi_cent%s%s",cent_Title[k],gTrgSetupTitle[j]),bin_cut,-1);
	  TH1F *htmp = (TH1F*)hMcTrkPtVsZdc[0][j][k]->ProjectionX(Form("htmp_cent%s%s",cent_Title[k],gTrgSetupTitle[j]),bin_cut,-1);
	  hTpcEffVsLumi[j][k]->Divide(htmp);
	}
    }
  TCanvas *c = new TCanvas("TpcEff_vs_Lumi","TpcEff_vs_Lumi",1100,700);
  c->Divide(2,2);
  for(int k=1; k<nCentBins; k++)
    {
      c->cd(k);
      SetPadMargin(gPad,0.13, 0.13, 0.03,0.1);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hTpcEffVsLumi[j][k]->SetMarkerStyle(20+j-1);
	  hTpcEffVsLumi[j][k]->SetLineColor(color[j-1]);
	  hTpcEffVsLumi[j][k]->SetMarkerColor(color[j-1]);
	  hTpcEffVsLumi[j][k]->GetYaxis()->SetRangeUser(0,1);
	  hTpcEffVsLumi[j][k]->GetXaxis()->SetRangeUser(0,120);
	  hTpcEffVsLumi[j][k]->GetXaxis()->SetTitle("ZdcRate (kHz)");
	  hTpcEffVsLumi[j][k]->GetYaxis()->SetTitle("Efficiency");
	  ScaleHistoTitle(hTpcEffVsLumi[j][k],0.06,1,0.05,0.06,0.9,0.05,62);
	  if(j==1) hTpcEffVsLumi[j][k]->Draw();
	  else     hTpcEffVsLumi[j][k]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency (p_{T,mc} > %1.0f GeV/c, %s%%)",pt_cut,cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/McTpcEff_vs_Lumi.pdf",run_type));

  // efficiency vs centrality
  list->Clear();
  TH1F *hTpcEffVsCent[gNTrgSetup];
  for(int j=1; j<gNTrgSetup; j++)
    {
      int bin_cut = hMcTrkPtVsCent[1][j]->GetYaxis()->FindFixBin(pt_cut);
      hTpcEffVsCent[j] = (TH1F*)hMcTrkPtVsCent[1][j]->ProjectionX(Form("hTpcEffVsCent%s",gTrgSetupTitle[j]),bin_cut,-1);
      TH1F *htmp = (TH1F*)hMcTrkPtVsCent[0][j]->ProjectionX(Form("htmp%s",gTrgSetupTitle[j]),bin_cut,-1);
      hTpcEffVsCent[j]->Divide(htmp);
      hTpcEffVsCent[j]->SetMarkerStyle(20+j-1);
      hTpcEffVsCent[j]->SetLineColor(color[j-1]);
      list->Add(hTpcEffVsCent[j]);
      for(int bin=1; bin<=hTpcEffVsCent[j]->GetNbinsX(); bin++)
	{
	  hTpcEffVsCent[j]->GetXaxis()->SetBinLabel(bin,Form("%d-%d%%",80-bin*10,90-bin*10));
	}
      hTpcEffVsCent[j]->GetXaxis()->SetLabelSize(0.05);
    }
  c = drawHistos(list,"TpcEff_vs_Cent",Form("TPC tracking efficiency (p_{T,mc} > %1.0f GeV/c);;Efficiency",pt_cut),false,2.0,3.8,true,0,1,kFALSE,kTRUE,legName_trg,true,"|#eta_{MC}| < 0.5",0.5,0.7,0.2,0.45,kTRUE,0.04,0.04,false,1,false,true); 
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/McTpcEff_vs_Cent.pdf",run_type));

  // other efficiencies
  const char *trkEffTitle[4] = {"Tpc tracking","MTD matching","muon PID","MTD trigger"};
  TCanvas *c = new TCanvas("OtherEff_vs_Cent","OtherEff_vs_Cent",1100,700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.4,0.2,0.7,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int i=1; i<5; i++)
    {
      c->cd(i);
      SetPadMargin(gPad,0.13, 0.13, 0.03,0.1);
      for(int k=1; k<nCentBins; k++)
	{
	  if(i==2)  hMcTrkPtEff[i][0][k]->GetYaxis()->SetRangeUser(0,1);
	  else      hMcTrkPtEff[i][0][k]->GetYaxis()->SetRangeUser(0,1);
	  hMcTrkPtEff[i][0][k]->SetMarkerStyle(19+k);
	  hMcTrkPtEff[i][0][k]->SetMarkerColor(k);
	  hMcTrkPtEff[i][0][k]->SetTitle(";p_{T}^{mc} (GeV/c);Efficiency");
	  ScaleHistoTitle(hMcTrkPtEff[i][0][k],0.06,1,0.05,0.06,0.9,0.05,62);
	  if(k==1) hMcTrkPtEff[i][0][k]->Draw();
	  else     hMcTrkPtEff[i][0][k]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("Single muon efficiency of %s (|#eta_{#mu}| < 0.5)",trkEffTitle[i-1]),0.06);
	  t1->Draw();
	  if(i==1) leg->AddEntry(hMcTrkPtEff[i][j][k],Form("%s%%",cent_Name[k]),"P");
	}
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/McTrkEff_vs_Cent.pdf",run_type));
}

//================================================
void effVsEta(const int savePlot = 1)
{
  TList *list = new TList;

  // eta dependence
  const int nEtaCuts = 3;
  const double eta_cut[nEtaCuts] = {1.5,0.8,0.5};
  THnSparseF *hnTrkInfo[2];
  hnTrkInfo[0] = (THnSparseF*)f->Get("mhMcTrkInfo_di_mu");
  hnTrkInfo[1] = (THnSparseF*)f->Get("mhMcTrkInfoTpc_di_mu");

  TH2F *hMcTrkEtaVsPt = (TH2F*)hnTrkInfo[0]->Projection(1,0);
  hMcTrkEtaVsPt->GetXaxis()->SetRangeUser(0,20);
  c = draw2D(hMcTrkEtaVsPt,"#eta vs p_{T} for input MC tracks");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/McTrkEtaVsPt.pdf",run_type));

  TH1F *hTrkPtEtaCut[2][nEtaCuts];
  TH1F *hTrkPtEffEtaCut[nEtaCuts];
  for(int i=0; i<2; i++)
    {
      hnTrkInfo[i]->Sumw2();
      hnTrkInfo[i]->GetAxis(4)->SetRange(13,16);
      for(int j=0; j<nEtaCuts; j++)
	{
	  hnTrkInfo[i]->GetAxis(1)->SetRangeUser(-1*eta_cut[j]+1e-4, eta_cut[j]-1e-4);
	  hTrkPtEtaCut[i][j] = (TH1F*)hnTrkInfo[i]->Projection(0);
	  hTrkPtEtaCut[i][j]->SetName(Form("hTrkPt_%d_etaCut%d",i,j));
	}
    }
  TString legName[nEtaCuts];
  for(int j=0; j<nEtaCuts; j++)
    {
      hTrkPtEffEtaCut[j] = (TH1F*)hTrkPtEtaCut[1][j]->Clone(Form("hTrkPtEff_etaCut%d",j));
      hTrkPtEffEtaCut[j]->Divide(hTrkPtEtaCut[0][j]);
      list->Add(hTrkPtEffEtaCut[j]);
      legName[j] = Form("|#eta_{MC}| < %1.1f",eta_cut[j]);
    }
  c = drawHistos(list,"TrkPtEff_EtaCut","TPC tracking efficiency for single muons;p_{T} (GeV/c);Efficiency",true,0,20,true,0,1,kFALSE,kTRUE,legName,kTRUE,"Au+Au 0-20%",0.5,0.7,0.25,0.5,kTRUE);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/McTrkPtEff_EtaCuts.pdf",run_type)); 
}

//================================================
void efficiency(TString inName, const bool savePlot = 1, const bool saveHisto = 1)
{
  const int* nCentBins      = nCentBins_npart; 
  const int* centBins_low   = centBins_low_npart;
  const int* centBins_high  = centBins_high_npart;
  const char** cent_Name    = cent_Name_npart;
  const char** cent_Title   = cent_Title_npart;
  const int kNCent          = nCentBins[0];

  TFile *fin = 0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"update");
  else fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"read");
  TList *list = new TList;

  // TPC tracking efficiency
  TH3F *hMcTrkPtEtaPhi[nCentBins][4][3];
  TH1F *hMcTrkPt[nCentBins][4][3];
  TH2F *hMcTrkEtaPhi[nCentBins][4][3];
  TH1F *hMcTrkEta[nCentBins][4][3];
  TH1F *hMcTrkPhi[nCentBins][4][3];
  TH1F *hMcTrkPtInEtaPhi[nCentBins][4][2][2][12];

  TH2F *hMcTrkEtaPhiEff[nCentBins][4][3];
  TH1F *hMcTrkEtaEff[nCentBins][4][3];
  TH1F *hMcTrkPhiEff[nCentBins][4][3];
  TH1F *hMcTrkPtInEtaPhiEff[nCentBins][4][3][2][12];

  double eta_bounds[3] = {-0.5,0,0.5};
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<4; i++)
	    {
	      hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      int low_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(-0.5+1e-4);
	      int high_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(0.5-1e-4);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),low_ybin,high_ybin,0,-1);
	      
	      int low_xbin = hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->FindFixBin(1.5);
	      hMcTrkEta[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionY(Form("McTrkEta%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),low_xbin,-1,0,-1);
	      hMcTrkEtaEff[k][i][j] = (TH1F*)hMcTrkEta[k][i][j]->Clone(Form("%s_Eff",hMcTrkEta[k][i][j]->GetName()));
	      hMcTrkEtaEff[k][i][j]->Divide(hMcTrkEta[k][0][j]);

	      hMcTrkPhi[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionZ(Form("McTrkPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),low_xbin,-1,low_ybin,high_ybin);
	      hMcTrkPhiEff[k][i][j] = (TH1F*)hMcTrkPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPhi[k][i][j]->GetName()));
	      hMcTrkPhiEff[k][i][j]->Divide(hMcTrkPhi[k][0][j]);
	     
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRangeUser(1.5,20);
	      hMcTrkEtaPhi[k][i][j] = (TH2F*)hMcTrkPtEtaPhi[k][i][j]->Project3D("zy");
	      hMcTrkEtaPhi[k][i][j]->SetName(Form("McTrkEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      hMcTrkEtaPhiEff[k][i][j] = (TH2F*)hMcTrkEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkEtaPhi[k][i][j]->GetName()));
	      hMcTrkEtaPhiEff[k][i][j]->Divide(hMcTrkEtaPhi[k][0][j]);
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRange(0,-1);
	    }
	  for(int i=0; i<2; i++)
	    {
	      for(int ieta=0; ieta<2; ieta++)
	      	{
	      	  int low_etabin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(eta_bounds[ieta]+1e-4);
	      	  int high_etabin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(eta_bounds[ieta+1]-1e-4);
	      	  for(int iphi=0; iphi<12; iphi++)
	      	    {
	      	      int low_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12+pi/6*iphi + 1e-4);
	      	      int high_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12+pi/6*(iphi+1) - 1e-4);
	      	      hMcTrkPtInEtaPhi[k][i][j][ieta][iphi] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_eta%d_phi%d_cent%d",det_name[i],charge_name[j],ieta,iphi,k),low_etabin,high_etabin,low_phibin,high_phibin);
	      	      //hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Rebin(2);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi] = (TH1F*)hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Clone(Form("%s_Eff",hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->GetName()));
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->Divide(hMcTrkPtInEtaPhi[k][0][j][ieta][iphi]);
	      	    }
	      	}
	    }
	}
    }

  TString legName[4] = {"MC input","Reconstructed in TPC","Match to MTD","Total efficiency"};
  const char *eff_type[4] = {"MC input","reconstructed in TPC","matched to MTD","fulfill PID/trigger"};
  TString legName_charge[3] = {"All","Positive","Negative"};
  for(int k=0; k<nCentBins; k++)
    {
      // efficiency vs eta vs phi
      TCanvas *c = new TCanvas(Form("EffVsEtaVsPhi_%s",cent_Name[k]),Form("EffVsEtaVsPhi_%s",cent_Name[k]),1100,700);
      c->Divide(2,2);
      for(int i=1; i<4; i++)
	{
	  hMcTrkEtaPhiEff[k][i][0]->GetXaxis()->SetRangeUser(-1,1);
	  hMcTrkEtaPhiEff[k][i][0]->GetZaxis()->SetRangeUser(0,1);
	  c->cd(i);
	  ScaleHistoTitle(hMcTrkEtaPhiEff[k][i][0],0.05,0.9,0.04,0.05,0.8,0.04,62);
	  hMcTrkEtaPhiEff[k][i][0]->SetTitle(";#eta;#varphi");
	  hMcTrkEtaPhiEff[k][i][0]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("Efficiency of MC muons %s (p_{T,mc}>1.5 GeV/c, %s%%)",eff_type[i],cent_Name[k]),0.045);
	  t1->Draw();
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaPhiEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaPhiEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      // efficiency vs eta
      c = new TCanvas(Form("McTrkEffEta_Cent%d",k),Form("McTrkEffEta_Cent%d",k),800,600);
      double xmin = 0.15, xmax = 0.3, ymin = 0.7, ymax = 0.88;
      int counter = 0;
      for(int i=1; i<4; i++)
	{
	  leg = new TLegend(xmin+counter*0.2,ymin,xmax+counter*0.2,ymax);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.035);
	  if(i==1) leg->SetHeader("TPC tracking");
	  if(i==2) leg->SetHeader("MTD matched");
	  if(i==3) leg->SetHeader("Muon PID/trigger");
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkEtaEff[k][i][j]->SetMarkerStyle(19+i);
	      else     
		{
		  hMcTrkEtaEff[k][i][j]->SetMarkerStyle(23+i);
		  if(i==3) hMcTrkEtaEff[k][i][j]->SetMarkerStyle(32);
		}
	      hMcTrkEtaEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->SetLineColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->GetXaxis()->SetRangeUser(-0.6,0.6);
	      hMcTrkEtaEff[k][i][j]->GetYaxis()->SetRangeUser(0,1.4);
	      hMcTrkEtaEff[k][i][j]->SetTitle(";#eta;Efficiency");
	      if(j==0 && i==1) hMcTrkEtaEff[k][i][j]->Draw("P");
	      else hMcTrkEtaEff[k][i][j]->Draw("sames");
	      leg->AddEntry(hMcTrkEtaEff[k][i][j],legName_charge[j].Data(),"P");
	    }
	  leg->Draw();
	  counter++;
	}
      t1 = GetTitleText("Efficiency of MC muons (p_{T,mc} > 1.5 GeV/c)");
      if(nCentBins>1)   t1 = GetTitleText(Form("Efficiency of MC muon track (p_{T,mc} > 1.5 GeV/c, %s%%)",cent_Name[k]));
      t1->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}
 
      // efficiency vs phi
      TCanvas *c = new TCanvas(Form("EffVsPhi_%s",cent_Name[k]),Form("EffVsPhi_%s",cent_Name[k]),1100,700);
      c->Divide(2,2);

      c->cd(1);
      leg = new TLegend(0.15,ymin,0.3,ymax);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
      int counter = 0;
      for(int i=1; i<4; i++)
	{
	  TH1F *htmp = (TH1F*)hMcTrkPhiEff[k][i][0]->Clone(Form("%s_clone",hMcTrkPhiEff[k][i][0]->GetName()));
	  htmp->SetLineColor(color[counter]);
	  htmp->GetYaxis()->SetRangeUser(0,1.4);
	  htmp->SetTitle(";#varphi_{mc};Efficiency");
	  if(i==1) htmp->DrawCopy("HIST");
	  else     htmp->DrawCopy("sames HIST");
	  leg->AddEntry(htmp,legName[counter+1].Data(),"L");
	  counter++;
	}
      leg->Draw();
      t1 = GetTitleText("Efficiency of MC muons (p_{T,mc} > 1.5 GeV/c)",0.045);
      t1->Draw();
      for(int i=1; i<4; i++)
	{
	  c->cd(i+1);
	  leg = new TLegend(xmin,ymin,xmax,ymax);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkPhiEff[k][i][j]->SetMarkerStyle(19+i);
	      else     
		{
		  hMcTrkPhiEff[k][i][j]->SetMarkerStyle(23+i);
		  if(i==4) hMcTrkPhiEff[k][i][j]->SetMarkerStyle(32);
		}
	      hMcTrkPhiEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkPhiEff[k][i][j]->SetLineColor(color[j]);
	      hMcTrkPhiEff[k][i][j]->GetYaxis()->SetRangeUser(0,1.4);
	      hMcTrkPhiEff[k][i][j]->SetTitle(";#varphi_{mc};Efficiency");
	      ScaleHistoTitle(hMcTrkPhiEff[k][i][j],0.05,0.9,0.04,0.05,0.8,0.04,62);
	      if(j==0 && i==1) hMcTrkPhiEff[k][i][j]->Draw("HIST");
	      else hMcTrkPhiEff[k][i][j]->Draw("sames HIST");
	      leg->AddEntry(hMcTrkPhiEff[k][i][j],legName_charge[j].Data(),"L");
	    }
	  leg->Draw();
	  t1 = GetTitleText(Form("Efficiency of MC muons %s (p_{T,mc} > 1.5 GeV/c)",eff_type[i]),0.045);
	  t1->Draw();
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPhiEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPhiEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}


      // efficiency vs pt
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<4; i++)
	    {
	      if(k==nCentBins-1) hMcTrkPt[k][i][j]->Rebin(2);
	      scaleHisto(hMcTrkPt[k][i][j], 1, 1, true, false, false);
	    }
	}
      list->Clear();
      for(int i=0; i<4; i++)
	{
	  list->Add(hMcTrkPt[k][i][1]);
	}
      char *centrality = Form(" (%s%%)",cent_Name[k]);
      if(year==2013) centrality = "";
      c = sysCompare(list,Form("McTrkEffVsPt_Cent%d",k),Form("p_{T} distribution of MC muons %s",centrality),"Efficiency of single muon;p_{T}^{mc} (GeV/c);Efficiency",kTRUE,0,11,kFALSE,0.1,10,kTRUE,0,1.1,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.5,0.7,0.6,0.85,kTRUE);
      for(int i=0; i<4; i++)
	{
	  hMcTrkPt[k][i][2]->SetMarkerStyle(hMcTrkPt[k][i][1]->GetMarkerStyle()+4);
	  hMcTrkPt[k][i][2]->SetMarkerColor(hMcTrkPt[k][i][1]->GetMarkerColor());
	  hMcTrkPt[k][i][2]->SetLineColor(hMcTrkPt[k][i][1]->GetLineColor());
	  c->cd(1);
	  hMcTrkPt[k][i][2]->Draw("samesPE");
	}
      for(int i=1; i<4; i++)
	{
	  TH1F *htmpEff = (TH1F*)hMcTrkPt[k][i][2]->Clone(Form("%s_clone",hMcTrkPt[k][i][2]->GetName()));
	  htmpEff->Divide(hMcTrkPt[k][0][2]);
	  c->cd(2);
	  htmpEff->Draw("samesPE");
	}
      leg = new TLegend(0.5,0.75,0.7,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hMcTrkPt[k][0][1],"positive","P");
      leg->AddEntry(hMcTrkPt[k][0][2],"negative","P");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}
    }

  list->Clear();
  TH1F *hMcTrkPtEff[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hMcTrkPtEff[k] = (TH1F*)hMcTrkPt[k][1][0]->Clone(Form("%s_Eff",hMcTrkPt[k][1][0]->GetName()));
      hMcTrkPtEff[k]->Divide(hMcTrkPt[k][0][0]);
    }

  TString legName_cent[nCentBins-1];
  for(int k=0; k<nCentBins-1; k++)
    {
      TH1F *htmp = (TH1F*)hMcTrkPtEff[k+1]->Clone(Form("%s_clone",hMcTrkPtEff[k+1]->GetName()));
      list->Add(htmp);
      legName_cent[k] = Form("%s%%",cent_Name[k+1]);
    }
  c = drawHistos(list,Form("TrkEff"),Form("%s: TPC tracking efficiency of muon tracks;p_{T,mc} (GeV/c);Efficiency",run_type),kTRUE,0,11,kTRUE,0,1.1,kFALSE,kTRUE,legName_cent,kTRUE,"",0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtEff_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtEff_CentBins.png",run_type,run_cfg_name.Data()));
    }

  // track momentum resolution
  TH2F *hResVsTruePt[nCentBins];
  TH1F *hTrkPtRes[nCentBins];
  TH1F *hTrkPtShift[nCentBins];
  TF1 *funcRes[nCentBins];
  list->Clear();
  for(int k=0; k<nCentBins; k++)
    {
      hResVsTruePt[k] = (TH2F*) fin->Get(Form("pTrkRes_vs_TruePt_%s",cent_Title[k]));
      hResVsTruePt[k]->GetXaxis()->SetRangeUser(0,10);
      if(k<4) hResVsTruePt[k]->RebinX(5);
      else    hResVsTruePt[k]->RebinX(10);
      hTrkPtRes[k] = (TH1F*)hResVsTruePt[k]->ProjectionX(Form("pTrkRes_%s",cent_Title[k]));
      hTrkPtRes[k]->Reset();
      hTrkPtShift[k] = (TH1F*)hResVsTruePt[k]->ProjectionX(Form("pTrkShift_%s",cent_Title[k]));
      hTrkPtShift[k]->Reset();

      TCanvas *c = new TCanvas(Form("FitTrkRes_cent%d",k),Form("TrkRes_cent%d",k),1200,800);
      c->Divide(5,8);
      for(int ibin=1; ibin<=hTrkPtRes[k]->GetNbinsX(); ibin++)
	{
	  TH1F *htmp = (TH1F*)hResVsTruePt[k]->ProjectionY(Form("TrkPt_bin%d_cent%d",ibin,k),ibin,ibin);
	  htmp->GetXaxis()->SetRangeUser(-0.2,0.2);
	  htmp->SetMarkerStyle(20);
	  TF1 *func = new TF1(Form("func_bin%d",ibin),"gaus",-0.15,0.15);
	  htmp->Fit(func,"IR0Q");
	  c->cd(ibin);
	  htmp->Draw("HIST");
	  func->SetLineColor(2);
	  func->Draw("same");
	  hTrkPtRes[k]->SetBinContent(ibin,func->GetParameter(2));
	  hTrkPtRes[k]->SetBinError(ibin,func->GetParError(2));
	  hTrkPtShift[k]->SetBinContent(ibin,func->GetParameter(1));
	  hTrkPtShift[k]->SetBinError(ibin,func->GetParError(1));
	}
    }

  TCanvas *c = new TCanvas("FitTrkPtRes","FitTrkPtRes",1100,700);
  c->Divide(3,2);
  TH1F *hFit = 0x0;
  for(int k=0; k<nCentBins; k++)
    {
      funcRes[k] = new TF1(Form("pTrkResFit_%s",cent_Title[k]),"sqrt([0]^2*x^2+[1]^2)",0.8,20);
      funcRes[k]->SetParNames("a","b");
      funcRes[k]->SetParameter(0.005,0.0);
      hFit = (TH1F*)hTrkPtRes[k]->Clone(Form("Fit_%s",hTrkPtRes[k]->GetName()));
      hFit->Fit(funcRes[k],"IR0");
      hFit->SetMarkerStyle(21);
      hFit->GetYaxis()->SetRangeUser(0,0.1);
      c->cd(k+1);
      hFit->SetTitle(";p_{T,mc} (GeV/c);#sigma(p_{T})/p_{T}");
      hFit->Draw("P");
      TPaveText *t1 = GetTitleText(Form("p_{T} resolution of primary tracks (%s%%)",cent_Name[k]),0.05);
      t1->Draw();
      funcRes[k]->SetLineColor(2);
      funcRes[k]->Draw("same");
    }
  c->cd(6);
  leg = new TLegend(0.2,0.4,0.4,0.7);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetHeader(run_type);
  leg->AddEntry(hFit,"Embedding data","P");
  leg->AddEntry(funcRes[0],"Fit: #sqrt{(a*p_{T})^{2}+b^{2}}","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sFitTrkPtRes_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sFitTrkPtRes_CentBins.png",run_type,run_cfg_name.Data()));
    }

  bool drawLegend = true;
  if(year==2013) drawLegend = false;
  list->Clear();
  for(int k=0; k<nCentBins-1; k++)
    {
      TH1F *htmp = (TH1F*)hTrkPtRes[k+1]->Clone(Form("%s_clone",hTrkPtRes[k+1]->GetName()));
      list->Add(htmp);
    }
  c = drawHistos(list,Form("TrkPtRes"),Form("p_{T} resolution of primary muon tracks;p_{T,true} (GeV/c);#sigma(p_{T})/p_{T}"),kTRUE,0,20,kTRUE,0,0.1,kFALSE,drawLegend,legName_cent,drawLegend,"",0.3,0.5,0.6,0.85,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtRes_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtRes_CentBins.png",run_type,run_cfg_name.Data()));
    }

  list->Clear();
  for(int k=0; k<nCentBins-1; k++)
    {
      TH1F *htmp = (TH1F*)hTrkPtShift[k+1]->Clone(Form("%s_clone",hTrkPtShift[k+1]->GetName()));
      htmp->Scale(100);
      list->Add(htmp);
    }
  c = drawHistos(list,Form("TrkPtShift"),Form("p_{T} shift of primary muon tracks;p_{T,true} (GeV/c);<#Deltap_{T}/p_{T}> (%%)"),kTRUE,0,11,kTRUE,-0.5,0.5,kFALSE,drawLegend,legName_cent,drawLegend,"",0.25,0.45,0.15,0.4,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtShift_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtShift_CentBins.png",run_type,run_cfg_name.Data()));
    }

  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  hTrkPtRes[k]->Write("",TObject::kOverwrite);
	  hTrkPtShift[k]->Write("",TObject::kOverwrite);
	  funcRes[k]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void TrkEff3D(TString inName, TString outPDFName, const bool savePlot = 0, const bool saveHisto = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);

  //----------------------------------------------------------------------------
  // Title page
  TPaveText *t1 = GetPaveText(0.28,0.7,0.5,0.7,0.07,62);
  t1->AddText("Efficiency of single muons from embedding");
  if(year==2013) t1->AddText("in pp 500 GeV from Run13");
  if(year==2014) t1->AddText("in Au+Au 200 GeV from Run14");
  t1->Draw();
  c1->Print(Form("PDF/%s(",outPDFName.Data()));

  //----------------------------------------------------------------------------
  // Get TPDF
  TPDF *pdf = 0;
  TSeqCollection *col = gROOT->GetListOfSpecials();
  for(Int_t i=0; i<col->GetEntries(); i++)
    {
      TObject *obj = (TObject*)col->At(i);
      if( obj->IsA()->InheritsFrom("TPDF") && strcmp(obj->GetName(),outPDFName.Data())==0 )
        pdf = dynamic_cast<TPDF*>obj;
    }
  if(!pdf) 
    {
      printf("No pointer to PDF file is available.\n");
      return;
    }
  pdf->Off();

  t1 = GetPaveText(0.28,0.7,0.5,0.7,0.05,62);
  t1->AddText("TPC tracking efficiency");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  // single tracks
  TFile *fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"read");
  TH3F *hMcTrkPtEtaPhi[nCentBins][nDet][3];
  TH1F *hMcTrkPt[nCentBins][nDet][3];
  TH1F *hMcTrkPtEff[nCentBins][nDet][3];
  TH2F *hMcTrkEtaPhi[nCentBins][nDet][3];
  TH2F *hMcTrkEtaPhiEff[nCentBins][nDet][3];
  TH1F *hMcTrkEta[nCentBins][nDet][3];
  TH1F *hMcTrkEtaEff[nCentBins][nDet][3];
  TH1F *hMcTrkPhi[nCentBins][nDet][3];
  TH1F *hMcTrkPhiEff[nCentBins][nDet][3];
  double eta_bounds[3] = {-0.5,0,0.5};
  int style[3] = {20,24,24};
  TString legName_charge[3] = {"#mu^{+}+#mu^{-}","#mu^{+}","#mu^{-}"};
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<nDet; i++)
	    {
	      hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      int low_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(-0.5+1e-4);
	      int high_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(0.5-1e-4);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_%d",det_name[i],charge_name[j],k),low_ybin,high_ybin,0,-1);
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);
	      hMcTrkPtEff[k][i][j]->SetMarkerStyle(style[j]);
	      hMcTrkPtEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkPtEff[k][i][j]->SetLineColor(color[j]);

	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRangeUser(1.5,20);
	      hMcTrkEtaPhi[k][i][j] = (TH2F*)hMcTrkPtEtaPhi[k][i][j]->Project3D("zy");
	      hMcTrkEtaPhi[k][i][j]->SetName(Form("McTrkEtaPhi%s%s_%d",det_name[i],charge_name[j],k));
	      hMcTrkEtaPhiEff[k][i][j] = (TH2F*)hMcTrkEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkEtaPhi[k][i][j]->GetName()));
	      hMcTrkEtaPhiEff[k][i][j]->Divide(hMcTrkEtaPhi[k][0][j]);
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRange(0,-1);

	      int low_xbin = hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->FindFixBin(1.5);
	      hMcTrkEta[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionY(Form("McTrkEta%s%s_%d",det_name[i],charge_name[j],k),low_xbin,-1,0,-1);
	      hMcTrkPhi[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionZ(Form("McTrkPhi%s%s_%d",det_name[i],charge_name[j],k),low_xbin,-1,low_ybin,high_ybin);
	      hMcTrkEtaEff[k][i][j] = (TH1F*)hMcTrkEta[k][i][j]->Clone(Form("%s_Eff",hMcTrkEta[k][i][j]->GetName()));
	      hMcTrkEtaEff[k][i][j]->Divide(hMcTrkEta[k][0][j]);
	      hMcTrkEtaEff[k][i][j]->SetMarkerStyle(style[j]);
	      hMcTrkEtaEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->SetLineColor(color[j]);

	      hMcTrkPhiEff[k][i][j] = (TH1F*)hMcTrkPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPhi[k][i][j]->GetName()));
	      hMcTrkPhiEff[k][i][j]->Divide(hMcTrkPhi[k][0][j]);
	      hMcTrkPhiEff[k][i][j]->SetMarkerStyle(style[j]);
	      hMcTrkPhiEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkPhiEff[k][i][j]->SetLineColor(color[j]);
	    }
	}
      c = new TCanvas(Form("McTrkPtEff_Cent%d",k),Form("McTrkPtEff_Cent%d",k),800,600);
      leg = new TLegend(0.3,0.65,0.5,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
      for(int j=0; j<3; j++)
	{
	  hMcTrkPtEff[k][1][j]->GetXaxis()->SetRangeUser(0,10);
	  hMcTrkPtEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.4);
	  hMcTrkPtEff[k][1][j]->SetTitle(";p_{T,mc} (GeV/c);Efficiency");
	  if(j==0) hMcTrkPtEff[k][1][j]->Draw("P");
	  else     hMcTrkPtEff[k][1][j]->Draw("Psames");
	  leg->AddEntry(hMcTrkPtEff[k][1][j],legName_charge[j],"P");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (|#eta_{mc}| < 0.5)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPtTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      c = new TCanvas(Form("McTrkEtaEff_Cent%d",k),Form("McTrkEtaEff_Cent%d",k),800,600);
      for(int j=0; j<3; j++)
	{
	  hMcTrkEtaEff[k][1][j]->GetXaxis()->SetRangeUser(-1,1);
	  hMcTrkEtaEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.4);
	  hMcTrkEtaEff[k][1][j]->SetTitle(";#eta_{mc};Efficiency");
	  if(j==0) hMcTrkEtaEff[k][1][j]->Draw("P");
	  else     hMcTrkEtaEff[k][1][j]->Draw("PEsames");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (p_{T,mc} > 1.5 GeV/c)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      c = new TCanvas(Form("McTrkPhiEff_Cent%d",k),Form("McTrkPhiEff_Cent%d",k),800,600);
      for(int j=0; j<3; j++)
      	{
      	  hMcTrkPhiEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.4);
      	  hMcTrkPhiEff[k][1][j]->SetTitle(";#varphi_{mc};Efficiency");
      	  if(j==0) hMcTrkPhiEff[k][1][j]->Draw("HIST");
      	  else     hMcTrkPhiEff[k][1][j]->Draw("sames HIST");
      	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (|#eta_{mc}|<0.5)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPhiTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkPhiTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

     c = new TCanvas(Form("McTrkEtaPhiEff_Cent%d",k),Form("McTrkEtaPhiEff_Cent%d",k),1100,700);
     c->Divide(2,2);
     for(int j=0; j<3; j++)
      	{
	  c->cd(j+1);
     	  hMcTrkEtaPhiEff[k][1][j]->GetXaxis()->SetRangeUser(-1,1);
      	  hMcTrkEtaPhiEff[k][1][j]->GetZaxis()->SetRangeUser(0,1);
      	  hMcTrkEtaPhiEff[k][1][j]->SetTitle(";#eta_{mc};#varphi_{mc}");
	  ScaleHistoTitle(hMcTrkEtaPhiEff[k][1][j],0.06,0.7,0.05,0.06,0.7,0.05,62);
	  hMcTrkEtaPhiEff[k][1][j]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency for %s",legName_charge[j].Data()),0.06);
	  t1->Draw();
      	}
     c->cd(4);
     TPaveText *t1 = GetPaveText(0.3,0.5,0.4,0.6,0.07,62);
     t1->AddText(Form("%s%%",cent_Name[k]));
     t1->Draw();
     PrintCanvasToPDF(c,pdf);
     if(savePlot)
       {
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaPhiTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/%sMcTrkEtaPhiTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
       }
    }

  TList *list = new TList;
  TString legName[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      list->Add(hMcTrkPtEff[k][1][0]);
      legName[k] = Form("%s%%",cent_Name[k]);
    }
  c = drawHistos(list,Form("McTrkEff"),Form("Efficiency of embedded tracks;p_{T} (GeV/c);Efficiency"),kTRUE,0,10,kTRUE,0,1.4,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.3,0.45,0.6,0.85,kTRUE);
  PaintCanvasToPDF(c,pdf);


  // in (phi,eta) bins
  const int index = 1; // 1 - tpc; 2 - mtd; 3 - final
  const int neta = 8;
  const int nphi = 180;
  const int nCanvas = nphi*neta/15;
  TCanvas *cTpcEff[nCanvas];
  for(int i=0; i<nCanvas; i++)
    {
      cTpcEff[i] = new TCanvas(Form("TpcEff_%d",i),Form("TpcEff_%d",i),1200,750);
      cTpcEff[i]->Divide(5,3);
    }
  TH3F *hMcTrkPtEtaPhiEff[nCentBins][3];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  hMcTrkPtEtaPhiEff[k][j] = (TH3F*)hMcTrkPtEtaPhi[k][index][j]->Clone(Form("%s_Eff",hMcTrkPtEtaPhi[k][index][j]->GetName()));
	  int nbinsy = hMcTrkPtEtaPhiEff[k][j]->GetNbinsY();
	  int nbinsz = hMcTrkPtEtaPhiEff[k][j]->GetNbinsZ();
	  hMcTrkPtEtaPhiEff[k][j]->RebinZ(nbinsz/nphi);
	  TH3F *h3tmp = (TH3F*)hMcTrkPtEtaPhi[k][0][j]->Clone(Form("%s_tmp",hMcTrkPtEtaPhi[k][0][j]->GetName()));
	  h3tmp->RebinZ(nbinsz/nphi);
	  hMcTrkPtEtaPhiEff[k][j]->Divide(h3tmp);
	}
    }

  for(int k=0; k<1; k++)
    {
      if(nCentBins>1)
	{
	  t1 = GetPaveText(0.28,0.7,0.5,0.7,0.05,62);
	  if(index==1) t1->AddText(Form("%s%%: TPC tracking",cent_Name[k]));
	  if(index==2) t1->AddText(Form("%s%%: TPC tracking + MTD matching",cent_Name[k]));
	  if(index==3) t1->AddText(Form("%s%%: total efficiency",cent_Name[k]));
	  c1->Clear();
	  t1->Draw();
	  PaintCanvasToPDF(c1,pdf);
	}

      leg = new TLegend(0.6,0.2,0.9,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      for(int ieta=1; ieta<=neta; ieta++)
	{
	  for(int iphi=1; iphi<=nphi; iphi++)
	    {
	      int ind = (ieta-1)*nphi + iphi - 1;
	      cTpcEff[ind/15]->cd(ind%15+1);
	      SetPadMargin(gPad,0.15,0.15,0.05,0.05);
	      for(int j=1; j<3; j++)
		{
		  TH1F *htmp = (TH1F*)hMcTrkPtEtaPhiEff[k][j]->ProjectionX(Form("%s_%d_%d",hMcTrkPtEtaPhiEff[k][j]->GetName(),ieta,iphi),ieta,ieta,iphi,iphi);
		  ScaleHistoTitle(htmp,0.06,1.1,0.05,0.06,1.1,0.05,62);
		  htmp->GetXaxis()->SetRangeUser(0,20);
		  htmp->GetYaxis()->SetRangeUser(0,1.2);
		  htmp->SetMarkerStyle(20+(j-1)*4);
		  if(index==1) htmp->SetTitle(";p_{T,mc} (GeV/c);Tpc Efficiency");
		  if(index==2) htmp->SetTitle(";p_{T,mc} (GeV/c);Tpc+Mtd Efficiency");
		  if(index==3) htmp->SetTitle(";p_{T,mc} (GeV/c);Total Efficiency");
		  if(j==1) htmp->Draw("PE");
		  else     htmp->Draw("PEsames");
		  TPaveText *t1 = GetPaveText(0.25,0.6,0.75,0.9,0.06);
		  t1->AddText(Form("%1.1f < #eta < %1.1f",hMcTrkPtEtaPhiEff[k][j]->GetYaxis()->GetBinLowEdge(ieta),hMcTrkPtEtaPhiEff[k][j]->GetYaxis()->GetBinUpEdge(ieta)));
		  t1->AddText(Form("%1.0f < #varphi < %1.0f",hMcTrkPtEtaPhiEff[k][j]->GetZaxis()->GetBinLowEdge(iphi)/pi*180,hMcTrkPtEtaPhiEff[k][j]->GetZaxis()->GetBinUpEdge(iphi)/pi*180));
		  t1->Draw();

		  if(ieta==1 && iphi==1)
		    {
		      if(j==1) leg->AddEntry(htmp,"#mu^{+}","P");
		      if(j==2) leg->AddEntry(htmp,"#mu^{-}","P");
		    }
		}
	      if(ind%15==0) leg->Draw();
	    }
	}
    }
  
  for(int i=0; i<nCanvas; i++)
    {
      PrintCanvasToPDF(cTpcEff[i],pdf);
    }

  // in (backleg,module)
  t1 = GetPaveText(0.28,0.7,0.5,0.7,0.05,62);
  t1->AddText(Form("%s%%: MTD efficiencies",cent_Name[0]));
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  TCanvas *cMtdEffBlMod[10];
  for(int i=0; i<10; i++)
    {
      cMtdEffBlMod[i] = new TCanvas(Form("MtdEff_%d",i),Form("MtdEff_%d",i),1200,750);
      cMtdEffBlMod[i]->Divide(5,3);
    }
  TH1F *hRcTrkMtdEff[30][5];
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  int bl = i+1; 
	  int mod = j+1;
	  hRcTrkMtdEff[i][j] = (TH1F*)fin->Get(Form("RecoTrkPt_MtdEff_BL%d_Mod%d",bl,mod));
	  cMtdEffBlMod[i/3]->cd((i%3)*5+j+1);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.05);
	  ScaleHistoTitle(hRcTrkMtdEff[i][j],0.06,1.1,0.05,0.06,1.1,0.05,62);
	  hRcTrkMtdEff[i][j]->GetXaxis()->SetRangeUser(0,20);
	  hRcTrkMtdEff[i][j]->GetYaxis()->SetRangeUser(0,1.2);
	  hRcTrkMtdEff[i][j]->SetMarkerStyle(20);
	  hRcTrkMtdEff[i][j]->SetTitle(";p_{T,rec} (GeV/c);Mtd Efficiency");
	  hRcTrkMtdEff[i][j]->Draw();
	  TPaveText *t1 = GetPaveText(0.25,0.6,0.75,0.9,0.06);
	  t1->AddText(Form("BL = %d, Mod = %d",bl,mod));
	  t1->Draw();
	}
    }
  for(int i=0; i<10; i++)
    {
      PrintCanvasToPDF(cMtdEffBlMod[i],pdf);
    }


  TH2F *hMcTrkPtEtaMtd[nCentBins][3][2][3];
  TH2F *hMcTrkPtEtaMtdEff[nCentBins][3][2];
  for(int k=0; k<0; k++)
    {
      TCanvas *cMtdEff = new TCanvas(Form("MtdMthEff_%d",k),Form("MtdMthEff_%d",k),1200,700);
      cMtdEff->Divide(5,4);
      for(int j=0; j<3; j++)
	{
	  int nbin1 = hMcTrkPtEtaPhi[k][1][j]->GetZaxis()->FindFixBin(3.77);
	  int nbin2 = hMcTrkPtEtaPhi[k][1][j]->GetZaxis()->FindFixBin(5.65);
	  double bins[4] = {0,nbin1,nbin2,hMcTrkPtEtaPhi[k][1][j]->GetNbinsZ()};
	  for(int m=0; m<2; m++)
	    {
	      for(int l=0; l<3; l++)
		{
		  hMcTrkPtEtaPhi[k][m+1][j]->GetZaxis()->SetRange(bins[l]+1,bins[l+1]);
		  hMcTrkPtEtaMtd[k][j][m][l] = (TH2F*)hMcTrkPtEtaPhi[k][m+1][j]->Project3D("yx");
		  hMcTrkPtEtaMtd[k][j][m][l]->SetName(Form("%s_%d",hMcTrkPtEtaPhi[k][m+1][j]->GetName(),l));
		}
	    }

	  hMcTrkPtEtaMtdEff[k][j][1] = (TH2F*)hMcTrkPtEtaMtd[k][j][1][0]->Clone(Form("McTrkPtEtaMtdEff%s_5tray_%s",charge_name[j],cent_Title[k]));
	  hMcTrkPtEtaMtdEff[k][j][1]->Add(hMcTrkPtEtaMtd[k][j][1][2]);
	  hMcTrkPtEtaMtdEff[k][j][1]->Rebin2D(5,1);
	  h2tmp = (TH2F*)hMcTrkPtEtaMtd[k][j][0][0]->Clone(Form("McTrkPtEtaMtdEff%s_%s_tmp2",charge_name[j],cent_Title[k]));
	  h2tmp->Add(hMcTrkPtEtaMtd[k][j][0][2]);
	  h2tmp->Rebin2D(5,1);
	  hMcTrkPtEtaMtdEff[k][j][1]->Divide(h2tmp);
	  draw2D(hMcTrkPtEtaMtdEff[k][j][1]);

	  int begin = hMcTrkPtEtaMtdEff[k][j][1]->GetYaxis()->FindFixBin(-0.6);
	  int end = hMcTrkPtEtaMtdEff[k][j][1]->GetYaxis()->FindFixBin(0.6);
	  for(int ieta=begin; ieta<=end; ieta++)
	    {
	      cMtdEff->cd(ieta-begin+1);
	      TH1F *htmp = (TH1F*)hMcTrkPtEtaMtdEff[k][j][1]->ProjectionX(Form("%s_%d",hMcTrkPtEtaMtdEff[k][j][1]->GetName(),ieta),ieta,ieta);
	      htmp->SetMarkerStyle(style[j]+4);
	      htmp->SetMarkerColor(color[j]);
	      htmp->SetLineColor(color[j]);
	      if(j==0) htmp->Draw();
	      else     htmp->Draw("sames");
	    }
	  
	  hMcTrkPtEtaMtdEff[k][j][0] = (TH2F*)hMcTrkPtEtaMtd[k][j][1][1]->Clone(Form("McTrkPtEtaMtdEff%s_3tray_%s",charge_name[j],cent_Title[k]));
	  hMcTrkPtEtaMtdEff[k][j][0]->Rebin2D(5,1);
	  TH2F *h2tmp = (TH2F*)hMcTrkPtEtaMtd[k][j][0][1]->Clone(Form("McTrkPtEtaMtdEff%s_%s_tmp",charge_name[j],cent_Title[k]));
	  h2tmp->Rebin2D(5,1);
	  hMcTrkPtEtaMtdEff[k][j][0]->Divide(h2tmp);
	  draw2D(hMcTrkPtEtaMtdEff[k][j][0]);

	  for(int ieta=begin; ieta<=end; ieta++)
	    {
	      cMtdEff->cd(ieta-begin+1);
	      TH1F *htmp = (TH1F*)hMcTrkPtEtaMtdEff[k][j][0]->ProjectionX(Form("%s_%d",hMcTrkPtEtaMtdEff[k][j][0]->GetName(),ieta),ieta,ieta);
	      htmp->SetMarkerStyle(style[j]);
	      htmp->SetMarkerColor(color[j]);
	      htmp->SetLineColor(color[j]);
	      htmp->Draw("sames");
	    }
	}
      PrintCanvasToPDF(cMtdEff,pdf);
    }

  pdf->On();
  pdf->Close();
}

