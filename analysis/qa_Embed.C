const int part_type = 0;
const char *part_name[4] = {"Jpsi","Upsilon1S","Upsilon2S","Upsilon3S"};
const char *part_text[4] = {"J/#psi","#Upsilon(1S)","#Upsilon(2S)","#Upsilon(3S)"};
const double part_mass[4] = {3.09,9.46, 10.023, 10.34};
const Bool_t iPico = 0;
const int year = 2014;
TString run_cfg_name;

TFile *fdata, *fmc;

//================================================
void qa_Embed()
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
      fmc   = TFile::Open(Form("./output/Run13.pp500.jpsi.%sMC.root",run_config),"read");
      fdata = TFile::Open(Form("./output/Run13.pp500.jpsi.%sPicoData.root",run_config),"read");
    }
  if(year==2014)
    {
      run_type = "Run14_AuAu200";
      fmc   = TFile::Open(Form("./output/Run14.AuAu200.%s.Embed.root",part_name[part_type]),"read");
      fdata = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.root"),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  compWithData();
  //makePDF();
}

//================================================
void compWithData(const int savePlot = 1)
{
  // TPC vz distribution
  TH1F *hTpcVz[2];
  hTpcVz[0] = (TH1F*)fdata->Get("mhTpcVzWithCut_di_mu");
  hTpcVz[1] = (TH1F*)fmc  ->Get("mhDataVtxZ_di_mu");
  for(int i=0; i<2; i++)
    {
      hTpcVz[i]->Sumw2();
      hTpcVz[i]->Rebin(2);
      hTpcVz[i]->Scale(1./hTpcVz[i]->Integral());
      hTpcVz[i]->SetMarkerStyle(20+1);
      hTpcVz[i]->SetLineColor(i+1);
      hTpcVz[i]->SetMarkerColor(i+1);
      hTpcVz[i]->GetXaxis()->SetRangeUser(-110,110);
    }
  c = draw1D(hTpcVz[1],"Run14_AuAu: z distribution of TPC vertex");
  hTpcVz[0]->Draw("samesHIST");
  TLegend *leg = new TLegend(0.7,0.7,0.85,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hTpcVz[0],"Data","L");
  leg->AddEntry(hTpcVz[1],"Embed","P");
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Embed/DataVsEmbed_TpcVzWithCut.pdf",run_type));

  // gRefMult
  TH1F *hgRefMult[2];
  hgRefMult[0] = (TH1F*)fdata->Get("mhgRefMultCorr_di_mu");
  hgRefMult[1] = (TH1F*)fmc  ->Get("mhgRefMultCorr_di_mu");
  for(int i=0; i<2; i++)
    {
      hgRefMult[i]->SetName(Form("hgRefMultCorr_%d",i));
      hgRefMult[i]->Sumw2();
      hgRefMult[i]->Rebin(2);
      hgRefMult[i]->Scale(1./hgRefMult[i]->Integral());
      hgRefMult[i]->SetMarkerStyle(20+1);
      hgRefMult[i]->SetLineColor(i+1);
      hgRefMult[i]->SetMarkerColor(i+1);
      hgRefMult[i]->GetXaxis()->SetRangeUser(0,800);
    }
  c = draw1D(hgRefMult[1],"Run14_AuAu: corrected global reference multiplicity distribution");
  hgRefMult[0]->Draw("samesHIST");
  TLegend *leg = new TLegend(0.7,0.7,0.85,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hgRefMult[0],"Data","L");
  leg->AddEntry(hgRefMult[1],"Embed","P");
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Embed/DataVsEmbed_gReflMultCorr.pdf",run_type));

  // centrality
  TH1F *hCentrality[2];
  hCentrality[0] = (TH1F*)fdata->Get("mhCentrality_di_mu");
  hCentrality[1] = (TH1F*)fmc  ->Get("mhCentrality_di_mu");
  for(int i=0; i<2; i++)
    {
      hCentrality[i]->SetName(Form("hCentrality_%d",i));
      hCentrality[i]->Sumw2();
      hCentrality[i]->Rebin(2);
      hCentrality[i]->Scale(1./hCentrality[i]->Integral());
      hCentrality[i]->SetMarkerStyle(20+1);
      hCentrality[i]->SetLineColor(i+1);
      hCentrality[i]->SetMarkerColor(i+1);
      hCentrality[i]->GetXaxis()->SetRangeUser(0,16);
    }
  c = draw1D(hCentrality[1],"Run14_AuAu: centrality distribution");
  gPad->SetLogy();
  hCentrality[0]->Draw("samesHIST");
  TLegend *leg = new TLegend(0.6,0.25,0.8,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hCentrality[0],"Data","L");
  leg->AddEntry(hCentrality[1],"Embed","P");
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Embed/DataVsEmbed_Centrality.pdf",run_type));
  
}

//================================================
void makePDF(char *outPDFName="")
{
  if(year==2013) outPDFName = "Run13_pp500_EmbedQA_Jpsi.pdf";
  if(year==2014) outPDFName = Form("Run14_AuAu200_EmbedQA_%s.pdf",part_name[part_type]);
  TDatime time;
  Int_t Year  = time.GetYear();
  Int_t month = time.GetMonth();
  Int_t day   = time.GetDay();

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);

  //----------------------------------------------------------------------------
  // Title page
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText(Form("Embedding QA for %s #rightarrow #mu^{+}#mu^{-}",part_text[part_type]));
  if(year==2013) t1->AddText("in pp 500 GeV from Run13");
  if(year==2014) t1->AddText("in Au+Au 200 GeV from Run14");
  t1->Draw();
  t1  = new TPaveText(0.28,0.25,0.7,0.4,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.05);
  t1->AddText(Form("%02d/%02d/%04d",month,day,Year));
  t1->AddText("Rongrong Ma");
  t1->Draw();
  c1->Print(Form("%s(",outPDFName));

  //----------------------------------------------------------------------------
  // Get TPDF
  TPDF *pdf = 0;
  TSeqCollection *col = gROOT->GetListOfSpecials();
  for(Int_t i=0; i<col->GetEntries(); i++)
    {
      TObject *obj = (TObject*)col->At(i);
      if( obj->IsA()->InheritsFrom("TPDF") && strcmp(obj->GetName(),outPDFName)==0 )
        pdf = dynamic_cast<TPDF*>obj;
    }
  if(!pdf) 
    {
      printf("No pointer to PDF file is available.\n");
      return;
    }
  pdf->Off();

  //----------------------------------------------------------------------------
  // analysis cuts
  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("Analysis cuts");
  
  t1  = new TPaveText(0.05,0.1,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);

  t1->AddText("Trigger: di-muon");

  TH1F *h1 = (TH1F*)fmc->Get("hAnalysisCuts");
  Int_t counter = (Int_t)(h1->GetBinContent(3)/1e4);
  Double_t vtxz = h1->GetBinContent(1)/counter;
  if(1)
    {
      t1->AddText(Form("Vertex cut: |vz_{TPC}| < %d cm, |vz_{TPC}-vz_{VPD}| < 3 cm, vr < 2 cm",100));
    }
  else
    {
      t1->AddText(Form("Vertex cut: none"));
    }

  t1->AddText(Form("Track selection:"));
  t1->AddText("    Primary tracks");
  //t1->AddText(Form("    p_{T} > %1.1f GeV/c",h1->GetBinContent(2)/counter));
  t1->AddText(Form("    p_{T} > 1 GeV/c"));
  t1->AddText(Form("    |#eta| < %1.1f",h1->GetBinContent(4)/counter));  
  t1->AddText(Form("    NHitsFit >= %1.0f",h1->GetBinContent(5)/counter));     
  t1->AddText(Form("    NHitsDedx >= %1.0f",h1->GetBinContent(6)/counter));
  t1->AddText(Form("    NHitsFit/NHitsPoss >= %1.2f",h1->GetBinContent(10)/counter));
  if(TMath::Abs(h1->GetBinContent(7)/counter)<10)
    t1->AddText(Form("    global dca <= %1.1f cm",h1->GetBinContent(7)/counter));  
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("PID cuts in data");
  
  t1  = new TPaveText(0.05,0.1,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  t1->AddText("Pion PID");
  if(year==2013) 
    t1->AddText("    0.016 < m^{2} < 0.021 (GeV/c^{2})^{2}");
  t1->AddText("    |n#sigma_{#pi}| < 2");
  t1->AddText("Muon PID");
  t1->AddText(Form("    %1.1f < n#sigma_{#pi} < %1.1f",h1->GetBinContent(8)/counter,h1->GetBinContent(9)/counter));
  t1->AddText(Form("    |#Deltaz| < %1.0f cm",h1->GetBinContent(14)/counter));
  t1->AddText(Form("    Under J/#psi mass peak of [3.0,3.2] GeV/c^{2}"));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  TList *list = new TList;

  TH1F *hNJpsi = (TH1F*)fmc->Get("hNEmbedJpsi_di_mu");
  TCanvas *c = draw1D(hNJpsi,"MC input: # of embedded J/psi per event",kFALSE,kFALSE);
  PaintCanvasToPDF(c,pdf);

  // MC input
  TCanvas *c = new TCanvas("mc_input","mc_input",800,600);
  c->Divide(2,2);
  THnSparseF *hJpsiUS_mc = (THnSparseF*)fmc->Get("hJpsiInfo_MC_di_mu");

  const char *mc_name[3] = {"invariant mass","p_{T}","rapidity"};
  for(int i=0; i<3; i++)
    {
      TH1F *hJpsi_true = (TH1F*)hJpsiUS_mc->Projection(i);
      hJpsi_true->SetName("hJpsi_true");
      hJpsi_true->SetTitle("");
      hJpsi_true->SetMinimum(0);
      if(i==0) 
	{
	  hJpsi_true->SetTitle(";M_{#mu#mu} (GeV/c^{2})");
	  if(part_type==0) hJpsi_true->GetXaxis()->SetRangeUser(0,4);
	  else hJpsi_true->GetXaxis()->SetRangeUser(9,12);
	}
      if(i==1)
	{
	  hJpsi_true->GetXaxis()->SetRangeUser(0,20);
	}
      if(i==2)
	{
	  hJpsi_true->GetXaxis()->SetRangeUser(-1.5,1.5);
	}
      c->cd(i+1);
      title = GetTitleText(Form("MC input: %s distribution",mc_name[i]),0.045);
      hJpsi_true->Draw("HIST");
      title->Draw();
    }
  PrintCanvasToPDF(c,pdf);

  // TPC
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("TPC");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  const int nHistos = 2;
  const char *trkinfo_title[3] = {"recontructed MC track in embedding","pion candidates in data","J/#psi muon in data"};
  TH2F *hTrkInfo[8][3];
  THnSparseF *hnTrk[4];
  hnTrk[0] = (THnSparseF*)fmc->Get("hTrkEtaPhi_MCreco_di_mu");
  if(year==2013) hnTrk[1] = (THnSparseF*)fdata->Get("hTrkEtaPhi_DataPion_di_mu");
  if(year==2014) hnTrk[1] = (THnSparseF*)fmc->Get("hTrkEtaPhi_DataPion_di_mu");
  if(nHistos>2)
    {
      hnTrk[2] = (THnSparseF*)fdata->Get("hTrkEtaPhi_DataMuonUL_di_mu");
      hnTrk[3] = (THnSparseF*)fdata->Get("hTrkEtaPhi_DataMuonLS_di_mu");
      hnTrk[2]->Add(hnTrk[3],-1);
    }

  // === eta vs phi
  c = new TCanvas("track_eta_vs_phi","track_eta_vs_phi",800,600);
  c->Divide(2,2);
  TH2F *hTrkEtaPhi[3];
  for(int j=0; j<nHistos; j++)
    {
      hTrkEtaPhi[j] = (TH2F*)hnTrk[j]->Projection(2,1);
      hTrkEtaPhi[j]->RebinY(4);
      hTrkEtaPhi[j]->SetTitle(";#eta;#varphi");
      c->cd(j+1);
      hTrkEtaPhi[j]->Draw("colz");
      title = GetTitleText(Form("#varphi vs #eta of %s",trkinfo_title[j]),0.045);
      title->Draw();
    }
  PrintCanvasToPDF(c,pdf);

  // === vs pt
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  if(i==0)
	    {
	      hTrkInfo[i][j] = (TH2F*)hnTrk[j]->Projection(1,0);
	      hTrkInfo[i][j]->SetName(Form("hTrkEtaVsPt_%d",j));
	    }
	  else if(i==1)
	    {
	      hnTrk[j]->GetAxis(1)->SetRangeUser(-1,-0.01);
	      hTrkInfo[i][j] = (TH2F*)hnTrk[j]->Projection(2,0);
	      hTrkInfo[i][j]->SetName(Form("hTrkPhiVsPt_NegEta_%d",j));
	      hnTrk[j]->GetAxis(1)->SetRange(0,-1);
	    }
	  else if(i==2)
	    {
	      hnTrk[j]->GetAxis(1)->SetRangeUser(0.01,1);
	      hTrkInfo[i][j] = (TH2F*)hnTrk[j]->Projection(2,0);
	      hTrkInfo[i][j]->SetName(Form("hTrkPhiVsPt_PosEta_%d",j));
	      hnTrk[j]->GetAxis(1)->SetRange(0,-1);
	    }
	  hTrkInfo[i][j]->Sumw2();
	}
    }
  const char *trkinfo_name[5] = {"Dca","NHitsFit","NHitsFrac","NHitsDedx","NSigmaPi"};
  for(int i=0; i<5; i++)
    {
      hTrkInfo[i+3][0] = (TH2F*)fmc->Get(Form("hTrk%s_MCreco_di_mu",trkinfo_name[i]));
      hTrkInfo[i+3][0]->Sumw2();
      if(year==2013) hTrkInfo[i+3][1] = (TH2F*)fdata->Get(Form("hTrk%s_DataPion_di_mu",trkinfo_name[i]));
      if(year==2014) hTrkInfo[i+3][1] = (TH2F*)fmc->Get(Form("hTrk%s_DataPion_di_mu",trkinfo_name[i]));
      hTrkInfo[i+3][1]->Sumw2();
      if(nHistos>2)
	{
	  hTrkInfo[i+3][2] = (TH2F*)fdata->Get(Form("hTrk%s_DataMuonUL_di_mu",trkinfo_name[i]));
	  hTrkInfo[i+3][2]->Sumw2();
	  TH2F *hLS = (TH2F*)fdata->Get(Form("hTrk%s_DataMuonLS_di_mu",trkinfo_name[i]));
	  hTrkInfo[i+3][2]->Add(hLS,-1);
	}
    }

  
  const char *trk_title[8] = {"#eta","#varphi","#varphi","Dca","NHitsFit","NHitsFrac","NHitsDedx","n#sigma_{#pi}"};
  const double pt_cuts[6] = {1.0,1.5,2.0,3.0,5.0,20.0};
  for(int i=0; i<8; i++)
    {
      c = new TCanvas(Form("track_%s_vs_pt",trk_title[i]),Form("track_%s_vs_pt",trk_title[i]),800,600);
      c->Divide(3,2);
      TLegend *leg = new TLegend(0.1,0.6,0.7,0.83);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      if(i==1) leg->SetHeader("-1 < #eta < 0");
      if(i==2) leg->SetHeader("0 < #eta < 1");
      for(int k=0; k<5; k++)
	{
	  c->cd(k+1);
	  for(int j=1; j>-1; j--)
	    {
	      int low_bin = hTrkInfo[i][j]->GetXaxis()->FindFixBin(pt_cuts[k]+0.01);
	      int hi_bin  = hTrkInfo[i][j]->GetXaxis()->FindFixBin(pt_cuts[k+1]-0.01);
	      TH1F *h1 = (TH1F*)hTrkInfo[i][j]->ProjectionY(Form("%s_pro%d",hTrkInfo[i][j]->GetName(),k),low_bin,hi_bin);
	      if(i==1 || i==2) 
		{
		  h1->Rebin(10);
		  h1->SetXTitle("#varphi");
		}
	      h1->Scale(1./h1->Integral());
	      if(i==3) h1->SetMaximum(1.5*h1->GetMaximum());
	      else     h1->SetMaximum(1.2*h1->GetMaximum());
	      h1->SetMinimum(0);
	      h1->SetTitle("");
	      if(j==0)
		{
		  h1->SetMarkerStyle(21);
		  h1->SetMarkerColor(2);
		  h1->SetLineColor(2);
		  h1->SetMarkerSize(0.7);
		  h1->SetLineColor(2);
		  h1->Draw("samesP");
		  if(k==0) leg->AddEntry(h1,"Reco MC #mu tracks","PE");
		}
	      else if(j==1)
		{
		  if(i==3) 
		    {
		      if(year==2013) h1->GetXaxis()->SetRangeUser(0,3);
		      if(year==2014) h1->GetXaxis()->SetRangeUser(0,1);
		    }
		  if(i==5) h1->GetXaxis()->SetRangeUser(0.4,1);
		  if(i==7) h1->GetXaxis()->SetRangeUser(-6,6);
		  h1->SetLineColor(1);
		  h1->Draw("HIST");
		  if(k==0) leg->AddEntry(h1,"#pi candidates in data","L");
		}
	      else if(j==2)
		{
		  if(i==3) 		    
		    {
		      if(year==2013) h1->GetXaxis()->SetRangeUser(0,3);
		      if(year==2014) h1->GetXaxis()->SetRangeUser(0,1);
		    }
		  if(i==5) h1->GetXaxis()->SetRangeUser(0.4,1);
		  if(i==7) h1->GetXaxis()->SetRangeUser(-6,6);
		  h1->SetMarkerStyle(20);
		  h1->SetMarkerColor(4);
		  h1->SetLineColor(4);
		  h1->Draw("P");
		  if(k==0) leg->AddEntry(h1,"#mu candidates in data","PE");
		}
	    }
	  title = GetTitleText(Form("%s distribution (%1.1f < p_{T} < %1.1f)",trk_title[i],pt_cuts[k],pt_cuts[k+1]),0.05);
	  title->Draw();
	}
      c->cd(6);
      leg->Draw();
      PrintCanvasToPDF(c,pdf);
    }
  
  // === tracking efficency
  c = new TCanvas("track_pt_eff","track_pt_eff",800,600);
  c->Divide(2,2);
  const TString legName[3] = {"MC tracks","Matched to TPC","Matched to TPC+MTD"};
  THnSparseF *hnTrkInfo[3];
  hnTrkInfo[0] = (THnSparseF*)fmc->Get("mhMcTrkInfo_di_mu");
  hnTrkInfo[1] = (THnSparseF*)fmc->Get("mhMcTrkInfoTpc_di_mu");
  hnTrkInfo[2] = (THnSparseF*)fmc->Get("mhMcTrkInfoMtd_di_mu");
  TH1F *hMcTrkPt[3];
  c->cd(1);
  TLegend *leg = new TLegend(0.5,0.6,0.7,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(Int_t i=0; i<3; i++)
    {
      hnTrkInfo[i]->GetAxis(1)->SetRangeUser(-0.5+1e-4,0.5-1e-4);
      hMcTrkPt[i] = (TH1F*)hnTrkInfo[i]->Projection(0);
      hMcTrkPt[i]->SetName(Form("hMcTrkPt_%d",i));
      scaleHisto( hMcTrkPt[i], 1, 1, kTRUE);
      hMcTrkPt[i]->GetXaxis()->SetTitleOffset(1.1);
      hMcTrkPt[i]->Sumw2();
      hMcTrkPt[i]->SetMarkerStyle(20);
      hMcTrkPt[i]->SetMarkerColor(color[i]);
      hMcTrkPt[i]->SetLineColor(color[i]);
      if(part_type==0) hMcTrkPt[i]->GetXaxis()->SetRangeUser(0,10);
      hMcTrkPt[i]->SetTitle(";p_{T,true} (GeV/c);dN/dp_{T}");
      leg->AddEntry(hMcTrkPt[i],legName[i].Data(),"P");
      if(i==0) hMcTrkPt[i]->Draw();
      else     hMcTrkPt[i]->Draw("samesP");
      hnTrkInfo[i]->GetAxis(1)->SetRange(0,-1);
    }
  leg->Draw();
  title = GetTitleText("p_{T} distribution of single muons (|#eta_{true}|<0.5)",0.045);
  title->Draw();

  c->cd(2);
  gPad->SetGridy();
  const TString legName2[2] = {"Matched to TPC","Matched to TPC+MTD"};
  TLegend *leg = new TLegend(0.15,0.65,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  TH1F *hMcTrkPtRatio[3];
  for(Int_t i=1; i<3; i++)
    {
      hMcTrkPtRatio[i] = (TH1F*)hMcTrkPt[i]->Clone(Form("%s_ratio",hMcTrkPt[i]->GetName()));
      hMcTrkPtRatio[i]->Divide(hMcTrkPt[0]);
      hMcTrkPtRatio[i]->SetMarkerColor(color[i]);
      hMcTrkPtRatio[i]->SetLineColor(color[i]);
      if(part_type==0) hMcTrkPtRatio[i]->GetXaxis()->SetRangeUser(0,10);
      hMcTrkPtRatio[i]->GetYaxis()->SetRangeUser(0,1.3);
      hMcTrkPtRatio[i]->SetTitle(";p_{T,true} (GeV/c);Tracking efficiency");
      leg->AddEntry(hMcTrkPtRatio[i],legName2[i-1].Data(),"P");
      if(i==1) hMcTrkPtRatio[i]->Draw();
      else     hMcTrkPtRatio[i]->Draw("samesP");
    }
  leg->Draw();
  title = GetTitleText("Tracking efficiency of single muons (|#eta_{true}|<0.5)",0.045);
  title->Draw();

  c->cd(3);
  gPad->SetGridy();
  TH1F *hMcTrkPhi[3], *hMcTrkPhiRatio[3];
  TLegend *leg = new TLegend(0.5,0.65,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("p_{T} > 2 GeV/c");
  for(Int_t i=0; i<3; i++)
    {
      hnTrkInfo[i]->GetAxis(0)->SetRangeUser(2+1e-4,100);
      hnTrkInfo[i]->GetAxis(1)->SetRangeUser(-0.5+1e-4,0.5-1e-4);
      hMcTrkPhi[i] = (TH1F*)hnTrkInfo[i]->Projection(2);
      hMcTrkPhi[i]->SetName(Form("hMcTrkPhi_%d",i));
      hMcTrkPhiRatio[i] = (TH1F*)hMcTrkPhi[i]->Clone(Form("%s_ratio",hMcTrkPhi[i]->GetName()));
      hMcTrkPhiRatio[i]->Divide(hMcTrkPhi[0]);
      hMcTrkPhiRatio[i]->SetMarkerColor(color[i]);
      hMcTrkPhiRatio[i]->SetLineColor(color[i]);
      hMcTrkPhiRatio[i]->GetXaxis()->SetRangeUser(0,10);
      hMcTrkPhiRatio[i]->GetYaxis()->SetRangeUser(0,1.3);
      hMcTrkPhiRatio[i]->SetTitle(";#varphi;Tracking efficiency");
      if(i>0)
	{
	  leg->AddEntry(hMcTrkPhiRatio[i],legName2[i-1].Data(),"L");
	  if(i==1) hMcTrkPhiRatio[i]->Draw("HIST");
	  else     hMcTrkPhiRatio[i]->Draw("samesHIST");
	}
      hnTrkInfo[i]->GetAxis(0)->SetRange(0,-1);
      hnTrkInfo[i]->GetAxis(1)->SetRange(0,-1);
    }
  leg->Draw();
  title = GetTitleText("Tracking efficiency of single muons (|#eta_{true}|<0.5)",0.045);
  title->Draw();

  c->cd(4);
  gPad->SetGridy();
  TH1F *hMcTrkEta[3], *hMcTrkEtaRatio[3];
  TLegend *leg = new TLegend(0.5,0.65,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("p_{T} > 2 GeV/c");
  for(Int_t i=0; i<3; i++)
    {
      hnTrkInfo[i]->GetAxis(0)->SetRangeUser(2+1e-4,100);
      hMcTrkEta[i] = (TH1F*)hnTrkInfo[i]->Projection(1);
      hMcTrkEta[i]->SetName(Form("hMcTrkPhi_%d",i));
      hMcTrkEtaRatio[i] = (TH1F*)hMcTrkEta[i]->Clone(Form("%s_ratio",hMcTrkEta[i]->GetName()));
      hMcTrkEtaRatio[i]->Divide(hMcTrkEta[0]);
      hMcTrkEtaRatio[i]->SetMarkerColor(color[i]);
      hMcTrkEtaRatio[i]->SetLineColor(color[i]);
      hMcTrkEtaRatio[i]->GetYaxis()->SetRangeUser(0,1.1);
      hMcTrkEtaRatio[i]->GetXaxis()->SetRangeUser(-1,1);
      hMcTrkEtaRatio[i]->SetTitle(";#eta;Tracking efficiency");
      if(i>0)
	{
	  leg->AddEntry(hMcTrkEtaRatio[i],legName2[i-1].Data(),"L");
	  if(i==1) hMcTrkEtaRatio[i]->Draw("HIST");
	  else     hMcTrkEtaRatio[i]->Draw("samesHIST");
	}
      hnTrkInfo[i]->GetAxis(0)->SetRange(0,-1);
    }
  leg->Draw();
  title = GetTitleText("Tracking efficiency of single muons",0.045);
  title->Draw();
  PrintCanvasToPDF(c,pdf);

  // === track pt resolution
  c = new TCanvas("track_pt_res","track_pt_res",800,600);
  c->Divide(2,2);
  THnSparseF *hnpTrkRes = (THnSparseF*)fmc->Get("mhpTrkPtRes_di_mu");
  TH2F *hpTrkRes = (TH2F*)hnpTrkRes->Projection(0,2);
  hpTrkRes->Sumw2();
  hpTrkRes->SetName("hpTrkPtRes");
  hpTrkRes->GetXaxis()->SetTitleOffset(1.1);
  hpTrkRes->GetXaxis()->SetRangeUser(0,15);
  c->cd(1);
  gPad->SetLogz();
  hpTrkRes->SetTitle("");
  hpTrkRes->Draw("colz");
  title = GetTitleText("p_{T} resolution of pimary tracks (embedding muons)",0.045);
  title->Draw();

  THnSparseF *hngTrkRes = (THnSparseF*)fmc->Get("mhgTrkPtRes_di_mu");
  TH2F *hgTrkRes = (TH2F*)hngTrkRes->Projection(0,2);
  hgTrkRes->Sumw2();
  hgTrkRes->SetName("hgTrkPtRes");
  hgTrkRes->GetXaxis()->SetTitleOffset(1.1);
  hgTrkRes->GetXaxis()->SetRangeUser(0,15);
  c->cd(2);
  gPad->SetLogz();
  hgTrkRes->SetTitle("");
  hgTrkRes->Draw("colz");
  title = GetTitleText("p_{T} resolution of global tracks (embedding muons)",0.045);
  title->Draw();

  TObjArray pSlices;
  hpTrkRes->FitSlicesY(0, 0, -1, 0, "QNR", &pSlices);
  TH1D *hpRes = (TH1D*)pSlices[2]->Clone("hpRes");
  hpRes->SetYTitle("#sigma_{p_{T}}/p_{T}");
  hpRes->SetMarkerStyle(21); 
  hpRes->SetMarkerColor(2);
  hpRes->SetLineColor(2);
  hpRes->GetXaxis()->SetRangeUser(0,15);
  TObjArray gSlices;
  hgTrkRes->FitSlicesY(0, 0, -1, 0, "QNR", &gSlices);
  TH1D *hgRes = (TH1D*)gSlices[2]->Clone("hgRes");
  hgRes->SetTitle(";p_{T,true} (GeV/c);#sigma_{p_{T}}/p_{T}");
  hgRes->SetMarkerStyle(21);
  hgRes->GetXaxis()->SetRangeUser(0,15);
  hgRes->GetYaxis()->SetRangeUser(0,0.3);
  c->cd(3);
  gPad->SetGridy();
  hgRes->Draw("");
  hpRes->Draw("samesP");
  title = GetTitleText("p_{T} resolution of embedded muon tracks",0.045);
  title->Draw();
  TLegend *leg = new TLegend(0.15,0.6,0.35,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hpRes,"Primary tracks","P");
  leg->AddEntry(hgRes,"Global tracks","P");
  leg->Draw();
  TH1D *hgShift = (TH1D*)gSlices[1]->Clone("hgShift");
  hgShift->SetTitle(";p_{T,true} (GeV/c);<#Deltap_{T}/p_{T}>");
  hgShift->SetMarkerStyle(21);
  hgShift->GetXaxis()->SetRangeUser(0,15);
  hgShift->GetYaxis()->SetRangeUser(-0.01,0.01);
  TH1D *hpShift = (TH1D*)pSlices[1]->Clone("hpShift");
  hpShift->SetMarkerStyle(21); 
  hpShift->SetMarkerColor(2);
  hpShift->SetLineColor(2);
  c->cd(4);
  gPad->SetGridy();
  hgShift->Draw("");
  hpShift->Draw("samesP");
  PrintCanvasToPDF(c,pdf);

  c = new TCanvas("fit_jpsi","fit_jpsi",800,600);
  c->Divide(2,2);
  THnSparseF *hnJpsiTpc = (THnSparseF*)fmc->Get("hJpsiInfo_Tpc_di_mu_w");
  TH2F *hJpsiMassVsPtTpc = (TH2F*)hnJpsiTpc->Projection(0,1);
  hJpsiMassVsPtTpc->SetName("hJpsiMassVsPt_TPC");
  c->cd(1);
  TH1F *hJpsiPtTpc = (TH1F*)hJpsiMassVsPtTpc->ProjectionX("hJpsiPt_TPC");
  gPad->SetLogy();
  hJpsiPtTpc->SetTitle(";p_{T} (GeV/c)");
  hJpsiPtTpc->SetMarkerStyle(21);
  hJpsiPtTpc->Draw();
  title = GetTitleText("p_{T} of J/#psi in TPC with weighting",0.06);
  title->Draw();
  c->cd(2);
  TH1F *hJpsiMassTpc = (TH1F*)hJpsiMassVsPtTpc->ProjectionY("hJpsiMass_TPC");
  hJpsiMassTpc->GetXaxis()->SetRangeUser(part_mass[part_type]-1,part_mass[part_type]+1);
  TF1 *func = new TF1("func","gaus",part_mass[part_type]-0.1,part_mass[part_type]+0.1);
  hJpsiMassTpc->Fit(func,"IR0Q");
  hJpsiMassTpc->SetTitle(";M_{#mu#mu} (GeV/c^{2})");
  hJpsiMassTpc->Draw("HIST");
  title = GetTitleText("Reconstructed J/#psi peak in TPC",0.06);
  title->Draw();
  func->SetLineColor(2);
  func->Draw("sames");
  c->cd(3);
  TH1D *hJpsiRes = (TH1D*)hJpsiPtTpc->Clone("hJpsiRes");
  hJpsiRes->Reset();
  TH1D *hJpsiPeak = (TH1D*)hJpsiPtTpc->Clone("hJpsiPeak");
  hJpsiPeak->Reset();
  for(int bin=1; bin<=hJpsiPtTpc->GetNbinsX(); bin++)
    {
      TH1F *hMass = (TH1F*)hJpsiMassVsPtTpc->ProjectionY(Form("hMass_%d",bin),bin,bin);
      TF1 *func = new TF1(Form("func_%d",bin),"gaus",part_mass[part_type]-0.1,part_mass[part_type]+0.1);
      hMass->Fit(func,"IR0Q");
      hJpsiRes->SetBinContent(bin, func->GetParameter(2));
      hJpsiRes->SetBinError(bin, func->GetParError(2));
      hJpsiPeak->SetBinContent(bin, func->GetParameter(1));
      hJpsiPeak->SetBinError(bin, func->GetParError(1));
    }

  hJpsiRes->SetTitle(";p_{T} (GeV/c);#sigma");
  hJpsiRes->SetMarkerStyle(20); 
  hJpsiRes->SetMarkerColor(2);
  hJpsiRes->SetLineColor(2);
  hJpsiRes->GetXaxis()->SetRangeUser(0,15);
  hJpsiRes->Draw();
  title = GetTitleText("J/#psi resolution in TPC",0.06);
  title->Draw();
  c->cd(4);

  hJpsiPeak->SetTitle(";p_{T} (GeV/c);mass");
  hJpsiPeak->SetMarkerStyle(21); 
  hJpsiPeak->SetMarkerColor(2);
  hJpsiPeak->SetLineColor(2);
  hJpsiPeak->GetXaxis()->SetRangeUser(0,15);
  hJpsiPeak->GetYaxis()->SetRangeUser(3.085,3.1);
  hJpsiPeak->Draw();
  title = GetTitleText("J/#psi peak in TPC",0.06);
  title->Draw();
  PrintCanvasToPDF(c,pdf); 
  

  // MTD
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("MTD");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  c = new TCanvas("MTD_dy","MTD_dy",800,600);
  c->Divide(2,2);
  c->cd(1);
  gPad->SetLogz();
  TH2F *hDyVsRcTrkPt = (TH2F*)fmc->Get("hTrkDyVsPt_MCreco_di_mu");
  hDyVsRcTrkPt->GetXaxis()->SetTitleOffset(1.1);
  hDyVsRcTrkPt->GetXaxis()->SetRangeUser(0,10);
  hDyVsRcTrkPt->SetTitle("");
  hDyVsRcTrkPt->Draw("colz");
  title = GetTitleText("#Deltay vs p_{T} of reconstructed MC tracks (primary)",0.05);
  title->Draw();

  c->cd(2);
  TH1F *hDy = (TH1F*)hDyVsRcTrkPt->ProjectionY("hDeltay");
  hDy->GetXaxis()->SetRangeUser(-50,50);
  hDy->Draw("HIST");
  title = GetTitleText("#Deltay of reconstructed MC tracks (primary)",0.05);
  title->Draw();

  TObjArray dySlices;
  hDyVsRcTrkPt->FitSlicesY(0, 0, -1, 0, "QNR", &dySlices);
  c->cd(3);
  TH1D *hdyRes = (TH1D*)dySlices[2]->Clone("hdyRes");
  hdyRes->SetTitle(";p_{T} (GeV/c);#sigma(#Deltay) [cm]");
  hdyRes->SetMarkerStyle(25); 
  hdyRes->GetXaxis()->SetRangeUser(0,15);
  hdyRes->GetYaxis()->SetRangeUser(0,20);
  hdyRes->Draw("");
  title = GetTitleText("#Deltay resolution vs p_{T}",0.06);
  title->Draw();
  TF1 *func= new TF1("func_dy","sqrt(([0]/x/x)+[1])",1,15);
  func->SetParameters(100,1);
  hdyRes->Fit(func,"IR0Q");
  func->SetLineColor(2);
  func->Draw("sames");
  TLegend *leg = new TLegend(0.4,0.5,0.7,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(func,"Fit: #sqrt{[0]/x^{2}+[1]}","L");
  leg->Draw();
  c->cd(4);
  TH1D *hdyShift = (TH1D*)dySlices[1]->Clone("hdyShift");
  hdyShift->SetTitle(";p_{T} (GeV/c);<#Deltay> [cm]");
  hdyShift->SetMarkerStyle(25); 
  hdyShift->GetXaxis()->SetRangeUser(0,15);
  hdyShift->GetYaxis()->SetRangeUser(-1,1);
  hdyShift->Draw("");
  title = GetTitleText("#Deltay shift vs p_{T}",0.06);
  title->Draw();
  PrintCanvasToPDF(c,pdf);

  c = new TCanvas("MTD_dz","MTD_dz",800,600);
  c->Divide(2,2);
  c->cd(1);
  gPad->SetLogz();
  TH2F *hDzVsRcTrkPt = (TH2F*)fmc->Get("hTrkDzVsPt_MCreco_di_mu");
  hDzVsRcTrkPt->GetXaxis()->SetTitleOffset(1.1);
  hDzVsRcTrkPt->GetXaxis()->SetRangeUser(0,10);
  hDzVsRcTrkPt->SetTitle("");
  hDzVsRcTrkPt->Draw("colz");
  title = GetTitleText("#Deltaz vs p_{T} of reconstructed MC tracks (primary)",0.05);
  title->Draw();

  c->cd(2);
  TH1F *hDz = (TH1F*)hDzVsRcTrkPt->ProjectionY("hDeltaz");
  hDz->GetXaxis()->SetRangeUser(-50,50);
  hDz->Draw("HIST");
  title = GetTitleText("#Deltaz of reconstructed MC tracks (primary)",0.05);
  title->Draw();

  TObjArray dzSlices;
  hDzVsRcTrkPt->FitSlicesY(0, 0, -1, 0, "QNR", &dzSlices);
  c->cd(3);
  TH1D *hdzRes = (TH1D*)dzSlices[2]->Clone("hdyRes");
  hdzRes->SetTitle(";p_{T} (GeV/c);#sigma(#Deltaz) [cm]");
  hdzRes->SetMarkerStyle(25); 
  hdzRes->GetXaxis()->SetRangeUser(0,15);
  hdzRes->GetYaxis()->SetRangeUser(0,20);
  hdzRes->Draw("");
  title = GetTitleText("#Deltaz resolution vs p_{T}",0.06);
  title->Draw();
  TF1 *func= new TF1("func_dz","sqrt(([0]/x/x)+[1])",1,15);
  func->SetParameters(100,1);
  hdzRes->Fit(func,"IR0Q");
  func->SetLineColor(2);
  func->Draw("sames");
  TLegend *leg = new TLegend(0.4,0.5,0.7,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(func,"Fit: #sqrt{[0]/x^{2}+[1]}","L");
  leg->Draw();
  c->cd(4);
  TH1D *hdzShift = (TH1D*)dzSlices[1]->Clone("hdyShift");
  hdzShift->SetTitle(";p_{T} (GeV/c);<#Deltay> [cm]");
  hdzShift->SetMarkerStyle(25); 
  hdzShift->GetXaxis()->SetRangeUser(0,15);
  hdzShift->GetYaxis()->SetRangeUser(-1,1);
  hdzShift->Draw("");
  title = GetTitleText("#Deltaz shift vs p_{T}",0.06);
  title->Draw();
  PrintCanvasToPDF(c,pdf);

  pdf->On();
  pdf->Close();
}
