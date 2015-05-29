const char *run_config = "EmbedQA.";
const Bool_t iPico = 1;
const int year = 2013;
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
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  makePDF();
}


//================================================
void makePDF(char *outPDFName="")
{
  if(year==2013) outPDFName = "EmbedQA_Run13_pp500_Jpsi.pdf";
  TDatime time;
  Int_t year  = time.GetYear();
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
  t1->AddText("Embedding QA for J/#psi #rightarrow #mu^{+}#mu^{-}");
  t1->AddText("in pp 500 GeV from Run13");
  t1->Draw();
  t1  = new TPaveText(0.28,0.25,0.7,0.4,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.05);
  t1->AddText(Form("%02d/%02d/%04d",month,day,year));
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
  if(TMath::Abs(vtxz)<1e3)
    {
      t1->AddText(Form("Vertex cut: |vtx_z| < %1.0f cm",vtxz));
    }
  else
    {
      t1->AddText(Form("Vertex cut: none"));
    }

  t1->AddText(Form("Track selection:"));
  t1->AddText("    Primary tracks");
  t1->AddText(Form("    p_{T} > %1.1f GeV/c",h1->GetBinContent(2)/counter));
  t1->AddText(Form("    |#eta| < %1.1f",h1->GetBinContent(4)/counter));  
  t1->AddText(Form("    NHitsFit >= %1.0f",h1->GetBinContent(5)/counter));     
  t1->AddText(Form("    NHitsDedx >= %1.0f",h1->GetBinContent(6)/counter));
  t1->AddText(Form("    NHitsFit/NHitsPoss > %1.2f",h1->GetBinContent(10)/counter));
  if(TMath::Abs(h1->GetBinContent(7)/counter)<10)
    t1->AddText(Form("    global dca < %1.1f cm",h1->GetBinContent(7)/counter));  
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
  t1->AddText("    0.016 < m^{2} < 0.021 (GeV/c^{2})^{2}");
  t1->AddText("    |n#sigma_{#pi}| < 4");
  t1->AddText("Muon PID");
  t1->AddText(Form("    %1.1f < n#sigma_{#pi} < %1.1f",h1->GetBinContent(8)/counter,h1->GetBinContent(9)/counter));
  t1->AddText(Form("    |#Deltaz| < %1.0f cm",h1->GetBinContent(14)/counter));
  t1->AddText(Form("    Under J/#psi mass peak of [3.0,3.2] GeV/c^{2}"));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  TList *list = new TList;

  // MC input
  TCanvas *c = new TCanvas("mc_input","mc_input",800,600);
  c->Divide(2,2);
  THnSparseF *hJpsiUS_mc = (THnSparseF*)fmc->Get("hJpsiInfo_di_mu");
  hJpsiUS_mc->GetAxis(7)->SetRange(4,4);

  const char *mc_name[4] = {"invariant mass","p_{T}","#eta","#varphi"};
  for(int i=0; i<4; i++)
    {
      TH1F *hJpsi_true = (TH1F*)hJpsiUS_mc->Projection(i);
      hJpsi_true->SetName("hJpsi_true");
      hJpsi_true->SetTitle("");
      hJpsi_true->SetMinimum(0);
      if(i==0) 
	{
	  hJpsi_true->SetTitle(";M_{#mu#mu} (GeV/c^{2})");
	  hJpsi_true->GetXaxis()->SetRangeUser(0,4);
	}
      if(i==1)
	{
	  hJpsi_true->GetXaxis()->SetRangeUser(0,12);
	}
      if(i==2)
	{
	  hJpsi_true->GetXaxis()->SetRangeUser(-1.5,1.5);
	}
      c->cd(i+1);
      title = GetTitleText(Form("MC input: %s distribution",mc_name[i]),0.045);
      hJpsi_true->Draw();
      title->Draw();
    }
  hJpsiUS_mc->GetAxis(7)->SetRange(0,-1);
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

  const char *trkinfo_title[3] = {"recontructed MC track in embedding","pion candidates in data","J/#psi muon in data"};
  TH2F *hTrkInfo[8][3];
  THnSparseF *hnTrk[4];
  hnTrk[0] = (THnSparseF*)fmc->Get("hTrkEtaPhi_MCreco_di_mu");
  hnTrk[1] = (THnSparseF*)fdata->Get("hTrkEtaPhi_DataPion_di_mu");
  hnTrk[2] = (THnSparseF*)fdata->Get("hTrkEtaPhi_DataMuonUL_di_mu");
  hnTrk[3] = (THnSparseF*)fdata->Get("hTrkEtaPhi_DataMuonLS_di_mu");
  hnTrk[2]->Add(hnTrk[3],-1);

  // === eta vs phi
  c = new TCanvas("track_eta_vs_phi","track_eta_vs_phi",800,600);
  c->Divide(2,2);
  TH2F *hTrkEtaPhi[3];
  for(int j=0; j<3; j++)
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
      for(int j=0; j<3; j++)
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
	}
    }
  const char *trkinfo_name[5] = {"Dca","NHitsFit","NHitsFrac","NHitsDedx","NSigmaPi"};
  for(int i=0; i<5; i++)
    {
      hTrkInfo[i+3][0] = (TH2F*)fmc->Get(Form("hTrk%s_MCreco_di_mu",trkinfo_name[i]));
      hTrkInfo[i+3][1] = (TH2F*)fdata->Get(Form("hTrk%s_DataPion_di_mu",trkinfo_name[i]));
      hTrkInfo[i+3][2] = (TH2F*)fdata->Get(Form("hTrk%s_DataMuonUL_di_mu",trkinfo_name[i]));
      TH2F *hLS = (TH2F*)fdata->Get(Form("hTrk%s_DataMuonLS_di_mu",trkinfo_name[i]));
      hTrkInfo[i+3][2]->Add(hLS,-1);
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
	      hTrkInfo[i][j]->Sumw2();
	      TH1F *h1 = (TH1F*)hTrkInfo[i][j]->ProjectionY(Form("%s_pro%d",hTrkInfo[i][j]->GetName(),k),low_bin,hi_bin);
	      if(i==1 || i==2) 
		{
		  h1->Rebin(10);
		  h1->SetXTitle("#varphi");
		}
	      h1->Scale(1./h1->Integral());
	      if(i==3) h1->SetMaximum(2.2*h1->GetMaximum());
	      else     h1->SetMaximum(1.2*h1->GetMaximum());
	      h1->SetMinimum(0);
	      h1->SetTitle("");
	      if(j==0)
		{
		  h1->SetMarkerStyle(21);
		  h1->SetMarkerColor(2);
		  h1->SetMarkerSize(0.7);
		  h1->SetLineColor(2);
		  h1->Draw("samesP");
		  if(k==0) leg->AddEntry(h1,"Reconstructed MC tracks","PE");
		}
	      else if(j==1)
		{
		  if(i==3) h1->GetXaxis()->SetRangeUser(0,3);
		  if(i==5) h1->GetXaxis()->SetRangeUser(0.4,1);
		  if(i==7) h1->GetXaxis()->SetRangeUser(-6,6);
		  h1->SetLineColor(1);
		  h1->Draw("HIST");
		  if(k==0) leg->AddEntry(h1,"#pi candidates in data","L");
		}
	      else if(j==2)
		{
		  h1->SetMarkerStyle(20);
		  h1->SetMarkerColor(4);
		  h1->SetLineColor(4);
		  h1->Draw("sames P");
		}
	    }
	  title = GetTitleText(Form("%s distribution (%1.1f < p_{T} < %1.1f)",trk_title[i],pt_cuts[k],pt_cuts[k+1]),0.05);
	  title->Draw();
	}
      c->cd(6);
      leg->Draw();
      PrintCanvasToPDF(c,pdf);
    }

  // === track pt resolution
  c = new TCanvas("track_pt_res","track_pt_res",800,600);
  c->Divide(2,2);
  TH2F *hpTrkRes = (TH2F*)fmc->Get("hpTrkPtRes_di_mu");
  hpTrkRes->GetXaxis()->SetTitleOffset(1.1);
  hpTrkRes->GetXaxis()->SetRangeUser(0,10);
  c->cd(1);
  hpTrkRes->SetTitle("");
  hpTrkRes->Draw("colz");
  title = GetTitleText("p_{T} resolution of pimary tracks (embedding muons)",0.045);
  title->Draw();

  TH2F *hgTrkRes = (TH2F*)fmc->Get("hgTrkPtRes_di_mu");
  hgTrkRes->GetXaxis()->SetTitleOffset(1.1);
  hgTrkRes->GetXaxis()->SetRangeUser(0,10);
  c->cd(2);
  hgTrkRes->SetTitle("");
  hgTrkRes->Draw("colz");
  title = GetTitleText("p_{T} resolution of global tracks (embedding muons)",0.045);
  title->Draw();

  TObjArray pSlices;
  hpTrkRes->FitSlicesY(0, 0, -1, 0, "QNR", &pSlices);
  TH1D *hpRes = (TH1D*)pSlices[2];
  hpRes->SetYTitle("#sigma_{p_{T}}/p_{T}");
  hpRes->SetMarkerStyle(21); 
  hpRes->SetMarkerColor(2);
  hpRes->SetLineColor(2);
  hpRes->GetXaxis()->SetRangeUser(0,10);
  TObjArray gSlices;
  hgTrkRes->FitSlicesY(0, 0, -1, 0, "QNR", &gSlices);
  TH1D *hgRes = (TH1D*)gSlices[2];
  hgRes->SetTitle(";p_{T,true} (GeV/c);#sigma_{p_{T}}/p_{T}");
  hgRes->SetMarkerStyle(21);
  hgRes->GetXaxis()->SetRangeUser(0,10);
  hgRes->GetYaxis()->SetRangeUser(0,0.3);
  c->cd(3);
  gPad->SetGridy();
  hgRes->Draw("P");
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
  PrintCanvasToPDF(c,pdf);

  c = new TCanvas("fit_jpsi","fit_jpsi",800,600);
  hJpsiUS_mc->GetAxis(7)->SetRange(5,5);
  TH1F *hReco = (TH1F*)hJpsiUS_mc->Projection(0);
  hReco->SetName("hRecoSignal_TPC");
  hReco->GetXaxis()->SetRangeUser(2.4,4);
  TF1 *func = new TF1("func","gaus",3,3.2);
  hReco->Fit(func,"IR0");
  c = draw1D(hReco,Form("Reconstructed J/#psi peak using primary tracks in TPC;M_{#mu#mu} (GeV/c^{2})"),kFALSE,kFALSE);
  func->SetLineColor(2);
  func->Draw("sames");
  PrintCanvasToPDF(c,pdf);

  // === tracking efficency
  c = new TCanvas("track_pt_eff","track_pt_eff",800,600);
  c->Divide(2,2);
  const TString legName[3] = {"MC tracks","Matched to TPC","Matched to TPC+MTD"};
  TH2F *hTrkPhiPt[3];
  hTrkPhiPt[0] = (TH2F*)fmc->Get("hMcTrkPhiVsPt_di_mu");
  hTrkPhiPt[1] = (TH2F*)fmc->Get("hMcTrkPhiVsPtTpc_di_mu");
  hTrkPhiPt[2] = (TH2F*)fmc->Get("hMcTrkPhiVsPtMtd_di_mu");
  TH1F *hMcTrkPt[3];
  c->cd(1);
  TLegend *leg = new TLegend(0.5,0.6,0.7,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(Int_t i=0; i<3; i++)
    {
      hMcTrkPt[i] = (TH1F*)hTrkPhiPt[i]->ProjectionX(Form("hMcTrkPt_%d",i));
      scaleHisto( hMcTrkPt[i], 1, 1, kTRUE);
      hMcTrkPt[i]->GetXaxis()->SetTitleOffset(1.1);
      hMcTrkPt[i]->Sumw2();
      hMcTrkPt[i]->SetMarkerStyle(20);
      hMcTrkPt[i]->SetMarkerColor(color[i]);
      hMcTrkPt[i]->SetLineColor(color[i]);
      hMcTrkPt[i]->GetXaxis()->SetRangeUser(0,10);
      hMcTrkPt[i]->SetTitle(";p_{T,true} (GeV/c);dN/dp_{T}");
      leg->AddEntry(hMcTrkPt[i],legName[i].Data(),"P");
      if(i==0) hMcTrkPt[i]->Draw();
      else     hMcTrkPt[i]->Draw("samesP");
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
      hMcTrkPtRatio[i]->GetXaxis()->SetRangeUser(0,10);
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
      hTrkPhiPt[i]->GetXaxis()->SetRangeUser(2,100);
      hMcTrkPhi[i] = (TH1F*)hTrkPhiPt[i]->ProjectionY(Form("hMcTrkPhi_%d",i));
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
    }
  leg->Draw();
  title = GetTitleText("Tracking efficiency of single muons (|#eta_{true}|<0.5)",0.045);
  title->Draw();

  c->cd(4);
  gPad->SetGridy();
  TH2F *hTrkEtaPt[3];
  hTrkEtaPt[0] = (TH2F*)fmc->Get("hMcTrkEtaVsPt_di_mu");
  hTrkEtaPt[1] = (TH2F*)fmc->Get("hMcTrkEtaVsPtTpc_di_mu");
  hTrkEtaPt[2] = (TH2F*)fmc->Get("hMcTrkEtaVsPtMtd_di_mu");
  TH1F *hMcTrkEta[3], *hMcTrkEtaRatio[3];
  TLegend *leg = new TLegend(0.5,0.65,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("p_{T} > 2 GeV/c");
  for(Int_t i=0; i<3; i++)
    {
      hTrkEtaPt[i]->GetXaxis()->SetRangeUser(2,100);
      hMcTrkEta[i] = (TH1F*)hTrkEtaPt[i]->ProjectionY(Form("hMcTrkPhi_%d",i));
      hMcTrkEtaRatio[i] = (TH1F*)hMcTrkEta[i]->Clone(Form("%s_ratio",hMcTrkEta[i]->GetName()));
      hMcTrkEtaRatio[i]->Divide(hMcTrkEta[0]);
      hMcTrkEtaRatio[i]->SetMarkerColor(color[i]);
      hMcTrkEtaRatio[i]->SetLineColor(color[i]);
      hMcTrkEtaRatio[i]->GetYaxis()->SetRangeUser(0,1.1);
      hMcTrkEtaRatio[i]->SetTitle(";#eta;Tracking efficiency");
      if(i>0)
	{
	  leg->AddEntry(hMcTrkEtaRatio[i],legName2[i-1].Data(),"L");
	  if(i==1) hMcTrkEtaRatio[i]->Draw("HIST");
	  else     hMcTrkEtaRatio[i]->Draw("samesHIST");
	}
    }
  leg->Draw();
  title = GetTitleText("Tracking efficiency of single muons",0.045);
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

  TH2F *hDyVsRcTrkPt = (TH2F*)fmc->Get("hDyVsRcTrkPt_MCreco_true_di_mu");
  hDyVsRcTrkPt->GetXaxis()->SetTitleOffset(1.1);
  hDyVsRcTrkPt->GetXaxis()->SetRangeUser(0,10);
  c = draw2D(hDyVsRcTrkPt,"#Deltay vs p_{T} of reconstructed MC tracks (primary)");
  PaintCanvasToPDF(c,pdf);

  TH2F *hDzVsRcTrkPt = (TH2F*)fmc->Get("hDzVsRcTrkPt_MCreco_true_di_mu");
  hDzVsRcTrkPt->GetXaxis()->SetTitleOffset(1.1);
  hDzVsRcTrkPt->GetXaxis()->SetRangeUser(0,10);
  hDzVsRcTrkPt->GetYaxis()->SetRangeUser(-100,100);
  c = draw2D(hDzVsRcTrkPt,"#Deltaz vs p_{T} of reconstructed MC tracks (primary)");
  PaintCanvasToPDF(c,pdf);

  TH2F *hMcTof = (TH2F*)fmc->Get("hMcTof_di_mu");
  hMcTof->GetXaxis()->SetTitleOffset(1.1);
  hMcTof->GetXaxis()->SetRangeUser(0,10);
  c = draw2D(hMcTof,"Time-of-flight for MTD hits calculated by GEANT");
  PaintCanvasToPDF(c,pdf);

  TH2F *hMcDeltaTof = (TH2F*)fmc->Get("hMcDeltaTof_di_mu");
  hMcDeltaTof->GetXaxis()->SetTitleOffset(1.1);
  hMcDeltaTof->GetXaxis()->SetRangeUser(0,10);
  hMcDeltaTof->GetYaxis()->SetRangeUser(-5,5);
  c = draw2D(hMcDeltaTof,"#Deltatof vs p_{T} of reconstructed MC tracks (primary)");
  PaintCanvasToPDF(c,pdf);

  pdf->On();
  pdf->Close();
}
