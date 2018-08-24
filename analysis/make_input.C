//================================================
void make_input()
{
  gStyle->SetOptStat(0);
  //Run13pp500();
  Run14AuAu200();
  //Run14AuAu200_2();
  //Run16AuAu200();
}


//================================================
void Run16AuAu200(const int saveHisto = 1)
{
  run_type = "Run16_AuAu200";

  // MTD response efficiency
  TF1  *funcRespEffCos[30][5];
  TF1  *funcRespEffEmb[30][5];
  TFile *fRespEff = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");
  TFile *fEmb2014 = TFile::Open("Rootfiles/Run14_AuAu200.Input.root","read");
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  int bl = i + 1;
	  int mod = j + 1; 
	  if(bl>9 && bl<23) 
	    {
	      funcRespEffCos[i][j] = (TF1*)fRespEff->Get(Form("Cosmic_FitRespEff_BL%d_Mod%d",bl,mod));
	    }
	  else                  
	    {
	      funcRespEffCos[i][j] = (TF1*)fRespEff->Get(Form("Cosmic_TempRespEff_BL%d_Mod%d",bl,mod));
	    }
	  funcRespEffCos[i][j]->SetName(Form("Cosmic_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	  funcRespEffEmb[i][j] = (TF1*)fEmb2014->Get(Form("Embed_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	}
    }


  // trigger efficiency
  TFile *fTrig =  TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type));
  TF1 *hMuonTrigEff = (TF1*)fTrig->Get(Form("%s_gTacDiffEffFinal_prod",run_type));
  TF1 *funcTrigElecEff = (TF1*)fTrig->Get(Form("%s_TrigElecEff_FitFunc",run_type));

  TFile *fTrigSys =  TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type));
  TH1F *hTrigSys = (TH1F*)fTrigSys->Get(Form("%s_Ups1SEffVsPt_Sys_TacDiffEff",run_type));

  TH1F *hplot = new TH1F("hplot",Form("%s: MTD trigger efficiency;p_{T} (GeV/c);Eff",run_type),100,0,10);
  hplot->GetYaxis()->SetRangeUser(0,1);
  TCanvas *c = new TCanvas("MtdTrigEff","MtdTrigEff",800,600);
  hplot->GetYaxis()->SetRangeUser(0.5,1);
  hplot->DrawCopy();
  hMuonTrigEff->Draw("sames");
  funcTrigElecEff->SetLineColor(4);
  funcTrigElecEff->Draw("sames");
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hMuonTrigEff,"Trigger efficiency","L");
  leg->AddEntry(funcTrigElecEff, "Trigger electronics eff.", "L");
  leg->Draw();


  
  // save
  if(saveHisto)
    {
      TFile *fout =  TFile::Open(Form("Rootfiles/%s.Input.root",run_type),"recreate");
      hMuonTrigEff->Write("MtdTrigEff_FitFunc");
      funcTrigElecEff->Write("TrigElecEff_FitFunc");
      hTrigSys->Write("MtdTrigEff_Ups1SSysVsPt");
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      funcRespEffEmb[i][j]->Write();
	      funcRespEffCos[i][j]->Write();
	    }
	}
    } 
}

//================================================
void Run14AuAu200(const int savePlot = 1, const int saveHisto = 1)
{
  TFile *fData = TFile::Open("output/Run14_AuAu200.jpsi.root","read");
  TFile *fEmb = TFile::Open(Form("./output/Run14_AuAu200.Embed.Jpsi.root"),"read");

  // compare vertex
  TH1F *hVertex[2];
  hVertex[0] = (TH1F*)fEmb->Get("mhDataVtxZ_di_mu");
  //hVertex[0]->Rebin(2);
  hVertex[1] = (TH1F*)fData->Get("mhTpcVzWithCut_di_mu");
  TList *list = new TList;
  const TString legName[2] = {"Embedding","Data"};
  for(int i=0; i<2; i++)
    {
      hVertex[i]->Sumw2();
      hVertex[i]->Rebin(2);
      hVertex[i]->Scale(1./hVertex[i]->Integral());
      list->Add(hVertex[i]);
    }
  c = drawHistos(list,"Compare_TpcVz","Z distribution of TPC vertex",kFALSE,0,800,kFALSE,0,0,kFALSE,kTRUE,legName,kTRUE,"",0.65,0.8,0.7,0.85,kFALSE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/TpcVz_EmbedVsData.pdf",run_type));
  if(gSaveAN)  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_TpcVz_EmbedVsData.pdf"));

  TH1F *hRatio = (TH1F*) hVertex[1]->Clone(Form("VtxWeight"));
  hRatio->Divide(hVertex[0]);
  hRatio->SetTitle(";vz (cm);Ratio=Data/Embed");
  hRatio->SetMarkerStyle(21);
  hRatio->GetXaxis()->SetRangeUser(-100,100);
  c = draw1D(hRatio);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/TpcVz_EmbedOverData.pdf",run_type));
  if(gSaveAN)  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_TpcVz_EmbedOverData.pdf"));

  // TPC acceptance loss
  const int nHistos = 2;
  const char *name[nHistos] = {"McMuon","DataPion"};
  const char *title[nHistos] = {"recontructed MC track in embedding","pion candidates in data"};
  THnSparseF *hnTrk[nHistos];
  hnTrk[0] = (THnSparseF*)fEmb->Get("hTrkEtaPhi_MCreco_di_mu");
  hnTrk[1] = (THnSparseF*)fEmb->Get("hTrkEtaPhi_DataPion_di_mu");

  // === eta vs phi
  TH2F *hTrkEtaPhi[nHistos];
  for(int j=0; j<nHistos; j++)
    {
      hTrkEtaPhi[j] = (TH2F*)hnTrk[j]->Projection(2,1);
      hTrkEtaPhi[j]->RebinY(2);
      hTrkEtaPhi[j]->SetTitle(";#eta;#varphi");
      c = draw2D(hTrkEtaPhi[j], Form("#varphi vs #eta of %s",title[j]), 0.04, false);
      if(savePlot)  
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/TrkPhivsEta_%s.pdf",run_type,name[j]));
    }

  const double eta_cuts[3] = {-1.0,0.2,1.0};
  const double pt_cuts[6] = {1.0,1.5,2.0,3.0,5.0,20.0};
  TH1F *hTrkPhi[2][nHistos][5];
  TCanvas *cPhi[2];
  for(int i=0; i<2; i++)
    {
      cPhi[i] = new TCanvas(Form("track_phi_%d",i),Form("track_phi_%d",i),1000,600);
      cPhi[i]->Divide(3,2);
      TLegend *leg = new TLegend(0.1,0.6,0.7,0.83);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(Form("%1.1f < #eta < %1.1f",eta_cuts[i], eta_cuts[i+1]));
      for(int j=0; j<nHistos; j++)
	{
	  hnTrk[j]->GetAxis(1)->SetRangeUser(eta_cuts[i]+0.01, eta_cuts[i+1]-0.01);
	  for(int k=0; k<5; k++)
	    {
	      hnTrk[j]->GetAxis(0)->SetRangeUser(pt_cuts[k]+0.01, pt_cuts[k+1]-0.01);
	      hTrkPhi[i][j][k] = (TH1F*)hnTrk[j]->Projection(2);
	      hTrkPhi[i][j][k]->SetName(Form("hTrkPhi_%s_%d_%d",name[j],i,k));
	      hTrkPhi[i][j][k]->SetMarkerStyle(21);
	      hTrkPhi[i][j][k]->SetMarkerColor(2-j);
	      hTrkPhi[i][j][k]->SetLineColor(2-j);
	      hTrkPhi[i][j][k]->Rebin(10);
	      hTrkPhi[i][j][k]->Scale(1./hTrkPhi[i][j][k]->Integral());
	      hTrkPhi[i][j][k]->GetYaxis()->SetRangeUser(0,0.04);
	      hTrkPhi[i][j][k]->SetTitle(";#varphi");
	      cPhi[i]->cd(k+1);
	      if(j==0) 
		{
		  hTrkPhi[i][j][k]->Draw("P");
		  t = GetTitleText(Form("#varphi distribution (%1.1f < p_{T} < %1.1f)",pt_cuts[k],pt_cuts[k+1]),0.05);
		  t->Draw();
		}
	      else     hTrkPhi[i][j][k]->Draw("samesHIST");
	      if(k==0)
		{
		  if(j==0) leg->AddEntry(hTrkPhi[i][j][k],"Reco MC #mu tracks","PE");
		  if(j==1) leg->AddEntry(hTrkPhi[i][j][k],"#pi candidates in data","L");
		}
	      hnTrk[j]->GetAxis(0)->SetRange(0,-1);
	    }
	  hnTrk[j]->GetAxis(1)->SetRange(0,-1);
	}
      cPhi[i]->cd(6);
      leg->Draw();
    }

  TH2F *hTpcCorr = new TH2F("hTpcCorr","Correction factor of TPC inefficiency for #eta < 0.2;p_{T} (GeV/c);#varphi",5,pt_cuts,36,0,2*pi);
  TH1F *hTrkPhiCorr[5];
  for(int k=0; k<5; k++)
    {
      hTrkPhiCorr[k] = (TH1F*)hTrkPhi[0][0][k]->Clone(Form("%s_corr",hTrkPhi[0][0][k]->GetName()));
      for(int bin=1; bin<=hTrkPhiCorr[k]->GetNbinsX(); bin++)
	{
	  if(bin<=31 || bin==36) hTpcCorr->SetBinContent(k+1,bin,1);
	  else
	    {
	      double eff = hTrkPhi[0][1][k]->GetBinContent(bin)/hTrkPhi[0][0][k]->GetBinContent(bin);
	      hTrkPhiCorr[k]->SetBinContent(bin,hTrkPhiCorr[k]->GetBinContent(bin)*eff);
	      hTpcCorr->SetBinContent(k+1,bin,eff);
	      if(k==4) hTpcCorr->SetBinContent(k+1,bin,hTpcCorr->GetBinContent(4,bin));
	    }
	}
      cPhi[0]->cd(k+1);
      hTrkPhiCorr[k]->SetMarkerStyle(25);
      hTrkPhiCorr[k]->SetMarkerColor(4);
      hTrkPhiCorr[k]->Draw("samesP");
    }
  TLegend *leg = new TLegend(0.1,0.5,0.7,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hTrkPhiCorr[0], "Reco MC #mu tracks after correction","P");
  cPhi[0]->cd(6);
  leg->Draw();
  c = draw2D(hTpcCorr,"",0.04,false);
  if(savePlot)  
    {
      cPhi[0]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/CompTrkPhi_NegEta.pdf",run_type));
      cPhi[1]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/CompTrkPhi_PosEta.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbTrkEff/TpcAccCorrFactor.pdf",run_type));
    }
  if(gSaveAN)   
    {
      cPhi[0]->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_CompTrkPhi_NegEta.pdf"));
      cPhi[1]->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_CompTrkPhi_PosEta.pdf"));
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_TpcAccCorr.pdf"));
    }


  // input Jpsi shape
  TH1F *hInPutJpsiPt[4];
  TFile *fInJpsi = TFile::Open("Rootfiles/models.root","read");
  for(int i=0; i<4; i++)
    {
      if(i<3) hInPutJpsiPt[i] = (TH1F*)fInJpsi->Get(Form("TBW_JpsiYield_AuAu200_cent%s",cent_Title_pt[i+1]));
      else    hInPutJpsiPt[i] = (TH1F*)fInJpsi->Get(Form("TBW_JpsiYield_AuAu200_cent%s",cent_Title_pt[3]));
      hInPutJpsiPt[i]->SetName(Form("hInputJpsiShape_Cent%d",i));
      hInPutJpsiPt[i]->Scale(1./hInPutJpsiPt[i]->GetBinContent(hInPutJpsiPt[i]->FindFixBin(1)));
      c = draw1D(hInPutJpsiPt[i],"",true,true);
    }

  // MTD response efficiency
  TFile *fMtd =  TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");
  TF1 *funcRespEffEmb[30][5];
  TF1 *funcRespEffCos[30][5];
  for(int i=0; i<30; i++)
    {
      int bl = i+1;
      for(int j=0; j<5; j++)
	{
	  funcRespEffEmb[i][j] = (TF1*)fMtd->Get(Form("Embed_FitRespEff_BL%d_Mod%d",i+1,j+1));
	  funcRespEffEmb[i][j]->SetName(Form("Embed_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	  if(bl>9 && bl<23)
	    {
	      funcRespEffCos[i][j] = (TF1*)fMtd->Get(Form("Cosmic_FitRespEff_BL%d_Mod%d",i+1,j+1));
	    }
	  else
	    {
	      funcRespEffCos[i][j] = (TF1*)fMtd->Get(Form("Cosmic_TempRespEff_BL%d_Mod%d",i+1,j+1));
	    }
	  funcRespEffCos[i][j]->SetName(Form("Cosmic_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	}
    }


  TH1F *hplot = new TH1F("hplot",Form("%s: MTD response efficiency from embedding;p_{T} (GeV/c);Eff",run_type),100,0,10);
  hplot->GetYaxis()->SetRangeUser(0,1);
  TCanvas *cRespEff[6];
  for(int i=0; i<6; i++)
    {
      cRespEff[i] = new TCanvas(Form("cRespEff_BL%d-%d",i*5+1, i*5+5), Form("cRespEff_BL%d-%d",i*5+1, i*5+5), 1100, 700);
      cRespEff[i]->Divide(5,5);
    }
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  cRespEff[i/5]->cd(i%5*5+j+1);
	  hplot->DrawCopy();
	  funcRespEffEmb[i][j]->SetLineColor(1);
	  funcRespEffEmb[i][j]->Draw("sames");
	  funcRespEffCos[i][j]->SetLineColor(2);
	  funcRespEffCos[i][j]->Draw("sames");
	}
    }

  // trigger efficiency
  TFile *fTrig =  TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type));
  TF1 *hMuonTrigEff = (TF1*)fTrig->Get("Run14_AuAu200_Muon_TacDiffEff");

  // trigger electronics efficiency
  TFile *fTrigElec =  TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type));
  TF1 *funcTrigElecEff = (TF1*)fTrigElec->Get("Run14_AuAu200_TrigElecEff_FitFunc");

  // dTOF efficiency
  TFile *fdtof = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type));
  TF1 *funcDtofCut = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_FitFunc");
  funcDtofCut->SetRange(1.3,20);


  TCanvas *c = new TCanvas("MtdTrigEff","MtdTrigEff",800,600);
  hplot->GetYaxis()->SetRangeUser(0.5,1);
  hplot->DrawCopy();
  hMuonTrigEff->Draw("sames");
  funcTrigElecEff->SetLineColor(4);
  funcTrigElecEff->Draw("sames");
  funcDtofCut->SetLineColor(6);
  funcDtofCut->Draw("sames");
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hMuonTrigEff,"Trigger efficiency","L");
  leg->AddEntry(funcTrigElecEff, "Trigger electronics eff.", "L");
  leg->AddEntry(funcDtofCut, "#Deltatof cut eff.", "L");
  leg->Draw();
 
  // save
  if(saveHisto)
    {
      TFile *fout =  TFile::Open("Rootfiles/Run14_AuAu200.Input.root","recreate");
      hRatio->Write();
      for(int i=0; i<4; i++)
	{
	  hInPutJpsiPt[i]->Write();
	}
      hMuonTrigEff->Write("MtdTrigEff_FitFunc");
      funcTrigElecEff->Write("TrigElecEff_FitFunc");
      funcDtofCut->Write("DtofEff0.75_FitFunc");
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      funcRespEffEmb[i][j]->Write();
	      funcRespEffCos[i][j]->Write();
	    }
	}
      hTpcCorr->Write("",TObject::kOverwrite);
    } 
}


//================================================
void Run14AuAu200_2(const int saveHisto = 1)
{
  // MTD response efficiency
  TFile *fMtd =  TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");
  cout << fMtd->GetName() << endl;
  TF1 *funcRespEffCos[30][5];
  for(int i=0; i<30; i++)
    {
      int bl = i+1;
      for(int j=0; j<5; j++)
	{
	  if(bl>9 && bl<23)
	    {
	      funcRespEffCos[i][j] = (TF1*)fMtd->Get(Form("Cosmic_FitRespEff_BL%d_Mod%d",i+1,j+1));
	    }
	  else
	    {
	      funcRespEffCos[i][j] = (TF1*)fMtd->Get(Form("Cosmic_TempRespEff_BL%d_Mod%d",i+1,j+1));
	    }
	  funcRespEffCos[i][j]->SetName(Form("Cosmic_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	}
    }


  
  // save
  if(saveHisto)
    {
      TFile *fout =  TFile::Open("Rootfiles/Run14_AuAu200.Input.root","update");
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      funcRespEffCos[i][j]->Write("",TObject::kOverwrite);
	    }
	}
    } 
}

//================================================
void Run13pp500(const int saveHisto = 1)
{
  const char* run_type = "Run13_pp500";

  // MTD response efficiency
  TFile *fMtd =  TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");
  TF1 *funcRespEffEmb[30][5];
  TF1 *funcRespEffCos[30][5];
  for(int i=0; i<30; i++)
    {
      int bl = i+1;
      for(int j=0; j<5; j++)
	{
	  funcRespEffEmb[i][j] = (TF1*)fMtd->Get(Form("Embed_FitRespEff_BL%d_Mod%d",i+1,j+1));
	  funcRespEffEmb[i][j]->SetName(Form("Embed_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	  funcRespEffCos[i][j] = (TF1*)fMtd->Get(Form("Cosmic_TempRespEff_BL%d_Mod%d",i+1,j+1));
	  funcRespEffCos[i][j]->SetName(Form("Cosmic_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	}
    }


  TH1F *hplot = new TH1F("hplot",Form("%s: MTD response efficiency from embedding;p_{T} (GeV/c);Eff",run_type),100,0,10);
  hplot->GetYaxis()->SetRangeUser(0,1);
  TCanvas *cRespEff[6];
  for(int i=0; i<6; i++)
    {
      cRespEff[i] = new TCanvas(Form("cRespEff_BL%d-%d",i*5+1, i*5+5), Form("cRespEff_BL%d-%d",i*5+1, i*5+5), 1100, 700);
      cRespEff[i]->Divide(5,5);
    }
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  cRespEff[i/5]->cd(i%5*5+j+1);
	  hplot->DrawCopy();
	  funcRespEffEmb[i][j]->SetLineColor(1);
	  funcRespEffEmb[i][j]->Draw("sames");
	  funcRespEffCos[i][j]->SetLineColor(2);
	  funcRespEffCos[i][j]->Draw("sames");
	}
    }

  TFile *fRespSys = TFile::Open("Rootfiles/Run13_pp500.Sys.MtdRespEff.root", "read");
  TH1F *hRespSys = (TH1F*)fRespSys->Get("Jpsi_hMtdRespEffSysAll");
  hRespSys->SetName("MtdRespEff_Sys");

  // trigger efficiency
  TFile *fTrig =  TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type));
  TF1 *hMuonTrigEff = (TF1*)fTrig->Get("Run13_pp500_gTacDiffEffFinalFit_BinCount_prod_Run15_pp200");
  TF1 *funcTrigElecEff = (TF1*)fTrig->Get("Run13_pp500_TrigElecEff_FitFunc");

  TCanvas *c = new TCanvas("MtdTrigEff","MtdTrigEff",800,600);
  hplot->GetYaxis()->SetRangeUser(0.5,1);
  hplot->DrawCopy();
  hMuonTrigEff->Draw("sames");
  funcTrigElecEff->SetLineColor(4);
  funcTrigElecEff->Draw("sames");
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hMuonTrigEff,"Trigger efficiency","L");
  leg->AddEntry(funcTrigElecEff, "Trigger electronics eff.", "L");
  leg->Draw();
 
  // save
  if(saveHisto)
    {
      TFile *fout =  TFile::Open("Rootfiles/Run13_pp500.Input.root","recreate");
      hMuonTrigEff->Write("MtdTrigEff_FitFunc");
      funcTrigElecEff->Write("TrigElecEff_FitFunc");
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      funcRespEffEmb[i][j]->Write();
	      funcRespEffCos[i][j]->Write();
	    }
	}
      hRespSys->Write();
    } 
}
