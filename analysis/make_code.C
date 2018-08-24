//================================================
void make_code()
{
  gStyle->SetOptStat(0);
  //data();
  //lumi();
  //embed();

  //JpsiWidth();
  //JpsiPtShape();
  //MtdAcc();

  publication();
  //efficiency();
}


//================================================
void efficiency(const int saveHisto = 1)
{
  // Online timing cut efficiency
  TFile *fTrigEff = TFile::Open("Rootfiles/Run14_AuAu200.Sys.MtdTrigEff.root","read");
  TF1 *fucnTrig = (TF1*)fTrigEff->Get("Run14_AuAu200_Muon_TacDiffEff");
  TF1 *fucnTrigUp = (TF1*)fTrigEff->Get("Run14_AuAu200_Muon_TacDiffEff_Sysup");
  TF1 *fucnTrigDown = (TF1*)fTrigEff->Get("Run14_AuAu200_Muon_TacDiffEff_Sysdown");

  // PID efficiency
  TFile *fpid = TFile::Open("Rootfiles/Run14_AuAu200.EmbTrkEff.root","read");
  TH1F *hTrkPtMtdMth = (TH1F*)fpid->Get("McTrkPt_MtdMth_cent0080");
  TH1F *hTrkPtPid = (TH1F*)fpid->Get("McTrkPt_MuonPid_cent0080");

  // dTof efficiency
  TFile *fdtof = TFile::Open("Rootfiles/Run14_AuAu200.DtofEff.root","read");
  TF1 *fucnDtof = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_FitFunc");
  TF1 *fucnDtofUp = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_Sysup");
  TF1 *fucnDtofDown = (TF1*)fdtof->Get("TagAndProbe_Muon_Dtof0.75Eff_Sysdown");

  // nsigmaPi, dy and dz efficiency
  TFile *femb = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root","read");
  TH1F *hTrkPt[4];
  hTrkPt[0] = (TH1F*)femb->Get("mhTrkPtDis_Tpc_di_mu");
  hTrkPt[1] = (TH1F*)femb->Get("mhTrkPtDis_NsigmaPi_di_mu");
  hTrkPt[2] = (TH1F*)femb->Get("mhTrkPtDis_MtdMth_di_mu");
  hTrkPt[3] = (TH1F*)femb->Get("mhTrkPtDis_DyDz_di_mu");

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/JpsiAna.input.root"),"update");
      fucnTrig->Write("",TObject::kOverwrite);
      fucnTrigUp->Write("",TObject::kOverwrite);
      fucnTrigDown->Write("",TObject::kOverwrite);
      hTrkPtMtdMth->Write("",TObject::kOverwrite);
      hTrkPtPid->Write("",TObject::kOverwrite);
      fucnDtof->Write("",TObject::kOverwrite);
      fucnDtofUp->Write("",TObject::kOverwrite);
      fucnDtofDown->Write("",TObject::kOverwrite);
      fout->Close();

      TFile *fout2 = TFile::Open(Form("paper_code/Rootfiles/Run14_AuAu200.embed.jpsi.root"),"update");
      for(int i=0; i<4; i++)
	{
	  hTrkPt[i]->Write("",TObject::kOverwrite);
	}
      fout2->Close();
    }
}

//================================================
void publication(const int saveHisto = 1)
{
  TFile *frun12 = TFile::Open("Rootfiles/jpsi_xsec_pp200_run12.root","read");
  TGraphAsymmErrors *gRun12Sys = (TGraphAsymmErrors*)frun12->Get("gJpsiXsecCombSys"); 
  gRun12Sys->SetName("STAR_2012_xsec_sys");
  TGraphAsymmErrors *gRun12 = (TGraphAsymmErrors*)frun12->Get("gJpsiXsecCombAsy");
  gRun12->SetName("STAR_2012_xsec");

  TFile *fjpsi = TFile::Open(Form("Rootfiles/2016sQM/jpsi_xsec_pp200_run12.root"),"read");
  TGraphAsymmErrors *gPhenixSys = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix_Systematics");
  gPhenixSys->SetName("PHENIX_xsec_sys");
  TGraphAsymmErrors *gPhenix = (TGraphAsymmErrors*)fjpsi->Get("gYieldVsPt_pp_Phenix");
  gPhenix->SetName("PHENIX_xsec");

  TFile *fpub = TFile::Open(Form("Rootfiles/Paper/Publication.Jpsi.200GeV.root"),"read");
  const char* star_cent[4] = {"0020","2040","4060","0060"};
  TGraphAsymmErrors *gRaaLowPt[4];
  TGraphAsymmErrors *gRaaLowPtSys[4];
  TGraphAsymmErrors *gRaaHighPt[4];
  TGraphAsymmErrors *gRaaHighPtSys[4];
  for(int k=0; k<4; k++)
    {
      gRaaLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",star_cent[k]));
      gRaaLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",star_cent[k]));
      gRaaHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",star_cent[k]));
      gRaaHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",star_cent[k]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/JpsiAna.input.root",run_type),"update");
      gRun12->Write("",TObject::kOverwrite);
      gRun12Sys->Write("",TObject::kOverwrite);

      gPhenix->Write("",TObject::kOverwrite);
      gPhenixSys->Write("",TObject::kOverwrite);

      for(int k=0; k<4; k++)
	{
	  gRaaLowPt[k]->Write("",TObject::kOverwrite);
	  gRaaLowPtSys[k]->Write("",TObject::kOverwrite);
	  gRaaHighPt[k]->Write("",TObject::kOverwrite);
	  gRaaHighPtSys[k]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void embed(const int saveHisto = 1)
{
  TFile *fin = TFile::Open("./output/Run14_AuAu200.Embed.Jpsi.root","read");
  THnSparseF *hMcTrkMc = (THnSparseF*)fin->Get("mhMcTrkPtEff_MC_di_mu");
  THnSparseF *hMcTrkTpc = (THnSparseF*)fin->Get("mhMcTrkPtEff_Tpc_di_mu");

  THnSparseF *hJpsiMatch = (THnSparseF*)fin->Get("mhJpsiMatch_di_mu");
  THnSparseF *hJpsiMc = (THnSparseF*)fin->Get("hJpsiInfo_MC_di_mu");
  THnSparseF *hJpsiMtd = (THnSparseF*)fin->Get("hJpsiInfo_TrigUnit_di_mu");

 if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/%s.embed.jpsi.root",run_type),"update");
      hMcTrkMc->Write("",TObject::kOverwrite);
      hMcTrkTpc->Write("",TObject::kOverwrite);
      hJpsiMatch->Write("",TObject::kOverwrite);
      hJpsiMc->Write("",TObject::kOverwrite);
      hJpsiMtd->Write("hJpsiInfo_Mtd_di_mu",TObject::kOverwrite);
      fout->Close();
    }
}

//================================================
void lumi(const int saveHisto = 1)
{
  TFile *frf = TFile::Open("./output/Run14_AuAu200.RejectFactor.root","read");
  TH1F *hEvtAll = (TH1F*)frf->Get("hEvtAll");
  TH1F *hEvtAcc = (TH1F*)frf->Get("hEvtAcc");

  const char* lumiName[2] = {"prod_low", "prod_high"};
  TFile *fMB[2];
  THnSparseF *hn[2][2];
  for(int i=0; i<2; i++)
    {
      fMB[i] = TFile::Open(Form("./output/Run14_AuAu200.MB.VtxEff.%s.root",lumiName[i]),"read");
      hn[i][0] = (THnSparseF*)fMB[i]->Get("mhMbEvtEff");
      hn[i][1] = (THnSparseF*)fMB[i]->Get("mhMbEvtEffWeight");
      for(int j=0; j<2; j++)
	{
	  hn[i][j]->SetName(Form("%s_%s",hn[i][j]->GetName(),lumiName[i]));
	}
    }

 if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/%s.data.jpsi.root",run_type),"update");
      hEvtAll->Write("RF_hEvtAll",TObject::kOverwrite);
      hEvtAcc->Write("RF_hEvtAcc",TObject::kOverwrite);

      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hn[i][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void MtdAcc(const int saveHisto = 1)
{
  const int nRunRange = 7;
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type),"read");
  TH1F *hAccLoss[nRunRange];
  for(int i=0; i<nRunRange; i++)
    {
      hAccLoss[i] = (TH1F*)fAcc->Get(Form("hAccepLoss_RunRange%d",i+1));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/JpsiAna.input.root",run_type),"update");
      for(int i=0; i<nRunRange; i++)
	{
	  hAccLoss[i]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void JpsiPtShape(const int saveHisto = 1)
{
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");  
  TH1F *hInputJpsi[4];
  hInputJpsi[0] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  hInputJpsi[1] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0020");
  hInputJpsi[2] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent2040");
  hInputJpsi[3] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent4060");

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/JpsiAna.input.root",run_type),"update");
      for(int i=0; i<4; i++)
	{
	  hInputJpsi[i]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void JpsiWidth(const int saveHisto = 1)
{
  // get the Jpsi width from embedding
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"read");
  TH1F *hEmbJpsiWidthVsPt[nCentBins_pt];
  for(int i=0; i<nCentBins_pt; i++)
    {
      hEmbJpsiWidthVsPt[i] = (TH1F*)fscan->Get(Form("SmearEmb_JpsiWidth_cent%s_def",cent_Title_pt[i]));
    }

  TH1F *hEmbJpsiWidthVsCent[nPtBins_npart];
  for(int i=0; i<nPtBins_npart; i++)
    {
      hEmbJpsiWidthVsCent[i] = (TH1F*)fscan->Get(Form("SmearEmb_JpsiWidthIntegr_Pt%s_def",pt_Name_npart[i]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/JpsiAna.input.root",run_type),"update");
      for(int i=0; i<nCentBins_pt; i++)
	{
	  hEmbJpsiWidthVsPt[i]->Write("",TObject::kOverwrite);
	}
      for(int i=0; i<nPtBins_npart; i++)
	{
	  hEmbJpsiWidthVsCent[i]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void data(const int saveHisto = 1)
{
  TFile *fin = TFile::Open(Form("./output/%s.jpsi.%sroot",run_type,run_config),"read");
  TH1F *hStat = (TH1F*)fin->Get("hEventStat");
  TH1F *hEvtRun = (TH1F*)fin->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fin->Get("mhEvtRunAcc_di_mu");
  TH1F *hEvtRunAccW = (TH1F*)fin->Get("mhEvtRunAccWeight_di_mu");

  // same event
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[2][3] = {0x0};
  char name[512];
  for(int w=0; w<2; w++) // event weights
    { 
      for(Int_t j=0; j<3; j++) // pair type
	{ 
	  if(w==0) sprintf(name,"m%s_%s",hName[j],trigName[kTrigType]);
	  else     sprintf(name,"m%sWeight_%s",hName[j],trigName[kTrigType]);
	  hnInvMass[w][j] = (THnSparseF*)fin->Get(name);
	}
    }

  // mixed event
  TFile *fmix = TFile::Open(Form("Output/%s.Mix.%spt%1.1f.pt%1.1f.root",run_type,run_config,pt1_cut,pt2_cut),"read");
  printf("INFO: using Shuai's mixed events: %s\n",fmix->GetName());
  TH3D *hMixMmumuvsPtCen[3];
  hMixMmumuvsPtCen[0] = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
  hMixMmumuvsPtCen[1] = (TH3D*)fmix->Get("hMixLPosMmumuvsPtCen");
  hMixMmumuvsPtCen[2] = (TH3D*)fmix->Get("hMixLNegMmumuvsPtCen");

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("paper_code/Rootfiles/%s.data.jpsi.root",run_type),"update");
      hStat->Write("",TObject::kOverwrite);
      hEvtRun->Write("",TObject::kOverwrite);
      hEvtRunAcc->Write("",TObject::kOverwrite);
      hEvtRunAccW->Write("",TObject::kOverwrite);

      for(int w=1; w<2; w++)
	{ 
	  for(Int_t j=0; j<3; j++)
	    { 
	      hnInvMass[w][j]->Write("",TObject::kOverwrite);
	    }
	}
      
      for(Int_t j=0; j<3; j++)
	{
	  hMixMmumuvsPtCen[j]->Write("",TObject::kOverwrite);
	}
      fout->Close();
	      
    }
}
