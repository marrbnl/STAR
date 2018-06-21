//================================================
void make_code()
{
  gStyle->SetOptStat(0);
  //data();
  //JpsiWidth();
  lumi();

  //embed();
  //input();
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
