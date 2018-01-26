
//================================================
void make_AnaInfo()
{
  Run14_AuAu200();
}

//================================================
void Run14_AuAu200(const int saveHisto = 1)
{
  // ---------------------------------------------
  // Jpsi efficiency
  TFile *feff = TFile::Open(Form("Rootfiles/Run14_AuAu200.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  const int nHistos = 6;
  const char *trkEffType[nHistos] = {"MC","Tpc","MtdMth","MuonPid","MtdTrig","TrigUnit"};
  TH1F *hJpsiPt[nHistos][nCentBins_pt];
  TH1F *hJpsiPtEff[nHistos][nCentBins_pt];
  for(int i=0; i<nHistos; i++)
    {
      for(int k=0; k<nCentBins_pt; k++)
	{
	  hJpsiPt[i][k] = (TH1F*)feff->Get(Form("hJpsiPt_%s_cent%s",trkEffType[i],cent_Title_pt[k]));
	  int index = i-1;
	  if(i==0) index = 0;
	  hJpsiPtEff[i][k] = DivideTH1ForEff(hJpsiPt[i][k],hJpsiPt[index][k],Form("hJpsiPtEff_%s_cent%s",trkEffType[i],cent_Title_pt[k]));
	}
    }

  // ---------------------------------------------
  // Jpsi efficiency for uncertainty studies
  //
  TH1F *hJpsiEffVsPt[11];
  hJpsiEffVsPt[0] = (TH1F*)feff->Get(Form("JpsiEffVsPt_cent%s_final",cent_Title_pt[0]));
  hJpsiEffVsPt[0]->SetName(Form("JpsiEffVsPt_cent%s_default",cent_Title_pt[0]));
  const TString sysName[10] = {"dcaUp","dcaDown","NHitsUp","NDedxUp",
			       "dzUp","dzDown","dyUp","dyDown","nSigPiUp","nSigPiDown"};

  TFile *fsys[2];
  fsys[0] = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.TpcTracking.root"),"read");
  fsys[1] = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.MuonPid.root"),"read");
  TH1F *hJpsiTpcEffVsPt[2];
  for(int i=0; i<2; i++)
    {
      hJpsiTpcEffVsPt[i] = (TH1F*)fsys[i]->Get(Form("Sysdefault_JpsiEffVsPtEmb_cent%s",cent_Title_pt[0]));
      hJpsiTpcEffVsPt[i]->SetName(Form("%s_%d",hJpsiTpcEffVsPt[i]->GetName(),i));
    }
  c = draw1D(hJpsiEffVsPt[0]);
  for(int s=0; s<10; s++)
    {
      int index = 0;
      if(s>=4) index = 1; 
      hJpsiEffVsPt[s+1] = (TH1F*)fsys[index]->Get(Form("Sys%s_JpsiEffVsPtEmb_cent%s",sysName[s].Data(),cent_Title_pt[0]));
      hJpsiEffVsPt[s+1]->SetName(Form("JpsiEffVsPt_cent%s_%s",cent_Title_pt[0],sysName[s].Data()));
      hJpsiEffVsPt[s+1]->Divide(hJpsiTpcEffVsPt[index]);
      hJpsiEffVsPt[s+1]->Multiply(hJpsiEffVsPt[0]);
      hJpsiEffVsPt[s+1]->SetMarkerStyle(21+s);
      hJpsiEffVsPt[s+1]->Draw("sames");
    }


  // ---------------------------------------------
  // Fully corrected Jpsi invariant yield
  TFile *fxsec = TFile::Open("Rootfiles/Run14_AuAu200.JpsiXsec.pt1.5.pt1.3.root","read");
  TH1F *hJpsiYieldVsPt[nCentBins_pt];
  TH1F *hJpsiRaaVsPt[nCentBins_pt];
  TF1 *funcJpsiYield[nCentBins_pt];
  for(int i=0; i<nCentBins_pt; i++)
    {
      hJpsiYieldVsPt[i] = (TH1F*)fxsec->Get(Form("Jpsi_InvYieldVsPt_cent%s",cent_Title_pt[i]));
      hJpsiRaaVsPt[i] = (TH1F*)fxsec->Get(Form("Jpsi_RaaVsPt_cent%s",cent_Title_pt[i]));
      funcJpsiYield[i] = (TF1*)fxsec->Get(Form("func_JpsiYieldVsPt_cent%s",cent_Title_pt[i]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.AnaInfo.root","update");
      for(int s=0; s<11; s++)
	{
	  hJpsiEffVsPt[s]->Write("", TObject::kOverwrite);
	}
      for(int k=0; k<nCentBins_pt; k++)
	{
	  for(int i=0; i<nHistos; i++)
	    {
	      if(i<4) hJpsiPt[i][k]->Write(Form("Emb_JpsiPt_%s_cent%s",trkEffType[i],cent_Title_pt[k]), TObject::kOverwrite);
	      if(i>0 && i<4) hJpsiPtEff[i][k]->Write(Form("Emb_JpsiEffVsPt_%s_cent%s",trkEffType[i],cent_Title_pt[k]), TObject::kOverwrite);
	    }
	}
      for(int i=0; i<nCentBins_pt; i++)
	{
	  hJpsiYieldVsPt[i]->Write(Form("JpsiInvYieldVsPt_cent%s",cent_Title_pt[i]), TObject::kOverwrite);
	  funcJpsiYield[i]->Write(Form("JpsiInvYieldVsPtFit_cent%s",cent_Title_pt[i]), TObject::kOverwrite);
	}
    }
}
