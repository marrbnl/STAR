const int year = YEAR;
const char *trkEffType[6] = {"MC","Tpc","MtdMth","Fake","MuonPid","MtdTrig"};
const char *weight_name[2] = {"","_w"};

TFile *f;

//================================================
void make_Embed()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      f = TFile::Open(Form("./output/Run13.pp500.jpsi.Embed.root"),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Run14.AuAu200.Jpsi.Embed.root"),"read");
    }
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("[i] # of events: %4.4e\n",hStat->GetBinContent(3));

  //makeTracks();
  makeJets();
}

//================================================
void makeJets(const bool saveHisto = 1)
{
  THnSparseF *hnJpsiInfo[6][2];
  TH1F *hJpsiInvMass[6][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPt[6][gNTrgSetup][nCentBins][2];
  TH2F *hJpsiMassVsPt[6][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiRapdity[6][gNTrgSetup][nCentBins][2];

  for(int i=0; i<6; i++)
    {
      for(int w=0; w<2; w++)
	{
	  hnJpsiInfo[i][w] = (THnSparseF*)f->Get(Form("hJpsiInfo_%s_di_mu%s",trkEffType[i],weight_name[w]));
	  if(i>0)
	    {
	      hnJpsiInfo[i][w]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	      hnJpsiInfo[i][w]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	    }
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      if(j>0) hnJpsiInfo[i][w]->GetAxis(6)->SetRange(j,j);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnJpsiInfo[i][w]->GetAxis(5)->SetRange(centBins_low[k],centBins_high[k]);

		  hJpsiInvMass[i][j][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(0);
		  hJpsiInvMass[i][j][k][w]->SetName(Form("hJpsiInvMass_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiInvMass[i][j][k][w]->SetTitle("");
		  hJpsiInvMass[i][j][k][w]->Sumw2();

		  hJpsiPt[i][j][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(1);
		  hJpsiPt[i][j][k][w]->SetName(Form("hJpsiPt_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiPt[i][j][k][w]->SetTitle("");
		  hJpsiPt[i][j][k][w]->SetBinContent(hJpsiPt[i][j][k][w]->GetNbinsX()+1,0); // reset overflow bin
		  hJpsiPt[i][j][k][w]->Sumw2();

		  hJpsiMassVsPt[i][j][k][w] = (TH2F*)hnJpsiInfo[i][w]->Projection(0,1);
		  hJpsiMassVsPt[i][j][k][w]->SetName(Form("hJpsiMassVsPt_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiMassVsPt[i][j][k][w]->SetTitle("");
		  hJpsiMassVsPt[i][j][k][w]->Sumw2();

		  hJpsiRapdity[i][j][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(2);
		  hJpsiRapdity[i][j][k][w]->SetName(Form("hJpsiRapdity_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiRapdity[i][j][k][w]->SetTitle("");
		  hJpsiRapdity[i][j][k][w]->Sumw2();

		  hnJpsiInfo[i][w]->GetAxis(5)->SetRange(0,-1);
		}
	      hnJpsiInfo[i][w]->GetAxis(6)->SetRange(0,-1);
	    }
	  hnJpsiInfo[i][w]->GetAxis(3)->SetRange(0,-1);
	  hnJpsiInfo[i][w]->GetAxis(4)->SetRange(0,-1);
	}
    }
  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"recreate");
      for(int i=0; i<6; i++)
	{
	  for(int w=0; w<2; w++)
	    {
	      for(int j=0; j<gNTrgSetup; j++)
		{
		  for(int k=0; k<nCentBins; k++)
		    {
		      hJpsiMassVsPt[i][j][k][w]->Write();
		      hJpsiInvMass[i][j][k][w]->Write();
		      hJpsiPt[i][j][k][w]->Write();
		      hJpsiRapdity[i][j][k][w]->Write();
		    }
		}
	    }
	}
    }
}

//================================================
void makeTracks(const bool saveHisto = 1)
{
  // tracking efficiency
  THnSparseF *hnMcTrkPt[6];
  TH1F *hMcTrkPt[6][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsZdc[6][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsCent[6][gNTrgSetup];

  for(int i=0; i<6; i++)
    {
      hnMcTrkPt[i] = (THnSparseF*)f->Get(Form("mhMcTrkPtEff_%s_di_mu",trkEffType[i]));
      hnMcTrkPt[i]->GetAxis(0)->SetRangeUser(0,20);
      hnMcTrkPt[i]->GetAxis(1)->SetRangeUser(-0.45,0.45);

      for(int j=0; j<gNTrgSetup; j++)
	{
	  if(j>0) hnMcTrkPt[i]->GetAxis(4)->SetRange(j,j);
	  hMcTrkPtVsCent[i][j] = (TH2F*)hnMcTrkPt[i]->Projection(0,2);
	  hMcTrkPtVsCent[i][j]->SetName(Form("hMcTrkPtVsCent_%s%s",trkEffType[i],gTrgSetupTitle[j]));
	  hMcTrkPtVsCent[i][j]->SetTitle("");

	  for(int k=0; k<nCentBins; k++)
	    {
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(centBins_low[k],centBins_high[k]);
	      hMcTrkPt[i][j][k] = (TH1F*)hnMcTrkPt[i]->Projection(0);
	      hMcTrkPt[i][j][k]->SetName(Form("hMcTrkPt_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      hMcTrkPt[i][j][k]->SetTitle("");

	      hMcTrkPtVsZdc[i][j][k] = (TH2F*)hnMcTrkPt[i]->Projection(0,3);
	      hMcTrkPtVsZdc[i][j][k]->SetName(Form("hMcTrkPtVsZdcRate_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      hMcTrkPtVsZdc[i][j][k]->SetTitle("");
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(0,-1);
	    }
	  hnMcTrkPt[i]->GetAxis(4)->SetRange(0,-1);
	}
    }

  TH1F *hMcTrkPt2[6][gNTrgSetup][6];
  for(int i=0; i<6; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  if(j>0) hnMcTrkPt[i]->GetAxis(4)->SetRange(j,j);
	  for(int k=0; k<6; k++)
	    {
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(5+k*2,6+k*2);
	      hMcTrkPt2[i][j][k] = (TH1F*)hnMcTrkPt[i]->Projection(0);
	      hMcTrkPt2[i][j][k]->SetName(Form("hMcTrkPt_%s_cent%d%d%s",trkEffType[i],50-k*10,60-k*10,gTrgSetupTitle[j]));
	      hMcTrkPt2[i][j][k]->SetTitle("");
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(0,-1);
	    }
	  hnMcTrkPt[i]->GetAxis(4)->SetRange(0,-1);
	}
    }

  // momentum resolution for primary tracks
  THnSparseF *hnTrkPtRes = (THnSparseF*)f->Get("mhpTrkPtRes_di_mu");
  TH2F *hResVsRecoPt[nCentBins];
  TH2F *hResVsTruePt[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hnTrkPtRes->GetAxis(3)->SetRange(centBins_low[k],centBins_high[k]);
      hResVsRecoPt[k] = (TH2F*)hnTrkPtRes->Projection(0,1);
      hResVsRecoPt[k]->SetName(Form("PrimTrkRes_vs_RecoPt_cent%s",cent_Title[k]));

      hResVsTruePt[k] = (TH2F*)hnTrkPtRes->Projection(0,2);
      hResVsTruePt[k]->SetName(Form("PrimTrkRes_vs_TruePt_cent%s",cent_Title[k]));

      hnTrkPtRes->GetAxis(3)->SetRange(0,-1);
    }

  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  // single tracks
	  hResVsRecoPt[k]->Sumw2();
	  hResVsRecoPt[k]->Write();

	  hResVsTruePt[k]->Sumw2();
	  hResVsTruePt[k]->Write();
	}

      for(int i=0; i<6; i++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      hMcTrkPtVsCent[i][j]->Sumw2();
	      hMcTrkPtVsCent[i][j]->Write();
	      for(int k=0; k<nCentBins; k++)
		{
		  hMcTrkPt[i][j][k]->Sumw2();
		  hMcTrkPt[i][j][k]->Write();
		  hMcTrkPtVsZdc[i][j][k]->Sumw2();
		  hMcTrkPtVsZdc[i][j][k]->Write();
		}
	      for(int k=0; k<6; k++)
		{
		  hMcTrkPt2[i][j][k]->Sumw2();
		  hMcTrkPt2[i][j][k]->Write();
		}
	    }
	}
    }
}

//================================================
void makeHistos(TString outName, const bool save = 1)
{
  gStyle->SetOptStat(1);

  printf("+++ Single track efficiency +++\n");
  THnSparseF *hMcTrkInfo[5], *hRcTrkInfo[4];
  TH1F *hMcTrkPt[nCentBins][5][3], *hRcTrkPt[nCentBins][4][3];
  TH1F *hMcTrkPtEff[nCentBins][5][3], *hRcTrkPtEff[nCentBins][4][3];
  TH3F *hMcTrkPtEtaPhi[nCentBins][5][3], *hRcTrkPtEtaPhi[nCentBins][4][3];
  TH3F *hMcTrkPtEtaPhiEff[nCentBins][5][3], *hRcTrkPtEtaPhiEff[nCentBins][4][3];
  for(int i=0; i<5; i++)
    {
      hMcTrkInfo[i] = (THnSparseF*)f->Get(Form("mhMcTrkInfo%s_di_mu",det_name[i]));
      hMcTrkInfo[i]->SetTitle(Form("%s_jpsi",hMcTrkInfo[i]->GetName()));
      hMcTrkInfo[i]->Sumw2();

      for(int k=0; k<nCentBins; k++)
	{
	  if(nCentBins>1) hMcTrkInfo[i]->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkInfo[i]->GetAxis(3)->SetRange(1,3);
	      if(j==1) hMcTrkInfo[i]->GetAxis(3)->SetRange(3,3);
	      if(j==2) hMcTrkInfo[i]->GetAxis(3)->SetRange(1,1);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkInfo[i]->Projection(0);
	      hMcTrkPt[k][i][j]->Sumw2();
	      hMcTrkPt[k][i][j]->SetName(Form("McTrkPt%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);

	      TH3F *h3 = (TH3F*)hMcTrkInfo[i]->Projection(0,1,2);
	      h3->SetName(Form("McTrkPtEtaPhi%d%d%d",i,j,k));
	      hMcTrkPtEtaPhi[k][i][j] = new TH3F(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),h3->GetTitle(),200,0,20,40,-2,2,360,pi/12,2*pi+pi/12);
	      for(int binx=1; binx<=h3->GetNbinsX(); binx++)
		{
		  for(int biny=1; biny<=h3->GetNbinsY(); biny++)
		    {
		      for(int binz=1; binz<=h3->GetNbinsZ(); binz++)
			{
			  int new_binz = binz-15;
			  if(binz<=15) new_binz = binz + 360 - 15;
			  hMcTrkPtEtaPhi[k][i][j]->SetBinContent(binx,biny,new_binz,h3->GetBinContent(binx,biny,binz));
			}
		    }
		}
	      hMcTrkPtEtaPhi[k][i][j]->Sumw2();
	      hMcTrkPtEtaPhiEff[k][i][j] = (TH3F*)hMcTrkPtEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPtEtaPhi[k][i][j]->GetName()));
	      hMcTrkPtEtaPhiEff[k][i][j]->Divide(hMcTrkPtEtaPhi[k][0][j]);
	      hMcTrkInfo[i]->GetAxis(3)->SetRange(0,-1);
	    }
	  hMcTrkInfo[i]->GetAxis(4)->SetRange(0,-1);
	}
    }

  for(int i=0; i<4; i++)
    {
      hRcTrkInfo[i] = (THnSparseF*)f->Get(Form("mhRcTrkInfo%s_di_mu",det_name[i+1]));
      hRcTrkInfo[i]->SetTitle(Form("%s_jpsi",hRcTrkInfo[i]->GetName()));
      hRcTrkInfo[i]->Sumw2();

      for(int k=0; k<nCentBins; k++)
	{
	  if(nCentBins>1) hRcTrkInfo[i]->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hRcTrkInfo[i]->GetAxis(3)->SetRange(1,3);
	      if(j==1) hRcTrkInfo[i]->GetAxis(3)->SetRange(3,3);
	      if(j==2) hRcTrkInfo[i]->GetAxis(3)->SetRange(1,1);
	      hRcTrkPt[k][i][j] = (TH1F*)hRcTrkInfo[i]->Projection(0);
	      hRcTrkPt[k][i][j]->Sumw2();
	      hRcTrkPt[k][i][j]->SetName(Form("RcTrkPt%s%s_%s",det_name[i+1],charge_name[j],cent_Title[k]));
	      hRcTrkPtEff[k][i][j] = (TH1F*)hRcTrkPt[k][i][j]->Clone(Form("%s_Eff",hRcTrkPt[k][i][j]->GetName()));
	      hRcTrkPtEff[k][i][j]->Divide(hRcTrkPt[k][0][j]);

	      TH3F *h3 = (TH3F*)hRcTrkInfo[i]->Projection(0,1,2);
	      h3->SetName(Form("McTrkPtEtaPhi%d%d%d",i,j,k));
	      TString name_tmp = Form("RcTrkPtEtaPhi%s%s_%s",det_name[i+1],charge_name[j],cent_Title[k]);
	      hRcTrkPtEtaPhi[k][i][j] = new TH3F(name_tmp.Data(),h3->GetTitle(),200,0,20,40,-2,2,360,pi/12,2*pi+pi/12);
	      for(int binx=1; binx<=h3->GetNbinsX(); binx++)
		{
		  for(int biny=1; biny<=h3->GetNbinsY(); biny++)
		    {
		      for(int binz=1; binz<=h3->GetNbinsZ(); binz++)
			{
			  int new_binz = binz-15;
			  if(binz<=15) new_binz = binz + 360 - 15;
			  hRcTrkPtEtaPhi[k][i][j]->SetBinContent(binx,biny,new_binz,h3->GetBinContent(binx,biny,binz));
			}
		    }
		}
	      hRcTrkPtEtaPhi[k][i][j]->Sumw2();
	      hRcTrkPtEtaPhiEff[k][i][j] = (TH3F*)hRcTrkPtEtaPhi[k][i][j]->Clone(Form("%s_Eff",hRcTrkPtEtaPhi[k][i][j]->GetName()));
	      hRcTrkPtEtaPhiEff[k][i][j]->Divide(hRcTrkPtEtaPhi[k][0][j]);
	      hRcTrkInfo[i]->GetAxis(3)->SetRange(0,-1);
	    }
	  hRcTrkInfo[i]->GetAxis(4)->SetRange(0,-1);
	}
    }
  
  if(save)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outName.Data()),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  // single tracks
	  hResVsRecoPt[k]->Write();
	  hResVsTruePt[k]->Write();
	  for(int i=0; i<5; i++)
	    {
	      for(int j=0; j<3; j++)
		{
		  hMcTrkPt[k][i][j]->Write();
		  hMcTrkPtEff[k][i][j]->Write();
		  hMcTrkPtEtaPhi[k][i][j]->Write();
		  hMcTrkPtEtaPhiEff[k][i][j]->Write();
		}
	    }

	  for(int i=0; i<4; i++)
	    {
	      for(int j=0; j<3; j++)
		{
		  hRcTrkPt[k][i][j]->Write();
		  hRcTrkPtEff[k][i][j]->Write();
		  hRcTrkPtEtaPhi[k][i][j]->Write();
		  hRcTrkPtEtaPhiEff[k][i][j]->Write();
		}
	    }
	}
    }
}
