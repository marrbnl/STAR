const int year = YEAR;
TFile *f;

//================================================
void make_EmbTrkEff()
{
  gStyle->SetOptStat(1);

  if(year==2013)
    {
      f = TFile::Open(Form("./output/Run13.pp500.jpsi.Embed.root"),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Run14_AuAu200.Embed.Jpsi.root"),"read");
    }
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("[i] # of events: %4.4e\n",hStat->GetBinContent(3));

  makeHistos();
}

//================================================
void makeHistos(const int saveHistos = 1)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  printf("+++ Single track resolution +++\n");
  const int nTypes = 2;
  const char *trackType[nTypes] = {"p","g"};
  THnSparseF *hnTrkPtRes[nTypes];
  const int nCentBins_tmp = 11;
  const int centBins_low_tmp[nCentBins_tmp]  = {1, 13,9, 5,1,15,13,11,9, 7,5};
  const int centBins_high_tmp[nCentBins_tmp] = {16,16,12,8,4,16,14,12,10,8,6};
  const char *cent_Name_tmp[nCentBins_tmp] = {"0-80","0-20","20-40","40-60","60-80","0-10","10-20","20-30","30-40","40-50","50-60"};
  const char *cent_Title_tmp[nCentBins_tmp] = {"0080","0020","2040","4060","6080","0010","1020","2030","3040","4050","5060"};
  TH2F *hResVsRecoPt[nTypes][nCentBins_tmp];
  TH2F *hResVsTruePt[nTypes][nCentBins_tmp];
  for(int i=0; i<2; i++)
    {
      hnTrkPtRes[i]= (THnSparseF*)f->Get(Form("mh%sTrkPtRes_di_mu",trackType[i]));
      for(int k=0; k<nCentBins_tmp; k++)
	{
	  hnTrkPtRes[i]->GetAxis(3)->SetRange(centBins_low_tmp[k],centBins_high_tmp[k]);
	  hResVsRecoPt[i][k] = (TH2F*)hnTrkPtRes[i]->Projection(0,1);
	  hResVsRecoPt[i][k]->SetName(Form("%sTrkRes_vs_RecoPt_%s",trackType[i],cent_Title_tmp[k]));
	  hResVsRecoPt[i][k]->SetTitle();
	  
	  hResVsTruePt[i][k] = (TH2F*)hnTrkPtRes[i]->Projection(0,2);
	  hResVsTruePt[i][k]->SetName(Form("%sTrkRes_vs_TruePt_%s",trackType[i],cent_Title_tmp[k]));
	  hResVsTruePt[i][k]->SetTitle();
      
	  hnTrkPtRes[i]->GetAxis(3)->SetRange(0,-1);
	}
    }

  printf("+++ Single track efficiency +++\n");
  const int nDet = 4;
  const char *det_name[nDet] = {"","Tpc","Mtd","Final"};
  const char *charge_name[3] = {"","_pos","_neg"};
  const char *charge_title[3] = {" "," positive "," negative "};
  THnSparseF *hMcTrkInfo[nDet];
  TH1F *hMcTrkPt[nCentBins][nDet][3];
  TH1F *hMcTrkPtEff[nCentBins][nDet][3];
  TH3F *hMcTrkPtEtaPhi[nCentBins][nDet][3];
  TH3F *hMcTrkPtEtaPhiEff[nCentBins][nDet][3];
  const int nbinsX = 23, nbinsY = 8, nbinsZ = 360;
  const double xbinsX[nbinsX+1] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6,7,8,9,10,15,20};
  const double xbinsY[nbinsY+1] = {-0.7,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,0.7};
  double xbinsZ[nbinsZ+1];
  for(int i=0; i<=nbinsZ; i++)
    {
      xbinsZ[i] = pi/12 + 2*pi/360*i;
    }
  
  for(int i=0; i<nDet; i++)
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
	      hMcTrkPt[k][i][j]->SetTitle();
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);

	      TH3F *h3 = (TH3F*)hMcTrkInfo[i]->Projection(0,1,2);
	      h3->SetName(Form("McTrkPtEtaPhi%d%d%d",i,j,k));
	      hMcTrkPtEtaPhi[k][i][j] = new TH3F(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),"",nbinsX,xbinsX,nbinsY,xbinsY,nbinsZ,xbinsZ);
	      hMcTrkPtEtaPhi[k][i][j]->Sumw2();
	      for(int binx=1; binx<=nbinsX; binx++)
		{
		  for(int biny=1; biny<=nbinsY; biny++)
		    {
		      for(int binz=1; binz<=nbinsZ; binz++)
			{
			  int new_binz = binz+15;
			  if(binz>360) new_binz = binz - 360;
			  int new_binx_low = h3->GetXaxis()->FindFixBin(xbinsX[binx-1]+1e-4);
			  int new_binx_hig = h3->GetXaxis()->FindFixBin(xbinsX[binx]-1e-4);
			  int new_biny_low = h3->GetYaxis()->FindFixBin(xbinsY[biny-1]+1e-4);
			  int new_biny_hig = h3->GetYaxis()->FindFixBin(xbinsY[biny]-1e-4);
			  double value = 0;
			  for(int new_binx=new_binx_low; new_binx<=new_binx_hig; new_binx++)
			    {
			      for(int new_biny=new_biny_low; new_biny<=new_biny_hig; new_biny++)
				{
				  value += h3->GetBinContent(new_binx, new_biny, new_binz);
				}
			    }
			  hMcTrkPtEtaPhi[k][i][j]->SetBinContent(binx,biny,binz,value);
			  hMcTrkPtEtaPhi[k][i][j]->SetBinError(binx,biny,binz,sqrt(value));
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

  printf("+++ MTD PID+trigger efficiency +++\n");
  TH1F *hMtdEff[30][5];
  TH3F *hMtdInPtBlMod = (TH3F*)f->Get("mhTrkPtBlModMtd_di_mu");
  TH3F *hMtdOutPtBlMod = (TH3F*)f->Get("mhTrkPtBlModTrig_di_mu");
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  int bl = i+1; 
	  int mod = j+1;
	  hMtdEff[i][j] = (TH1F*)hMtdOutPtBlMod->ProjectionX(Form("RecoTrkPt_MtdEff_BL%d_Mod%d",bl,mod),bl,bl,mod,mod);
	  hMtdEff[i][j]->Rebin(5);
	  hMtdEff[i][j]->Sumw2();
	  TH1F *htmp = (TH1F*)hMtdInPtBlMod->ProjectionX(Form("Temp_BL%d_Mod%d",bl,mod),bl,bl,mod,mod);
	  htmp->Rebin(5);
	  hMtdEff[i][j]->Divide(htmp);
	  hMtdEff[i][j]->SetTitle();
	}
    }

  printf("+++ Tracking efficiency vs. Zdc +++\n");
  const int nEffType = 6;
  const char *trkEffType[nEffType] = {"MC","Tpc","MtdMth","Fake","MuonPid","MtdTrig"};
  THnSparseF *hnMcTrkPt[nEffType];
  TH2F *hMcTrkPtVsCent[nEffType][gNTrgSetup];
  TH1F *hMcTrkPtVsPt[nEffType][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsZdc[nEffType][gNTrgSetup][nCentBins];
  for(int i=0; i<nEffType; i++)
    {
      hnMcTrkPt[i] = (THnSparseF*)f->Get(Form("mhMcTrkPtEff_%s_di_mu",trkEffType[i]));
      hnMcTrkPt[i]->GetAxis(0)->SetRangeUser(0,20);
      hnMcTrkPt[i]->GetAxis(1)->SetRangeUser(-0.5,0.5);

      for(int j=0; j<gNTrgSetup; j++)
	{
	  if(j>0) hnMcTrkPt[i]->GetAxis(4)->SetRange(j,j);
	  hMcTrkPtVsCent[i][j] = (TH2F*)hnMcTrkPt[i]->Projection(0,2);
	  hMcTrkPtVsCent[i][j]->SetName(Form("McTrkPtVsCent_%s%s",trkEffType[i],gTrgSetupTitle[j]));
	  hMcTrkPtVsCent[i][j]->SetTitle("");

	  for(int k=0; k<nCentBins; k++)
	    {
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(centBins_low[k],centBins_high[k]);
	      hMcTrkPtVsPt[i][j][k] = (TH1F*)hnMcTrkPt[i]->Projection(0);
	      hMcTrkPtVsPt[i][j][k]->SetName(Form("McTrkPt_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      hMcTrkPtVsPt[i][j][k]->SetTitle("");

	      hMcTrkPtVsZdc[i][j][k] = (TH2F*)hnMcTrkPt[i]->Projection(0,3);
	      hMcTrkPtVsZdc[i][j][k]->SetName(Form("McTrkPtVsZdc_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      hMcTrkPtVsZdc[i][j][k]->SetTitle("");
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(0,-1);
	    }
	  hnMcTrkPt[i]->GetAxis(4)->SetRange(0,-1);
	}
    }
  
  // efficiency vs. centrality
  const int kNCent = nCentBins_npart[0];
  TH1F *hMcTrkPtVsNpart[nEffType][gNTrgSetup][kNCent];
  for(int i=0; i<nEffType; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  if(j>0) hnMcTrkPt[i]->GetAxis(4)->SetRange(j,j);
	  for(int k=0; k<kNCent; k++)
	    {
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(centBins_low_npart[k],centBins_high_npart[k]);
	      hMcTrkPtVsNpart[i][j][k] = (TH1F*)hnMcTrkPt[i]->Projection(0);
	      hMcTrkPtVsNpart[i][j][k]->SetName(Form("McTrkPt_%s_cent%s%s",trkEffType[i],cent_Title_npart[k],gTrgSetupTitle[j]));
	      hMcTrkPtVsNpart[i][j][k]->SetTitle("");
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(0,-1);
	    }
	  hnMcTrkPt[i]->GetAxis(4)->SetRange(0,-1);
	}
    }

  TH1F *hMcTrkPtInZdc[2][gNZdcRate][kNCent];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<gNZdcRate; j++)
	{
	  hnMcTrkPt[i]->GetAxis(3)->SetRange(j+1,j+1);
	  for(int k=0; k<kNCent; k++)
	    {
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(centBins_low_npart[k],centBins_high_npart[k]);
	      hMcTrkPtInZdc[i][j][k] = (TH1F*)hnMcTrkPt[i]->Projection(0);
	      hMcTrkPtInZdc[i][j][k]->SetName(Form("McTrkPt_%s_cent%s_Zdc%d-%d",trkEffType[i],cent_Title_npart[k],j*10,j*10+10));
	      hMcTrkPtInZdc[i][j][k]->SetTitle("");
	      hnMcTrkPt[i]->GetAxis(2)->SetRange(0,-1);
	    }
	  hnMcTrkPt[i]->GetAxis(3)->SetRange(0,-1);
	}
    }

  if(saveHistos)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.EmbTrkEff.root",run_type),"recreate");
      for(int k=0; k<nCentBins_tmp; k++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hResVsRecoPt[i][k]->Write();
	      hResVsTruePt[i][k]->Write();
	    }
	}
      for(int k=0; k<nCentBins; k++)
	{
	  for(int i=0; i<nDet; i++)
	    {
	      for(int j=0; j<3; j++)
		{
		  hMcTrkPt[k][i][j]->Write();
		  hMcTrkPtEff[k][i][j]->Write();
		  hMcTrkPtEtaPhi[k][i][j]->Write();
		  hMcTrkPtEtaPhiEff[k][i][j]->Write();
		}
	    }
	}

      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      hMtdEff[i][j]->Write();
	    }
	}

     for(int i=0; i<nEffType; i++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      hMcTrkPtVsCent[i][j]->Sumw2();
	      hMcTrkPtVsCent[i][j]->Write();
	      for(int k=0; k<nCentBins; k++)
		{
		  hMcTrkPtVsPt[i][j][k]->Sumw2();
		  hMcTrkPtVsPt[i][j][k]->Write();
		  hMcTrkPtVsZdc[i][j][k]->Sumw2();
		  hMcTrkPtVsZdc[i][j][k]->Write();
		}
	    }
	}

     for(int i=0; i<nEffType; i++)
       {
	 for(int j=0; j<gNTrgSetup; j++)
	   {
	     for(int k=0; k<kNCent; k++)
	       {
		 hMcTrkPtVsNpart[i][j][k]->Write();
	       }
	   }
       }

     for(int i=0; i<2; i++)
       {
	 for(int j=0; j<gNZdcRate; j++)
	   {
	     for(int k=0; k<kNCent; k++)
	       {
		 hMcTrkPtInZdc[i][j][k]->Write();
	       }
	   }
       }

    }
}
