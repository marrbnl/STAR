TFile *f = 0x0;
const int year = YEAR;
TString run_cfg_name;

//================================================
void make_JpsiYield()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  f = TFile::Open(Form("./output/%s.jpsi.%sroot",run_type,run_config),"read");
  run_cfg_name = Form("%s",run_config);

  if(f)
    {
      TH1F *hStat = (TH1F*)f->Get("hEventStat");
      printf("%s\n",f->GetName());
      printf("all         events: %4.4e\n",hStat->GetBinContent(1));
      printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
      printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
    }

  make_histo_pt();
  //make_histo_npart();
}

//================================================
void make_histo_pt()
{
  const int nPtBins        = nPtBins_pt;
  const double* ptBinsLow  = ptBins_low_pt;
  const double* ptBinsHigh = ptBins_high_pt;
  const char** ptName      = pt_Name_pt;
  const int nCentBins      = nCentBins_pt; 
  const int* centBinsLow   = centBins_low_pt;
  const int* centBinsHigh  = centBins_high_pt;
  const char** centName    = cent_Name_pt;

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[2][3] = {0x0};
  TH1F *hInvMass[2][5][nCentBins][nPtBins][3] = {0x0};
  
  // same event
  char name[512];
  for(int w=0; w<2; w++) // event weights
    { 
      for(Int_t j=0; j<3; j++) // pair type
	{ 
	  if(w==0) sprintf(name,"m%s_%s",hName[j],trigName[kTrigType]);
	  else     sprintf(name,"m%sWeight_%s",hName[j],trigName[kTrigType]);
	  hnInvMass[w][j] = (THnSparseF*)f->Get(name);
	  if(!hnInvMass[w][j]) continue; 
	  hnInvMass[w][j]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[w][j]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
	  for(Int_t i=0; i<nPtBins; i++) // pt bins
	    {
	      hnInvMass[w][j]->GetAxis(1)->SetRangeUser(ptBinsLow[i]+0.01,ptBinsHigh[i]-0.01);
	      for(int k=0; k<nCentBins; k++) // centrality bins
		{
		  hnInvMass[w][j]->GetAxis(4)->SetRange(centBinsLow[k],centBinsHigh[k]);
		  for(int t=0; t<gNTrgSetup; t++) // trigger setup
		    {
		      if(t>0) hnInvMass[w][j]->GetAxis(5)->SetRange(t,t);
		      hInvMass[w][t][k][i][j] = (TH1F*)hnInvMass[w][j]->Projection(0);
		      hInvMass[w][t][k][i][j]->SetName(Form("%d_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_P%d",w,hName[j],trigName[kTrigType],i,k,t));
		      hInvMass[w][t][k][i][j]->Sumw2();
		      hnInvMass[w][j]->GetAxis(5)->SetRange(0,-1);
		    }
		  hnInvMass[w][j]->GetAxis(4)->SetRange(0,-1);
		}
	      hnInvMass[w][j]->GetAxis(1)->SetRange(0,-1);
	    }
	  hnInvMass[w][j]->GetAxis(2)->SetRange(0,-1);
	  hnInvMass[w][j]->GetAxis(3)->SetRange(0,-1);
	}
    }

  for(int w=0; w<2; w++) 
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      for(Int_t i=0; i<nPtBins; i++)
		{
		  if(hInvMass[w][t][k][i][1])
		    hInvMass[w][t][k][i][1]->Add(hInvMass[w][t][k][i][2]);
		}
	    }
	}
    }

  // mixed event
  TFile *fmix = 0;
  if(year==2014) 
    {
      char *mixName = Form("%s.Mix.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
      fmix = TFile::Open(Form("Output/%s",mixName),"read");

      cout << "Mix file: " << fmix->GetName() << endl;
      TH1F *hMixInvMass[nCentBins][nPtBins][3];
      printf("INFO: using Shuai's mixed events\n");
      TH3D *hMixMmumuvsPtCen[3];
      hMixMmumuvsPtCen[0] = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
      hMixMmumuvsPtCen[1] = (TH3D*)fmix->Get("hMixLPosMmumuvsPtCen");
      hMixMmumuvsPtCen[2] = (TH3D*)fmix->Get("hMixLNegMmumuvsPtCen");
      for(Int_t j=0; j<3; j++)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      int ybin_min = hMixMmumuvsPtCen[j]->GetYaxis()->FindFixBin(ptBinsLow[i]+1e-4);
	      int ybin_max = hMixMmumuvsPtCen[j]->GetYaxis()->FindFixBin(ptBinsHigh[i]-1e-4);
	      for(int k=0; k<nCentBins; k++)
		{
		  TH1F *htmp = (TH1F*)hMixMmumuvsPtCen[j]->ProjectionZ(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_tmp",hName[j],trigName[kTrigType],i,k),centBinsLow[k],centBinsHigh[k],ybin_min,ybin_max);
		  hMixInvMass[k][i][j] = new TH1F(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",hName[j],trigName[kTrigType],i,k),htmp->GetTitle(),1400,0,14);
		  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
		    {
		      hMixInvMass[k][i][j]->SetBinContent(bin,htmp->GetBinContent(bin));
		      hMixInvMass[k][i][j]->SetBinError(bin,htmp->GetBinError(bin));
		    }
		}
	    }
	}
      for(int k=0; k<nCentBins; k++)
        {
          for(Int_t i=0; i<nPtBins; i++)
            {
              hMixInvMass[k][i][1]->Add(hMixInvMass[k][i][2]);
            }
        }
    }

  TString fileName = Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TFile *fout = NULL;
  if(nCentBins_pt==5) fout = TFile::Open(fileName,"recreate");
  else  fout = TFile::Open(fileName,"update");
  const char* pair_name[2] = {"UL","LS"};
  for(int w=0; w<2; w++) 
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      for(Int_t i=0; i<nPtBins; i++)
		{
		  if(hInvMass[w][t][k][i][0])
		    {
		      for(int j=0; j<2; j++)
			{
			  hInvMass[w][t][k][i][j]->SetTitle("");
			  hInvMass[w][t][k][i][j]->Write(Form("InvMass_%s_pt%s_cent%s%s%s",pair_name[j],ptName[i],centName[k],gWeightName[w],gTrgSetupName[t]),TObject::kOverwrite);
			}
		    }
		}
	    }
	}
    }
  
  
  if(fmix)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hMixInvMass[k][i][j]->SetTitle("");
		  hMixInvMass[k][i][j]->Write(Form("Mix_InvMass_%s_pt%s_cent%s",pair_name[j],ptName[i],centName[k]),TObject::kOverwrite);
		}
	    }
	}
    }
  fout->Close();
}


//================================================
void make_histo_npart()
{
  const int nPtBins         = nPtBins_npart;
  const double* ptBinsLow   = ptBins_low_npart;
  const double* ptBinsHigh  = ptBins_high_npart;
  const char** ptName       = pt_Name_npart;
  const int* nCentBins      = nCentBins_npart; 
  const int* centBinsLow    = centBins_low_npart;
  const int* centBinsHigh   = centBins_high_npart;
  const char** centName     = cent_Name_npart;
  const int kNCent          = nCentBins[0];

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[2][3] = {0x0};
  TH1F *hInvMass[2][5][kNCent+1][nPtBins][3] = {0x0};
  
  // same event
  char name[512];
  for(int w=0; w<2; w++) // event weights
    { 
      for(Int_t j=0; j<3; j++) // pair type
	{ 
	  if(w==0) sprintf(name,"m%s_%s",hName[j],trigName[kTrigType]);
	  else     sprintf(name,"m%sWeight_%s",hName[j],trigName[kTrigType]);
	  hnInvMass[w][j] = (THnSparseF*)f->Get(name);
	  if(!hnInvMass[w][j]) continue; 
	  hnInvMass[w][j]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[w][j]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
	  for(Int_t i=0; i<nPtBins; i++) // pt bins
	    {
	      hnInvMass[w][j]->GetAxis(1)->SetRangeUser(ptBinsLow[i]+0.01,ptBinsHigh[i]-0.01);
	      for(int k=0; k<nCentBins[i]+1; k++) // centrality bins
		{
		  if(k<nCentBins[i]) hnInvMass[w][j]->GetAxis(4)->SetRange(centBinsLow[i*kNCent+k],centBinsHigh[i*kNCent+k]);
		  else hnInvMass[w][j]->GetAxis(4)->SetRange(1,16);
		  for(int t=0; t<gNTrgSetup; t++) // trigger setup
		    {
		      if(t>0) hnInvMass[w][j]->GetAxis(5)->SetRange(t,t);
		      hInvMass[w][t][k][i][j] = (TH1F*)hnInvMass[w][j]->Projection(0);
		      hInvMass[w][t][k][i][j]->SetName(Form("%d_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_P%d",w,hName[j],trigName[kTrigType],i,k,t));
		      hInvMass[w][t][k][i][j]->Sumw2();
		      hnInvMass[w][j]->GetAxis(5)->SetRange(0,-1);
		    }
		  hnInvMass[w][j]->GetAxis(4)->SetRange(0,-1);
		}
	      hnInvMass[w][j]->GetAxis(1)->SetRange(0,-1);
	    }
	  hnInvMass[w][j]->GetAxis(2)->SetRange(0,-1);
	  hnInvMass[w][j]->GetAxis(3)->SetRange(0,-1);
	}
    }

  for(int w=0; w<2; w++) 
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      for(int k=0; k<nCentBins[i]+1; k++)
		{
		  if(hInvMass[w][t][k][i][1])
		    hInvMass[w][t][k][i][1]->Add(hInvMass[w][t][k][i][2]);
		}
	    }
	}
    }

  // mixed event
  TFile *fmix = 0;
  if(year==2014) 
    {
      char *mixName = Form("%s.Mix.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
      fmix = TFile::Open(Form("Output/%s",mixName),"read");

      cout << "Mix file: " << fmix->GetName() << endl;
      TH1F *hMixInvMass[kNCent+1][nPtBins][3];
      printf("INFO: using Shuai's mixed events\n");
      TH3D *hMixMmumuvsPtCen[3];
      hMixMmumuvsPtCen[0] = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
      hMixMmumuvsPtCen[1] = (TH3D*)fmix->Get("hMixLPosMmumuvsPtCen");
      hMixMmumuvsPtCen[2] = (TH3D*)fmix->Get("hMixLNegMmumuvsPtCen");
      for(Int_t j=0; j<3; j++)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      int ybin_min = hMixMmumuvsPtCen[j]->GetYaxis()->FindFixBin(ptBinsLow[i]+1e-4);
	      int ybin_max = hMixMmumuvsPtCen[j]->GetYaxis()->FindFixBin(ptBinsHigh[i]-1e-4);
	      for(int k=0; k<nCentBins[i]+1; k++)
		{
		  TH1F *htmp = 0x0;
		  if(k<nCentBins[i]) htmp = (TH1F*)hMixMmumuvsPtCen[j]->ProjectionZ(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_tmp",hName[j],trigName[kTrigType],i,k),centBinsLow[i*kNCent+k],centBinsHigh[i*kNCent+k],ybin_min,ybin_max);
		  else htmp = (TH1F*)hMixMmumuvsPtCen[j]->ProjectionZ(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_tmp",hName[j],trigName[kTrigType],i,k),1,16,ybin_min,ybin_max);
		  hMixInvMass[k][i][j] = new TH1F(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",hName[j],trigName[kTrigType],i,k),htmp->GetTitle(),1400,0,14);
		  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
		    {
		      hMixInvMass[k][i][j]->SetBinContent(bin,htmp->GetBinContent(bin));
		      hMixInvMass[k][i][j]->SetBinError(bin,htmp->GetBinError(bin));
		    }
		}
	    }
	}
      for(Int_t i=0; i<nPtBins; i++)
	{
	  for(int k=0; k<nCentBins[i]+1; k++)
	    {
              hMixInvMass[k][i][1]->Add(hMixInvMass[k][i][2]);
            }
        }
    }

  TString fileName = Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TFile *fout = TFile::Open(fileName,"update");
  const char* pair_name[2] = {"UL","LS"};
  for(int w=0; w<2; w++) 
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      for(int k=0; k<nCentBins[i]+1; k++)
		{
		  if(hInvMass[w][t][k][i][0])
		    {
		      for(int j=0; j<2; j++)
			{
			  hInvMass[w][t][k][i][j]->SetTitle("");
			  if(k<nCentBins[i]) hInvMass[w][t][k][i][j]->Write(Form("InvMass_%s_pt%s_cent%s%s%s",pair_name[j],ptName[i],centName[i*kNCent+k],gWeightName[w],gTrgSetupName[t]),TObject::kOverwrite);
			  else hInvMass[w][t][k][i][j]->Write(Form("InvMass_%s_pt%s_cent0-80%s%s",pair_name[j],ptName[i],gWeightName[w],gTrgSetupName[t]),TObject::kOverwrite);
			}
		    }
		}
	    }
	}
    }
  
  
  if(fmix)
    {
      for(Int_t i=0; i<nPtBins; i++)
	{
	  for(int k=0; k<nCentBins[i]+1; k++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hMixInvMass[k][i][j]->SetTitle("");
		  if(k<nCentBins[i]) hMixInvMass[k][i][j]->Write(Form("Mix_InvMass_%s_pt%s_cent%s",pair_name[j],ptName[i],centName[i*kNCent+k]),TObject::kOverwrite);
		  else hMixInvMass[k][i][j]->Write(Form("Mix_InvMass_%s_pt%s_cent0-80",pair_name[j],ptName[i]),TObject::kOverwrite);
		}
	    }
	}
    }
  fout->Close();
}

