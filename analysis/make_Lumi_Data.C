#include "TStyle.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cmath>
using namespace std;


void makeHistoData(const int savePlot = 0, const int saveHisto = 0);

//================================================
void make_Lumi_Data()
{

  makeHistoData(1, 1);
}

//================================================
void makeHistoData(const int savePlot, const int saveHisto)
{
  int year = 2014;
  if(year!=2014)
    {
      printf("[e] Not suitable for %d\n",year);
      return;
    }
  const char *run_type = "Run14_AuAu200";
  const char *trigName = "VPD-ZDC-novtx-mon";
  const int max_vz = 100;
  const int max_dz = 3;
  const int max_vr = 2;

  const char *wName[2] = {"","_w"};
  const char *setupName[5] = {"all","prod","prod_low","prod_mid","prod_high"};
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  const int nCentBins = 17;
  const int centBins_low[nCentBins]  = {5,  13, 9,  5, 15, 11, 5,  13, 11, 9,  7, 5, 1,  1, 3, 1, 1};
  const int centBins_high[nCentBins] = {16, 16, 12, 8, 16, 14, 10, 14, 12, 10, 8, 6, 16, 4, 4, 2, 1};
  const char *cent_Name[nCentBins] = {"0-60","0-20","20-40","40-60","0-10","10-30","30-60","10-20","20-30","30-40","40-50","50-60","0-80","60-80","60-70","70-80","75-80"};
  const char *cent_Title[nCentBins] = {"0060","0020","2040","4060","0010","1030","3060","1020","2030","3040","4050","5060","0080","6080","6070","7080","7580"};
  TH1F *hNEvents[5][2];
  TH1F *hTpcVz[5][nCentBins];
  TH1F *hDiffVz[5][nCentBins];
  TH1F *hTpcVr[5][nCentBins];
  TH1F *hNEventsAll[5][nCentBins][2];
  TH1F *hNEventsVr[5][nCentBins][2];
  TH1F *hNEventsVz[5][nCentBins][2];
  TH1F *hNEventsAcc[5][nCentBins][2];
  TH1F *hCent[5][2];

  TH2F *hTpcVzVsRun[nCentBins];
  TH2F *hTpcVrVsRun[nCentBins];
  TH2F *hDiffVzVsRun[nCentBins];
  TH2F *hCentVsRun[2];

  TFile *fMB = TFile::Open(Form("./output/Run14_AuAu200.MB.VtxEff.root"),"read");
  THnSparseF *hn[2];
  hn[0] = (THnSparseF*)fMB->Get("mhMbEvtEff");
  hn[1] = (THnSparseF*)fMB->Get("mhMbEvtEffWeight");

  // Assign the prod_high to MB
  TFile *fhigh = TFile::Open(Form("./output/Run14_AuAu200.MB.VtxEff.prod_high.root"),"read");
  THnSparseF *hn_high[2];
  hn_high[0] = (THnSparseF*)fhigh->Get("mhMbEvtEff");
  hn_high[1] = (THnSparseF*)fhigh->Get("mhMbEvtEffWeight");
  for(int i=0; i<2; i++)
    {
      hn_high[i]->SetName(Form("%s_high",hn_high[i]->GetName()));
      const int nAxis = 6;
      int axisNbins[nAxis];
      for(int iaxis=0; iaxis<nAxis; iaxis++)
	{
	  axisNbins[iaxis] = hn_high[i]->GetAxis(iaxis)->GetNbins();
	}

      // the index correspond to bin number starting from 1
      TH1F *hRunStat = (TH1F*)hn_high[i]->Projection(0);
      hRunStat->SetName(Form("hRunStat%s_prod_high",wName[i]));
      for(int ibin=1; ibin<=axisNbins[0]; ibin++)
	{
	  if(hRunStat->GetBinContent(ibin)<=0) continue;
	  cout << "Haha: " << ibin << endl;
	  for(int jbin=1; jbin<=axisNbins[1]; jbin++)
	    {
	      for(int kbin=1; kbin<=axisNbins[2]; kbin++)
		{
		  for(int lbin=1; lbin<=axisNbins[3]; lbin++)
		    {
		      for(int mbin=1; mbin<=axisNbins[4]; mbin++)
			{
			  for(int nbin=1; nbin<=axisNbins[5]; nbin++)
			    {
			      int idx[6] = {ibin, jbin, kbin, lbin, mbin, nbin};
			      hn[i]->SetBinContent(idx, hn_high[i]->GetBinContent(idx));
			      hn[i]->SetBinError(idx, hn_high[i]->GetBinError(idx));
			    }
			}
		    }
		}
	    }
	}
    }
  return;
  

  for(int i=0; i<2; i++)
    {
      hNEvents[0][i] = (TH1F*)hn[i]->Projection(0);
      hNEvents[0][i]->SetName(Form("NEvents_%s%s",setupName[0],wName[i]));

      for(int j=0; j<nCentBins; j++)
	{
	  hn[i]->GetAxis(5)->SetRange(centBins_low[j],centBins_high[j]);
	  hNEventsAll[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsAll[0][j][i]->SetName(Form("NEvents_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	  if(i==0)
	    {
	      hTpcVz[0][j] = (TH1F*)hn[i]->Projection(1);
	      hTpcVz[0][j]->SetName(Form("TpcVz_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      hTpcVzVsRun[j] = (TH2F*)hn[i]->Projection(1,0);
	      hTpcVzVsRun[j]->SetName(Form("TpcVzVsRun_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      //draw2D(hTpcVzVsRun[j]);

	      hTpcVr[0][j] = (TH1F*)hn[i]->Projection(3);
	      hTpcVr[0][j]->SetName(Form("TpcVr_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      hTpcVrVsRun[j] = (TH2F*)hn[i]->Projection(3,0);
	      hTpcVrVsRun[j]->SetName(Form("TpcVrVsRun_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	      hDiffVz[0][j] = (TH1F*)hn[i]->Projection(2);
	      hDiffVz[0][j]->SetName(Form("TpcVpdDz_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      hDiffVzVsRun[j] = (TH2F*)hn[i]->Projection(2,0);
	      hDiffVzVsRun[j]->SetName(Form("DiffVzVsRun_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	    }

	  hn[i]->GetAxis(3)->SetRangeUser(-1,max_vr-1e-4);
	  hNEventsVr[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsVr[0][j][i]->SetName(Form("NEvents_VrCut_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	  hn[i]->GetAxis(1)->SetRangeUser(-1*max_vz+1e-4,max_vz-1e-4);
	  hNEventsVz[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsVz[0][j][i]->SetName(Form("NEvents_VrVzCut_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	  hn[i]->GetAxis(2)->SetRangeUser(-1*max_dz+1e-4,max_dz-1e-4);
	  hNEventsAcc[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsAcc[0][j][i]->SetName(Form("NEvents_VrVzDzCut_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	  hn[i]->GetAxis(5)->SetRange(0,-1);
	  if(j==0)
	    {
	      hCent[0][i] = (TH1F*)hn[i]->Projection(5);
	      hCent[0][i]->SetName(Form("Cent_%s%s",setupName[0],wName[i]));
	      hCentVsRun[i] = (TH2F*)hn[i]->Projection(5,0);
	      hCentVsRun[i]->SetName(Form("CentVsRun_%s%s",setupName[0],wName[i]));
	    }

	  hn[i]->GetAxis(1)->SetRange(0,-1);
	  hn[i]->GetAxis(2)->SetRange(0,-1);
	  hn[i]->GetAxis(3)->SetRange(0,-1);

	}
    }


  for(int k=0; k<4; k++)
    {
      for(int i=0; i<2; i++)
	{
	  hNEvents[k+1][i] = (TH1F*)hNEvents[0][i]->Clone(Form("NEvents_%s%s",setupName[k+1],wName[i]));
	  hNEvents[k+1][i]->Reset();
	  for(int j=0; j<nCentBins; j++)
	    {
	      if(i==0)
		{
		  hTpcVz[k+1][j] = (TH1F*)hTpcVz[0][j]->Clone(Form("TpcVz_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  hTpcVz[k+1][j]->Reset();

		  hTpcVr[k+1][j] = (TH1F*)hTpcVr[0][j]->Clone(Form("TpcVr_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  hTpcVr[k+1][j]->Reset();

		  hDiffVz[k+1][j] = (TH1F*)hDiffVz[0][j]->Clone(Form("TpcVpdDz_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  hDiffVz[k+1][j]->Reset();
		}
	      hNEventsAll[k+1][j][i] = (TH1F*)hNEventsAll[0][j][i]->Clone(Form("NEvents_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsAll[k+1][j][i]->Reset();

	      hNEventsVr[k+1][j][i] = (TH1F*)hNEventsVr[0][j][i]->Clone(Form("NEvents_VrCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsVr[k+1][j][i]->Reset();

	      hNEventsVz[k+1][j][i] = (TH1F*)hNEventsVz[0][j][i]->Clone(Form("NEvents_VrVzCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsVz[k+1][j][i]->Reset();

	      hNEventsAcc[k+1][j][i] = (TH1F*)hNEventsAcc[0][j][i]->Clone(Form("NEvents_VrVzDzCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsAcc[k+1][j][i]->Reset();

	      if(j==0)
		{
		  hCent[k+1][i] = (TH1F*)hCent[0][i]->Clone(Form("Cent_%s%s",setupName[k+1],wName[i]));
		  hCent[k+1][i]->Reset();
		}
	    }
	}
    }

  for(int k=0; k<4; k++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s_2014.list",run_type,trgSetupName[k]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int bin = hNEvents[0][0]->FindFixBin(runnumber);

	  for(int i=0; i<2; i++)
	    {
	      hNEvents[k+1][i]->SetBinContent(bin,hNEvents[0][i]->GetBinContent(bin));
	      hNEvents[k+1][i]->SetBinError(bin,hNEvents[0][i]->GetBinError(bin));
	      for(int j=0; j<nCentBins; j++)
		{
		  if(i==0)
		    {
		      TH1F *htmp = (TH1F*)hTpcVzVsRun[j]->ProjectionY(Form("%s_%d",hTpcVzVsRun[j]->GetName(),bin),bin,bin);
		      hTpcVz[k+1][j]->Add(htmp);

		      htmp = (TH1F*)hTpcVrVsRun[j]->ProjectionY(Form("%s_%d",hTpcVrVsRun[j]->GetName(),bin),bin,bin);
		      hTpcVr[k+1][j]->Add(htmp);

		      htmp = (TH1F*)hDiffVzVsRun[j]->ProjectionY(Form("%s_%d",hDiffVzVsRun[j]->GetName(),bin),bin,bin);
		      hDiffVz[k+1][j]->Add(htmp);
		    }
	  
		  hNEventsAll[k+1][j][i]->SetBinContent(bin,hNEventsAll[0][j][i]->GetBinContent(bin));
		  hNEventsAll[k+1][j][i]->SetBinError(bin,hNEventsAll[0][j][i]->GetBinError(bin));

		  hNEventsVr[k+1][j][i]->SetBinContent(bin,hNEventsVr[0][j][i]->GetBinContent(bin));
		  hNEventsVr[k+1][j][i]->SetBinError(bin,hNEventsVr[0][j][i]->GetBinError(bin));

		  hNEventsVz[k+1][j][i]->SetBinContent(bin,hNEventsVz[0][j][i]->GetBinContent(bin));
		  hNEventsVz[k+1][j][i]->SetBinError(bin,hNEventsVz[0][j][i]->GetBinError(bin));

		  hNEventsAcc[k+1][j][i]->SetBinContent(bin,hNEventsAcc[0][j][i]->GetBinContent(bin));
		  hNEventsAcc[k+1][j][i]->SetBinError(bin,hNEventsAcc[0][j][i]->GetBinError(bin));
		  if(j==0)
		    {
		      TH1F *htmp = (TH1F*)hCentVsRun[i]->ProjectionY(Form("%s_%d",hCentVsRun[i]->GetName(),bin),bin,bin);
		      hCent[k+1][i]->Add(htmp);
		    }
		}
	    }
	}
    }


  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
      for(int k=0; k<5; k++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hNEvents[k][i]->SetTitle("");
	      hNEvents[k][i]->Write("",TObject::kOverwrite);
	      for(int j=0; j<nCentBins; j++)
		{
		  hNEventsAll[k][j][i]->SetTitle("");
		  hNEventsVr[k][j][i]->SetTitle("");
		  hNEventsVz[k][j][i]->SetTitle("");
		  hNEventsAcc[k][j][i]->SetTitle("");
		  hNEventsAll[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsVr[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsVz[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsAcc[k][j][i]->Write("",TObject::kOverwrite);
		  if(j==0)
		    {
		      hCent[k][i]->SetTitle("");
		      hCent[k][i]->Write("",TObject::kOverwrite);
		    }
		  if(i==0)
		    {
		      hTpcVz[k][j]->SetTitle("");
		      hTpcVr[k][j]->SetTitle("");
		      hDiffVz[k][j]->SetTitle("");
		      hTpcVz[k][j]->Write("",TObject::kOverwrite);
		      hTpcVr[k][j]->Write("",TObject::kOverwrite);
		      hDiffVz[k][j]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}
