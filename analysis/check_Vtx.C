const int year = YEAR;

//================================================
void check_Vtx()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //makeMbVz();
  anaMbVz();
  
  //makeDmVz();
  //anaDmVz();

  //JpsiEffVsVz();
  //compDmMb();
  //prod_mid();
}

//================================================
void prod_mid(const int mode = 1, const int savePlot = 0)
{
  if(mode==0)
    {
      TFile *fin = TFile::Open("output/Run14_AuAu200.Study.Vtx.root", "update");
      THnSparseF* hnVertex = (THnSparseF*)fin->Get("mhnVertex_di_mu");
      hnVertex->GetAxis(4)->SetRange(3,3);
      hnVertex->GetAxis(3)->SetRange(1,16);
      
      TH2F *hTpcVzVsRun = (TH2F*)hnVertex->Projection(1,0);
      hTpcVzVsRun->SetName("hTpcVzVsRun_di_mu");
      hTpcVzVsRun->Write("",TObject::kOverwrite);
    }
  else
    {
      // TFile *fin = TFile::Open("output/Run14_AuAu200.Study.Vtx.root", "read");
      // TH2F* hTpcVzVsRun = (TH2F*)fin->Get("hTpcVzVsRun_di_mu");

      TFile *fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root", "read");
      TH2F* hTpcVzVsRun = (TH2F*)fin->Get("hTpcVzVsRun_cent0080__prod_mid");
      hTpcVzVsRun->GetYaxis()->SetRangeUser(-40,40);
      c = draw2D(hTpcVzVsRun);
      TProfile *hpro = (TProfile*)hTpcVzVsRun->ProfileX("hpro");
      hpro->SetMarkerStyle(20);
      hpro->Draw("sames");

      hTpcVzVsRun->GetYaxis()->SetRangeUser(-200, 200);
      TCanvas *cDay[5];
      for(int i=0; i<5; i++)
	{
	  cDay[i] = new TCanvas(Form("cTpcVz_%d",i), Form("cTpcVz_%d",i), 1100, 700);
	  cDay[i]->Divide(4,4);	}
      TH1F *hTpcVz[75];
      for(int i=0; i<75; i++)
	{
	  int bin1 = hTpcVzVsRun->GetXaxis()->FindFixBin(15000000+(94+i)*1000);
	  int bin2 = hTpcVzVsRun->GetXaxis()->FindFixBin(15000000+(95+i)*1000);
	  hTpcVz[i] = (TH1F*)hTpcVzVsRun->ProjectionY(Form("hTpcVz_day%d",94+i), bin1, bin2);
	  cDay[i/16]->cd(i%16+1);
	  hTpcVz[i]->Draw("HIST");
	  TPaveText *t1 = GetTitleText(Form("Day %d",94+i), 0.075);
	  t1->Draw();
	}
      for(int i=0; i<5; i++)
	{
	  if(savePlot) cDay[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_TpcVz_prod_mid_%d.pdf",run_type,i));
	}
    }
}

//================================================
void JpsiEffVsVz(const int savePlot = 0)
{
  const double jpsi_rapidity = 0.5;
  const int nEffType = 4;
  const char *trkEffType[nEffType] = {"MC","Tpc","MtdMth","TrigUnit"};
  const char *trkEffTitle[nEffType] = {"MC","TPC","MTD","Total"};
  THnSparseF *hnJpsiInfo[nEffType];
  TH2F *hJpsiPtVsTpcVz[nEffType];
  TFile *fin = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root","read");
  for(int i=0; i<nEffType; i++)
    {
      hnJpsiInfo[i] = (THnSparseF*)fin->Get(Form("hJpsiInfo_%s_di_mu",trkEffType[i]));
      hnJpsiInfo[i]->GetAxis(0)->SetRangeUser(2.8, 3.4); // cut on jpsi mass
      hnJpsiInfo[i]->GetAxis(2)->SetRangeUser(-1*jpsi_rapidity+0.01, jpsi_rapidity-0.01); // cut on jpsi rapidity
      hnJpsiInfo[i]->GetAxis(5)->SetRange(1, 16); // 0-80% centrality
      if(i>0)
	{
	  hnJpsiInfo[i]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	  hnJpsiInfo[i]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	}
      hJpsiPtVsTpcVz[i] = (TH2F*)hnJpsiInfo[i]->Projection(1, 7);
      hJpsiPtVsTpcVz[i]->SetName(Form("hJpsiPtVsTpcVz_%s",trkEffType[i]));
      hJpsiPtVsTpcVz[i]->RebinX(5);
    }
  TH2F *hJpsiEffVsTpcVz[nEffType-1];
  for(int i=0; i<nEffType-1; i++)
    {
      hJpsiEffVsTpcVz[i] = (TH2F*)hJpsiPtVsTpcVz[i+1]->Clone(Form("hJpsiEffVsTpcVz_%s",trkEffType[i]));
      hJpsiEffVsTpcVz[i]->GetYaxis()->SetRangeUser(0,10);
      hJpsiEffVsTpcVz[i]->Divide(hJpsiPtVsTpcVz[0]);
      c = draw2D(hJpsiEffVsTpcVz[i],Form("Run14_AuAu200: J/psi efficiency vs. vertex z (%s, |#eta| < 0.5, 0-80%%);vz_{TPC} (cm);Efficiency",trkEffTitle[i+1]),0.04,false);
      gPad->SetRightMargin(0.13);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_JpsiEffVsTpcVz_%s.pdf",run_type,trkEffTitle[i+1]));
    }

  // calculate inclusive J/psi efficiency vs. Vz using true J/psi 
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  TH1F *hModel = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  TH1F *hWeight = (TH1F*)hJpsiEffVsTpcVz[0]->ProjectionY("hWeight");
  hWeight->Reset("AC");
  for(int bin=1; bin<=hWeight->GetNbinsX(); bin++)
    {
      double low_pt  = hWeight->GetXaxis()->GetBinLowEdge(bin);
      double high_pt = hWeight->GetXaxis()->GetBinUpEdge(bin);
      double jpsi_yield = hModel->Integral(hModel->FindFixBin(low_pt+1e-5),hModel->FindFixBin(high_pt-1e-5));
      hWeight->SetBinContent(bin, jpsi_yield);
    }
  hWeight->Scale(1./hWeight->Integral());
  
  TList *list = new TList();
  TH1F *hJpsiEffInTpcVz[nEffType-1];
  for(int i=0; i<nEffType-1; i++)
    {
      hJpsiEffInTpcVz[i] = (TH1F*)hJpsiEffVsTpcVz[i]->ProjectionX(Form("hJpsiEffInTpcVz_%s",trkEffTitle[i+1]));
      hJpsiEffInTpcVz[i]->Reset("AC");
      for(int xbin=1; xbin<=hJpsiEffInTpcVz[i]->GetNbinsX(); xbin++)
	{
	  double eff = 0, err = 0;
	  for (int ybin=1; ybin<=hJpsiEffVsTpcVz[i]->GetNbinsY(); ybin++)
	    {
	      eff += hJpsiEffVsTpcVz[i]->GetBinContent(xbin, ybin) * hWeight->GetBinContent(ybin);
	      err += pow(hJpsiEffVsTpcVz[i]->GetBinError(xbin, ybin) * hWeight->GetBinContent(ybin), 2);
	    }
	  err = sqrt(err);
	  hJpsiEffInTpcVz[i]->SetBinContent(xbin, eff);
	  hJpsiEffInTpcVz[i]->SetBinError(xbin, err);
	}
      hJpsiEffInTpcVz[i]->SetMarkerStyle(20);
      hJpsiEffInTpcVz[i]->SetMarkerColor(pow(2, i));
      hJpsiEffInTpcVz[i]->SetLineColor(pow(2, i));
      if(i==0) hJpsiEffInTpcVz[i]->GetYaxis()->SetRangeUser(0.05, 0.12);
      if(i==1) hJpsiEffInTpcVz[i]->GetYaxis()->SetRangeUser(0.003, 0.012);
      if(i==2) hJpsiEffInTpcVz[i]->GetYaxis()->SetRangeUser(0, 0.005);
      c = draw1D(hJpsiEffInTpcVz[i], Form("Run14_AuAu200: J/psi efficiency (%s, p_{T} > 0, |#eta| < 0.5, 0-80%%);vz_{TPC} (cm);Efficiency",trkEffTitle[i+1]));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_JpsiEffInTpcVz_%s.pdf",run_type,trkEffTitle[i+1]));
    }

  // calculate efficiency for different vz distribution
  TFile *f = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  const int nData = 2;
  const char* dataName[nData] = {"dimuon", "MB"};
  const int gDataType = 0;
  TH1F *hVtxDis[gNTrgSetup][nData];
  for(int d=0; d<nData; d++)
    {
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  if(d==0)
	    {
	      hVtxDis[k+1][d] = (TH1F*)f->Get(Form("hTpcVz_cent%s%s_%s",cent_Title_pt[0],gTrgSetupTitle[k+1],dataName[d]));
	    }
	  if(d==1)
	    {
	      for(int r=0; r<8; r++)
		{
		  TH1F *htmp = (TH1F*)f->Get(Form("hTpcVz_cent%s_%s_run%d_%s",cent_Title_pt[0],gTrgSetupTitle[k+1],r,dataName[d]));
		  if(r==0) hVtxDis[k+1][d] = (TH1F*)htmp->Clone(Form("hTpcVz_cent%s%s_%s",cent_Title_pt[0],gTrgSetupTitle[k+1],dataName[d]));
		  else     hVtxDis[k+1][d]->Add(htmp);
		}
	    }
	  
	  if(k==0) hVtxDis[0][d] = (TH1F*)hVtxDis[k+1][d]->Clone(Form("hTpcVz_cent%s%s_%s",cent_Title_pt[0],gTrgSetupTitle[0],dataName[d]));
	  else     hVtxDis[0][d]->Add(hVtxDis[k+1][d]);
	  hVtxDis[k+1][d]->Scale(1./hVtxDis[k+1][d]->Integral(hVtxDis[k+1][d]->FindFixBin(-100+0.1),hVtxDis[k+1][d]->FindFixBin(100-0.1)));
	}
      hVtxDis[0][d]->Scale(1./hVtxDis[0][d]->Integral(hVtxDis[0][d]->FindFixBin(-100+0.1),hVtxDis[0][d]->FindFixBin(100-0.1)));
    }
  TString legName[5];
  for(int k=0; k<5; k++)
    {
      list->Add(hVtxDis[k][gDataType]);
      legName[k] = Form("%s%s",run_type,gTrgSetupTitle[k]);
    }
  c = drawHistos(list,Form("TpcVz_%s",dataName[gDataType]),Form("%s: TPC v_{z} distribution;v_{z} [cm]",dataName[gDataType]),false,0,0.35,true,0,0.018,false,kTRUE,legName,kTRUE,"",0.15,0.3,0.68,0.88,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_TpcVzComp_%s.pdf",run_type,dataName[gDataType]));

  double jpsiEff[gNTrgSetup], jpsiEff_err[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      jpsiEff[k] = 0;
      jpsiEff_err[k] = 0;
      int nbins = hJpsiEffInTpcVz[2]->GetNbinsX();
      double total = 0;
      for(int xbin=1; xbin<=nbins; xbin++)
	{
	  if(fabs(hJpsiEffInTpcVz[2]->GetBinCenter(xbin))>100) continue;
	  double xmin = hJpsiEffInTpcVz[2]->GetXaxis()->GetBinLowEdge(xbin);
	  double xmax = hJpsiEffInTpcVz[2]->GetXaxis()->GetBinUpEdge(xbin);
	  int bin1 = hVtxDis[k][gDataType]->FindFixBin(xmin+0.1);
	  int bin2 = hVtxDis[k][gDataType]->FindFixBin(xmax-0.1);
	  double weight = hVtxDis[k][gDataType]->Integral(bin1, bin2);
	  jpsiEff[k] += hJpsiEffInTpcVz[2]->GetBinContent(xbin) * weight;
	  jpsiEff_err[k] += pow(hJpsiEffInTpcVz[2]->GetBinError(xbin)*weight, 2);
	  total += weight;
	}
      cout << total << endl;
      jpsiEff_err[k] = sqrt(jpsiEff_err[k]);
      printf("[i] Vertex dependent efficiency for %s%s: %4.4f +/- %4.4f %%\n",run_type,gTrgSetupTitle[k],jpsiEff[k]*100,jpsiEff_err[k]*100);
    }

}

//================================================
void compDmMb(const int savePlot = 0)
{
  TFile *fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");

  const int nData = 2;
  const char* dataName[nData] = {"dimuon", "MB"};

  const int nHistos = 2;
  const char* histoName[nHistos] = {"DiffVz", "TpcVz"};
  const char* histoTitle[nHistos] = {"TPC-VPD v_{z}", "TPC v_{z}"};
  const char* legTitle[nHistos] = {"|#Deltav_{z}| < 3 cm", "|v_{z}| < 100 cm"};

  TH1F *hVtxDis[4][nHistos][nData];
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  for(int d=0; d<nData; d++)
	    {
	      if(d==0)
		{
		  hVtxDis[k][j][d] = (TH1F*)fin->Get(Form("h%s_cent%s%s_%s",histoName[j],cent_Title_pt[0],gTrgSetupTitle[k+1],dataName[d]));
		}
	      if(d==1)
		{
		  for(int r=0; r<8; r++)
		    {
		      TH1F *htmp = (TH1F*)fin->Get(Form("h%s_cent%s_%s_run%d_%s",histoName[j],cent_Title_pt[0],gTrgSetupTitle[k+1],r,dataName[d]));
		      if(r==0) hVtxDis[k][j][d] = (TH1F*)htmp->Clone(Form("h%s_cent%s%s_%s",histoName[j],cent_Title_pt[0],gTrgSetupTitle[k+1],dataName[d]));
		      else     hVtxDis[k][j][d]->Add(htmp);
		    }
		}
	      hVtxDis[k][j][d]->Scale(1./hVtxDis[k][j][d]->GetBinContent(hVtxDis[k][j][d]->FindBin(0)));
	    }
	}
    }

  TCanvas *cData[nHistos];
  TString legNameData[nData] = {"dimuon", "VPD-ZDC-no-vtx"};
  for(int j=0; j<nHistos; j++)
    {
      cData[j] = new TCanvas(Form("cData_%s",histoName[j]), Form("cData_%s",histoName[j]), 1100, 700);
      cData[j]->Divide(2, 2);
      TLegend *leg =  0x0;
      if(j==0) leg = new TLegend(0.15, 0.6, 0.3, 0.85);
      if(j==1) leg = new TLegend(0.6, 0.6, 0.8, 0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  cData[j]->cd(k+1);
	  if(j==0) gPad->SetLogy();
	  for(int d=0; d<nData; d++)
	    {
	      TH1F *hplot = (TH1F*)hVtxDis[k][j][d]->Clone(Form("%s_clone",hVtxDis[k][j][d]->GetName()));
	      hplot->SetMarkerStyle(20+2*d);
	      hplot->SetMarkerColor(color[d]);
	      hplot->SetLineColor(color[d]);
	      hplot->SetMarkerSize(1.0);
	      hplot->SetTitle("");
	      if(j==0) 
		{
		  hplot->GetYaxis()->SetRangeUser(1e-7,100);
		}
	      if(j==1) 
		{
		  // hplot->GetXaxis()->SetRangeUser(0,1.5);
		  hplot->GetYaxis()->SetRangeUser(0, 1.3);
		}
	      if(d==0) hplot->DrawCopy("P");
	      else     hplot->DrawCopy("samesP");
	      if(k==0) leg->AddEntry(hplot, legNameData[d].Data(), "P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s%s (0-80%%)",run_type,gTrgSetupTitle[k+1]), 0.055);
	  t1->Draw();
	}
      cData[j]->cd(1);
      leg->Draw();
      if(savePlot) cData[j]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_%sVsData.pdf",run_type,histoName[j]));
    } 
}


//================================================
void anaDmVz(const int savePlot = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  TFile *fin = TFile::Open("Rootfiles/Run14_AuAu200.CheckVtx.root","read");
  const int nHistos = 2;
  const char* histoName[nHistos]  = {"TpcVz", "DiffVz"};

  // check different run periods
  const int nRunRanges = 7;
  const char* runRangeName[nRunRanges] = {"Day 074-084", "Day 084-085", "Day 086-100", "15100100-15100103", "Day 100-121", "Day 122-162", "Day 163-167"};
  TH1F *hVtxDis[4][nHistos][nRunRanges];
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  for(int r=0; r<nRunRanges; r++)
	    {
	      hVtxDis[k][j][r] = (TH1F*)fin->Get(Form("DM_h%sDis_run%d_cent%s%s",histoName[j],r,cent_Title_pt[0],gTrgSetupTitle[k+1]));
	    }
	}
    }

  TCanvas *cVtx[nHistos][2];
  for(int j=0; j<nHistos; j++)
    {
      for(int i=0; i<2; i++)
	{
	  cVtx[j][i] = new TCanvas(Form("c%s_%d",histoName[j],i), Form("c%s_%d",histoName[j],i), 1100, 700);
	  cVtx[j][i]->Divide(2,2);
	}
      for(int r=0; r<nRunRanges; r++)
	{
	  cVtx[j][r/4]->cd(r%4+1);
	  if(j==1) gPad->SetLogy();
	  TLegend *leg =  new TLegend(0.55, 0.6, 0.75, 0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  for(int k=0; k<4; k++)
	    {
	      if(hVtxDis[k][j][r]->Integral()>0) hVtxDis[k][j][r]->Scale(1./hVtxDis[k][j][r]->Integral());
	      hVtxDis[k][j][r]->SetMarkerStyle(20+k);
	      hVtxDis[k][j][r]->SetMarkerColor(color[k]);
	      hVtxDis[k][j][r]->SetLineColor(color[k]);
	      if(j==0) 
		{
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-5, 0.02);
		}
	      else if(j==1) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(-10, 10);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-6, 10);
		}


	      if(k==0) hVtxDis[k][j][r]->Draw();
	      else     hVtxDis[k][j][r]->Draw("sames");
	      if(hVtxDis[k][j][r]->GetEntries()>0)
		{
		  int bin1 = -1, bin2 = -1;
		  if(j==0)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(-100+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(100-1e-4);
		    }
		  else if(j==1)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(-3+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(3-1e-4);
		    }

		  double eff = hVtxDis[k][j][r]->Integral(bin1, bin2)/hVtxDis[k][j][r]->Integral(0,-1);
		  leg->AddEntry(hVtxDis[k][j][r], Form("%s: f = %4.1f%%",gTrgSetupTitle[k+1], eff*100), "P");
		}
	    }
	  TPaveText *t1 = GetTitleText(Form("%s MB trigger (%s)",run_type,runRangeName[r]),0.05);
	  t1->Draw();
	  leg->Draw();
	}

     for(int i=0; i<2; i++)
	{
	  if(savePlot) cVtx[j][i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_DM_%sDis_run%d.pdf",run_type,histoName[j],i));
	}
    }
}

//================================================
void makeDmVz(const int saveHisto = 1)
{  
  // dimuon trigger
  TFile *fdm = TFile::Open("output/Run14_AuAu200.Study.Vtx.root", "read");
  THnSparseF *hnVtx = (THnSparseF*)fdm->Get("mhnVertex_di_mu");
  hnVtx->GetAxis(3)->SetRange(1, 16); // 0-80%;
  const int icent = 0;

  TH1F *hTpcVz[gNTrgSetup];
  TH1F *hDiffVz[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hnVtx->GetAxis(4)->SetRange(k, k);
      hTpcVz[k] = (TH1F*)hnVtx->Projection(1);
      hTpcVz[k]->Sumw2();
      hTpcVz[k]->SetName(Form("DM_hTpcVzDis_cent%s%s",cent_Title_pt[icent],gTrgSetupTitle[k]));

      hDiffVz[k] = (TH1F*)hnVtx->Projection(2);
      hDiffVz[k]->Sumw2();
      hDiffVz[k]->SetName(Form("DM_hDiffVzDis_cent%s%s",cent_Title_pt[icent],gTrgSetupTitle[k]));
    }
  hnVtx->GetAxis(4)->SetRange(0,-1);


  const int nHistos = 2;
  const char* histoName[nHistos]  = {"TpcVz", "DiffVz"};
  const int nRunRanges = 7;
  const int runRangeLow[nRunRanges]  = {15074000, 15084500, 15088500, 15100100, 15100105, 15121500, 15162500};
  const int runRangeHigh[nRunRanges] = {15084500, 15085500, 15100050, 15100103, 15121000, 15162500, 15167500};
  TH2F *hVtxDisVsRun[gNTrgSetup][nHistos];
  TH1F *hVtxDis[gNTrgSetup][nHistos][nRunRanges];
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hnVtx->GetAxis(4)->SetRange(k,k);
      for(int j=0; j<nHistos; j++)
	{
	  hVtxDisVsRun[k][j] = (TH2F*)hnVtx->Projection(j+1, 0);
	  hVtxDisVsRun[k][j]->SetName(Form("DM_h%sDisVsRun_cent%s%s",histoName[j],cent_Title_pt[icent],gTrgSetupTitle[k]));
	  for(int r=0; r<nRunRanges; r++)
	    {
	      int bin1 = hVtxDisVsRun[k][j]->GetXaxis()->FindBin(runRangeLow[r]);
	      int bin2 = hVtxDisVsRun[k][j]->GetXaxis()->FindBin(runRangeHigh[r]);
	      hVtxDis[k][j][r] = (TH1F*)hVtxDisVsRun[k][j]->ProjectionY(Form("DM_h%sDis_run%d_cent%s%s",histoName[j],r,cent_Title_pt[icent],gTrgSetupTitle[k]), bin1, bin2);
	    }
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.CheckVtx.root","update");
      for(int k=0; k<gNTrgSetup; k++)
	{
	  hTpcVz[k]->SetTitle("");
	  hTpcVz[k]->Write("",TObject::kOverwrite);

	  hDiffVz[k]->SetTitle("");
	  hDiffVz[k]->Write("",TObject::kOverwrite);
	  for(int j=0; j<nHistos; j++)
	    {
	      hVtxDisVsRun[k][j]->SetTitle("");
	      hVtxDisVsRun[k][j]->Write("",TObject::kOverwrite);
	      for(int r=0; r<nRunRanges; r++)
		{
		  hVtxDis[k][j][r]->SetTitle("");
		  hVtxDis[k][j][r]->Write("",TObject::kOverwrite);
		}
	    }
	}
      fout->Close();
    }
}


//================================================
void anaMbVz(const int savePlot = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  const int weight = 0;
  const char* gWeightName[2] = {"","Weight"};

  TFile *fin = TFile::Open("Rootfiles/Run14_AuAu200.CheckVtx.root","read");
  const int nHistos = 5;
  const char* histoName[nHistos]  = {"Cent", "TpcVr", "DiffVz", "TpcVz", "Acc"};
  const char* histoTitle[nHistos] = {"gRefMult", "TPC v_{r}", "TPC-VPD v_{z}", "TPC v_{z}", "centrality"};
  const char* cutTitle[nHistos]   = {"No cut", "0-80%", "0-80%, |v_{r}| < 2 cm", "0-80%, |v_{r}| < 2 cm, |#Deltav_{z}| < 3 cm", "0-80%, vtx cuts"};

  // vertex distribution in refMult
  TH1F *hVtxDisInRefMult[gNTrgSetup][3][2];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hVtxDisInRefMult[k][j][i] = (TH1F*)fin->Get(Form("MB_h%sDis_RefMult%d%s%s",histoName[j+1],i,gTrgSetupTitle[k],gWeightName[weight]));
	      hVtxDisInRefMult[k][j][i]->Scale(1./hVtxDisInRefMult[k][j][i]->Integral());
	      hVtxDisInRefMult[k][j][i]->SetLineColor(i+1);
	    }
	}
    }
  TCanvas *cVtxDisInRefMult[3];
  for(int j=0; j<3; j++)
    {
      cVtxDisInRefMult[j] = new TCanvas(Form("cVtxDisInRefMult_%s",histoName[j+1]), Form("cVtxDisInRefMult_%s",histoName[j+1]), 1100, 700);
      cVtxDisInRefMult[j]->Divide(2,2); 
      for(int k=1; k<gNTrgSetup; k++)
	{
	  cVtxDisInRefMult[j]->cd(k);
	  if(j<2) gPad->SetLogy();
	  hVtxDisInRefMult[k][j][0]->DrawCopy("HIST");
	  hVtxDisInRefMult[k][j][1]->DrawCopy("samesHIST");
	}
    }

  TCanvas *cVtxDisVsLumiInRefMult[3];
  for(int j=0; j<3; j++)
    {
      cVtxDisVsLumiInRefMult[j] = new TCanvas(Form("cVtxDisVsLumiInRefMult_%s",histoName[j+1]), Form("cVtxDisVsLumiInRefMult_%s",histoName[j+1]), 800, 600);
      if(j<2) gPad->SetLogy();
      for(int k=1; k<gNTrgSetup; k++)
	{
	  hVtxDisInRefMult[k][j][0]->SetLineColor(color[k-1]);
	  if(k==1) hVtxDisInRefMult[k][j][0]->DrawCopy("HIST");
	  else     hVtxDisInRefMult[k][j][0]->DrawCopy("samesHIST");
	}
    }
  
  TH2F *hVtxDisVsRun[gNTrgSetup][nHistos];
  TH1F *hNevent[gNTrgSetup][nHistos];
  TH1F *hCutEff[gNTrgSetup][nHistos];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  hNevent[k][j] = (TH1F*)fin->Get(Form("MB_h%sVsRun%s_Nevts%s",histoName[j],gTrgSetupTitle[k],gWeightName[weight]));
	}
    }
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  if(j<nHistos-1)
	    {
	      hCutEff[k][j] = (TH1F*)hNevent[k][j+1]->Clone(Form("h%sVsRun%s_Eff",histoName[j],gTrgSetupTitle[k]));
	      hCutEff[k][j]->Divide(hNevent[k][j]);
	    }
	  else
	    {
	      hCutEff[k][j] = (TH1F*)hNevent[k][j]->Clone(Form("h%sVsRun%s_Eff",histoName[j],gTrgSetupTitle[k]));
	      hCutEff[k][j]->Divide(hNevent[k][0]);
	    }
	}
    }

  // check efficiencies
  TCanvas *cEff[nHistos];
  for(int j=0; j<nHistos; j++)
    {
      cEff[j] = new TCanvas(Form("cEff_%s",histoName[j]), Form("cEff_%s",histoName[j]), 800, 600);
      //cEff[j]->Divide(2,2);
      for(int k=1; k<gNTrgSetup; k++)
	{
	  hCutEff[k][j]->SetMarkerSize(1.5);
	  hCutEff[k][j]->SetMarkerStyle(19+k+k/4);
	  hCutEff[k][j]->SetMarkerColor(color[k-1]);
	  hCutEff[k][j]->SetLineColor(color[k-1]);
	  hCutEff[k][j]->GetYaxis()->SetRangeUser(0,1);
	  if(k==1) hCutEff[k][j]->Draw("P");
	  else     hCutEff[k][j]->Draw("samesP");
	}
    }

  // check different run periods
  //const char* legTitle[nHistos] = { "No cut", "0-80%", "0-80%, |v_{r}| < 2 cm", "0-80%, |v_{r}| < 2 cm, |#Deltav_{z}| < 3 cm"};
  const int nRunRanges = 7;
  const char* runRangeName[nRunRanges] = {"Day 074-084", "Day 084-085", "Day 086-100", "15100100-15100103", "Day 100-121", "Day 122-162", "Day 163-167"};
  TH1F *hVtxDis[4][nHistos][nRunRanges];
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  for(int r=0; r<nRunRanges; r++)
	    {
	      hVtxDis[k][j][r] = (TH1F*)fin->Get(Form("MB_h%sDis%s_run%d%s",histoName[j],gTrgSetupTitle[k+1],r,gWeightName[weight]));
	    }
	}
    }

  TCanvas *cVtx[nHistos][2];
  for(int j=0; j<nHistos; j++)
    {
      for(int i=0; i<2; i++)
	{
	  cVtx[j][i] = new TCanvas(Form("c%s_%d",histoName[j],i), Form("c%s_%d",histoName[j],i), 1100, 700);
	  cVtx[j][i]->Divide(2,2);
	}
      for(int r=0; r<nRunRanges; r++)
	{
	  cVtx[j][r/4]->cd(r%4+1);
	  if(j<3) gPad->SetLogy();
	  TLegend *leg =  new TLegend(0.55, 0.6, 0.75, 0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(cutTitle[j]);
	  for(int k=0; k<4; k++)
	    {
	      if(hVtxDis[k][j][r]->Integral()>0) hVtxDis[k][j][r]->Scale(1./hVtxDis[k][j][r]->Integral());
	      hVtxDis[k][j][r]->SetMarkerStyle(20+k);
	      hVtxDis[k][j][r]->SetMarkerColor(color[k]);
	      hVtxDis[k][j][r]->SetLineColor(color[k]);
	      if(j==0) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(0, 700);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-5, 0.5);
		}
	      else if(j==1) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(0, 5);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-6, 10);
		}
	      else if(j==2) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(-10, 10);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-5, 10);
		}
	      else if(j==3) 
		{
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-5, 0.02);
		}
	      else if(j==4) 
		{
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-5, 0.15);
		}

	      if(k==0) hVtxDis[k][j][r]->Draw();
	      else     hVtxDis[k][j][r]->Draw("sames");
	      if(hVtxDis[k][j][r]->GetEntries()>0)
		{
		  int bin1 = -1, bin2 = -1;
		  if(j==3)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(-100+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(100-1e-4);
		    }
		  else if(j==2)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(-3+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(3-1e-4);
		    }
		  else if(j==1)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(0+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(2-1e-4);
		    }
		  else if(j==0)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(10+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(1000-1e-4);
		    }
		  else if(j==4)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(0.5+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(16.5-1e-4);
		    }

		  double eff = hVtxDis[k][j][r]->Integral(bin1, bin2)/hVtxDis[k][j][r]->Integral(0,-1);
		  leg->AddEntry(hVtxDis[k][j][r], Form("%s: f = %4.1f%%",gTrgSetupTitle[k+1], eff*100), "P");
		}
	    }
	  TPaveText *t1 = GetTitleText(Form("%s MB trigger (%s)",run_type,runRangeName[r]),0.05);
	  t1->Draw();
	  leg->Draw();
	}

     for(int i=0; i<2; i++)
	{
	  if(savePlot) cVtx[j][i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Vtx_MB_%sDis_run%d.pdf",run_type,histoName[j],i));
	}
    }
}

//================================================
void makeMbVz(const int saveHisto = 1)
{
  const int weight = 0;
  const char* gWeightName[2] = {"","Weight"};

  const int nHistos = 5;
  const char* histoName[nHistos] = {"Cent", "TpcVr", "DiffVz", "TpcVz", "Acc"};

  // MB data
  TFile *fMB = TFile::Open(Form("./output/Run14_AuAu200.MB.VtxEff.root"),"read");
  THnSparseF *hn = (THnSparseF*)fMB->Get(Form("mhMbEvtEff%s",gWeightName[weight]));
  
  // vertex distributions for different gRefMultCorr ranges
  TH2F *hVtxDisVsRunInRefMult[gNTrgSetup][3][2];
  for(int i=0; i<2; i++)
    {
      if(i==0) hn->GetAxis(4)->SetRangeUser(10,1000);
      else     hn->GetAxis(4)->SetRangeUser(0, 10);
      for(int j=0; j<3; j++)
	{
	  hVtxDisVsRunInRefMult[0][j][i] = (TH2F*)hn->Projection(3-j, 0);
	  hVtxDisVsRunInRefMult[0][j][i]->SetTitle("");
	  hVtxDisVsRunInRefMult[0][j][i]->Sumw2();
	  hVtxDisVsRunInRefMult[0][j][i]->SetName(Form("MB_h%sDisVsRun_RefMult%d%s%s",histoName[j+1],i,gTrgSetupTitle[0],gWeightName[weight]));
	}
    }
  hn->GetAxis(4)->SetRange(0,-1);

  // 2D vs. run
  TH2F *hVtxDisVsRun[gNTrgSetup][nHistos];
  hVtxDisVsRun[0][0] = (TH2F*)hn->Projection(4, 0);
  hVtxDisVsRun[0][0]->RebinY(5);
  
  hn->GetAxis(5)->SetRange(1,16); // + 0-80% cut
  hVtxDisVsRun[0][1] = (TH2F*)hn->Projection(3, 0);

  hn->GetAxis(3)->SetRange(1,2); // apply vr < 2 cm cut
  hVtxDisVsRun[0][2] = (TH2F*)hn->Projection(2, 0);

  hn->GetAxis(2)->SetRangeUser(-3+1e-4,3-1e-4); // + |dz| < 3 cm cut
  hVtxDisVsRun[0][3] = (TH2F*)hn->Projection(1, 0);

  hn->GetAxis(1)->SetRangeUser(-100+1e-4,100-1e-4); // + |vz| < 100 cm cut
  hVtxDisVsRun[0][4] = (TH2F*)hn->Projection(5, 0);

  for(int j=0; j<nHistos; j++)
    {
      hVtxDisVsRun[0][j]->SetTitle("");
      hVtxDisVsRun[0][j]->Sumw2();
      hVtxDisVsRun[0][j]->SetName(Form("MB_h%sVsRun%s%s",histoName[j],gTrgSetupTitle[0],gWeightName[weight]));
      printf("[i] %s has %1.0f entries\n",hVtxDisVsRun[0][j]->GetName(),hVtxDisVsRun[0][j]->Integral(0,-1));
    }
  for(int j=1; j<=5; j++) hn->GetAxis(j)->SetRange(0,-1);

  // different luminosity configurations
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  hVtxDisVsRun[k+1][j] = (TH2F*)hVtxDisVsRun[0][j]->Clone(Form("MB_h%sVsRun%s%s",histoName[j],gTrgSetupTitle[k+1],gWeightName[weight]));
	  hVtxDisVsRun[k+1][j]->Reset();
	}

      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hVtxDisVsRunInRefMult[k+1][j][i] = (TH2F*)hVtxDisVsRunInRefMult[0][j][i]->Clone(Form("MB_h%sDisVsRun_RefMult%d%s%s",histoName[j+1],i,gTrgSetupTitle[k+1],gWeightName[weight]));
	      hVtxDisVsRunInRefMult[k+1][j][i]->Reset();
	    }
	}
    }

  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  for(int k=0; k<4; k++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/Run14_AuAu200/AuAu_200_%s_2014.list",trgSetupName[k]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  if(runnumber==15121062 || runnumber==15119042 || runnumber==15102024 || runnumber==15121060) continue;
	  int xbin = hVtxDisVsRun[0][0]->GetXaxis()->FindFixBin(runnumber);

	  for(int j=0; j<nHistos; j++)
	    {
	      int nybins = hVtxDisVsRun[0][j]->GetNbinsY();
	      for(int ybin=0; ybin<=nybins+1; ybin++)
		{
		  hVtxDisVsRun[k+1][j]->SetBinContent(xbin, ybin, hVtxDisVsRun[0][j]->GetBinContent(xbin, ybin));
		  hVtxDisVsRun[k+1][j]->SetBinError(xbin, ybin, hVtxDisVsRun[0][j]->GetBinError(xbin, ybin));
		}
	    }

	  for(int j=0; j<3; j++)
	    {
	      for(int i=0; i<2; i++)
		{
		  int nybins = hVtxDisVsRunInRefMult[0][j][i]->GetNbinsY();
		  for(int ybin=0; ybin<=nybins+1; ybin++)
		    {
		      hVtxDisVsRunInRefMult[k+1][j][i]->SetBinContent(xbin, ybin, hVtxDisVsRunInRefMult[0][j][i]->GetBinContent(xbin, ybin));
		      hVtxDisVsRunInRefMult[k+1][j][i]->SetBinError(xbin, ybin, hVtxDisVsRunInRefMult[0][j][i]->GetBinError(xbin, ybin));
		    }
		}
	    }
	}
    }

  // efficiencies
  TH1F *hNevent[gNTrgSetup][nHistos];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  hNevent[k][j] = (TH1F*)hVtxDisVsRun[k][j]->ProjectionX(Form("MB_h%sVsRun%s_Nevts%s",histoName[j],gTrgSetupTitle[k],gWeightName[weight]),0,-1);
	}
    }

  TH1F *hVtxDisInRefMult[gNTrgSetup][3][2];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hVtxDisInRefMult[k][j][i] = (TH1F*)hVtxDisVsRunInRefMult[k][j][i]->ProjectionY(Form("MB_h%sDis_RefMult%d%s%s",histoName[j+1],i,gTrgSetupTitle[k],gWeightName[weight]));
	    }
	}
    }

  const int nRunRanges = 7;
  const int runRangeLow[nRunRanges]  = {15074000, 15084500, 15088500, 15100100, 15100105, 15121500, 15162500};
  const int runRangeHigh[nRunRanges] = {15084500, 15085500, 15100050, 15100103, 15121000, 15162500, 15167500};
  TH1F *hVtxDis[gNTrgSetup][nHistos][nRunRanges];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  for(int r=0; r<nRunRanges; r++)
	    {
	      int bin1 = hVtxDisVsRun[k][j]->GetXaxis()->FindBin(runRangeLow[r]);
	      int bin2 = hVtxDisVsRun[k][j]->GetXaxis()->FindBin(runRangeHigh[r]);
	      hVtxDis[k][j][r] = (TH1F*)hVtxDisVsRun[k][j]->ProjectionY(Form("MB_h%sDis%s_run%d%s",histoName[j],gTrgSetupTitle[k],r,gWeightName[weight]), bin1, bin2);
	      //printf("[i] %s has underflow bin %1.0f\n",hVtxDis[k][j][r]->GetName(),hVtxDis[k][j][r]->GetBinContent(0));
	    }
	}
    }


  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.CheckVtx.root","update");
      for(int k=0; k<gNTrgSetup; k++)
	{
	  for(int j=0; j<nHistos; j++)
	    {
	      hVtxDisVsRun[k][j]->Write("",TObject::kOverwrite);
	      hNevent[k][j]->Write("",TObject::kOverwrite);
	      // hCutEff[k][j]->Write("",TObject::kOverwrite);
	      for(int r=0; r<nRunRanges; r++)
	      	{
	      	  hVtxDis[k][j][r]->Write("",TObject::kOverwrite);
	      	}
	    }

	  for(int j=0; j<3; j++)
	    {
	      for(int i=0; i<2; i++)
		{
		  hVtxDisInRefMult[k][j][i]->Write("",TObject::kOverwrite);
		}
	    }
	}
      fout->Close();
    }
}
