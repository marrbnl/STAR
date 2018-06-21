//================================================
void comp_All()
{
  gStyle->SetOptStat(0);
  All();
  //SignalShape();
}

//================================================
void All(const int compCount = 1, const int compEff = 1, const int compRef = 1, const int savePlot = 1)
{
  const char* dataType[2] = {"New","Old"};
  // Jpsi efficiency vs. pT
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** ptName       = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;
  const int kNCent = nCentBins_npart[0];

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TList *list = new TList;

  if(compCount)
    {
      TFile *fSig[2];
      fSig[0] = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
      fSig[1] = TFile::Open(Form("Rootfiles/old.%s.JpsiYield.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
      TH1F *hJpsiCounts[2][7];
      TCanvas *c = new TCanvas("comp_JpsiCounts", "comp_JpsiCounts", 1100, 500);
      c->Divide(4,2);
      TLegend *leg = new TLegend(0.3,0.3,0.6,0.6);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(62);
      leg->SetTextSize(0.08);
      for(int k=0; k<7; k++)
	{
	  for(int j=0; j<2; j++)
	    {
	      if(k<5) hJpsiCounts[j][k] = (TH1F*)fSig[j]->Get(Form("Jpsi_FitYield_cent%s_weight",cent_Title[k]));
	      else    hJpsiCounts[j][k] = (TH1F*)fSig[j]->Get(Form("Jpsi_FitYield_pt%s_weight",pt_Name_npart[k-5]));
	      hJpsiCounts[j][k]->SetName(Form("%s_%d",hJpsiCounts[j][k]->GetName(),j));
	      hJpsiCounts[j][k]->SetMarkerStyle(21+j*4);
	      hJpsiCounts[j][k]->SetMarkerColor(j+1);
	      hJpsiCounts[j][k]->SetLineColor(j+1);

	      c->cd(k+1);
	      gPad->SetLogy();
	      if(k==2 || k==3) hJpsiCounts[j][k]->GetXaxis()->SetRangeUser(0.5,10);
	      if(k==4) hJpsiCounts[j][k]->GetXaxis()->SetRangeUser(0.5,6);
	      hJpsiCounts[j][k]->SetTitle(";;Counts");
	      if(k<5) hJpsiCounts[j][k]->SetXTitle("p_{T} [GeV/c]");
	      if(j==0) hJpsiCounts[j][k]->Draw();
	      else     hJpsiCounts[j][k]->Draw("sames");
	      if(k==0) leg->AddEntry(hJpsiCounts[j][k], dataType[j], "P");
	    }
	  if(k<5) TPaveText *t1 = GetTitleText(Form("%s: %s%%",run_type,cent_Name[k]),0.06);
	  else    TPaveText *t1 = GetTitleText(Form("%s: %s",run_type,pt_Title_npart[k-5]),0.06);
	  t1->Draw();
	}
      c->cd(8);
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_CompAll/Compare_JpsiCounts.pdf",run_type));
    }
  if(compEff)
    {
      const char *trkEffType[6] = {"MC","Tpc","MtdMth","MuonPid","MtdTrig","TrigUnit"};
      TFile *fEff[2];
      fEff[0] = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
      fEff[1] = TFile::Open(Form("Rootfiles/old.%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");

      const int nHistos = 6;
      TH1F *hJpsiPt[2][nHistos][nCentBins];
      TH1F *hJpsiPtEffs[2][nHistos][nCentBins];
      for(int j=0; j<2; j++)
	{
	  for(int i=0; i<nHistos; i++)
	    {
	      for(int k=0; k<nCentBins; k++)
		{
		  hJpsiPt[j][i][k] = (TH1F*)fEff[j]->Get(Form("hJpsiPt_%s_cent%s",trkEffType[i],cent_Title[k]));
		  hJpsiPt[j][i][k]->SetName(Form("%s_file%d",hJpsiPt[j][i][k]->GetName(),j));
		  hJpsiPt[j][i][k]->Rebin(4);
		  int index = i-1;
		  if(i==0) index = 0;
		  hJpsiPtEffs[j][i][k] = DivideTH1ForEff(hJpsiPt[j][i][k],hJpsiPt[j][index][k],Form("hJpsiPtEff_%s_cent%s_file%d",trkEffType[i],cent_Title[k],j));
		}
	    }
	}

      // various efficiency
      const int kcent = 0;
      for(int i=1; i<nHistos; i++)
	{
	  hJpsiPtEffs[0][i][kcent]->Divide(hJpsiPtEffs[1][i][kcent]);
	  list->Add(hJpsiPtEffs[0][i][kcent]);
	}
      TString legName2[5] = {"TPC tracking + p_{T,#mu} cut","MTD acceptance & response","Muon PID","MTD triggering","Trigger unit"};
      c = drawHistos(list,"JpsiEff_AllEffs",Form("%s: efficiencies for J/#psi ;p_{T} (GeV/c);New/Old",run_type),true,0,15,true,0.8,1.2,false,kTRUE,legName2,true,Form("%s%%",cent_Name[kcent]),0.2,0.4,0.63,0.88,kTRUE,0.04,0.035);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_CompAll/Compare_JpsiEff_AllTypes.pdf",run_type));
      list->Clear();
    }

  if(compRef)
    {
      TFile *fpp[2];
      fpp[0] = TFile::Open(Form("Rootfiles/Paper.%s.Jpsi.root",run_type),"read");
      fpp[1] = TFile::Open(Form("Rootfiles/Comb2.Paper.%s.Jpsi.root",run_type),"read");
      TGraphAsymmErrors	*hppJpsiVsPt[2];
      TGraphAsymmErrors	*hppJpsiVsCent[2];
      double x, y, x1, y1;
      for(int j=0; j<2; j++)
	{
	  hppJpsiVsPt[j]   = (TGraphAsymmErrors*)fpp[j]->Get("hpp200JpsiVsPtFinalSys");
	  hppJpsiVsPt[j]->SetName(Form("%s_%d",hppJpsiVsPt[j]->GetName(),j));
	  hppJpsiVsPt[j]->SetMarkerStyle(21+j*4);
	  hppJpsiVsPt[j]->SetMarkerColor(j+1);
	  hppJpsiVsPt[j]->SetLineColor(j+1);
	  offset_x(hppJpsiVsPt[j], 0.1*j);
	  hppJpsiVsCent[j] = (TGraphAsymmErrors*)fpp[j]->Get("hpp200JpsiVsCentFinalSys");
	  hppJpsiVsCent[j]->SetName(Form("%s_%d",hppJpsiVsCent[j]->GetName(),j));
	  hppJpsiVsCent[j]->SetMarkerStyle(21+j*4);
	  hppJpsiVsCent[j]->SetMarkerColor(j+1);
	  hppJpsiVsCent[j]->SetLineColor(j+1);
	  offset_x(hppJpsiVsCent[j], 0.1*j);
	}

      for(int j=0; j<2; j++)
	{
	  for(int ipoint=0; ipoint<hppJpsiVsPt[j]->GetN(); ipoint++)
	    {
	      hppJpsiVsPt[1]->GetPoint(ipoint, x, y);
	      hppJpsiVsPt[j]->GetPoint(ipoint, x1, y1);
	      hppJpsiVsPt[j]->SetPoint(ipoint, x1, y1/y);
	      hppJpsiVsPt[j]->SetPointError(ipoint, hppJpsiVsPt[j]->GetErrorXlow(ipoint), hppJpsiVsPt[j]->GetErrorXhigh(ipoint),
					    hppJpsiVsPt[j]->GetErrorYlow(ipoint)/y, hppJpsiVsPt[j]->GetErrorYhigh(ipoint)/y);
	    }

	  for(int ipoint=0; ipoint<hppJpsiVsCent[j]->GetN(); ipoint++)
	    {
	      hppJpsiVsCent[1]->GetPoint(ipoint, x, y);
	      hppJpsiVsCent[j]->GetPoint(ipoint, x1, y1);
	      hppJpsiVsCent[j]->SetPoint(ipoint, x1, y1/y);
	      hppJpsiVsCent[j]->SetPointError(ipoint, hppJpsiVsCent[j]->GetErrorXlow(ipoint), hppJpsiVsCent[j]->GetErrorXhigh(ipoint),
					    hppJpsiVsCent[j]->GetErrorYlow(ipoint)/y, hppJpsiVsCent[j]->GetErrorYhigh(ipoint)/y);
	    }
	}
      hppJpsiVsPt[1]->GetYaxis()->SetRangeUser(0.5,1.5);
      c = drawGraph(hppJpsiVsPt[1],"Ratio of pp reference;p_{T} [GeV/c];Ratio");
      hppJpsiVsPt[0]->Draw("samesPEZ");
      TLegend *leg = new TLegend(0.6,0.7,0.8,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(62);
      leg->SetTextSize(0.035);
      leg->AddEntry(hppJpsiVsPt[0], "(STAR+PHENIX)/STAR", "P");
      leg->AddEntry(hppJpsiVsPt[1], "STAR/STAR", "P");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_CompAll/Compare_ppRefVsPt.pdf",run_type));

      hppJpsiVsCent[1]->GetXaxis()->SetRangeUser(-0.5,6.5);
      hppJpsiVsCent[1]->GetYaxis()->SetRangeUser(0.5,1.5);
      c = drawGraph(hppJpsiVsCent[1],"Ratio of pp reference;;Ratio");
      hppJpsiVsCent[0]->Draw("samesPEZ");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_CompAll/Compare_ppRefVsCent.pdf",run_type));
    }
}


//================================================
void SignalShape(const int savePlot = 1)
{
  const int icent = 0;
  const int nfiles = 2;
  const char *filename[nfiles] = {"Pico.Run13.pp500.jpsi.pt1.5.pt1.0.yield.root","Pico.Run13.pp500.jpsi.VtxCut.pt1.5.pt1.0.yield.root"};
  const TString legName[nfiles] = {"W/o vtx cut","W/ vtx cut"};
  const char *saveTitle = "Run13_VtxCut";

  TFile *file[nfiles];
  TH1F *hInvMass[nfiles][nPtBins];
  TH1F *hMean[nfiles];
  TH1F *hSigma[nfiles];  
  TH1F *hYield[nfiles];

  for(int i=0; i<nfiles; i++)
    {
      file[i] = TFile::Open(Form("Rootfiles/%s",filename[i]),"read");

      hMean[i] = (TH1F*)file[i]->Get(Form("Jpsi_FitMean_cent%s",cent_Title[icent]));
      hMean[i]->SetName(Form("%s_%d",hMean[i]->GetName(),i));

      hSigma[i] = (TH1F*)file[i]->Get(Form("Jpsi_FitSigma_cent%s",cent_Title[icent]));
      hSigma[i]->SetName(Form("%s_%d",hSigma[i]->GetName(),i));

      hYield[i] = (TH1F*)file[i]->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]));
      hYield[i]->SetName(Form("%s_%d",hYield[i]->GetName(),i));

      for(int ipt=0; ipt<nPtBins; ipt++)
	{
	  hInvMass[i][ipt] = (TH1F*)file[i]->Get(Form("Jpsi_Signal_cent%s_pt%s_save",cent_Title[icent],pt_Name[ipt]));
	  hInvMass[i][ipt]->SetName(Form("%s_%d",hInvMass[i][ipt]->GetName(),i));
	}
    }

  // invariant mass distribution
  TCanvas *c = new TCanvas("InvMass","InvMass",1200,700);
  c->Divide(nPtBins/2+nPtBins%2,2);
  TLegend *leg = new TLegend(0.6,0.65,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int ipt=0; ipt<nPtBins; ipt++)
    {
      c->cd(ipt+1);
      SetPadMargin(gPad,0.13, 0.13, 0.03, 0.1);
      for(int i=0; i<nfiles; i++)
	{
	  int bin = hInvMass[i][ipt]->FindFixBin(3.09);
	  hInvMass[i][ipt]->Scale(1./hInvMass[i][ipt]->GetBinContent(bin));
	  hInvMass[i][ipt]->SetMaximum(1.5);
	  hInvMass[i][ipt]->SetMarkerStyle(21+i*4);
	  hInvMass[i][ipt]->SetMarkerColor(color[i]);
	  hInvMass[i][ipt]->SetLineColor(color[i]);
	  hInvMass[i][ipt]->SetYTitle("a.u.");
	  ScaleHistoTitle(hInvMass[i][ipt],0.05,1,0.035,0.05,1.4,0.035,62);
	  if(i==0) hInvMass[i][ipt]->Draw();
	  else     hInvMass[i][ipt]->Draw("sames");
	  if(ipt==0) leg->AddEntry(hInvMass[i][ipt],legName[i].Data(),"PL");
	}
      TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T,#mu#mu} < %1.0f GeV/c",ptBins_low[ipt],ptBins_high[ipt]),0.06);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/comp_JpsiSig/%s_InvMass_cent%s.pdf",run_type,saveTitle,cent_Title[icent]));

  // mean & sigma
  TList *list = new TList;
  for(int i=0; i<nfiles; i++)
    {
      hMean[i]->SetMarkerStyle(21+i*4);
      hMean[i]->SetMarkerColor(color[i]);
      hMean[i]->SetLineColor(color[i]);
      list->Add(hMean[i]);
    }
  c = drawHistos(list,"mean","Mean of J/#Psi mass peak",false,0,10,false,2.5,3.5,false,true,legName,true,"",0.5,0.7,0.3,0.5,true,0.04,0.04,false,1,false,false);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/comp_JpsiSig/%s_JpsiMean_cent%s.pdf",run_type,saveTitle,cent_Title[icent]));

  list->Clear();
  for(int i=0; i<nfiles; i++)
    {
      hSigma[i]->SetMarkerStyle(21+i*4);
      hSigma[i]->SetMarkerColor(color[i]);
      hSigma[i]->SetLineColor(color[i]);
      list->Add(hSigma[i]);
    }
  c = drawHistos(list,"sigma","Width of J/#Psi mass peak",false,0,10,false,2.5,3.5,false,true,legName,true,"",0.55,0.75,0.2,0.4,true,0.04,0.04,false,1,false,false);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/comp_JpsiSig/%s_JpsiSigma_cent%s.pdf",run_type,saveTitle,cent_Title[icent]));

  list->Clear();
  for(int i=0; i<nfiles; i++)
    {
      hYield[i]->Scale(1./hYield[i]->Integral());
      hYield[i]->SetMarkerStyle(21+i*4);
      hYield[i]->SetMarkerColor(color[i]);
      hYield[i]->SetLineColor(color[i]);
      list->Add(hYield[i]);
    }
  c = drawHistos(list,"yield","Raw distribution of J/#Psi signal;p_{T} (GeV/c);Prob.",false,0,10,true,0,0.4,false,true,legName,true,"",0.55,0.75,0.6,0.8,true,0.04,0.04,false,1,false,false);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/comp_JpsiSig/%s_RawCounts_cent%s.pdf",run_type,saveTitle,cent_Title[icent]));
}
