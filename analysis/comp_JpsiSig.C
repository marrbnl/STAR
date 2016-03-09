//================================================
void comp_JpsiSig()
{
  gStyle->SetOptStat(0);
  comparison();
}


//================================================
void comparison(const int savePlot = 1)
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
