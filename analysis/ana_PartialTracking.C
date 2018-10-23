TFile *f;
const int year = YEAR;


//================================================
void ana_PartialTracking()
{
  gStyle->SetOptStat(0);
  TString fileName;

  if(year==2013)
    {
    }
  else if(year==2014)
    {
      fileName = Form("Run14_AuAu200.Study.PartialTrk.root");
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");

  pairEff();
  //efficiency();
}

//================================================
void pairEff(const Int_t savePlot = 0)
{
  TFile *fstudy = TFile::Open("output/Run14_AuAu200.Study.PartialTrk.root", "read");
  TH3F *hStat = (TH3F*)fstudy->Get("mhEvtCountParTrk_di_mu");
  TH1F *hReject = (TH1F*)hStat->ProjectionZ("hReject", 1, 4, 1, 16);
  c = draw1D(hReject);
  TH1F *hCentAll[gNTrgSetup];
  TH1F *hCentAcc[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k==0)
	{
	  hCentAll[k] = (TH1F*)hStat->ProjectionY(Form("hCentAll%s",gTrgSetupTitle[k]),1,4,1,3);
	  hCentAcc[k] = (TH1F*)hStat->ProjectionY(Form("hCentAcc%s",gTrgSetupTitle[k]),1,4,3,3);
	}
      else
	{
	  hCentAll[k] = (TH1F*)hStat->ProjectionY(Form("hCentAll%s",gTrgSetupTitle[k]),k,k,1,3);
	  hCentAcc[k] = (TH1F*)hStat->ProjectionY(Form("hCentAcc%s",gTrgSetupTitle[k]),k,k,3,3);
	}
      hCentAcc[k]->Sumw2();
      hCentAcc[k]->Divide(hCentAll[k]);
    }
  TCanvas *cCent = new TCanvas("cCent", "cCent", 800, 600);
  TLegend *leg = new TLegend(0.15,0.65,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(int k=0; k<gNTrgSetup; k++)
    {
      hCentAcc[k]->SetTitle(";centrality;");
      hCentAcc[k]->SetMarkerStyle(20+k);
      hCentAcc[k]->SetMarkerColor(gColor[k]);
      hCentAcc[k]->SetLineColor(gColor[k]);
      if(k==0) hCentAcc[k]->Draw("PE");
      else     hCentAcc[k]->Draw("samesPE");
      leg->AddEntry(hCentAcc[k], gLegNameTrg[k].Data(), "PL");
    }
  TPaveText *t1 = GetTitleText("Retention fraction");
  t1->Draw();
  leg->Draw();

  // pair efficiency  
  const int nHistos = 2;
  const char* histoName[nHistos] = {"TrigMth", "TrigMuon"};
  TH1F *hPairPtAll[nHistos][gNTrgSetup];
  TH1F *hPairPtAcc[nHistos][gNTrgSetup];
  TH1F *hPairPtEff[nHistos][gNTrgSetup];
  THnSparseF *hnInvMass = (THnSparseF*)fstudy->Get("mhInvMassParTrkEff_di_mu");
  hnInvMass->GetAxis(0)->SetRange(1, 1);  // like-sign
  hnInvMass->GetAxis(1)->SetRangeUser(3.0, 3.2); // pair mass
  hnInvMass->GetAxis(3)->SetRange(1, 16); // centrality
  hnInvMass->GetAxis(7)->SetRange(2, 2); // fire trigger
  for(int i=0; i<nHistos; i++)
    {
      if(i==1) hnInvMass->GetAxis(6)->SetRange(2, 2);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hnInvMass->GetAxis(4)->SetRange(k, k);
	  hPairPtAll[i][k] = (TH1F*)hnInvMass->Projection(2);
	  hPairPtAll[i][k]->Sumw2();
	  hPairPtAll[i][k]->Rebin(2);
	  hPairPtAll[i][k]->SetName(Form("hPairPtAll_%s%s",histoName[i],gTrgSetupTitle[k]));
      
	  hnInvMass->GetAxis(5)->SetRange(3, 3);
	  hPairPtAcc[i][k] = (TH1F*)hnInvMass->Projection(2);
	  hPairPtAcc[i][k]->Sumw2();
	  hPairPtAcc[i][k]->Rebin(2);
	  hPairPtAcc[i][k]->SetName(Form("hPairPtAcc_%s%s",histoName[i],gTrgSetupTitle[k]));
	  hnInvMass->GetAxis(5)->SetRange(0, -1);
	  hnInvMass->GetAxis(4)->SetRange(0, -1);
      
	  hPairPtEff[i][k] = (TH1F*)hPairPtAcc[i][k]->Clone(Form("hPairPtEff_%s%s",histoName[i],gTrgSetupTitle[k]));
	  hPairPtEff[i][k]->Divide(hPairPtAll[i][k]);
	}
    }
    // show the efficiency
  const char* histoTitle[nHistos] = {"pair of matched tracks that fire trigger", 
				     "pair of muon candidates that fire trigger"};
  TCanvas *cEff[nHistos];
  TF1 *funcEffRatio[nHistos][gNTrgSetup];
  for(int i=0; i<nHistos; i++)
    {
      cEff[i] = new TCanvas(Form("cEff_%s",histoName[i]), Form("cEff_%s",histoName[i]), 1100, 500);
      cEff[i]->Divide(2,1);

      cEff[i]->cd(1);
      TLegend *leg = new TLegend(0.15,0.65,0.35,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  hPairPtEff[i][k]->SetTitle(";p_{T} (GeV/c);");
	  hPairPtEff[i][k]->SetMarkerStyle(20+k);
	  hPairPtEff[i][k]->SetMarkerColor(gColor[k]);
	  hPairPtEff[i][k]->SetLineColor(gColor[k]);
	  if(i==0) hPairPtEff[i][k]->GetYaxis()->SetRangeUser(0.8,1.2);
	  if(i==1) hPairPtEff[i][k]->GetYaxis()->SetRangeUser(0.8, 1.2);

	  if(k==0) hPairPtEff[i][k]->Draw("PE");
	  else     hPairPtEff[i][k]->Draw("samesPE");
	  leg->AddEntry(hPairPtEff[i][k], gLegNameTrg[k].Data(), "PL");
	}
      TPaveText *t1 = GetTitleText(Form("%s",histoTitle[i]));
      t1->Draw();
      leg->Draw();

      cEff[i]->cd(2);
      leg = new TLegend(0.15,0.15,0.8,0.25);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetNColumns(2);
      for(int k=1; k<gNTrgSetup; k++)
	{
	  TH1F *hRatio = (TH1F*)hPairPtEff[i][k]->Clone(Form("%s_ratio",hPairPtEff[i][k]->GetName()));
	  hRatio->Divide(hPairPtEff[i][0]);
	  funcEffRatio[i][k] = new TF1(Form("func%sEffRatio%s",histoName[i],gTrgSetupTitle[k]), "pol0", 1, 10);
	  hRatio->Fit(funcEffRatio[i][k],"0RQ");
	  hRatio->SetTitle(";p_{T} (GeV/c);Ratio to inclusive");
	  hRatio->GetYaxis()->SetRangeUser(0.8,1.2);
	  if(k==1) hRatio->Draw("PE");
	  else     hRatio->Draw("samesPE");
	  leg->AddEntry(hRatio, Form("%4.3f +/- %4.3f",funcEffRatio[i][k]->GetParameter(0),funcEffRatio[i][k]->GetParError(0)), "PL");
	  funcEffRatio[i][k]->SetLineColor(gColor[k]);
	  funcEffRatio[i][k]->SetLineStyle(2);
	  funcEffRatio[i][k]->Draw("sames");
	}
      t1 = GetTitleText(Form("Ratio to inclusive sample"));
      t1->Draw();
      leg->Draw();
      if(savePlot)
	cEff[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_ParTrkEff_%s.pdf",run_type.Data(),histoName[i]));
    }
  //KEY: TH3F	mhEvtCountParTrk_di_mu;1	TrigSetup vs. centrality vs. ShouldHaveReject
  //KEY: THnSparseT<TArrayF>	mhInvMassParTrkEff_di_mu;1	Type vs. M_{#mu#mu} vs. p_{T,#mumu} vs. centrality vs. trigSetup vs. RejectEvent vs. isMuon vs. isTrig
}

//================================================
void efficiency(const Int_t save = 0)
{
  TH1F *hEventCount = (TH1F*)f->Get("mPartialTrkEvents");
  for(Int_t ibin=1; ibin<=6; ibin++)
    {
      printf("%s: %d\n",hEventCount->GetXaxis()->GetBinLabel(ibin),hEventCount->GetBinContent(ibin));
    }
  c = draw1D(hEventCount,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Event_Count.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Event_Count.png",run_type));
    }

  const int nbins = 19;
  const double xbins[nbins+1] = {0,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,4.0,5.0,6.0,8,9,10};
  TH2F *hPartialTrkMuonPt = (TH2F*)f->Get("mPartialTrkMuonPt");
  TH1F *hMuonPt[2];
  TList *list = new TList();
  for(int i=0; i<2; i++)
    {
      TH1F *htmp = (TH1F*)hPartialTrkMuonPt->ProjectionX(Form("hMuonPt_%d_tmp",i),i+1,i+1);
      hMuonPt[i] = (TH1F*)htmp->Rebin(nbins,Form("hMuonPt_%d",i),xbins);
      scaleHisto(hMuonPt[i], 1, 1,kTRUE,kFALSE, kFALSE);
      list->Add(hMuonPt[i]);
    }
  TString legName[2] = {"All events with >=2 muons","Filtered events with >=2 muons"};
  c = drawHistos(list,"MuonPt","p_{T} distribution of muon candidates;p_{T} (GeV/c)",kFALSE,0,5,kFALSE,0.1,10,kTRUE,kTRUE,legName,kTRUE,"Au+Au @ 200 GeV",0.45,0.65,0.55,0.75);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/MuonPt.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/MuonPt.png",run_type));
    }

  TGraphAsymmErrors *ratio = new TGraphAsymmErrors();
  ratio->SetName("ratio");
  ratio->Divide(hMuonPt[1],hMuonPt[0],"cl=0.683 b(1,1) mode");
  ratio->SetMarkerStyle(21);
  ratio->GetYaxis()->SetRangeUser(0.98,1.01);
  c = drawGraph((TGraph*)ratio,"Efficiency of partial tracking;p_{T} (GeV/c);efficiency");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Efficiency_vs_pt.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Efficiency_vs_pt.png",run_type));
    }

}
