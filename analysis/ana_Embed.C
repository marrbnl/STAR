TFile *f;
const char *signal_name[2] = {"J/#psi","#Upsilon(1S)"};
Int_t trk_index = 0;
Int_t signal_index = 0;
const char *run_config = "Embed.Upsilon.PrimaryTrk";

//================================================
void ana_Embed()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.15); 

  TString cut_name = run_config;
  if(cut_name.Contains("Global"))
    trk_index = 1;
  if(cut_name.Contains("Upsilon"))
    signal_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("# of di-muon events: %d\n",hStat->GetBinContent(7));
  
  MCtruth();
  //embedded();
  //efficiency();
}


//================================================
void embedded(const Int_t save = 1)
{
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  const char *name[4] = {"Hybrid","Data","MCone","MCtwo"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3][4];
  Double_t pt1_cut = 1.0, pt2_cut = 1.0;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      for(Int_t i=0; i<4; i++)
	{
	  if(i==0) hnInvMass[j]->GetAxis(7)->SetRange(1,3);
	  else     hnInvMass[j]->GetAxis(7)->SetRange(i,i);
	  hInvMass[j][i] = (TH1F*)hnInvMass[j]->Projection(0);
	  hInvMass[j][i]->SetName(Form("%s_%s_InvMass_%s",hName[j],trigName[kTrigType],name[i]));
	  hInvMass[j][i]->Sumw2();
	  hInvMass[j][i]->Rebin(10);
	  hnInvMass[j]->GetAxis(7)->SetRange(0,-1);
	}
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }
  TList *list = new TList;
  TString legName[4] = {"Data+MC","Data","One track matched to MC track","Two tracks matched to MC track"};
  for(Int_t i=0; i<4; i++)
    {
      list->Add(hInvMass[0][i]);
    }
  c = drawHistos(list,"Hybrid_InvMass",Form("Au+Au %s embedding: invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kFALSE,kTRUE,legName,kTRUE,Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut),0.4,0.6,0.5,0.8,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.Hybrid_InvMass_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));
  return;

  TH1F *hSignal = (TH1F*)hInvMass[0][0]->Clone("hSignal");
  TH1F *hBkg    = (TH1F*)hInvMass[1][0]->Clone("hBkg");
  hBkg->Add(hInvMass[2][0]);
  c = draw1D(hSignal);
  hBkg->Draw("sames");
  TH1F *hDiff = (TH1F*)hSignal->Clone("hDiff");
  hDiff->Add(hBkg,-1);
  c = draw1D(hDiff);
  
}

//================================================
void efficiency(const Int_t save = 1)
{
  // tracking efficiency
  const char *trackType[5] = {"Hybrid","Data","MCreco","MCtruth","Reco"};
  TH1F *hTrkPt[5];
  TH1F *hMthTrkPt[5];
  for(Int_t i=0; i<4; i++)
    {
      hTrkPt[i] = (TH1F*)f->Get(Form("hTrkPt_%s_%s",trackType[i],trigName[kTrigType]));
      hMthTrkPt[i] = (TH1F*)f->Get(Form("hMthTrkPt_%s_%s",trackType[i],trigName[kTrigType]));
      scaleHisto( hTrkPt[i], 1, 1, kTRUE);
      scaleHisto( hMthTrkPt[i], 1, 1, kTRUE);
    }

  TList *list = new TList;
  for(Int_t i=0; i<2; i++)
    {
      hTrkPt[3-i]->SetMaximum(50.*hTrkPt[3-i]->GetMaximum());
      list->Add(hTrkPt[3-i]);
    }
  const TString legName[2] = {"MC tracks","Matched to reconstructed tracks"};
  c = drawHistos(list,"MC_trkPt",Form("Au+Au %s embedding: p_{T} distribution of MC muons;p_{T,true} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName,kTRUE,"",0.15,0.25,0.7,0.88);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_muon_pt.png",run_config));

  TH1F *hEff = (TH1F*)hTrkPt[2]->Clone("hTrkEff");
  hEff->Divide(hTrkPt[3]);
  hEff->SetMarkerStyle(21);
  hEff->SetMarkerColor(1);
  hEff->SetLineColor(1);
  hEff->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hEff,Form("Au+Au %s embedding: tracking efficiency of MC muons;p_{T,true} (GeV/c);tracking efficiency",trigName[kTrigType]));
  SetPadMargin(gPad,0.11,0.11);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_muon_eff.png",run_config));

  // Matching efficiency
  TString legName2[3] = {"Data+embedding","Data","Matched to MC tracks"};
  list->Clear();
  hTrkPt[4] = (TH1F*)hTrkPt[0]->Clone(Form("hTrkPt_%s_%s",trackType[4],trigName[kTrigType]));
  hTrkPt[4]->Add(hTrkPt[1],-1);
  for(Int_t ibin=1; ibin<=hTrkPt[4]->GetNbinsX(); ibin++)
    {
      hTrkPt[4]->SetBinError(ibin,TMath::Sqrt(hTrkPt[4]->GetBinContent(ibin)));
    }
  list->Add(hTrkPt[0]);
  list->Add(hTrkPt[1]);
  list->Add(hTrkPt[4]);
  c = drawHistos(list,"Reco_trkPt",Form("Au+Au %s embedding: p_{T} distribution of tracks;p_{T,reco} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName2,kTRUE,"",0.35,0.5,0.6,0.8);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.Trk_pt.png",run_config));

  list->Clear();
  hMthTrkPt[4] = (TH1F*)hMthTrkPt[0]->Clone(Form("hMthTrkPt_%s_%s",trackType[4],trigName[kTrigType]));
  hMthTrkPt[4]->Add(hMthTrkPt[1],-1);
  for(Int_t ibin=1; ibin<=hMthTrkPt[4]->GetNbinsX(); ibin++)
    {
      hMthTrkPt[4]->SetBinError(ibin,TMath::Sqrt(hMthTrkPt[4]->GetBinContent(ibin)));
    }
  list->Add(hMthTrkPt[0]);
  list->Add(hMthTrkPt[1]);
  list->Add(hMthTrkPt[4]);
  c = drawHistos(list,"Reco_MthTrkPt",Form("Au+Au %s embedding: p_{T} distribution of tracks matched to MTD hits;p_{T,reco} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName2,kTRUE,"",0.35,0.5,0.6,0.8);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MthTrk_pt.png",run_config));

  list->Clear();
  TString legName3[2] = {"Data","Matched to MC tracks"};
  TH1F *hMthEff[2];
  hMthEff[0] = (TH1F*)hMthTrkPt[1]->Clone("MthEff_Data");
  hMthEff[0]->Divide(hTrkPt[1]);
  hMthEff[1] = (TH1F*)hMthTrkPt[4]->Clone("MthEff_McReco");
  hMthEff[1]->Divide(hTrkPt[4]);
  hMthEff[0]->GetYaxis()->SetNdivisions(505);
  hMthEff[0]->GetYaxis()->SetTitleOffset(1.1);
  list->Add(hMthEff[0]);
  list->Add(hMthEff[1]);
  c = drawHistos(list,"MthEff",Form("Au+Au %s embedding: track matching efficiency to MTD hits;p_{T,reco} (GeV/c);efficiency",trigName[kTrigType]),kTRUE,0,20,kTRUE,0,0.5,kFALSE,kTRUE,legName3,kTRUE,"",0.35,0.5,0.6,0.8);
  SetPadMargin(gPad,0.11,0.11);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MthEff_vs_pt.png",run_config));
}

//================================================
void MCtruth(const Int_t save = 0)
{
  TList *list = new TList;
  // upsilon mass peak
  const char *name[3] = {"MCtruth","MC_TPC","MC_MTD"};
  const char *variable[4] = {"Mass","Pt","Eta","Phi"};
  TH1F *hJpsi[3][4], *hJpsiPt[3], *hJpsiEta[3], *hJpsiPhi[3];
  THnSparseF *hn = (THnSparseF*)f->Get(Form("hJpsiInfo_%s",trigName[kTrigType]));

  hn->GetAxis(7)->SetRange(4,4);
  TH2F *hDaugPt = (TH2F*)hn->Projection(4,5);
  c = draw2D(hDaugPt,Form("MC truth: p_{T} correlation of daughter muons from %s",signal_name[signal_index]));
  hn->GetAxis(7)->SetRange(0,-1);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_truth_DaugPtCorr.png",run_config));

  for(Int_t i=0; i<3; i++)
    {
      hn->GetAxis(7)->SetRange(i+4,i+4);
      for(Int_t j=0; j<4; j++)
	{
	  hJpsi[i][j] = (TH1F*)hn->Projection(j);
	  hJpsi[i][j]->SetName(Form("hJpsi%s_%s_%s",variable[j],name[i],trigName[kTrigType]));
	}
      hn->GetAxis(7)->SetRange(0,-1);
    }

  const char *title[4] = {"","p_{T}","#eta","#varphi"};
  TString legName2[3] = {"MC truth","Two muons matched to TPC tracks","Two muons matched to TPC+MTD"};
  for(Int_t j=0; j<4; j++)
    {
      hJpsi[0][j]->GetYaxis()->SetRangeUser(0.1,hJpsi[0][j]->GetMaximum()*1.2);
      if(j==0) printf("# of embedded signals: %d\n",hJpsi[0][j]->Integral());
      if(j==0) 
	{
	  hJpsi[0][j]->GetXaxis()->SetRangeUser(9,10);
	  c = draw1D(hJpsi[0][j],Form("Au+Au %s embedding: invariant mass distribution of di-muon pairs;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kTRUE,kFALSE);
	}
      else
	{
	  if(j==1) hJpsi[0][j]->GetXaxis()->SetRangeUser(0,20);
	  c = draw1D(hJpsi[0][j],Form("Au+Au %s embedding: %s distribution of embedded %s",trigName[kTrigType],title[j],signal_name[signal_index]),kFALSE,kFALSE);
	}
      if(j==2)
	{
	  printf("# of input signals within |eta|<0.5: %d\n",hJpsi[0][j]->Integral(hJpsi[0][j]->FindFixBin(-0.75),hJpsi[0][j]->FindFixBin(0.75)));
	}
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_Input_%s.png",run_config,variable[j]));
    }
  return;

  for(Int_t j=0; j<4; j++)
    {
      list->Clear();
      for(Int_t i=0; i<3; i++)
	{
	  list->Add(hJpsi[i][j]);
	  printf("# of: %d\n",hJpsi[i][j]->Integral());
	}
      if(j==0) c = drawHistos(list,Form("MC_signal_%s",variable[j]),Form("Au+Au %s embedding: invariant mass distribution of di-muon pairs;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kTRUE,0,12,kFALSE,0.1,10,kTRUE,kTRUE,legName2,kTRUE,"",0.15,0.25,0.7,0.88,kFALSE);
      else if(j==1) c = sysCompare(list,Form("MC_signal_%s",variable[j]),Form("Au+Au %s: %s distribution of embedded %s",trigName[kTrigType],title[j],signal_name[signal_index]),Form("Efficiency as a function of %s",title[j]),kTRUE,0,15,kFALSE,0.1,10,kTRUE,0,0.5,kFALSE,kTRUE,legName2,kTRUE,"",0.15,0.25,0.35,0.55,kFALSE);
      else if(j==2) c = sysCompare(list,Form("MC_signal_%s",variable[j]),Form("Au+Au %s: %s distribution of embedded %s",trigName[kTrigType],title[j],signal_name[signal_index]),Form("Efficiency as a function of %s",title[j]),kFALSE,0,20,kFALSE,0.1,10,kTRUE,0,0.5,kFALSE,kTRUE,legName2,kTRUE,"",0.15,0.25,0.76,0.89,kFALSE);
      else if(j==3) c = sysCompare(list,Form("MC_signal_%s",variable[j]),Form("Au+Au %s: %s distribution of embedded %s",trigName[kTrigType],title[j],signal_name[signal_index]),Form("Efficiency as a function of %s",title[j]),kFALSE,0,20,kFALSE,0.1,10,kTRUE,0,0.5,kFALSE,kTRUE,legName2,kTRUE,"",0.15,0.25,0.35,0.55,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_reco_%s.png",run_config,variable[j]));
    }


  TH1F *hReco = (TH1F*)hJpsi[1][0]->Clone("hRecoSignal_TPC");
  hReco->SetLineColor(1);
  hReco->GetXaxis()->SetRangeUser(7,12);
  TF1 *func = new TF1("func","gaus",9,10);
  hReco->Fit(func,"IR0");
  c = draw1D(hReco,Form("Au+Au %s embedding: reconstructed true di-muon pairs;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kFALSE,kFALSE);
  func->SetLineColor(2);
  func->Draw("sames");
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.Fit_MC_reco_mass.png",run_config));
}
