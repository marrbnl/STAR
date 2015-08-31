TFile *f;
const char *signal_name[2] = {"J/#psi","#Upsilon(1S)"};
Int_t trk_index = 0;
Int_t signal_index = 0;

const char *run_config = "EmbedQA.";
const int year = 2013;
TString run_cfg_name;

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

  if(year==2013)
    {
      run_type = "Run13_pp500";
      f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sMC.root",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sMC.root",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("# of di-muon events: %d\n",hStat->GetBinContent(7));
  
  MCtruth();
  //embedded();
  //efficiency();
  //fakeRate();
}


//================================================
void embedded(const Int_t save = 0)
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
void efficiency(const Int_t save = 0)
{
  // tracking efficiency
  TH1F *hMcTrkPt[3];
  hMcTrkPt[0] = (TH1F*)f->Get(Form("hMcTrkPt_%s",trigName[kTrigType]));
  hMcTrkPt[1] = (TH1F*)f->Get(Form("hMcTrkPtTpc_%s",trigName[kTrigType]));
  hMcTrkPt[2] = (TH1F*)f->Get(Form("hMcTrkPtMtd_%s",trigName[kTrigType]));

  TList *list = new TList;
  for(Int_t i=0; i<3; i++)
    {
      scaleHisto( hMcTrkPt[i], 1, 1, kTRUE);
      hMcTrkPt[i]->SetMaximum(10*hMcTrkPt[i]->GetMaximum());
      list->Add(hMcTrkPt[i]);
    }
  const TString legName[3] = {"MC tracks","Matched to TPC","Matched to TPC+MTD"};
  c = drawHistos(list,"MC_trkPt",Form("Au+Au %s embedding: p_{T} distribution of MC muons (|#eta|<0.5);p_{T,true} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName,kTRUE,Form("%s #rightarrow #mu^{+}#mu^{-}",signal_name[signal_index]),0.55,0.75,0.65,0.88);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_muon_pt.png",run_config));

  list->Clear();
  TH1F *hMcTrkPtRatio[3];
  for(Int_t i=1; i<3; i++)
    {
      hMcTrkPtRatio[i] = (TH1F*)hMcTrkPt[i]->Clone(Form("%s_ratio",hMcTrkPt[i]->GetName()));
      hMcTrkPtRatio[i]->Divide(hMcTrkPt[0]);
      hMcTrkPtRatio[i]->SetMarkerColor(color[i]);
      hMcTrkPtRatio[i]->SetLineColor(color[i]);
      list->Add(hMcTrkPtRatio[i]);
    }
  const TString legName2[2] = {"TPC","TPC+MTD"};
  c = drawHistos(list,"MC_trkPt_eff",Form("Au+Au %s embedding: efficiency of MC muons (|#eta|<0.5);p_{T,true} (GeV/c);efficiency",trigName[kTrigType]),kTRUE,0,20,kTRUE,0.,1.5,kFALSE,kTRUE,legName2,kTRUE,Form("%s #rightarrow #mu^{+}#mu^{-}",signal_name[signal_index]),0.55,0.75,0.65,0.88,kTRUE,0.04,0.04,kFALSE,1,kTRUE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_muon_eff_vs_pt.png",run_config));

  // muon vs hadron
  const char *trackType[3] = {"Hybrid","Data","MCreco"};
  TH1F *hTrkPt[3];
  TH1F *hMthTrkPt[3];
  for(Int_t i=0; i<3; i++)
    {
      hTrkPt[i] = (TH1F*)f->Get(Form("hRcTrkPtTpc_%s_%s",trackType[i],trigName[kTrigType]));
      hMthTrkPt[i] = (TH1F*)f->Get(Form("hRcTrkPtMtd_%s_%s",trackType[i],trigName[kTrigType]));
      scaleHisto( hTrkPt[i], 1, 1, kTRUE);
      scaleHisto( hMthTrkPt[i], 1, 1, kTRUE);
    }

  list->Clear();
  TString legName3[3] = {"Data+embedding","Data (hadrons)","Embedded (muons)"};
  for(Int_t i=0; i<3; i++)
    {
      list->Add(hTrkPt[i]);
    }
  c = drawHistos(list,"Reco_trkPt",Form("Au+Au %s embedding: p_{T} distribution of primary tracks;p_{T,reco} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName3,kTRUE,"|#eta| < 0.8",0.5,0.6,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.Trk_reco_pt.png",run_config));

  list->Clear();
  for(Int_t i=0; i<3; i++)
    {
      list->Add(hMthTrkPt[i]);
    }
  c = drawHistos(list,"Reco_MthTrkPt",Form("Au+Au %s embedding: p_{T} of primary tracks matched to MTD hits;p_{T,reco} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName3,kTRUE,"|#eta| < 0.8",0.5,0.6,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MthTrk_reco_pt.png",run_config));

  list->Clear();
  TString legName4[2] = {"Data (hadrons)","Embedded (muons)"};
  TH1F *hMthEff[2];
  hMthEff[0] = (TH1F*)hMthTrkPt[1]->Clone("MthEff_Data");
  hMthEff[0]->Divide(hTrkPt[1]);
  hMthEff[1] = (TH1F*)hMthTrkPt[2]->Clone("MthEff_McReco");
  hMthEff[1]->Divide(hTrkPt[2]);
  hMthEff[0]->GetYaxis()->SetNdivisions(505);
  hMthEff[0]->GetYaxis()->SetTitleOffset(1.1);
  list->Add(hMthEff[0]);
  list->Add(hMthEff[1]);
  c = drawHistos(list,"MthEff",Form("Au+Au %s embedding: track matching efficiency to MTD hits;p_{T,reco} (GeV/c);efficiency",trigName[kTrigType]),kTRUE,0,20,kTRUE,0,0.5,kFALSE,kTRUE,legName4,kTRUE,"|#eta| < 0.8",0.5,0.6,0.65,0.85);
  SetPadMargin(gPad,0.11,0.11);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MthEff_vs_reco_pt.png",run_config));
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
	  printf("# of input signals within |eta|<0.5: %d\n",hJpsi[0][j]->Integral(hJpsi[0][j]->FindFixBin(-0.5+0.01),hJpsi[0][j]->FindFixBin(0.5-0.01)));
	}
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.MC_Input_%s.png",run_config,variable[j]));
    }

  for(Int_t j=0; j<4; j++)
    {
      list->Clear();
      for(Int_t i=0; i<3; i++)
	{
	  list->Add(hJpsi[i][j]);
	  if(j==2) printf("# of signals within |eta|<0.5: %d\n",hJpsi[i][j]->Integral(hJpsi[i][j]->FindFixBin(-0.5+0.01),hJpsi[i][j]->FindFixBin(0.5-0.01)));
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


//================================================
void fakeRate(const Int_t save = 1)
{
  const char *trackType[5] = {"Hybrid","Data","MCreco","MCreco_true","MCreco_fake"};
  TH1F *hMthTrkPt[5];
  TH2F *hDyVsPt[5];
  TH1F *hDy[5];
  TH2F *hDzVsPt[5];
  TH1F *hDz[5];
  for(Int_t i=0; i<5; i++)
    {
      hMthTrkPt[i] = (TH1F*)f->Get(Form("hRcTrkPtMtd_%s_%s",trackType[i],trigName[kTrigType]));
      scaleHisto( hMthTrkPt[i], 1, 1, kTRUE);

      hDyVsPt[i] = (TH2F*)f->Get(Form("hDyVsTrkPtMtd_%s_%s",trackType[i],trigName[kTrigType]));
      hDzVsPt[i] = (TH2F*)f->Get(Form("hDzVsTrkPtMtd_%s_%s",trackType[i],trigName[kTrigType]));

      hDy[i] = (TH1F*)hDyVsPt[i]->ProjectionY(Form("%s_proy",hDyVsPt[i]->GetName()));  
      hDz[i] = (TH1F*)hDzVsPt[i]->ProjectionY(Form("%s_proy",hDzVsPt[i]->GetName()));
    }

  // fake rate
  TList *list = new TList;
  for(Int_t i=2; i<5; i++)
    {
      hMthTrkPt[i]->SetMinimum(0.1);
      list->Add(hMthTrkPt[i]);
    }  
  
  const TString legName[3] = {"Sum","Matched to MTD hits from embedding","Matched to MTD hits from data"};
  c = drawHistos(list,"MC_trkPt",Form("Au+Au %s embedding: p_{T} of TPC muon tracks matched MTD hits;p_{T,rec} (GeV/c);dN/dp_{T}",trigName[kTrigType]),kTRUE,0,20,kFALSE,-200,300,kTRUE,kTRUE,legName,kTRUE,"",0.25,0.55,0.15,0.35);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.Rec_muon_pt_matched.png",run_config));

  list->Clear();
  TH1F *hTrkPtRatio[2];
  for(Int_t i=0; i<2; i++)
    {
      hTrkPtRatio[i] = (TH1F*)hMthTrkPt[i+3]->Clone(Form("%s_ratio",hMthTrkPt[i]->GetName()));
      hTrkPtRatio[i]->Divide(hMthTrkPt[2]);
      hTrkPtRatio[i]->SetMarkerColor(color[i+1]);
      hTrkPtRatio[i]->SetLineColor(color[i+1]);
      list->Add(hTrkPtRatio[i]);
    }
  const TString legName2[2] = {"Matched to MTD hits from embedding","Matched to MTD hits from data"};
  c = drawHistos(list,"MC_trkPt_eff",Form("Au+Au %s embedding: fake rate of MTD matching;p_{T,rec} (GeV/c);efficiency",trigName[kTrigType]),kTRUE,0,20,kTRUE,0.01,1.5,kTRUE,kTRUE,legName2,kTRUE,"",0.3,0.5,0.5,0.7,kTRUE,0.04,0.04,kFALSE,1,kTRUE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.fake_rate_vs_pt.png",run_config));

  c = draw2D(hDyVsPt[1],"Data: #Deltay vs p_{T} of charged hadrons (primary);p_{T,rec} (GeV/c)");
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.dy_vs_pt_hadron.png",run_config));

  c = draw2D(hDyVsPt[3],"Embedding: #Deltay vs p_{T} of charged muons (primary);p_{T,rec} (GeV/c)");
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.dy_vs_pt_muon.png",run_config));

  list->Clear();
  list->Add(hDy[3]);
  list->Add(hDy[1]);
  const TString legName3[2] = {"Embedding: muons","Data: hadrons"};
  c = drawHistos(list,"dy_dis",Form("Au+Au %s embedding: #Deltay distribution;#Deltay (cm)",trigName[kTrigType]),kTRUE,-50,50,kFALSE,0.01,1.5,kFALSE,kTRUE,legName3,kTRUE,"",0.6,0.7,0.5,0.7,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.dy_hadron_vs_muon.png",run_config));

  // c = draw2D(hDzVsPt[1],"Data: #Deltaz vs p_{T} of charged hadrons (primary);p_{T,rec} (GeV/c)");
  // if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.dz_vs_pt_hadron.png",run_config));

  // c = draw2D(hDzVsPt[3],"Embedding: #Deltaz vs p_{T} of charged muons (primary);p_{T,rec} (GeV/c)");
  // if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.dz_vs_pt_muon.png",run_config));

  // list->Clear();
  // list->Add(hDz[3]);
  // list->Add(hDz[1]);
  // c = drawHistos(list,"dz_dis",Form("Au+Au %s embedding: #Deltaz distribution;#Deltaz (cm)",trigName[kTrigType]),kTRUE,-50,50,kFALSE,0.01,1.5,kFALSE,kTRUE,legName3,kTRUE,"",0.6,0.7,0.5,0.7,kFALSE);
  // if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_Embed/%s.dz_hadron_vs_muon.png",run_config));
}

