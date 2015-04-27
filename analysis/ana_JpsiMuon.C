const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;

const char *run_config = "TrigHit.";
const Bool_t iPico = 1;
const int year = 2013;
TString run_cfg_name;

TFile *f;

//================================================
void ana_JpsiMuon()
{
  gStyle->SetOptStat(0);

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  DeltaTof();
  //DeltaY();
  //DeltaZ();
}

//================================================
void DeltaTof(const Int_t save = 1)
{
  TList *list = new TList;

  TH2F *hDtofVsMod = (TH2F*)f->Get(Form("mhDeltaTof_%s",trigName[kTrigType]));
  hDtofVsMod->GetYaxis()->SetRangeUser(-5,5);
  c = draw2D(hDtofVsMod,Form("%s: #Deltatof of tracks matched to MTD;p_{T} (GeV/c);#Deltatof (ns)",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_vs_Mod_MthTrk.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_vs_Mod_MthTrk.pdf",run_type,run_cfg_name.Data()));
    }

  TH2F *hDtofVsPt[4];
  hDtofVsPt[0] = (TH2F*)f->Get(Form("mhDtofVsPt_%s",trigName[kTrigType]));

  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhMuonDzDy_%s",trigName[kTrigType]));
  hn->GetAxis(0)->SetRangeUser(low_mass+0.001,high_mass-0.001);
  for(int i=0; i<2; i++)
    {
      hn->GetAxis(1)->SetRange(i+1,i+1);
      hDtofVsPt[i+2] = (TH2F*)hn->Projection(5,2);
      hDtofVsPt[i+2]->SetName(Form("hDtofVsPt_di-mu_%d",i+2));
      cout << "# of di-muon pairs: " << hDtofVsPt[i+2]->GetEntries()/2 << endl;
    }
  hDtofVsPt[1] = (TH2F*)hDtofVsPt[2]->Clone("hDtofVsPt_US-LS_di-mu");
  hDtofVsPt[1]->Add(hDtofVsPt[3],-1);
  hDtofVsPt[1]->SetTitle(Form("%s: #Deltatof of muon tracks (US-LS, %1.1f<M_{#mu#mu}<%1.1f)",trigName[kTrigType],low_mass,high_mass));

  // Dtof vs pt
  hDtofVsPt[0]->SetTitle(Form("%s: #Deltatof of tracks matched to MTD;p_{T} (GeV/c);#Deltatof (ns)",trigName[kTrigType]));
  const char *hName[2] = {"MthTrk","MuonTrk"};
  for(int i=0; i<2; i++)
    {
      c = draw2D(hDtofVsPt[i]);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtofVsPt_%s.png",run_type,run_cfg_name.Data(),hName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtofVsPt_%s.pdf",run_type,run_cfg_name.Data(),hName[i]));
	}
    }

  // Dtof in Pt bins
  list->Clear();
  TString legName2[4];
  double ptcuts[5] = {1,1.5,2,3,20};
  for(int i=0; i<4; i++)
    {
      int low = hDtofVsPt[1]->GetXaxis()->FindFixBin(ptcuts[i]+0.01);
      int hi  = hDtofVsPt[1]->GetXaxis()->FindFixBin(ptcuts[i+1]-0.01);
      TH1F *h1 = (TH1F*)hDtofVsPt[1]->ProjectionY(Form("hDtof_PtBin%d",i),low,hi);
      legName2[i] = Form("%1.1f < p_{T} < %1.1f",ptcuts[i],ptcuts[i+1]);
      list->Add(h1);
    }
  c = drawHistos(list,"DTofInPtBin",Form("%s: #Deltatof of muon tracks (US-LS, %1.1f<M_{#mu#mu}<%1.1f);#Deltatof (ns)",trigName[kTrigType],low_mass,high_mass),kTRUE,-4,5,kFALSE,1e-5,0.13,kFALSE,kTRUE,legName2,kTRUE,"",0.6,0.8,0.6,0.8,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_InPtBin_US-LS.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_InPtBin_US-LS.pdf",run_type,run_cfg_name.Data()));
    }

  // Dtof of US, LS, US-LS
  TH1F *hDeltaTof[4];
  for(Int_t i=0; i<4; i++)
    {
      hDeltaTof[i] = (TH1F*)hDtofVsPt[i]->ProjectionY(Form("hDtof_%d_clone",i));
      hDeltaTof[i]->Sumw2();
      hDeltaTof[i]->SetLineColor(color[i]);
      hDeltaTof[i]->SetMarkerColor(color[i]);
      hDeltaTof[i]->GetXaxis()->SetRangeUser(-5,5);
      hDeltaTof[i]->SetMaximum(500);
    }
  hDeltaTof[1]->SetMarkerStyle(21);
  c = draw1D(hDeltaTof[1],Form("%s: #Deltatof of muon tracks (%s)",trigName[kTrigType],run_type),kFALSE,kTRUE);
  hDeltaTof[2]->Draw("sames HIST");
  hDeltaTof[3]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.62,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("%1.1f < M_{#mu#mu} < %1.1f (GeV/c^{2})",low_mass,high_mass));
  leg->AddEntry(hDeltaTof[2],"Unlike-sign","L");
  leg->AddEntry(hDeltaTof[3],"Like-sign","L");
  leg->AddEntry(hDeltaTof[1],"US-LS","P");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_US_LS.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_US_LS.pdf",run_type,run_cfg_name.Data()));
    }

  // side-band
  TH1F *hDtofSide[2][2];
  for(int i=0; i<2; i++)
    {
      if(i==0) hn->GetAxis(0)->SetRangeUser(2.8+0.001,3.0-0.001);
      if(i==1) hn->GetAxis(0)->SetRangeUser(3.2+0.001,3.4-0.001);
      for(int j=0; j<2; j++)
	{
	  hn->GetAxis(1)->SetRange(j+1,j+1);
	  hDtofSide[i][j] = (TH1F*)hn->Projection(5);
	  hDtofSide[i][j]->Sumw2();
	  hDtofSide[i][j]->SetName(Form("hDtof_di-mu_%d_%d",i,j));
	}
    }
  TH1F *hSide = (TH1F*)hDtofSide[0][0]->Clone("hDtof_side_band");
  hSide->Add(hDtofSide[1][0]);
  hSide->Add(hDtofSide[0][1],-1);
  hSide->Add(hDtofSide[1][1],-1);
  hSide->Scale(0.5);
  TH1F *hDtofClone = (TH1F*)hDeltaTof[1]->Clone("hDtof_clone");
  TH1F *hDtofSBsub = (TH1F*)hDeltaTof[1]->Clone("hDtof_side_sub");
  hDtofSBsub->Add(hSide,-1);
  hDtofClone->SetMaximum(0.7*hDtofClone->GetMaximum());
  c = draw1D(hDtofClone,Form("%s: #Deltatof of muon tracks (%s)",trigName[kTrigType],run_type),kFALSE,kFALSE);
  hSide->Draw("sames HIST");
  hDtofSBsub->SetMarkerColor(4);
  hDtofSBsub->SetLineColor(4);
  leg = new TLegend(0.15,0.62,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hDtofClone,"US-LS","L");
  leg->AddEntry(hSide,"SideBand","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_SideBand.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_SideBand.pdf",run_type,run_cfg_name.Data()));
    }

  // fit side-band subtracted distribution
  hDtofSBsub->SetLineColor(1);
  c = draw1D(hDtofSBsub,Form("%s: #Deltatof of muon tracks (%s)",trigName[kTrigType],run_type),kFALSE,kTRUE);
  TLine *line = GetLine(-5,0,5,0);
  line->Draw();
  TPaveText *t1 = GetPaveText(0.2,0.4,0.6,0.8);
  t1->AddText("Side-band subtracted");
  t1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_SideBand_sub.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDtof_SideBand_sub.pdf",run_type,run_cfg_name.Data()));
    }

  // compare Dtof
  list->Clear();
  for(Int_t i=0; i<2; i++)
    {
      TH1F *htmp = (TH1F*)hDtofVsPt[i]->ProjectionY(Form("hDtof_%d",i));
      htmp->Sumw2();
      htmp->Scale(1./htmp->Integral());
      htmp->SetMaximum(10*htmp->GetMaximum());
      list->Add(htmp);
    }
  TString legName[2] = {"All matched tracks", "Muon (US-LS, 3.0<M_{#mu#mu}<3.2)"};
  c = drawHistos(list,"DeltaTof",Form("%s: #Deltatof of tracks;#Deltatof (ns)",trigName[kTrigType]),kFALSE,-3000,-1000,kFALSE,1e-5,0.13,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sCompareDtof.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sCompareDtof.pdf",run_type,run_cfg_name.Data()));
    }
  c = drawHistos(list,"DeltaTof_zoomin",Form("%s: #Deltatof of tracks;#Deltatof (ns)",trigName[kTrigType]),kTRUE,-4,5,kTRUE,1e-5,0.16,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sCompareDtof_zoomin.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sCompareDtof_zoomin.pdf",run_type,run_cfg_name.Data()));
    }
}
//================================================
void DeltaY(const Int_t save = 0)
{
  TH2F *hDyVsMod = (TH2F*)f->Get(Form("mhDeltaY_%s",trigName[kTrigType]));
  hDyVsMod->GetYaxis()->SetRangeUser(-100,100);
  c = draw2D(hDyVsMod,Form("%s: #Deltay of tracks matched to MTD;p_{T} (GeV/c);#Deltay (cm)",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDy_vs_Mod_MthTrk.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDy_vs_Mod_MthTrk.pdf",run_type,run_cfg_name.Data()));
    }

  TH2F *hDyVsPt[4];
  hDyVsPt[0] = (TH2F*)f->Get(Form("mhDyVsPt_%s",trigName[kTrigType]));

  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhMuonDzDy_%s",trigName[kTrigType]));
  hn->GetAxis(0)->SetRangeUser(low_mass+0.001,high_mass-0.001);
  for(int i=0; i<2; i++)
    {
      hn->GetAxis(1)->SetRange(i+1,i+1);
      hDyVsPt[i+2] = (TH2F*)hn->Projection(4,2);
      hDyVsPt[i+2]->SetName(Form("hDyVsPt_di-mu_%d",i+2));
      cout << "# of di-muon pairs: " << hDyVsPt[i+2]->GetEntries()/2 << endl;
    }
  hDyVsPt[1] = (TH2F*)hDyVsPt[2]->Clone("hDyVsPt_US-LS_di-mu");
  hDyVsPt[1]->Add(hDyVsPt[3],-1);
  hDyVsPt[1]->SetTitle(Form("%s: #Deltay of muon tracks (US-LS, %1.1f<M_{#mu#mu}<%1.1f)",trigName[kTrigType],low_mass,high_mass));
  hDyVsPt[0]->SetTitle(Form("%s: #Deltay of tracks matched to MTD;p_{T} (GeV/c);#Deltay (cm)",trigName[kTrigType]));

  TH1F *hDeltaY[4];
  for(Int_t i=0; i<4; i++)
    {
      hDeltaY[i] = (TH1F*)hDyVsPt[i]->ProjectionY(Form("hDy_%d_clone",i));
      hDeltaY[i]->Sumw2();
      hDeltaY[i]->SetLineColor(color[i]);
      hDeltaY[i]->SetMarkerColor(color[i]);
      hDeltaY[i]->GetXaxis()->SetRangeUser(-50,50);
      hDeltaY[i]->SetMaximum(250);
    }
  hDeltaY[1]->SetMarkerStyle(21);
  c = draw1D(hDeltaY[1],Form("%s: #Deltay of muon tracks (%s)",trigName[kTrigType],run_type),kFALSE,kTRUE);
  hDeltaY[2]->Draw("sames HIST");
  hDeltaY[3]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.62,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("%1.1f < M_{#mu#mu} < %1.1f (GeV/c^{2})",low_mass,high_mass));
  leg->AddEntry(hDeltaY[2],"Unlike-sign","L");
  leg->AddEntry(hDeltaY[3],"Like-sign","L");
  leg->AddEntry(hDeltaY[1],"US-LS","P");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDy_US_LS.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDy_US_LS.pdf",run_type,run_cfg_name.Data()));
    }

  TList *list = new TList;
  for(Int_t i=0; i<2; i++)
    {
      TH1F *htmp = (TH1F*)hDyVsPt[i]->ProjectionY(Form("hDy_%d",i));
      htmp->Sumw2();
      htmp->Scale(1./htmp->Integral());
      htmp->SetMaximum(4*htmp->GetMaximum());
      list->Add(htmp);
    }
  TString legName[4] = {"All matched tracks", "Muon (US-LS, 3.0<M_{#mu#mu}<3.2)","Unlike-sign pair","Like-sign pair"};
  c = drawHistos(list,"Deltay",Form("%s: #Deltay of tracks;#Deltay (cm)",trigName[kTrigType]),kTRUE,-100,100,kFALSE,1e-5,0.13,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sCompareDy.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sCompareDy.pdf",run_type,run_cfg_name.Data()));
    }

  // Compare with embedding and cosmic
  TFile *fmc, *fcosmic;
  if(year==2013)
    {
      fmc = TFile::Open("output/Run13.pp500.jpsi.EmbedQA.MC.root","read");
      fcosmic = TFile::Open("../Calibration/Xinjie/Run13Position.root","read");
    }
  TH2F *hDyPt[4];
  hDyPt[0] = (TH2F*)hDyVsPt[0]->Clone("hDyVsPt_Data_MthTrk");
  hDyPt[1] = (TH2F*)hDyVsPt[1]->Clone("hDyVsPt_Data_Muon");
  hDyPt[2] = (TH2F*)fmc->Get("hDyVsRcTrkPt_MCreco_true_di_mu");
  hDyPt[3] = (TH2F*)fcosmic->Get("hdYPt");
  const char *hName[4] = {"Data_MthTrk","Data_MuonTrk","Embed_Muon","Cosmic_Muon"};
  const char *hTitle[4] = {"all matched track in data","muon tracks in data","embedded muon tracks","cosmic rays"};
  TH1F *hDy[4];
  for(int i=0; i<4; i++)
    {
      hDyPt[i]->SetTitle(Form("#Deltay vs p_{T} of %s (Run13 pp500);p_{T} (GeV/c);#Deltay (cm)",hTitle[i]));
      hDyPt[i]->GetXaxis()->SetRangeUser(0,20);
      c = draw2D(hDyPt[i]);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDyVsPt_%s.png",run_type,run_cfg_name.Data(),hName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDyVsPt_%s.pdf",run_type,run_cfg_name.Data(),hName[i]));
	}

      hDy[i] = (TH1F*)hDyPt[i]->ProjectionY(Form("hDy_%s",hName[i]));
      hDy[i]->Sumw2();
      hDy[i]->Scale(1./hDy[i]->Integral());
      hDy[i]->SetLineColor(color[i]);
      hDy[i]->SetMarkerColor(color[i]);
      hDy[i]->SetMarkerStyle(21);
    }
  hDy[0]->GetXaxis()->SetRangeUser(-50,50);
  hDy[0]->GetYaxis()->SetRangeUser(0,0.12);
  c = draw1D(hDy[0],"#Deltay distribution of matched track-hit pairs");
  hDy[1]->Draw("samesP");
  hDy[2]->Draw("sames HIST");
  hDy[3]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.62,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hDy[0],"Data: all matched tracks","P");
  leg->AddEntry(hDy[1],"Data: J/psi muon (US-LS)","P");
  leg->AddEntry(hDy[2],"Embedded muon","L");
  leg->AddEntry(hDy[3],"Cosmic ray","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDy_data_MC_cosmic.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDy_data_MC_cosmic.pdf",run_type,run_cfg_name.Data()));
    }
}

//================================================
void DeltaZ(const Int_t save = 1)
{
  TH2F *hDzVsPt[2];
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhMuonDzDy_%s",trigName[kTrigType]));
  hn->GetAxis(0)->SetRangeUser(low_mass+0.001,high_mass-0.001);
  for(int i=0; i<2; i++)
    {
      hn->GetAxis(1)->SetRange(i+1,i+1);
      hDzVsPt[i] = (TH2F*)hn->Projection(3,2);
      hDzVsPt[i]->SetName(Form("hzyVsPt_di-mu_%d",i));
      cout << "# of di-muon pairs: " << hDzVsPt[i]->GetEntries()/2 << endl;
    }

  // Compare with embedding and cosmic
  TFile *fmc, *fcosmic;
  if(year==2013)
    {
      fmc = TFile::Open("output/Run13.pp500.jpsi.EmbedQA.MC.root","read");
      fcosmic = TFile::Open("../Calibration/Xinjie/Run13Position.root","read");
    }
  TH2F *hDzPt[3];
  const char *hName[4] = {"Data_MthTrk","Data_MuonTrk","Embed_Muon","Cosmic_Muon"};
  const char *hTitle[4] = {"all matched track in data","muon tracks in data","embedded muon tracks","cosmic rays"};
  hDzPt[0] = (TH2F*)f->Get(Form("mhDzVsPt_%s",trigName[kTrigType]));

  hDzPt[1] = (TH2F*)hDzVsPt[0]->Clone("hDzVsPt_US-LS_di-mu");
  hDzPt[1]->Add(hDzVsPt[1],-1);

  hDzPt[2] = (TH2F*)fmc->Get("hDzVsRcTrkPt_MCreco_true_di_mu");

  TH1F *hDz[4];
  for(int i=0; i<3; i++)
    {
      hDzPt[i]->SetTitle(Form("#Deltaz vs p_{T} of %s (Run13 pp500);p_{T} (GeV/c);#Deltaz (cm)",hTitle[i]));
      hDzPt[i]->GetXaxis()->SetRangeUser(0,20);
      c = draw2D(hDzPt[i]);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDzVsPt_%s.png",run_type,run_cfg_name.Data(),hName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDzVsPt_%s.pdf",run_type,run_cfg_name.Data(),hName[i]));
	}

      hDz[i] = (TH1F*)hDzPt[i]->ProjectionY(Form("hDz_%s",hName[i]));
      hDz[i]->Sumw2();
    }

  for(int i=0; i<30; i++)
    {
      TH1F *h = (TH1F*)fcosmic->Get(Form("hZvsTBL_%d",i));
      if(i==0) hDz[3] = (TH1F*)h->Clone(Form("hDz_%s",hName[2]));
      else     hDz[3]->Add(h);
    }

  for(int i=0; i<4; i++)
    {
      hDz[i]->Scale(1./hDz[i]->GetBinContent(hDz[i]->FindFixBin(0)));
      hDz[i]->SetLineColor(color[i]);
      hDz[i]->SetMarkerColor(color[i]);
      hDz[i]->SetMarkerStyle(21);
    }
  hDz[0]->GetXaxis()->SetRangeUser(-100,100);
  hDz[0]->GetYaxis()->SetRangeUser(0,1.4);
  c = draw1D(hDz[0],"#Deltaz distribution of matched track-hit pairs");
  hDz[1]->Draw("sames P");
  hDz[2]->Draw("sames HIST");
  hDz[3]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.62,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hDz[0],"Data: all matched tracks","P");
  leg->AddEntry(hDz[1],"Data: J/psi muon (US-LS)","P");
  leg->AddEntry(hDz[2],"Embedded muon","L");
  leg->AddEntry(hDz[3],"Cosmic ray","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDz_data_MC_cosmic.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%sDz_data_MC_cosmic.pdf",run_type,run_cfg_name.Data()));
    }
}
