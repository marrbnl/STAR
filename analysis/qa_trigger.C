const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;
TFile *f;

const char *run_config = "trigger.";
const Bool_t iPico = 1;
const int year = 2013;
TString run_cfg_name;

//================================================
void qa_trigger()
{
  gStyle->SetOptStat(0);

  if(year==2013)
    {
      run_type = "Run13_pp500";

      if(iPico)
	f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else
	f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";

      if(iPico)
	f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else
	f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }

  signal();
  //QTcorrelation();
  //MtdVpdTacDiff();
  //MuonPairTacDiff();
}

//================================================
void MuonPairTacDiff(const Int_t save = 1, const Int_t saveHisto = 0)
{
  TFile *fdata = TFile::Open(Form("./output/Run13.pp500.jpsi.muon.root"),"read");
  TH1F *hTacDiffPair[2];
  hTacDiffPair[0] = (TH1F*)f->Get("mhMtdFastTacPairDiff_di-muon");

  THnSparseF *hnTacDiff = (THnSparseF*)fdata->Get("mhMuonPairTacDiff_di_mu");
  hnTacDiff->GetAxis(0)->SetRangeUser(low_mass+0.001,high_mass-0.001);
  hnTacDiff->GetAxis(1)->SetRange(1,1);
  TH1F *hTacDiffPair[1] = (TH1F*)hnTacDiff->Projection(2);
  hTacDiffPair[1]->SetName("hTacDiffPair_di-muon_US");
  cout << "# if unlike-sign muon pairs: " << hTacDiffPair[1]->GetEntries() << endl;
  hnTacDiff->GetAxis(1)->SetRange(2,2);
  TH1F *hLS = (TH1F*)hnTacDiff->Projection(2);
  hLS->SetName("hTacDiffPair_di-muon_LS");
  cout << "# if like-sign muon pairs: " << hLS->GetEntries() << endl;
  hTacDiffPair[1]->Add(hLS,-1);

  TList *list = new TList;
  for(Int_t i=0; i<2; i++)
    {
      TH1F *htmp = (TH1F*)hTacDiffPair[i]->Clone(Form("%s_clone",hTacDiffPair[i]->GetName()));
      htmp->Rebin(8);
      htmp->Scale(1./htmp->Integral());
      list->Add(htmp);
    }
  TString legName[2] = {"Fastest two QT signals","Muon pairs (US-LS, 3.0<M_{#mu#mu}<3.2)"};
  c = drawHistos(list,"TacDiffMtdVpd",Form("%s: TAC different between muon pairs;#Deltatac",trigName[kTrigType]),kFALSE,-3000,-1000,kTRUE,1e-5,0.13,kFALSE,kTRUE,legName,kTRUE,"",0.3,0.4,0.6,0.8,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/JpsiMuon_PairTacDiff.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/JpsiMuon_PairTacDiff.pdf",run_type));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run13.pp500.muon.timing.root","recreate");
      hTacDiffPair[0]->SaveAs("TacDiff_FastestTwoQTsingal",TObject::kOverwrite);
      hTacDiffPair[1]->SaveAs("TacDiff_Jpsi_muonPair",TObject::kOverwrite);
    }
}

//================================================
void MtdVpdTacDiff(const Int_t save = 1)
{
  TFile *fdata = TFile::Open(Form("./output/Run13.pp500.jpsi.muon.root"),"read");
  TH1F *hTacDiffMtdVpd[3];
  hTacDiffMtdVpd[0] = (TH1F*)f->Get("mhMtdVpdTacDiff_di-muon");
  hTacDiffMtdVpd[1] = (TH1F*)f->Get("mhMtdVpdTacDiffMatched_di-muon");

  THnSparseF *hnTacDiff = (THnSparseF*)fdata->Get("mhMuonDiffVPD_di_mu");
  hnTacDiff->GetAxis(0)->SetRangeUser(low_mass+0.001,high_mass-0.001);
  hnTacDiff->GetAxis(1)->SetRange(1,1);
  TH1F *hTacDiffMtdVpd[2] = (TH1F*)hnTacDiff->Projection(2);
  hTacDiffMtdVpd[2]->SetName("mhMuonMtdVpdTacDiff_di-muon_US");
  cout << "# if unlike-sign muons: " << hTacDiffMtdVpd[2]->GetEntries() << endl;
  hnTacDiff->GetAxis(1)->SetRange(2,2);
  TH1F *hLS = (TH1F*)hnTacDiff->Projection(2);
  hLS->SetName("mhMuonMtdVpdTacDiff_di-muon_LS");
  cout << "# if like-sign muons: " << hLS->GetEntries() << endl;
  hTacDiffMtdVpd[2]->Add(hLS,-1);

  TList *list = new TList;
  for(Int_t i=0; i<3; i++)
    {
      hTacDiffMtdVpd[i]->Rebin(8);
      hTacDiffMtdVpd[i]->Scale(1./hTacDiffMtdVpd[i]->Integral());
      list->Add(hTacDiffMtdVpd[i]);
    }
  TString legName[3] = {"All MTD hits","Matched MTD hits","Muons (US-LS, 3.0<M_{#mu#mu}<3.2)"};
  c = drawHistos(list,"TacDiffMtdVpd",Form("%s: TAC different between MTD and VPD;tac_{MTD} - tac_{VPD}",trigName[kTrigType]),kTRUE,-3000,-1000,kTRUE,1e-5,0.08,kFALSE,kTRUE,legName,kTRUE,"",0.15,0.25,0.6,0.88,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/JpsiMuon_TacDiffMtdVpd.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/JpsiMuon_TacDiffMtdVpd.pdf",run_type));
    }
}


//================================================
void signal(const Int_t save = 0)
{
  gStyle->SetOptStat(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
  
  TH1F *hNMtdHits = (TH1F*)f->Get("hNMtdHits_di-muon");
  c = draw1D(hNMtdHits,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfMtdHits.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfMtdHits.pdf",run_type));
    }

  TH1F *hNQTsignals = (TH1F*)f->Get("hNQTsignals_di-muon");
  c = draw1D(hNQTsignals,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfQTsignals.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfQTsignals.pdf",run_type));
    }

  TH1F *hNMIXsignals = (TH1F*)f->Get("hNMIXsignals_di-muon");
  c = draw1D(hNMIXsignals,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfMT101signals.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfMT101signals.pdf",run_type));
    }

  TH1F *hNMuons = (TH1F*)f->Get("hNMuons_di-muon");
  c = draw1D(hNMuons,"",kTRUE,kFALSE);
  cout << "2 muons: " << hNMuons->GetBinContent(3)/hNMuons->GetEntries() << endl;
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfTF201signals.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfTF201signals.pdf",run_type));
    }

  TH1F *hNTrigMtdHits = (TH1F*)f->Get("hNTrigMtdHits_di-muon");
  c = draw1D(hNTrigMtdHits,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfTrigMtdHits.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfTrigMtdHits.pdf",run_type));
    }

  TH1F *hNTrigGoodMtdHits = (TH1F*)f->Get("hNTrigGoodMtdHits_di-muon");
  c = draw1D(hNTrigGoodMtdHits,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfTrigGoodMtdHits.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/NumberOfTrigGoodMtdHits.pdf",run_type));
    }
}

//================================================
void QTcorrelation(const Int_t save = 0)
{
  gStyle->SetOptStat(0);
  
  TH2F *hMtdHitVsQtSignal = (TH2F*)f->Get("hMtdHitVsQtSignal_di-muon");
  hMtdHitVsQtSignal->SetXTitle("QT signal");
  hMtdHitVsQtSignal->SetYTitle("MTD hit");
  c = draw2D(hMtdHitVsQtSignal,"");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/MtdHitVsQtSignal.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_trigger/MtdHitVsQtSignal.pdf",run_type));
    }
  cout << "MTD hit: " << hMtdHitVsQtSignal->Integral(1,1,2,33)/hMtdHitVsQtSignal->Integral(1,33,2,33) << endl;

  // TH2F *hGoodMtdHitVsQtSignal = (TH2F*)f->Get("hGoodMtdHitVsQtSignal_di-muon");
  // hGoodMtdHitVsQtSignal->SetXTitle("QT signal");
  // hGoodMtdHitVsQtSignal->SetYTitle("MTD hit");
  // c = draw2D(hGoodMtdHitVsQtSignal,"");
  // if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/GoodMtdHitVsQtSignal.png");
  // cout << "MTD hit: " << hGoodMtdHitVsQtSignal->Integral(1,1,2,33)/hGoodMtdHitVsQtSignal->Integral(1,33,2,33) << endl;
}
