
//================================================
void make_TMVA()
{
  gStyle->SetOptStat(0);

  //Run14_NJpsi();
  Run14_bkg();
}

//================================================
void Run14_bkg()
{
  TFile *fall = TFile::Open("output/Run14_AuAu200.jpsi.root ", "read");
  TH1F *hStatAll = (TH1F*)fall->Get("hEventStat");
  hStatAll->SetName("hStatAll");

  TFile *f = TFile::Open("output/Run14_AuAu200.PreCut.root", "read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("%s\n",f->GetName());
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  THnSparseF *hnInvMass = (THnSparseF*)f->Get("mhBkgLSPos_di_mu");
  hnInvMass->Add((THnSparseF*)f->Get("mhBkgLSNeg_di_mu"));
  hnInvMass->GetAxis(0)->SetRangeUser(3.0, 3.2);
  hnInvMass->GetAxis(2)->SetRangeUser(1.5+0.01,100);
  hnInvMass->GetAxis(3)->SetRangeUser(1.3+0.01,100);
  hnInvMass->GetAxis(4)->SetRangeUser(1, 16);
  
  TH1F *hPairPt = (TH1F*)hnInvMass->Projection(1);
  c = draw1D(hPairPt,"Pair p_{T} distribution");
  gPad->SetLogy();
  double lowPt[3] = {0, 2, 4};
  double upPt[3] = {2, 4, 8};
  for(int i=0; i<3; i++)
    {
      double yield = hPairPt->Integral(hPairPt->FindBin(lowPt[i]+0.1), hPairPt->FindBin(upPt[i]-0.1));
      printf("[i] %1.0f < pt < %1.0f: yield = %1.0f, all = %1.0f\n",lowPt[i],upPt[i],yield,yield*hStatAll->GetBinContent(3)/hStat->GetBinContent(3));
    }
}


//================================================
void Run14_NJpsi()
{
  const int nbins = nPtBins_pt -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low_pt[i+1];
  xbins[nbins] = ptBins_high_pt[nbins];

  // PID efficiency from embedding+data
  TFile *femb = TFile::Open("Rootfiles/Run14_AuAu200.EmbJpsiEff.pt1.5.pt1.3.root", "read");
  TH1F *h1tmp = (TH1F*)femb->Get("hJpsiPt_MtdMth_cent0080");
  TH1F *hMtdMth = (TH1F*)h1tmp->Rebin(nbins, "hJpsiPt_MtdMth", xbins);
  h1tmp = (TH1F*)femb->Get("hJpsiPt_MuonPid_cent0080");
  TH1F *hMuonPid = (TH1F*)h1tmp->Rebin(nbins, "hJpsiPt_MuonPid", xbins);
  TH1F *hPidEff = (TH1F*)hMuonPid->Clone("hJpsiPidEffVsPt");
  hPidEff->Divide(hMtdMth);
  c = draw1D(hPidEff);

  // raw yield after muon PID
  TFile *fyield = TFile::Open("Rootfiles/Run14_AuAu200.JpsiYield.pt1.5.pt1.3.Dca3cm.root", "read");
  TH1F *hyield = (TH1F*)fyield->Get("Jpsi_FitYield_cent0080");
  c = draw1D(hyield);

  for(int bin=1; bin<=hyield->GetNbinsX(); bin++)
    {
      printf("[i] pt = %2.1f, yield = %2.1f, yield/eff = %2.1f\n",hyield->GetBinCenter(bin), hyield->GetBinContent(bin), hyield->GetBinContent(bin)/hPidEff->GetBinContent(bin));
    }
}
