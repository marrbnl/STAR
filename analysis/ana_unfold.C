//------------------------------------------------
void ana_unfold()
{
  //gSystem->Load("~/Work/RooUnfold-1.1.1/libRooUnfold.so");
  Int_t nbin = 2;
  Double_t xbin[3] = {0, 1, 2};
  gRandom->SetSeed(0);

  TH1D *htruth = new TH1D("htruth","htruth;p_{T} (GeV/c);entries", nbin, xbin);
  htruth->SetBinContent(1, 6000);
  htruth->SetBinContent(2, 16000);
  htruth->Sumw2();
  htruth->Scale(1./htruth->Integral());
  htruth->SetLineColor(kBlack);
  htruth->SetLineWidth(1);
  c = draw1D(htruth);

  TH1D *hmeasured = new TH1D("hmeasured","hmeasured;p_{T} (GeV/c);entries", nbin, xbin);
  hmeasured->SetBinContent(1, 12000);
  hmeasured->SetBinContent(2, 10000);
  hmeasured->Sumw2();
  hmeasured->SetLineColor(kRed);
  hmeasured->SetMarkerColor(kRed);
  hmeasured->SetMarkerStyle(21);

  TH2D *hresponse = new TH2D("hresponse","hresponse;p_{T}^{measured};p_{T}^{truth};entries", nbin, xbin, nbin, xbin);
  hresponse->SetBinContent(1, 1, 4000000);
  hresponse->SetBinContent(1, 2, 12000000);
  hresponse->SetBinContent(2, 1, 6000000);
  hresponse->SetBinContent(2, 2, 8000000);
  hresponse->Sumw2();
  draw2D(hresponse);


  RooUnfoldResponse response (0x0,0x0,hresponse,"unfold_jet","unfold_jet");  
  RooUnfoldBayes unfoldBayes (&response, hmeasured, 4);
  TH1F *hUnfolded = unfoldBayes.Hreco(2);

  TList *list = new TList;
  list->Add(htruth);
  list->Add(hmeasured);
  list->Add(hUnfolded);

  TString legName1[] = {"Truth","Measured","Unfolded (Bayesian, n=4)"};
  c = drawHistos(list,"jet_pt_compare","Test unfolding",false,0,220,false,0,0,false,kTRUE,legName1,kTRUE,"MC",0.4,0.6,0.6,0.8);
}
