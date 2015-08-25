const char *run_config = "";
const int year = 2013;
TString run_cfg_name;
const double pt1_cut = 1.5;
const double pt2_cut = 1.0;
const Double_t low_mass = 2.9;
const Double_t high_mass = 3.3;

TFile *f;

//================================================
void ana_CrossSection()
{
  ana();
}

//================================================
void ana()
{
  const int nbins = 6;
  const double xbins[nbins+1] = {0,1,2,3,4,6,10};
  TH1F *hJpsiCounts = new TH1F("hJpsi_ptRaw","p_{T} distributin of J/#psi;p_{T} (GeV/c)",nbins,xbins);
  double count[nbins] = {158.094, 207.805, 264.733, 198.954, 90.856, 36.1364};
  double error[nbins] = {37.8405, 44.5218, 40.0907, 29.1041, 12.945, 4.07319};

  for(int bin=1; bin<=nbins; bin++)
    {
      hJpsiCounts->SetBinContent(bin,count[bin-1]);
      hJpsiCounts->SetBinError(bin,error[bin-1]);
    }
  scaleHisto(hJpsiCounts, 1, 1,kTRUE,kFALSE, kTRUE);
  draw1D(hJpsiCounts);

  TFile *fEff = TFile::Open(Form("Rootfiles/Run13.pp500.jpsi.Eff.root"),"read");
  TH1F *hEff = (TH1F*)fEff->Get("JpsiPt_MTDreco_pt11.5_rebin_eff");
  draw1D(hEff);

  TH1F *hXsec = (TH1F*)hJpsiCounts->Clone("hJpsi_Xsec");
  for(int bin=1; bin<=nbins; bin++)
    {
      double eff = hEff->GetBinContent(bin);
      double pt = hJpsiCounts->GetBinCenter(bin);
      double phase = 2*pi*1.6;
      double scale = eff*pt*phase;
      //cout << eff << "  " << pt << "  " << phase << endl;

      hXsec->SetBinContent(bin,hXsec->GetBinContent(bin)/scale);
      hXsec->SetBinError(bin,hXsec->GetBinError(bin)/scale);
    }

  double lumi = 7.678 * 1e3; // nb-1
  hXsec->Scale(1./lumi);

  draw1D(hXsec,"",kTRUE,kTRUE);

  TFile *fHT = TFile::Open("Rootfiles/Spectrum_in_bin.root","read");
  TGraphErrors	*gSpectrum = (TGraphErrors*)fHT->Get("gall");
  gSpectrum->Draw("sames PE");
}
