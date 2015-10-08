void ana_ups()
{
  ppdata();
}

//================================================
void ppdata()
{
  TFile *fin = TFile::Open("Rootfiles/ups.root","read");
  TTree *tree = (TTree*)fin->Get("tree");
  tree->Draw("ups_pt>>hups(100,0,10)","abs(ups_y)<0.5");
  TH1F *hups = (TH1F*)gDirectory->Get("hups");
  draw1D(hups);
  TFile *fout = TFile::Open("Rootfiles/Upsion.pp200.root","recreate");
  hups->Write("Upsilon_pp200_midRap");
}
