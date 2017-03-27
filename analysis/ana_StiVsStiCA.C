
//================================================
void ana_StiVsStiCA()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  Run14vs16();
}

//================================================
void Run14vs16()
{
  const int nHisto = 2;
  const char* year[nHisto] = {"14","16"};
  TFile *fdata[nHisto];
  fdata[0] = TFile::Open("Pico.Run14.AuAu200.jpsi.root","read");
  fdata[1] = TFile::Open("Pico.Run14.AuAu200.jpsi.root","read");
}
