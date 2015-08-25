//================================================
void moveHisto()
{
  TFile *fin = TFile::Open("output/bk.Run14.AuAu200.jpsi.EmbedQA.MC.root","read");
  TH2F *h2 = (TH2F*)fin->Get("hMcDeltaTof_di_mu");
  THnSparseF *hn = (THnSparseF*)fin->Get("mhMcTofQA_di_mu");
  TFile *fout = TFile::Open("output/Run14.AuAu200.jpsi.EmbedQA.MC.root","update");
  h2->Write("",TObject::kOverwrite);
  hn->Write("",TObject::kOverwrite);
  fout->Close();
  fin->Close();
}
