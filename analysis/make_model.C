
//================================================
void make_model()
{
  TBW();
}

//================================================
void TBW()
{
  const char *fname[4] = {"060","020","2040","4060"};
  const char *hname[4] = {"0060","0020","2040","4060"};

  TFile *fin[4];
  TH1F *hTBW[4];
  TH1F *hTBW2[4];
  for(int i=0; i<4; i++)
    {
      fin[i] = TFile::Open(Form("Rootfiles/Published/Zebo/TBW/dimuon/star_mumu_ee/TBW_Jpsi_STAR_mumu_ee_%s.root",fname[i]),"read");
      hTBW[i] = (TH1F*)fin[i]->Get("hFit5");
      hTBW[i]->SetName(Form("TBW_JpsiInvYield_AuAu200_cent%s",cent_Title[i]));
      hTBW2[i] = (TH1F*)hTBW[i]->Clone(Form("TBW_JpsiYield_AuAu200_cent%s",cent_Title[i]));
      for(int bin=1; bin<=hTBW2[i]->GetNbinsX(); bin++)
	{
	  double pt = hTBW2[i]->GetBinCenter(bin);
	  hTBW2[i]->SetBinContent(bin,hTBW2[i]->GetBinContent(bin)*pt);
	}
    }
  
  TFile *fout = TFile::Open("Rootfiles/models.root","update");
  for(int i=0; i<4; i++)
    {
      hTBW[i]->Write("",TObject::kOverwrite);
      hTBW2[i]->Write("",TObject::kOverwrite);
    }
}
