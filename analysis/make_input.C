//================================================
void make_input()
{
  Run14AuAu200();
}

//================================================
void Run14AuAu200()
{
  // input Jpsi shape
  TH1F *hInPutJpsiPt[3];
  TFile *fInJpsi = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");
  for(int i=0; i<3; i++)
    {
      TH1F *htmp = (TH1F*)fInJpsi->Get(Form("TBW_Jpsi_Yield_cent%s",cent_Title[i]));
      cout << htmp->GetNbinsX() << endl;
      for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
      	{
      	  htmp->SetBinError(bin,1e-10);
      	}
      TF1 *func = new TF1(Form("func_%d",i),"exp([0]+[1]*x+[2]*x*x)",5,10);
      htmp->Fit(func,"IR0");
      c = new TCanvas(Form("c%d",i),Form("c%d",i),800,600);
      gPad->SetLogy();
      func->SetRange(0,20);
      func->SetLineColor(2);
      func->Draw("");
      htmp->Draw("sames");

      hInPutJpsiPt[i] = new TH1F(Form("hInputJpsiShape_Cent%d",i),"",2000,0,20);
      for(int bin=1; bin<=hInPutJpsiPt[i]->GetNbinsX(); bin++)
       	{
	  double x = hInPutJpsiPt[i]->GetBinCenter(bin);
	  double jbin = htmp->FindFixBin(x);
	  if(x<10)
	    hInPutJpsiPt[i]->SetBinContent(bin,htmp->GetBinContent(jbin));
	  else
	    hInPutJpsiPt[i]->SetBinContent(bin,func->Eval(x));
	  hInPutJpsiPt[i]->SetBinError(bin,1e-20);
       	}
      hInPutJpsiPt[i]->Scale(1./hInPutJpsiPt[i]->GetBinContent(hInPutJpsiPt[i]->FindFixBin(1)));
      c = draw1D(hInPutJpsiPt[i],"",true,true);
    }

  // MTD response efficiency
  TFile *fMtd =  TFile::Open("Rootfiles/Run14.AuAu200.MtdResponseEff.root");
  TH1F *hMtdRespEffEmbed = (TH1F*)fMtd->Get("MtdResponseEff_embed");
  TF1 *func = new TF1("func","[0]*exp(-pow([1]/x,[2])-pow([3]/x/x,[4])-pow([5]/x/x/x,[6]))",1,10);
  func->SetParameters(1,0.1,1,0.1,1,0.1,1);
  hMtdRespEffEmbed->Fit(func,"IR0");
  hMtdRespEffEmbed->GetYaxis()->SetRangeUser(0,1.2);
  hMtdRespEffEmbed->SetMarkerStyle(20);
  c = draw1D(hMtdRespEffEmbed,"MTD response efficiency from embedding;p_{T} (GeV/c);Efficiency");
  func->SetLineColor(4);
  func->Draw("sames");
  func->SetNpx(1000);
  TH1F *hMtdRespEffEmb = (TH1F*)func->GetHistogram();
  hMtdRespEffEmb->SetName("hMtdRespEffEmbed");

  TF1 *fMtdRespEff[30][5];
  TH1F *hMtdRespEff[30][5];
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  fMtdRespEff[i][j] = (TF1*)fMtd->Get(Form("MtdRespEffvsPt_Bkl%d_Mod%d",i+1,j+1));
	  fMtdRespEff[i][j]->SetNpx(1000);
	  hMtdRespEff[i][j] = (TH1F*)fMtdRespEff[i][j]->GetHistogram();
	  hMtdRespEff[i][j]->SetName(Form("MtdRespEffCosmic_BL%d_Mod%d",i+1,j+1));
	}
    }

  // trigger efficiency
  TFile *fTrig =  TFile::Open("Rootfiles/Run14.AuAu200.MuonTrigEff.root");
  TH1F *hMuonTrigEff[4][3];
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<3; k++)
	{
	  hMuonTrigEff[i][k] = (TH1F*)fTrig->Get(Form("CombinedMuonTrigEff_cent%s%s",cent_Title[k+1],gTrgSetupName[i+1]));
	  hMuonTrigEff[i][k]->SetName(Form("MuonTrigEff_Cent%d_P%d",k,i));
	  hMuonTrigEff[i][k]->SetTitle("");
	}
    }

  // save
  TFile *fout =  TFile::Open("Rootfiles/Run14.AuAu200.Input.root","recreate");
  for(int i=0; i<3; i++)
    {
      hInPutJpsiPt[i]->Write();
    }
  hMtdRespEffEmb->Write();
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<3; k++)
	{
	  hMuonTrigEff[i][k]->Write();
	}
    }
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  fMtdRespEff[i][j]->Write();
	}
    }
  
}
