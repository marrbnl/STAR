void ana_Upsilon()
{
  //AuAu();
  //ppdata();
  //Xinjie();
  ptShape();
}


//================================================
void Xinjie(const int saveHisto = 0)
{
  TFile *fin = TFile::Open("Rootfiles/160606.all.root","read");
  const char *hName[3] = {"mhUpsInfoUL","mhUpsInfoLSPos","mhUpsInfoLSNeg"};
  THnSparseF *hnInvMass[3] = {0x0};
  TH1F *hInvMass[3] = {0x0};
  
  // same event
  char name[512];
  for(Int_t j=0; j<3; j++) // pair type
    { 
      hnInvMass[j] = (THnSparseF*)fin->Get(hName[j]);
      hnInvMass[j]->GetAxis(3)->SetRange(0,13);
      hnInvMass[j]->GetAxis(4)->SetRange(0,6);
      hnInvMass[j]->GetAxis(5)->SetRange(0,23);
      hnInvMass[j]->GetAxis(6)->SetRange(0,14);

      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("InvMassHft_jpsi_%d",j));
      hInvMass[j]->SetTitle("");
      //hInvMass[j]->Rebin(2);
      hInvMass[j]->Sumw2();
    }
  //hInvMass[1]->Add(hInvMass[2]);

  TH1F *hSeUL  = (TH1F*)hInvMass[0]->Clone("hSeUL");
  TH1F *hSeLS  = (TH1F*)hInvMass[1]->Clone("hSeLS");

  hSeUL->SetMarkerStyle(21);
  hSeUL->SetMarkerColor(2);
  hSeUL->SetLineColor(2);
  hSeUL->GetXaxis()->SetRangeUser(8,12);
  c = draw1D(hSeUL,Form("Dimuon events: invariant mass distribution"),false);
  hSeLS->Draw("sames HIST");

  TFile *fout = TFile::Open("Rootfiles/sQM2016.Upsilon.root","recreate");
  hInvMass[0]->Write("Upsilon_UL_cent0080");
  hInvMass[1]->Write("Upsilon_LS_Pos_cent0080");
  hInvMass[2]->Write("Upsilon_LS_Neg_cent0080");
  fout->Close();
}

//================================================
void AuAu(const int savePlot = 0)
{
  gStyle->SetOptStat(0);
  const int ipt = 0, icent = 0;

  TFile *fin = TFile::Open("output/Pico.Run14.AuAu200.jpsi.root","read");
  TH1F *hStat = (TH1F*)fin->Get("hEventStat");
  printf("+++ check this +++\n");
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3] = {0x0};
  TH1F *hInvMass[3] = {0x0};
  
  // same event
  char name[512];
  for(Int_t j=0; j<3; j++) // pair type
    { 
      hnInvMass[j] = (THnSparseF*)fin->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(3)->SetRangeUser(1.5+0.01,100);
      hnInvMass[j]->GetAxis(4)->SetRangeUser(1.5+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRange(centBins_low[icent],centBins_high[icent]);
      hnInvMass[j]->GetAxis(1)->SetRangeUser(ptBins_low[ipt]+0.01,ptBins_high[ipt]-0.01);

      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("InvMassHft_jpsi_%d",j));
      hInvMass[j]->Rebin(25);
      hInvMass[j]->Sumw2();
    }
  hInvMass[1]->Add(hInvMass[2]);

  TH1F *hSeUL  = (TH1F*)hInvMass[0]->Clone("hSeUL");
  TH1F *hSeLS  = (TH1F*)hInvMass[1]->Clone("hSeLS");
 
  hSeUL->SetMarkerStyle(21);
  hSeUL->SetMarkerColor(2);
  hSeUL->SetLineColor(2);
  hSeUL->GetXaxis()->SetRangeUser(6,14);
  c = draw1D(hSeUL,Form("Dimuon events: %1.1f < p_{T} < %1.1f GeV/c (%s%%)",ptBins_low[ipt],ptBins_high[ipt],cent_Name[icent]),true);
  hSeLS->Draw("sames HIST");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/AuAu200_UpsilonInvMass.pdf",run_type));

  TH1F *hSeDiff  = (TH1F*)hInvMass[0]->Clone("hSeDiff");
  hSeDiff->SetMarkerStyle(20);
  hSeDiff->GetXaxis()->SetRangeUser(6,14);
  hSeDiff->Add(hSeLS,-1); 
  c = draw1D(hSeDiff,Form("Dimuon events: %1.1f < p_{T} < %1.1f GeV/c (%s%%, US-LS)",ptBins_low[ipt],ptBins_high[ipt],cent_Name[icent]),false);
  TLine *line = GetLine(6,0,14,0);
  line->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Upsilon/AuAu200_UpsilonInvMass_USminusLS.pdf",run_type));

}

//================================================
void ppdata()
{
  TFile *fin = TFile::Open("Rootfiles/Upsilon/ups.root","read");
  TTree *tree = (TTree*)fin->Get("tree");
  tree->Draw("ups_pt>>hups(100,0,10)","abs(ups_y)<0.5");
  TH1F *hups = (TH1F*)gDirectory->Get("hups");
  draw1D(hups);
  TFile *fout = TFile::Open("Rootfiles/Upsion.pp200.root","recreate");
  hups->Write("Upsilon_pp200_midRap");
}


//================================================
void ptShape()
{
  TF1 *fBol = new TF1("Boltzmann","x/(exp(x/[0]+1))",0,10);
  fBol->SetParameter(0,1.11);
  fBol->SetNpx(1000);
  TH1F *hBol = fBol->GetHistogram();
  TCanvas *c = new TCanvas("Boltzmann","Boltzmann",800,600);
  fBol->Draw();

  TFile *fin = TFile::Open("Rootfiles/Upsilon/Upsion.pp200.root","update");
  TH1F *hups = (TH1F*)fin->Get("Upsilon_pp200_midRap");
  hups->Sumw2();
  hups->Scale(1./hups->Integral());
  c = draw1D(hups);
  hBol->Scale(1./hBol->Integral());
  hBol->Draw("sames");
  fBol->Write("Upsilon_pp200_midRap_fBoltzmann");
  hBol->Write("Upsilon_pp200_midRap_hBoltzmann");
  
}

