
//================================================
void ana_InvMass()
{
  gStyle->SetOptStat(0);
  InvMass();
}



//================================================
void InvMass(Int_t save = 0)
{
  TFile *f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));
  printf("single-muon events: %4.2e\n",hStat->GetBinContent(8));
  printf("e-muon      events: %4.2e\n",hStat->GetBinContent(9));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[kNtrig][3];
  TH1F *hInvM[kNtrig][3];
  for(Int_t i=0; i<kNtrig; i++)
    {
      for(Int_t j=0; j<3; j++)
	{
	  hnInvMass[i][j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[i]));
	  hInvM[i][j] = (TH1F*)hnInvMass[i][j]->Projection(0);
	  hInvM[i][j]->SetName(Form("%s_%s_InvMass",hName[j],trigName[i]));
	}
    }

  TH1F *hInvMass[3];
  for(Int_t j=0; j<3; j++)
    {
      for(Int_t i=0; i<1; i++)
	{
	  if(i==0) hInvMass[j] = (TH1F*)hInvM[i][j]->Clone(Form("%s_InvMass",hName[j]));
	  else     hInvMass[j]->Add(hInvM[i][j]);
	}
      //hInvMass[j]->Rebin(1);
    }

  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->GetXaxis()->SetRangeUser(2.5,3.5);
  hInvMass[0]->SetYTitle("Counts");

  TH1F *hDiff = (TH1F*)hInvMass[0]->Clone("InvMass_US_minus_LS");
  hDiff->Sumw2();
  hDiff->Add(hInvMass[1],-1);
  hDiff->SetYTitle("US-LS");
  // for(Int_t ibin=1; ibin<=hDiff->GetNbinsX(); ibin++)
  //   {
  //     printf("Bin %d: error = %f, US = %f, LS = %f\n",ibin,hDiff->GetBinError(ibin),hInvMass[0]->GetBinContent(ibin),hInvMass[1]->GetBinContent(ibin));
  //   }
  c = draw1D(hDiff,"Invariant mass of di-muon pairs (MTD trigger)");
  TLine *line = GetLine(2.5,1,3.5,1,1);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_AuAu200_Run14.png",run_config));


  c = draw1D(hInvMass[0],"Invariant mass of di-muon pairs (MTD trigger)");
  hInvMass[1]->Draw("HIST sames");
  TLegend *leg = new TLegend(0.6,0.6,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hInvMass[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Like sign (++)+(--)","L");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_AuAu200_Run14.png",run_config));
}
