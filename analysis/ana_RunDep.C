TFile *f;

//================================================
void ana_RunDep()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  MtdVpdTacDiff();
}

//================================================
void MtdVpdTacDiff(const int savePlot = 1)
{
  TFile *f = TFile::Open("output/Pico.Run14.AuAu200.jpsi.RunDep.root","read");
  THnSparseF *hn = (THnSparseF*)f->Get("mhRunDepMtdVpdDiffMuon_di_mu");
  hn->GetAxis(2)->SetRangeUser(2.5,100);

  const int nruns = 4;
  const int runs[nruns] = {15117051, 15137027, 15159044, 15164023};
  const double zdccorr[nruns] = {57.1, 57.1, 56.5, 50.4};
  const double mtdrate[nruns] = {753.68, 850.73, 854, 806.31};

  TCanvas *c = new TCanvas("TacDiff","TacDiff",1100,750);
  c->Divide(2,2);
  TH1F *hTacDiff[nruns];
  TF1 *func[nruns];
  for(int i=0; i<nruns; i++)
    {
      int bin = runs[i] - 15074000 + 1; 
      hn->GetAxis(0)->SetRange(bin,bin);
      hTacDiff[i] = (TH1F*)hn->Projection(1);
      hTacDiff[i]->SetName(Form("hTacDiff_%d",runs[i]));
      c->cd(i+1);
      hTacDiff[i]->GetXaxis()->SetRangeUser(750,850);
      hTacDiff[i]->SetTitle("");
      func[i] = new TF1(Form("func_%d",runs[i]),"gaus",790,810);
      hTacDiff[i]->Fit(func[i],"IR0");
      hTacDiff[i]->Draw();
      func[i]->SetLineColor(2);
      func[i]->Draw("sames");
      TPaveText *t1 = GetTitleText("MTD-VPD tac difference for muon candidates (p_{T}>2)",0.045);
      t1->Draw();
      TPaveText *t1 = GetPaveText(0.12,0.25,0.6,0.8,0.045,62);
      t1->SetTextAlign(11);
      t1->AddText(Form("Run %d",runs[i]));
      t1->AddText(Form("ZDCcorr %3.1fkHz",zdccorr[i]));
      t1->AddText(Form("Dimuon rate %4.1fHz",mtdrate[i]));
      t1->Draw();
    }

  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_RunDep/MtdVpdTacDiff.pdf",run_type));
}
