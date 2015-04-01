TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;
const Double_t pt1_cut = 1.5;
const Double_t pt2_cut = 1.0;

const char *run_config = "muon";
const Bool_t iPico = kTRUE;
TString run_cfg_name;

//================================================
void ana_PtBin()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.28);                
  gStyle->SetStatH(0.28);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  if(iPico)
    f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%s.root",run_config),"read");
  else
    f = TFile::Open(Form("./output/Run13.pp500.jpsi.%s.root",run_config),"read");

  if(run_config=="CutRanking")
    run_cfg_name = Form("%s.",run_config);
  else
    run_cfg_name = "";

  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  ptBins();
}

//================================================
void ptBins(Int_t save = 1)
{
  const Int_t nPtBins = 6;
  const Double_t ptBins[nPtBins+1] = {0,1,2,3,4,5,10};
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};

  THnSparseF *hnInvMass[3];
  TH1F *hJpsiPt[3];
  TH1F *hInvMass[3][nPtBins];
  TH1F *hSignal[nPtBins];

  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);

      for(Int_t i=0; i<nPtBins; i++)
	{
	  hnInvMass[j]->GetAxis(1)->SetRangeUser(ptBins[i]+0.01,ptBins[i+1]-0.01);
	  hInvMass[j][i] = (TH1F*)hnInvMass[j]->Projection(0);
	  hInvMass[j][i]->SetName(Form("%s_%s_InvMass_jpsi_%1.1f_%1.1f",hName[j],trigName[kTrigType],ptBins[i],ptBins[i+1]));
	  hInvMass[j][i]->Sumw2();
	  hInvMass[j][i]->Rebin(4);
	  hnInvMass[j]->GetAxis(1)->SetRange(0,-1);
	}
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }

  for(Int_t i=0; i<nPtBins; i++)
    {
      hInvMass[1][i]->Add(hInvMass[2][i]);
      hSignal[i] = (TH1F*)hInvMass[0][i]->Clone(Form("InvMass_US_minus_LS_jpsi_%1.1f_%1.1f",ptBins[i],ptBins[i+1]));
      hSignal[i]->Add(hInvMass[1][i],-1);
      hSignal[i]->SetYTitle("US-LS");
    }

  TCanvas *canvas[nPtBins/6+1];
  for(Int_t i=0; i<nPtBins; i++)
    {
      if(i%6==0) 
	{
	  canvas[i/6] = new TCanvas(Form("canvas_%d",i/6),Form("canvas_%d",i/6),1100,700);
	  canvas[i/6]->Divide(3,2);
	}
      canvas[i/6]->cd(i%6+1);

      hInvMass[0][i]->GetXaxis()->SetRangeUser(2.5,3.5);
      hInvMass[0][i]->SetMaximum(1.3*hInvMass[0][i]->GetMaximum());
      hInvMass[0][i]->SetMinimum(0);
      hInvMass[0][i]->SetMarkerColor(1);
      hInvMass[0][i]->SetLineColor(1);
      hInvMass[0][i]->SetTitle(";M_{#mu#mu} (GeV/c)^{2};Counts");

      hInvMass[1][i]->SetMarkerColor(4);
      hInvMass[1][i]->SetLineColor(4);

      hSignal[i]->SetMarkerStyle(21);
      hSignal[i]->SetMarkerColor(2);
      hSignal[i]->SetLineColor(2);

      hInvMass[0][i]->Draw("HIST");
      hInvMass[1][i]->Draw("sames HIST");
      hSignal[i]->Draw("sames");

      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T,J/#psi} < %1.1f GeV/c",ptBins[i],ptBins[i+1]),0.05);
      t1->Draw();
      
      Int_t low_bin = hInvMass[0][i]->GetXaxis()->FindFixBin(low_mass+0.001);
      Int_t high_bin = hInvMass[0][i]->GetXaxis()->FindFixBin(high_mass-0.001);
      Double_t nBackground = hInvMass[1][i]->Integral(low_bin,high_bin);
      Double_t nSignal = hInvMass[0][i]->Integral(low_bin,high_bin) - nBackground;
      t1 = GetPaveText(0.1,0.7,0.8,0.9,0.05);
      t1->AddText(Form("S/B = %1.0f/%1.0f = %1.2e:1",nSignal,nBackground,nSignal/nBackground));
      t1->Draw();
      if(i%6==0)
	{
	  t1 = GetPaveText(0.1,0.5,0.65,0.8,0.05);
	  t1->AddText(Form("%1.1f < M_{#mu#mu} < %1.1f",low_mass,high_mass));
	  t1->AddText(Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f",pt1_cut,pt2_cut));
	  t1->Draw();
	}
      if(i%6==5 && save)
	{
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_Jpsi_PtBins_pt1_%1.0f_pt2_%1.0f_Bin%d.pdf",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10,i/6));
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sInvMass_Jpsi_PtBins_pt1_%1.0f_pt2_%1.0f_Bin%d.png",run_type,run_cfg_name.Data(),pt1_cut*10,pt2_cut*10,i/6));
	}
    }


  // TString legName[3] = {"Unlike-sign","Like-sign","US-LS"};
  // c = drawHistos(list,"Jpsi_pt",Form("%s: p_{T} distribution of J/psi candidates;p_{T} (GeV/c);counts",trigName[kTrigType]),kFALSE,0,100,kFALSE,-2,2,kFALSE,kTRUE,legName,kTRUE,Form("%1.1f < M_{#mu#mu} < %1.1f GeV/c^{2}",low_mass,high_mass),0.55,0.75,0.55,0.8,kTRUE);
  // if(save) 
  //   {
  //     c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt1_cut,pt2_cut));
  //     c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_InvMass/%sJpsi_pt_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt1_cut,pt2_cut));
  //   }

}


//================================================
void upsilon(Int_t save = 1)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all di-muon events: %4.2e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.2e\n",hStat->GetBinContent(7));

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3];
  Double_t pt1_cut = 4.5, pt2_cut = 1.2;
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("%s_%s_InvMass",hName[j],trigName[kTrigType]));
      hInvMass[j]->Sumw2();
    }
  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->SetYTitle("Counts");

  TH1F *hDiff = (TH1F*)hInvMass[0]->Clone("InvMass_US_minus_LS");
  hDiff->Add(hInvMass[1],-1);
  hDiff->SetYTitle("US-LS");

  TLegend *leg = new TLegend(0.5,0.63,0.7,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass[0],"Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Like sign (++)+(--)","L");

  // Upsilon mass range
  TH1F *hInvMass_upsilon[2];
  for(Int_t i=0; i<2; i++)
    {
      hInvMass_upsilon[i] = (TH1F*)hInvMass[i]->Clone(Form("%s_upsilon",hInvMass[i]->GetName()));
      hInvMass_upsilon[i]->Rebin(10);
      hInvMass_upsilon[i]->GetXaxis()->SetRangeUser(8.5,11);
    }
  c = draw1D(hInvMass_upsilon[0],Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kTRUE);
  hInvMass_upsilon[1]->Draw("HIST sames");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_upsion_pt1_%1.1f_pt2_%1.1f.png",run_cfg_name.Data(),pt1_cut,pt2_cut));

  TH1F *hdiff_upsilon = (TH1F*)hInvMass_upsilon[0]->Clone(Form("%s_upsilon",hDiff->GetName()));
  hdiff_upsilon->Add(hInvMass_upsilon[1],-1);
  c = draw1D(hdiff_upsilon,Form("Au+Au %s: invariant mass of di-muon pairs%s;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kTRUE);
  gPad->SetGridy();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/InvMass/%s.InvMass_US-LS_upsilon_pt1_%1.1f_pt2_%1.1f.png",run_cfg_name.Data(),pt1_cut,pt2_cut));
  
}

