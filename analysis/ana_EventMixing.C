TFile *f;
Int_t trk_index = 0;
const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void ana_EventMixing()
{
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.EventMixing.%s.root",run_config),"read");
  //acceptance();
  subtraction();
}

//================================================
void subtraction(const Int_t save = 1)
{
  Double_t pt1_cut = 1.0, pt2_cut = 1.0;
  //Double_t pt1_cut = 2.0, pt2_cut = 1.2;

  TFile *fdata = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[3];
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)fdata->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("%s_%s_InvMass_WithCut",hName[j],trigName[kTrigType]));
      hInvMass[j]->Sumw2();
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }

  hInvMass[1]->Add(hInvMass[2]);
  hInvMass[0]->SetMarkerStyle(20);
  hInvMass[0]->SetMarkerColor(2);
  hInvMass[0]->SetLineColor(2);
  hInvMass[0]->SetYTitle("Counts");

  TH1F *hMixUS = (TH1F*)f->Get(Form("US_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt1_cut,pt2_cut,trigName[kTrigType]));
  TH1F *hMixLS = (TH1F*)f->Get(Form("LS_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt1_cut,pt2_cut,trigName[kTrigType]));
  Double_t scale = hInvMass[1]->Integral()/hMixLS->Integral();
  hMixLS->Scale(scale);
  TList *list = new TList;
  TString legName[2] = {"Data: like sign","Mixed event: like sign"};
  list->Add(hInvMass[1]);
  list->Add(hMixLS);
  c = sysCompare(list,Form("EventMixing_pt1_%1.1f_pt2_%1.1f",pt1_cut,pt2_cut),Form("Au+Au %s: invariant mass of dimuon pairs",trigName[kTrigType]),"(Mixed event)/Data",kFALSE,0,15,kFALSE,0.1,10,kTRUE,0,2,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.7,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_EventMixing/%s.Compare_LS_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  c = draw1D(hInvMass[0],Form("Au+Au %s: invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kTRUE);
  hInvMass[1]->Draw("HIST sames");
  hMixUS->SetLineColor(4);
  hMixUS->SetLineWidth(2);
  hMixUS->Scale(scale);
  hMixUS->Draw("sames HIST");
  leg = new TLegend(0.45,0.65,0.55,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass[0],"Data: Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Data: Like sign (++)+(--)","L");
  leg->AddEntry(hMixUS,"Mixed event: Like sign (++)+(--)","L");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_EventMixing/%s.Compare_US_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  hInvMass[0]->GetXaxis()->SetRangeUser(2.5,3.5);
  hInvMass[0]->SetName(Form("%s_zoom",hInvMass[0]->GetName()));
  c = draw1D(hInvMass[0],Form("Au+Au %s: invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2})",trigName[kTrigType]),kFALSE);
  hInvMass[1]->Draw("HIST sames");
  hMixUS->Draw("sames HIST");
  leg = new TLegend(0.15,0.65,0.3,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  leg->AddEntry(hInvMass[0],"Data: Unlike sign","PLE");
  leg->AddEntry(hInvMass[1],"Data: Like sign (++)+(--)","L");
  leg->AddEntry(hMixUS,"Mixed event: Like sign (++)+(--)","L");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_EventMixing/%s.Compare_US_jpsi_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));

  Double_t low_mass = 3.0, high_mass = 3.2;
  Int_t low_bin = hInvMass[0]->GetXaxis()->FindFixBin(low_mass);
  Int_t high_bin = hInvMass[0]->GetXaxis()->FindFixBin(high_mass);
  Double_t nBackground = hMixUS->Integral(low_bin,high_bin);
  Double_t nSignal = hInvMass[0]->Integral(low_bin,high_bin) - nBackground;
  TPaveText *signif = GetPaveText(0.15,0.45,0.12,0.25);
  signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
  signif->AddText(Form("S/B = %1.0f/%1.0f = %1.2e:1",nSignal,nBackground,nSignal/nBackground));

  TH1F *hdiff_jpsi = (TH1F*)hInvMass[0]->Clone(Form("US-LS_jpsi"));
  hdiff_jpsi->Add(hMixUS,-1);
  hdiff_jpsi->GetXaxis()->SetRangeUser(3.05,3.15);
  c = draw1D(hdiff_jpsi,Form("Au+Au %s: invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2});US_{Data}-US_{mixing}",trigName[kTrigType]),kFALSE);
  gPad->SetGridy();
  TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85,0.04);
  t1->AddText(Form("p_{T,1} > %1.1f, p_{T,2} > %1.1f GeV/c",pt1_cut,pt2_cut));
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_EventMixing/%s.Data-EventMixing_jpsi_pt1_%1.1f_pt2_%1.1f.png",run_config,pt1_cut,pt2_cut));
  hdiff_jpsi->SetName(Form("%s_rebin4",hdiff_jpsi->GetName()));
  hdiff_jpsi->Rebin(4);
  hdiff_jpsi->GetXaxis()->SetRangeUser(2.5,3.5);
  c = draw1D(hdiff_jpsi,Form("Au+Au %s: invariant mass of di-muon pairs;M_{#mu#mu} (GeV/c^{2});US_{Data}-US_{mixing}",trigName[kTrigType]),kFALSE);
  gPad->SetGridy();
  t1->Draw();
  signif->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_EventMixing/%s.Data-EventMixing_jpsi_pt1_%1.1f_pt2_%1.1f_rebin4.png",run_config,pt1_cut,pt2_cut));
}

//================================================
void acceptance(const Int_t save = 0)
{
  Double_t pt_cuts_1[6] = {1.0,1.2,1.5,2.0,2.5,3.0};
  Double_t pt_cuts_2[2] = {1.0,1.2};

  TH1F *hUS[6][2];
  TH1F *hLS[6][2];

  TList *list = new TList;
  TString legName[2] = {"Unlike sign","Like sign"};
  for(Int_t i=0; i<6; i++)
    {
      for(Int_t k=0; k<2; k++)
	{
	  hUS[i][k] = (TH1F*)f->Get(Form("US_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]));
	  hLS[i][k] = (TH1F*)f->Get(Form("LS_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]));
	  list->Clear();
	  list->Add(hUS[i][k]);
	  list->Add(hLS[i][k]);
	  c = sysCompare(list,Form("EventMixing_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[k]),Form("Au+Au %s: invariant mass of dimuon pairs",trigName[kTrigType]),"Acceptance effects (US/LS)",kFALSE,0,15,kFALSE,0.1,10,kFALSE,0,0.5,kTRUE,kTRUE,legName,kTRUE,"Mixed events",0.25,0.4,0.35,0.55,kFALSE);
	  c->cd(1);
	  TPaveText *t1 = GetPaveText(0.65,0.75,0.75,0.85);
	  t1->AddText(Form("p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c",pt_cuts_1[i],pt_cuts_2[k]));
	  t1->Draw();
	  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_EventMixing/EventMixing.%s.InvMass_pt1_%1.1f_pt2_%1.1f.png",run_config,pt_cuts_1[i],pt_cuts_2[k]));
	}
    }
}
