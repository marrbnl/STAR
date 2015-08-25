const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;
TFile *f;

const char *run_config = "mix.";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void ana_EventMixing()
{
  gStyle->SetOptStat(0);
  TString fileName;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) fileName = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      fileName = Form("Run13.pp500.jpsi.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) fileName = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      fileName = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");

  run_cfg_name = run_config;

  eventBinning();
  //makeHisto(fileName);
  //acceptance(fileName);
  //subtraction();
}

//================================================
void eventBinning(const Int_t save = 0)
{
  TH1F *hVtxZ = (TH1F*)f->Get(Form("mhVertexZ_%s",trigName[kTrigType]));
  c = draw1D(hVtxZ,"",kFALSE,kFALSE);
  for(int i=0; i<19; i++)
    {
      double value = -100 + i*10 + 10;
      double height = hVtxZ->GetBinContent(hVtxZ->FindFixBin(value));
      TLine *line = GetLine(value,0,value,height,2,1,1);
      line->Draw();
    }
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sVertexZ_Binning.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sVertexZ_Binning.png",run_type,run_cfg_name.Data()));
    }

  TH1F *hgRefMultCorr = (TH1F*)f->Get(Form("mhgRefMultCorr_%s",trigName[kTrigType]));
  c = draw1D(hgRefMultCorr,"",kTRUE,kFALSE);
  //double bounds[9] = {472,401,283,193,126,77,44,23,11};
  double bounds[16] = {11,16,23,32,44,59,77,99,126,157,193,235,283,338,401,472};
  for(int i=0; i<16; i++)
    {
      double value = bounds[i];
      double height = hgRefMultCorr->GetBinContent(hgRefMultCorr->FindFixBin(value));
      TLine *line = GetLine(value,0,value,height,2,1,1);
      line->Draw();
    }
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sgRefMult_Binning.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sgRefMult_Binning.png",run_type,run_cfg_name.Data()));
    }

  // centrality
  TH1F *hgRefMult = (TH1F*)f->Get(Form("mhgRefMult_%s",trigName[kTrigType]));
  hgRefMult->SetLineColor(2);
  c = draw1D(hgRefMult,"",kTRUE,kFALSE);
  hgRefMultCorr->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.3,0.4,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hgRefMult,"Raw gRefMult","L");
  leg->AddEntry(hgRefMultCorr,"Corrected gRefMult","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sgRefMult.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sgRefMult.png",run_type,run_cfg_name.Data()));
    }


  TH1F *hCentrality = (TH1F*)f->Get(Form("mhCentrality_%s",trigName[kTrigType]));
  c = draw1D(hCentrality,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sCentrality.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sCentrality.png",run_type,run_cfg_name.Data()));
    }

  // event plane
  TH2F *hExcVsIncRawEP = (TH2F*)f->Get(Form("mhExcVsIncRawEP_%s",trigName[kTrigType]));
  c = draw2D(hExcVsIncRawEP);

  TH1F *hEventPlane = (TH1F*)f->Get(Form("mhEventPlane_%s",trigName[kTrigType]));
  c = draw1D(hEventPlane,"",kFALSE,kFALSE);
}

//================================================
void acceptance(TString fileName, const Int_t save = 0)
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s",fileName.Data()),"read");

  Double_t pt_cuts_1[2] = {1.0,1.5};
  Double_t pt_cuts_2[2] = {1.0,1.0};

  TList *list = new TList;

  // charge dependence
  TString legName[2] = {"Like-sign: ++","Like-sign: --"};
  const char *title[2] = {"positive","negative"};
  TH2F *hLSpair[2][2];
  for(int i=0; i<2; i++)
    {
      hLSpair[0][i] = (TH2F*)fin->Get(Form("LSpos_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      hLSpair[1][i] = (TH2F*)fin->Get(Form("LSneg_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      for(int j=0; j<2; j++)
	{
	  hLSpair[j][i]->GetXaxis()->SetRangeUser(0,8);
	  c = draw2D(hLSpair[j][i],Form("Mixed event: %s like-sign pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",title[j],pt_cuts_1[i],pt_cuts_2[i]));
	  if(save) 
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_vs_pt_LS_%s_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),title[j],pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_vs_pt_LS_%s_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),title[j],pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	    }
	}

      list->Clear();
      for(int j=0; j<2; j++)
	{
	  hLSpair[j][i]->GetXaxis()->SetRangeUser(0,8);
	  TH1F *h1 = (TH1F*)hLSpair[j][i]->ProjectionX(Form("%s_prox",hLSpair[j][i]->GetName()));
	  h1->Rebin(2);
	  list->Add(h1);
	}
      c = sysCompare(list,Form("Charge_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[i]),Form("Mixed event: di-muon pairs"),"Charge dependece;M_{#mu#mu} (GeV/c^{2});--/++",kFALSE,0,15,kFALSE,0.1,10,kTRUE,0.7,1.3,kTRUE,kTRUE,legName,kTRUE,"Mixed events",0.25,0.4,0.35,0.55,kFALSE);
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.7,0.8,0.7,0.85,0.04,62);
      t1->AddText(Form("p_{T,1}>%1.1f GeV/c",pt_cuts_1[i]));
      t1->AddText(Form("p_{T,2}>%1.1f GeV/c",pt_cuts_2[i]));
      t1->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_LS_ChargeDep_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_LS_ChargeDep_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
    }

  // acceptance
  TString legName1[2] = {"Unlike-sign","Like-sign: (--)+(++)"};
  TH2F *hLS[2], *hUL[2];
  for(int i=0; i<2; i++)
    {
      hUL[i] = (TH2F*)fin->Get(Form("US_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      hUL[i]->Sumw2();
      hUL[i]->GetXaxis()->SetRangeUser(0,8);
      c = draw2D(hUL[i],Form("Mixed event: unlike-sign pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt_cuts_1[i],pt_cuts_2[i]));
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_vs_pt_UL_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_vs_pt_UL_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
      hLS[i] = (TH2F*)fin->Get(Form("LS_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      hLS[i]->Sumw2();
      hLS[i]->GetXaxis()->SetRangeUser(0,8);
      c = draw2D(hLS[i],Form("Mixed event: like-sign pairs (--)+(++) (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt_cuts_1[i],pt_cuts_2[i]));
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_vs_pt_LS_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_vs_pt_LS_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}

      list->Clear();
      h1 = (TH1F*)hUL[i]->ProjectionX(Form("%s_prox",hUL[i]->GetName()));
      h1->Rebin(2);
      list->Add(h1);
      h1 = (TH1F*)hLS[i]->ProjectionX(Form("%s_prox",hLS[i]->GetName()));
      h1->Rebin(2);
      list->Add(h1);
      c = sysCompare(list,Form("Acceptance_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[i]),Form("Mixed event: di-muon pairs"),"Acceptance effect;M_{#mu#mu} (GeV/c^{2});LS/UL",kFALSE,0,15,kFALSE,0.1,10,kTRUE,0.7,1.3,kTRUE,kTRUE,legName1,kTRUE,"Mixed events",0.25,0.4,0.35,0.55,kFALSE);
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.7,0.8,0.7,0.85,0.04,62);
      t1->AddText(Form("p_{T,1}>%1.1f GeV/c",pt_cuts_1[i]));
      t1->AddText(Form("p_{T,2}>%1.1f GeV/c",pt_cuts_2[i]));
      t1->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Acceptance_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Acceptance_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}

      c = sysCompare(list,Form("Acceptance_pt1_%1.1f_pt2_%1.1f_zoom",pt_cuts_1[i],pt_cuts_2[i]),Form("Mixed event: di-muon pairs"),"Acceptance effect;M_{#mu#mu} (GeV/c^{2});LS/UL",kTRUE,0,4,kFALSE,0.1,10,kTRUE,0.7,1.3,kFALSE,kTRUE,legName1,kTRUE,"Mixed events",0.25,0.4,0.25,0.45,kFALSE);
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.3,0.4,0.7,0.85,0.04,62);
      t1->AddText(Form("p_{T,1}>%1.1f GeV/c",pt_cuts_1[i]));
      t1->AddText(Form("p_{T,2}>%1.1f GeV/c",pt_cuts_2[i]));
      t1->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Acceptance_pt1_%1.0f_pt2_%1.0f_zoomin.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Acceptance_pt1_%1.0f_pt2_%1.0f_zoomin.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
    }

  // Geometry mean
  TString legName2[2] = {"N_{++}+N_{--}", "2#sqrt{N_{++}*N_{--}}"};
  TH2F *hLSGeom[2];
  for(int i=0; i<2; i++)
    {
      hLSGeom[i] = (TH2F*)fin->Get(Form("LS_Geom_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      
      list->Clear();
      TH1F *htmp = (TH1F*)hLS[i]->ProjectionX(Form("%s_prox2",hLS[i]->GetName()));
      TH1F *h1 = (TH1F*)htmp->Rebin(nSpecMBins,Form("htmp_rebin",htmp->GetName()),specM);
      scaleHisto(h1, 1, 1,kTRUE,kFALSE, kTRUE);
      list->Add(h1);
      TH1F *h1 = (TH1F*)hLSGeom[i]->ProjectionX(Form("%s_prox",hLSGeom[i]->GetName()));
      scaleHisto(h1, 1, 1,kTRUE,kFALSE, kTRUE);
      list->Add(h1);

      c = sysCompare(list,Form("GeometricalMean_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[i]),Form("Mixed event: like-sign di-muon pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt_cuts_1[i],pt_cuts_2[i]),"Geometric/Arithmetic;M_{#mu#mu} (GeV/c^{2});",kTRUE,0,4,kFALSE,0.1,10,kTRUE,0.99,1.01,kFALSE,kTRUE,legName2,kTRUE,"Like-sign",0.25,0.4,0.25,0.45,kFALSE);
     if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Geom_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Geom_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
    }
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
void makeHisto(TString fileName)
{
  Double_t pt_cuts_1[6] = {1.0,1.2,1.5,2.0,4.0,4.5};
  Double_t pt_cuts_2[2] = {1.0,1.2};

  TH2F *hUS[6][2];
  TH2F *hLS[6][2];
  TH2F *hLSpos[6][2];
  TH2F *hLSneg[6][2];

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];

  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
    }

  for(Int_t i=0; i<6; i++)
    {
      for(Int_t k=0; k<2; k++)
	{
	  for(Int_t j=0; j<3; j++)
	    {
	      hnInvMass[j]->GetAxis(2)->SetRangeUser(pt_cuts_1[i]+0.01,100);
	      hnInvMass[j]->GetAxis(3)->SetRangeUser(pt_cuts_2[k]+0.01,100);
	      TH2F *h2 = (TH2F*)hnInvMass[j]->Projection(1,0);
	      h2->SetName(Form("%d_%d_%d",i,j,k));
	      if(j==0) hUS[i][k] = (TH2F*)h2->Clone(Form("US_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]));
	      if(j==1) hLSpos[i][k] = (TH2F*)h2->Clone(Form("LSpos_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]));
	      if(j==2) hLSneg[i][k] = (TH2F*)h2->Clone(Form("LSneg_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]));
	    }

	  hLS[i][k] = (TH2F*)hLSpos[i][k]->Clone(Form("LS_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]));
	  hLS[i][k]->Add(hLSneg[i][k]);
	}
    }

  TH2F *hUSRebin[6][2];
  TH2F *hLSGeom[6][2];
  TH2F *hLSposRebin[6][2];
  TH2F *hLSnegRebin[6][2];

  for(Int_t i=0; i<6; i++)
    {
      for(Int_t k=0; k<2; k++)
	{
	  hUSRebin[i][k] = new TH2F(Form("%s_Rebin",hUS[i][k]->GetName()), hUS[i][k]->GetTitle(), nSpecMBins, specM, nSpecPtBins, specPt);
	  CopyTH2(hUS[i][k], hUSRebin[i][k]);

	  hLSposRebin[i][k] = new TH2F(Form("%s_Rebin",hLSpos[i][k]->GetName()), hLSpos[i][k]->GetTitle(), nSpecMBins, specM, nSpecPtBins, specPt);
	  CopyTH2(hLSpos[i][k], hLSposRebin[i][k]);

	  hLSnegRebin[i][k] = new TH2F(Form("%s_Rebin",hLSneg[i][k]->GetName()), hLSneg[i][k]->GetTitle(), nSpecMBins, specM, nSpecPtBins, specPt);
	  CopyTH2(hLSneg[i][k], hLSnegRebin[i][k]);

	  hLSGeom[i][k] = new TH2F(Form("LS_Geom_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]), hLS[i][k]->GetTitle(), nSpecMBins, specM, nSpecPtBins, specPt);
	  GeometricMean(hLSposRebin[i][k], hLSnegRebin[i][k], hLSGeom[i][k]);
	}
    }
	  

  TFile *fout = TFile::Open(Form("Rootfiles/%s",fileName.Data()),"recreate");
  for(Int_t i=0; i<6; i++)
    {
      for(Int_t k=0; k<2; k++)
	{
	  hUS[i][k]->Write();
	  hLS[i][k]->Write();
	  hLSpos[i][k]->Write();
	  hLSneg[i][k]->Write();

	  hUSRebin[i][k]->Write();
	  hLSGeom[i][k]->Write();
	  hLSposRebin[i][k]->Write();
	  hLSnegRebin[i][k]->Write();
	}
    }
}

