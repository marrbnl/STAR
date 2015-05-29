TFile *f;
const char *run_config = "";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void ana_Background()
{
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.18);                
  gStyle->SetStatH(0.15); 
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

  run_cfg_name = "SameEvent.";

  //makeHisto(fileName);
  //ana(fileName);
  compareToMixEvent(fileName);
}

//================================================
void compareToMixEvent(TString fileName, const int save = 0)
{
  gStyle->SetOptFit(111);
  TFile *fin[2];
  fin[0] = TFile::Open(Form("Rootfiles/%s",fileName.Data()),"read");
  fin[1] = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.MixEvent.root"),"read");

  Double_t pt_cuts_1[2] = {1.0,1.5};
  Double_t pt_cuts_2[2] = {1.0,1.0};

  TList *list = new TList;
  TString legName[2] = {"SameEvent","MixEvent"};

  TH2F *hLS[2][2], *hLSGeom[2][2];
  TH1F *hInvMassLS[2][2], *hInvMassLSGeom[2][2];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hLS[i][j] = (TH2F*)fin[j]->Get(Form("LS_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
	  hLS[i][j]->Sumw2();
	  hLS[i][j]->SetName(Form("%s_%s",hLS[i][j]->GetName(),legName[j].Data()));
	  c = draw2D(hLS[i][j],Form("%s: like-sign pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",legName[j].Data(), pt_cuts_1[i],pt_cuts_2[i]));
	  if(save) 
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%s.InvMass_vs_pt_LS_pt1_%1.0f_pt2_%1.0f.pdf",run_type,legName[j].Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%s.InvMass_vs_pt_LS_pt1_%1.0f_pt2_%1.0f.png",run_type,legName[j].Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	    }
	  hInvMassLS[i][j] = (TH1F*)hLS[i][j]->ProjectionX(Form("%s_prox",hLS[i][j]->GetName()));
	  hInvMassLS[i][j]->Rebin(10);

	  hLSGeom[i][j] = (TH2F*)fin[j]->Get(Form("LS_Geom_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
	  hLSGeom[i][j]->Sumw2();
	  hLSGeom[i][j]->SetName(Form("%s_%s",hLSGeom[i][j]->GetName(),legName[j].Data()));
	  c = draw2D(hLSGeom[i][j],Form("%s: like-sign pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",legName[j].Data(), pt_cuts_1[i],pt_cuts_2[i]));
	  if(save) 
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%s.InvMass_vs_pt_LS_Geom_pt1_%1.0f_pt2_%1.0f.pdf",run_type,legName[j].Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%s.InvMass_vs_pt_LS_Geom_pt1_%1.0f_pt2_%1.0f.png",run_type,legName[j].Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	    }
	  hInvMassLSGeom[i][j] = (TH1F*)hLSGeom[i][j]->ProjectionX(Form("%s_prox",hLSGeom[i][j]->GetName()));
	  scaleHisto(hInvMassLSGeom[i][j], 1, 1,kTRUE,kFALSE, kTRUE);
	}
      TH1F *hRatio = (TH1F*)hInvMassLS[i][0]->Clone(Form("%s_ratio",hInvMassLS[i][0]->GetName()));
      hRatio->GetXaxis()->SetRangeUser(0,8);
      hRatio->GetYaxis()->SetRangeUser(6e-4,2e-3);
      hRatio->SetMarkerStyle(25);
      hRatio->Divide(hInvMassLS[i][1]);
      TF1 *func = new TF1(Form("%s_func",hRatio->GetName()),"pol0",1,2);
      hRatio->Fit(func,"IR0");
      c = draw1D(hRatio,Form("Like-sign: SameEvent/MixEvent (N_{++}+N_{--})"));
      func->SetLineColor(2);
      func->Draw("same");
      TPaveText *t1 = GetPaveText(0.2,0.3,0.2,0.35,0.04,62);
      t1->AddText(Form("p_{T,1}>%1.1f GeV/c",pt_cuts_1[i]));
      t1->AddText(Form("p_{T,2}>%1.1f GeV/c",pt_cuts_2[i]));
      t1->Draw();
      hInvMassLS[i][1]->Scale(func->GetParameter(0));
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_SameToMix_pt1_%1.0f_pt2_%1.0f.pdf",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_SameToMix_pt1_%1.0f_pt2_%1.0f.png",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
      list->Clear();
      list->Add(hInvMassLS[i][0]);
      list->Add(hInvMassLS[i][1]);
      c = sysCompare(list,Form("LS_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[i]),Form("Like-sign di-muon pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt_cuts_1[i],pt_cuts_2[i]),"MixEvent/SameEvent",kTRUE,0,5,kFALSE,0.1,10,kTRUE,0.7,1.1,kFALSE,kTRUE,legName,kTRUE,"N_{++}+N_{--}",0.25,0.4,0.25,0.45,kFALSE);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_SameVsMix_pt1_%1.0f_pt2_%1.0f.pdf",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_SameVsMix_pt1_%1.0f_pt2_%1.0f.png",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}

      hRatio = (TH1F*)hInvMassLSGeom[i][0]->Clone(Form("%s_ratio",hInvMassLSGeom[i][0]->GetName()));
      hRatio->GetXaxis()->SetRangeUser(0,8);
      hRatio->GetYaxis()->SetRangeUser(6e-4,2e-3);
      hRatio->SetMarkerStyle(25);
      hRatio->Divide(hInvMassLSGeom[i][1]);
      TF1 *func = new TF1(Form("%s_func",hRatio->GetName()),"pol0",1,2);
      hRatio->Fit(func,"IR0");
      c = draw1D(hRatio,Form("Like-sign: SameEvent/MixEvent (2#sqrt{N_{++}*N_{--}})"));
      func->SetLineColor(2);
      func->Draw("same");
      TPaveText *t1 = GetPaveText(0.2,0.3,0.2,0.35,0.04,62);
      t1->AddText(Form("p_{T,1}>%1.1f GeV/c",pt_cuts_1[i]));
      t1->AddText(Form("p_{T,2}>%1.1f GeV/c",pt_cuts_2[i]));
      t1->Draw();
      hInvMassLSGeom[i][1]->Scale(func->GetParameter(0));
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_Geom_SameToMix_pt1_%1.0f_pt2_%1.0f.pdf",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_Geom_SameToMix_pt1_%1.0f_pt2_%1.0f.png",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
      list->Clear();
      list->Add(hInvMassLSGeom[i][0]);
      list->Add(hInvMassLSGeom[i][1]);
      c = sysCompare(list,Form("LS_Geom_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[i]),Form("Like-sign di-muon pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt_cuts_1[i],pt_cuts_2[i]),"MixEvent/SameEvent",kTRUE,0,5,kFALSE,0.1,10,kTRUE,0.7,1.1,kFALSE,kTRUE,legName,kTRUE,"2#sqrt{N_{++}*N_{--}}",0.25,0.4,0.25,0.45,kFALSE);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_Geom_SameVsMix_pt1_%1.0f_pt2_%1.0f.pdf",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/LS_Geom_SameVsMix_pt1_%1.0f_pt2_%1.0f.png",run_type,pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}

    }
  
}

//================================================
void ana(TString fileName, const Int_t save = 0)
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
	  TH1F *h1 = (TH1F*)hLSpair[j][i]->ProjectionX(Form("%s_prox",hLSpair[j][i]->GetName()));
	  h1->Rebin(10);
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

  // Geometry mean
  TString legName2[2] = {"N_{++}+N_{--}", "2#sqrt{N_{++}*N_{--}}"};
  TH2F *hLSGeom[2], *hLS[2];
  for(int i=0; i<2; i++)
    {
      hLSGeom[i] = (TH2F*)fin->Get(Form("LS_Geom_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      hLS[i] = (TH2F*)fin->Get(Form("LS_pt_vs_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[i],trigName[kTrigType]));
      
      list->Clear();
      TH1F *htmp = (TH1F*)hLS[i]->ProjectionX(Form("%s_prox2",hLS[i]->GetName()));
      TH1F *h1 = (TH1F*)htmp->Rebin(nSpecMBins,Form("htmp_rebin",htmp->GetName()),specM);
      scaleHisto(h1, 1, 1,kTRUE,kFALSE, kTRUE);
      list->Add(h1);
      TH1F *h1 = (TH1F*)hLSGeom[i]->ProjectionX(Form("%s_prox",hLSGeom[i]->GetName()));
      scaleHisto(h1, 1, 1,kTRUE,kFALSE, kTRUE);
      list->Add(h1);

      c = sysCompare(list,Form("GeometricalMean_pt1_%1.1f_pt2_%1.1f",pt_cuts_1[i],pt_cuts_2[i]),Form("Mixed event: like-sign di-muon pairs (p_{T,1}>%1.1f, p_{T,2}>%1.1f GeV/c)",pt_cuts_1[i],pt_cuts_2[i]),"Geometric/Arithmetic;M_{#mu#mu} (GeV/c^{2});",kTRUE,0,4,kFALSE,0.1,10,kTRUE,0.93,1.07,kFALSE,kTRUE,legName2,kTRUE,"Like-sign",0.25,0.4,0.25,0.45,kTRUE);
     if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Geom_pt1_%1.0f_pt2_%1.0f.pdf",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EventMixing/%sInvMass_Geom_pt1_%1.0f_pt2_%1.0f.png",run_type,run_cfg_name.Data(),pt_cuts_1[i]*10,pt_cuts_2[i]*10));
	}
    }
}



//================================================
void makeHisto(TString fileName)
{
  gStyle->SetOptStat(1);
  Double_t pt_cuts_1[6] = {1.0,1.2,1.5,2.0,4.0,4.5};
  Double_t pt_cuts_2[2] = {1.0,1.2};

  TH2F *hUS[6][2];
  TH2F *hLS[6][2];
  TH2F *hLSpos[6][2];
  TH2F *hLSneg[6][2];

  const char *hName[3] = {"mhJpsiInfo","mhBkgLSPos","mhBkgLSNeg"};
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
	      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt_cuts_1[i]+0.01,100);
	      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt_cuts_2[k]+0.01,100);
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
