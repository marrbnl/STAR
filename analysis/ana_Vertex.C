TFile *f;
const char *run_config = "default.";
const Bool_t iPico = 0;
const int year = 2013;
TString run_cfg_name;
const Double_t low_mass = 3.0;
const Double_t high_mass = 3.2;

//================================================
void ana_Vertex()
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
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("di-muon     events: %4.4e\n",hStat->GetBinContent(9));

  run_cfg_name = run_config;

  vertexSelection();
  //vzDiff();
  //ranking();
}


//================================================
void vertexSelection(const bool save = 0)
{
  TList *list = new TList;
  TH1F *hVzDiff[3];
  TString legName[3] = {"Default","Closest","MtdMth"};
  for(int i=0; i<3; i++)
    {
      hVzDiff[i] = (TH1F*)f->Get(Form("mhVzDiff%s_%s",legName[i].Data(),trigName[kTrigType]));
      hVzDiff[i]->Sumw2();
      hVzDiff[i]->Scale(1./hVzDiff[i]->GetEntries());
      list->Add(hVzDiff[i]);
    }
  c = drawHistos(list,"VzDiff","Vz difference between TPC and VPD",kFALSE,0,15,kTRUE,1e-6,0.5,kTRUE,kTRUE,legName,kTRUE,"Vertex selection",0.2,0.4,0.6,0.85,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/Compare_VzDiff.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/Compare_VzDiff.png",run_type));
    }
  c = drawHistos(list,"VzDiff_linear","Vz difference between TPC and VPD",kFALSE,0,15,kTRUE,1e-6,0.05,kFALSE,kTRUE,legName,kTRUE,"Vertex selection",0.2,0.4,0.6,0.85,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/Compare_VzDiff_lin.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/Compare_VzDiff_lin.png",run_type));
    }

  // VtxIndexVsNMtdTrk
  TH2F *hVtxIndexVsNMtdTrk = (TH2F*)f->Get(Form("mhVtxIndexVsNMtdTrk_%s",trigName[kTrigType]));
  hVtxIndexVsNMtdTrk->GetXaxis()->SetRangeUser(0,25);
  c = draw2D(hVtxIndexVsNMtdTrk,"");
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_vs_NMtdTrk.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_vs_NMtdTrk.png",run_type));
    }
  TH1F *hMtdIndex = (TH1F*)hVtxIndexVsNMtdTrk->ProjectionX("MtdIndex_1",3,-1);
  c = draw1D(hMtdIndex,"Index of TPC vertex with at least 2 MTD matches",kTRUE,kFALSE);
  TPaveText *t1 = GetPaveText(0.5,0.6,0.7,0.8);
  t1->AddText(Form("Index = 0: %3.2f%%",hMtdIndex->GetBinContent(1)/hMtdIndex->GetEntries()*100));
  t1->Draw();
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_With2MtdMth.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_With2MtdMth.png",run_type));
    }

  // VtxClosestVsMtdMth
  TH2F *hVtxClosestVsMtdMth = (TH2F*)f->Get(Form("mhVtxClosestVsMtdMth_%s",trigName[kTrigType]));
  hVtxClosestVsMtdMth->GetXaxis()->SetRangeUser(-1,25);
  hVtxClosestVsMtdMth->GetYaxis()->SetRangeUser(-1,25);
  c = draw2D(hVtxClosestVsMtdMth,"");
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/Index_ClosestVtx_vs_MtdMth.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/Index_ClosestVtx_vs_MtdMth.png",run_type));
    }

  TH1F *hClosestIndex = (TH1F*)hVtxClosestVsMtdMth->ProjectionY("ClosestIndex");
  c = draw1D(hClosestIndex,"Index of TPC vertex closest to VPD",kTRUE,kFALSE);
  TPaveText *t1 = GetPaveText(0.5,0.6,0.7,0.8);
  t1->AddText(Form("Index = 0: %3.2f%%",hClosestIndex->GetBinContent(2)/hClosestIndex->GetEntries()*100));
  t1->Draw();
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_CloseToVPD.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_CloseToVPD.png",run_type));
    }

  TH1F *hMtdIndex = (TH1F*)hVtxClosestVsMtdMth->ProjectionX("MtdIndex");
  c = draw1D(hMtdIndex,"Index of TPC vertex with MTD match",kTRUE,kFALSE);
  TPaveText *t1 = GetPaveText(0.5,0.6,0.7,0.8);
  t1->AddText(Form("Index = 0: %3.2f%%",hMtdIndex->GetBinContent(2)/hMtdIndex->Integral(2,50)*100));
  t1->Draw();
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_MtdMth.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_vertex/VtxIndex_MtdMth.png",run_type));
    }

  double nClosest = 0;
  for(int bin=2; bin<=50; bin ++)
    {
      nClosest += hVtxClosestVsMtdMth->GetBinContent(bin,bin);
    }
  cout << "Total = " << hMtdIndex->Integral(2,50) << endl;
  cout << "TPC = " << hMtdIndex->GetBinContent(2) << " : " << hMtdIndex->GetBinContent(2)/hMtdIndex->Integral(2,50)*100 << "%" << endl;
  cout << "VPD = " << nClosest << " : " << nClosest/hMtdIndex->Integral(2,50)*100 << "%" << endl;
}

//================================================
void vzDiff(const bool save = 0)
{
  const double pt_cut = 1.5;
  const int diffcut = 6;
  const double center = 0;

  THnSparseF *hnInvMass = (THnSparseF*)f->Get(Form("mhJpsiCutStudy_%s",trigName[kTrigType]));
  hnInvMass->GetAxis(1)->SetRangeUser(pt_cut,100);
  
  TH2F *hInvMassVsVz[2];
  TH1F *hProjInvMass[2][3];
  int low_bin, hi_bin;
  const char *type[2] = {"UL","LS"};
  for(int i=0; i<2; i++)
    {
      hnInvMass->GetAxis(2)->SetRange(i+1,i+1);
      hInvMassVsVz[i] = (TH2F*)hnInvMass->Projection(0,4);
      hInvMassVsVz[i]->Sumw2();
      hInvMassVsVz[i]->SetName(Form("hInvMassVsVz_%s",type[i]));
      for(int j=0; j<3; j++)
	{
	  if(j==0) 
	    {
	      low_bin = hnInvMass->GetAxis(4)->FindBin(center-diffcut+1e-5);
	      hi_bin = hnInvMass->GetAxis(4)->FindBin(center+diffcut-1e-5);
	    }
	  else if(j==1)
	    {
	      low_bin = 0;
	      hi_bin = hnInvMass->GetAxis(4)->FindBin(center-diffcut-1e-5);
	    }
	  else 
	    {
	      low_bin = hnInvMass->GetAxis(4)->FindBin(center+diffcut+1e-5);
	      hi_bin = hnInvMass->GetAxis(4)->GetNbins()+1;
	    }
	  hnInvMass->GetAxis(4)->SetRange(low_bin, hi_bin);
	  hProjInvMass[i][j] = (TH1F*)hnInvMass->Projection(0);
	  hProjInvMass[i][j]->Sumw2();
	  hProjInvMass[i][j]->Rebin(4);
	  hProjInvMass[i][j]->SetName(Form("hProjInvMass_%d_%d",i,j));
	}
      hnInvMass->GetAxis(4)->SetRange(0,-1);
    }
  hInvMassVsVz[0]->Add(hInvMassVsVz[1],-1);
  draw2D(hInvMassVsVz[0],"Invariant mass vs vz_{TPC-VPD};vz_{TPC-VPD} (cm)");

  for(int i=0; i<2; i++)
    {
      hProjInvMass[i][1]->Add(hProjInvMass[i][2]);
    }

  TH1F *hDiff[2];
  const char *title[2] = {"small","large"};
  for(int j=0; j<2; j++)
    {
      hProjInvMass[0][j]->SetMarkerStyle(20);
      hProjInvMass[0][j]->SetMarkerColor(4);
      hProjInvMass[0][j]->SetLineColor(4);
      if(j==0) c = draw1D(hProjInvMass[0][j],Form("Invariant mass distribution with |vz_{TPC}-vz_{VPD}| < %d cm;M_{#mu#mu} (GeV/c^{2});counts",diffcut),kFALSE,kTRUE);
      if(j==1) c = draw1D(hProjInvMass[0][j],Form("Invariant mass distribution with |vz_{TPC}-vz_{VPD}| > 6 cm;M_{#mu#mu} (GeV/c^{2});counts",diffcut),kFALSE,kTRUE);
      hProjInvMass[0][j]->Draw("HIST sames");
      hProjInvMass[1][j]->Draw("HIST sames");
      TLegend *leg = new TLegend(0.5,0.63,0.7,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("p_{T,1} > %1.1f GeV/c",pt_cut));
      leg->AddEntry(hProjInvMass[0][j],"Unlike sign","PLE");
      leg->AddEntry(hProjInvMass[1][j],"Like sign (++)+(--)","L");
      leg->Draw();
      Int_t low_bin = hProjInvMass[0][j]->GetXaxis()->FindFixBin(low_mass+0.001);
      Int_t high_bin = hProjInvMass[0][j]->GetXaxis()->FindFixBin(high_mass-0.001);
      Double_t nBackground = hProjInvMass[1][j]->Integral(low_bin,high_bin);
      Double_t nSignal = hProjInvMass[0][j]->Integral(low_bin,high_bin) - nBackground;
      TPaveText *signif = GetPaveText(0.5,0.7,0.4,0.55);
      signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
      signif->AddText(Form("S/B = %1.0f/%1.0f = %1.2e:1",nSignal,nBackground,nSignal/nBackground));
      signif->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/%sInvMass_VzDiff%s%dcm_pt1_%1.0f.pdf",run_type,run_cfg_name.Data(),title[j],diffcut,pt_cut*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/%sInvMass_VzDiff%s%dcm_pt1_%1.0f.png",run_type,run_cfg_name.Data(),title[j],diffcut,pt_cut*10));
	}
    }
}

//================================================
void ranking(const bool save = 0)
{
  const double pt_cut = 1.5;
  THnSparseF *hnInvMass = (THnSparseF*)f->Get(Form("mhJpsiCutStudy_%s",trigName[kTrigType]));
  hnInvMass->GetAxis(1)->SetRangeUser(pt_cut,100);

  TH1F *hProjInvMass[2][2];
  int low_bin, hi_bin;
  const char *type[2] = {"UL","LS"};
  for(int i=0; i<2; i++)
    {
      hnInvMass->GetAxis(2)->SetRange(i+1,i+1);
      for(int j=0; j<2; j++)
	{
	  hnInvMass->GetAxis(3)->SetRange(3-j*2,3-j*2);
	  hProjInvMass[i][j] = (TH1F*)hnInvMass->Projection(0);
	  hProjInvMass[i][j]->Sumw2();
	  hProjInvMass[i][j]->Rebin(4);
	  hProjInvMass[i][j]->SetName(Form("hProjInvMass_%d_%d",i,j));
	}
    }
  TH1F *hDiff[2];
  const char *title[2] = {"pos","neg"};
  for(int j=0; j<2; j++)
    {
      hProjInvMass[0][j]->SetMarkerStyle(20);
      hProjInvMass[0][j]->SetMarkerColor(4);
      hProjInvMass[0][j]->SetLineColor(4);
      if(j==0) c = draw1D(hProjInvMass[0][j],Form("Invariant mass distribution with positive ranking;M_{#mu#mu} (GeV/c^{2});counts"),kFALSE,kTRUE);
      if(j==1) c = draw1D(hProjInvMass[0][j],Form("Invariant mass distribution with negative ranking;M_{#mu#mu} (GeV/c^{2});counts"),kFALSE,kTRUE);
      hProjInvMass[0][j]->Draw("HIST sames");
      hProjInvMass[1][j]->Draw("HIST sames");
      cout << hProjInvMass[1][j]->GetEntries() << endl;
      TLegend *leg = new TLegend(0.5,0.63,0.7,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("p_{T,1} > %1.1f GeV/c",pt_cut));
      leg->AddEntry(hProjInvMass[0][j],"Unlike sign","PLE");
      leg->AddEntry(hProjInvMass[1][j],"Like sign (++)+(--)","L");
      leg->Draw();
      Int_t low_bin = hProjInvMass[0][j]->GetXaxis()->FindFixBin(low_mass+0.001);
      Int_t high_bin = hProjInvMass[0][j]->GetXaxis()->FindFixBin(high_mass-0.001);
      Double_t nBackground = hProjInvMass[1][j]->Integral(low_bin,high_bin);
      Double_t nSignal = hProjInvMass[0][j]->Integral(low_bin,high_bin) - nBackground;
      TPaveText *signif = GetPaveText(0.5,0.7,0.4,0.55);
      signif->AddText(Form("[%1.1f,%1.1f] GeV/c^{2}",low_mass,high_mass));
      signif->AddText(Form("S/B = %1.0f/%1.0f = %1.2e:1",nSignal,nBackground,nSignal/nBackground));
      signif->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/%sInvMass_%sRanking_pt1_%1.0f.pdf",run_type,run_cfg_name.Data(),title[j],pt_cut*10));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Vertex/%sInvMass_%sRanking_pt1_%1.0f.png",run_type,run_cfg_name.Data(),title[j],pt_cut*10));
	}
    }
}

