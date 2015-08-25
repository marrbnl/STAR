TFile *f;
const char *run_config = "VPDMB.";
const Bool_t iPico = 0;
const int year = 2013;
TString run_cfg_name;
const int nMultBins = 4;
const int low_tofMult[nMultBins] = {0,0,4,7};
const int high_tofMult[nMultBins] = {30,3,6,30};

//================================================
void ana_Mult()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
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
  printf("all events: %4.4e\n",hStat->GetBinContent(1));
  printf("MB  events: %4.4e\n",hStat->GetBinContent(2));

  run_cfg_name = run_config;

  //mult();
  extrapolation();
}

//================================================
void extrapolation(const int save = 1)
{
  THnSparseF *hn = (THnSparseF*)f->Get("mhEventMult_qa_mb");
  TH2F *hMultVsZdc[2];
  TString title[2] = {"ranking>0","all"};
  for(int i=0; i<2; i++)
    {
      if(i==0) hn->GetAxis(4)->SetRange(3,3);
      if(i==1) hn->GetAxis(4)->SetRange(1,3);
      
      hMultVsZdc[i] = (TH2F*)hn->Projection(2,0);
      hMultVsZdc[i]->SetName(Form("TofMult_vs_zdc_%d",i));
      hMultVsZdc[i]->SetTitle(Form("MB: TofMult vs zdc rate (%s); zdc rate (kHz)",title[i].Data()));
      hMultVsZdc[i]->GetXaxis()->SetRangeUser(0,500);
      hMultVsZdc[i]->GetYaxis()->SetRangeUser(0,30);
      c = draw2D(hMultVsZdc[i]);
      TProfile *hpro = (TProfile*)hMultVsZdc[i]->ProfileX(Form("%s_profile",hMultVsZdc[i]->GetName()));
      hpro->SetMarkerStyle(20);
      hpro->DrawCopy("sames P");

      TF1 *func = new TF1(Form("Fit_TofMult_vs_zdc_%d",i),"pol1",10,150);
      hpro->Fit(func,"IR0");
      func->SetLineColor(2);
      func->Draw("sames");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%sTofMult_vs_ZdcRate_%d.pdf",run_type,run_cfg_name.Data(),i));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%sTofMult_vs_ZdcRate_%d.png",run_type,run_cfg_name.Data(),i));
	}

      TH2F *hTmp = (TH2F*)hMultVsZdc[i]->Clone(Form("%s_tmp",hMultVsZdc[i]->GetName()));
      hTmp->GetXaxis()->SetRangeUser(0,150);
      hTmp->GetYaxis()->SetRangeUser(0,10);
      c = draw2D(hTmp,Form("MB: TofMult vs zdc rate (%s); zdc rate (kHz)",title[i].Data()));
      hpro->DrawCopy("sames P");
      func->Draw("sames");

      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%sTofMult_vs_ZdcRate_%d_zoomin.pdf",run_type,run_cfg_name.Data(),i));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%sTofMult_vs_ZdcRate_%d_zoomin.png",run_type,run_cfg_name.Data(),i));
	}
    }
}

//================================================
void mult(const int save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get("mhEventMult_qa_mb");
  
  TH2F *hMultVsZdc[2][3];
  TH1F *hMultRanking[2][3];
  const char *name[2] = {"RefMult","TofMult"};
  TString title[3] = {"ranking>0","ranking<0","all"};
  const char *hname[3] = {"PosRank","NegRank","All"};
  TList *list = new TList;
  for(int i=0; i<2; i++)
    {
      list->Clear();
      for(int j=0; j<3; j++)
	{
	  if(j==0) hn->GetAxis(4)->SetRange(3,3);
	  if(j==1) hn->GetAxis(4)->SetRange(1,1);
	  if(j==2) hn->GetAxis(4)->SetRange(1,3);
	  hMultVsZdc[i][j] = (TH2F*)hn->Projection(i+1,0);
	  hMultVsZdc[i][j]->SetName(Form("%s_vs_zdc_%d",name[i],j));
	  hMultVsZdc[i][j]->SetTitle(Form("MB: %s vs zdc rate (%s); zdc rate (kHz)",name[i],title[j].Data()));
	  hMultVsZdc[i][j]->GetXaxis()->SetRangeUser(0,500);
	  if(i==0) hMultVsZdc[i][j]->GetYaxis()->SetRangeUser(0,50);
	  if(i==1) hMultVsZdc[i][j]->GetYaxis()->SetRangeUser(0,30);
	  c = draw2D(hMultVsZdc[i][j]);
	  TProfile *hpro = (TProfile*)hMultVsZdc[i][j]->ProfileX(Form("%s_profile",hMultVsZdc[i][j]->GetName()));
	  hpro->SetMarkerStyle(20);
	  hpro->DrawCopy("sames P");
	  if(save) 
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%s%s_vs_zdc_with_%s.pdf",run_type,run_cfg_name.Data(),name[i],hname[j]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%s%s_vs_zdc_with_%s.png",run_type,run_cfg_name.Data(),name[i],hname[j]));
	    }

	  hMultRanking[i][j] = (TH1F*)hMultVsZdc[i][j]->ProjectionY(Form("%s_proY",hMultVsZdc[i][j]->GetName()));
	  cout << title[j] << " with " << hMultRanking[i][j]->GetEntries() << endl;
	  hMultRanking[i][j]->Scale(1./hMultRanking[i][j]->Integral());
	  hMultRanking[i][j]->SetMaximum(10*hMultRanking[i][j]->GetMaximum());
	  list->Add(hMultRanking[i][j]);
	}
      c = drawHistos(list,name[i],Form("MB: distribution of %s;%s",name[i],name[i]),kFALSE,0,15,kFALSE,1e-6,0.05,kTRUE,kTRUE,title,kTRUE,"Vertex ranking",0.6,0.8,0.6,0.85,kFALSE);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%s%s_vs_ranking.pdf",run_type,run_cfg_name.Data(),name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%s%s_vs_ranking.png",run_type,run_cfg_name.Data(),name[i]));
	}
    }

  hn->GetAxis(4)->SetRange(3,3);
  TH2F *hMultVsVz[2];
  for(int i=0; i<2; i++)
    {
      hMultVsVz[i] = (TH2F*)hn->Projection(i+1,3);
      hMultVsVz[i]->SetName(Form("%s_vs_vz",name[i]));
      hMultVsVz[i]->SetTitle(Form("MB: %s vs TPC vz (ranking>0)",name[i]));
      if(i==0) hMultVsVz[i]->GetYaxis()->SetRangeUser(0,50);
      if(i==1) hMultVsVz[i]->GetYaxis()->SetRangeUser(0,30);
      c = draw2D(hMultVsVz[i]);
      TProfile *hpro = (TProfile*)hMultVsVz[i]->ProfileX(Form("%s_profile",hMultVsVz[i]->GetName()));
      hpro->SetMarkerStyle(20);
      hpro->DrawCopy("sames P");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%s%s_vs_TpcVz_PosRank.pdf",run_type,run_cfg_name.Data(),name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Mult/%s%s_vs_TpcVz_PosRank.png",run_type,run_cfg_name.Data(),name[i]));
	}
    }

  // Zdc rate in different TofMult Bin  
  list->Clear();
  TString legName[nMultBins];
  hn->GetAxis(4)->SetRange(0,-1);
  TH1F *hZdcRate[nMultBins];
  for(int i=0; i<nMultBins; i++)
    {
      hn->GetAxis(2)->SetRange(low_tofMult[i]+1,high_tofMult[i]+1);
      hZdcRate[i] = (TH1F*)hn->Projection(0);
      hZdcRate[i]->SetName(Form("ZdcRate_TofMult_%d_%d",low_tofMult[i],high_tofMult[i]));
      hZdcRate[i]->Rebin(10);
      hZdcRate[i]->Scale(1./hZdcRate[i]->Integral(15,35));
      list->Add(hZdcRate[i]);
      legName[i] = Form("%d <= TofMult <= %d",low_tofMult[i],high_tofMult[i]);
    }
  c = drawHistos(list,"ZdcRate","VPDMB: distribution of ZDC rate;ZDC rate",kFALSE,0,15,kFALSE,1e-6,0.05,kTRUE,kTRUE,legName,kTRUE,"",0.6,0.8,0.6,0.85,kFALSE);
  
  
}
