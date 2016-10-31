const char *run_config = "";
const int year = YEAR;
TString run_cfg_name;
const int nDet = 4;
const char *det_name[nDet] = {"","Tpc","Mtd","Final"};
const char *charge_name[3] = {"","_pos","_neg"};
const char *charge_title[3] = {" "," positive "," negative "};

TFile *f;

//================================================
void ana_Embed()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
  TString fileName;
  TString outName;
  TString outPDF;

  const char *particle = "Upsilon";

  if(year==2013)
    {
      fileName = Form("Run13.pp500.jpsi.Embed.root");
    }
  else if(year==2014)
    {
      fileName = Form("Run14.AuAu200.%s.Embed.root",particle);
    }
  outName = Form("%s.SingleTrkEff.%s.root",run_type,particle);
  outPDF = Form("%s_Embed_SingleTrkEff.%s.pdf",run_type,particle);

  run_cfg_name = run_config;

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("# of events: %4.4e\n",hStat->GetBinContent(3));

  //makeHistos(outName);
  //efficiency(outName);
  TrkEff3D(outName, outPDF);
}

//================================================
void efficiency(TString inName, const bool savePlot = 0, const bool saveHisto = 0)
{
  TFile *fin = 0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"update");
  else fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"read");
  TList *list = new TList;

  // TPC tracking efficiency
  TH3F *hMcTrkPtEtaPhi[nCentBins][4][3];
  TH1F *hMcTrkPt[nCentBins][4][3];
  TH2F *hMcTrkEtaPhi[nCentBins][4][3];
  TH1F *hMcTrkEta[nCentBins][4][3];
  TH1F *hMcTrkPhi[nCentBins][4][3];
  TH1F *hMcTrkPtInEtaPhi[nCentBins][4][2][2][12];

  TH1F *hMcTrkPtEff[nCentBins][4][3];
  TH2F *hMcTrkEtaPhiEff[nCentBins][4][3];
  TH1F *hMcTrkEtaEff[nCentBins][4][3];
  TH1F *hMcTrkPhiEff[nCentBins][4][3];
  TH1F *hMcTrkPtInEtaPhiEff[nCentBins][4][3][2][12];

  double eta_bounds[3] = {-0.5,0,0.5};
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<4; i++)
	    {
	      hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      int low_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(-0.5+1e-4);
	      int high_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(0.5-1e-4);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),low_ybin,high_ybin,0,-1);
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);

	      int low_xbin = hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->FindFixBin(1.5);
	      hMcTrkEta[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionY(Form("McTrkEta%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),low_xbin,-1,0,-1);
	      hMcTrkEtaEff[k][i][j] = (TH1F*)hMcTrkEta[k][i][j]->Clone(Form("%s_Eff",hMcTrkEta[k][i][j]->GetName()));
	      hMcTrkEtaEff[k][i][j]->Divide(hMcTrkEta[k][0][j]);

	      hMcTrkPhi[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionZ(Form("McTrkPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),low_xbin,-1,low_ybin,high_ybin);
	      hMcTrkPhiEff[k][i][j] = (TH1F*)hMcTrkPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPhi[k][i][j]->GetName()));
	      hMcTrkPhiEff[k][i][j]->Divide(hMcTrkPhi[k][0][j]);
	     
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRangeUser(1.5,20);
	      hMcTrkEtaPhi[k][i][j] = (TH2F*)hMcTrkPtEtaPhi[k][i][j]->Project3D("zy");
	      hMcTrkEtaPhi[k][i][j]->SetName(Form("McTrkEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      hMcTrkEtaPhiEff[k][i][j] = (TH2F*)hMcTrkEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkEtaPhi[k][i][j]->GetName()));
	      hMcTrkEtaPhiEff[k][i][j]->Divide(hMcTrkEtaPhi[k][0][j]);
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRange(0,-1);
	    }
	  for(int i=0; i<2; i++)
	    {
	      for(int ieta=0; ieta<2; ieta++)
	      	{
	      	  int low_etabin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(eta_bounds[ieta]+1e-4);
	      	  int high_etabin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(eta_bounds[ieta+1]-1e-4);
	      	  for(int iphi=0; iphi<12; iphi++)
	      	    {
	      	      int low_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12+pi/6*iphi + 1e-4);
	      	      int high_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12+pi/6*(iphi+1) - 1e-4);
	      	      hMcTrkPtInEtaPhi[k][i][j][ieta][iphi] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_eta%d_phi%d_cent%d",det_name[i],charge_name[j],ieta,iphi,k),low_etabin,high_etabin,low_phibin,high_phibin);
	      	      hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Rebin(2);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi] = (TH1F*)hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Clone(Form("%s_Eff",hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->GetName()));
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->Divide(hMcTrkPtInEtaPhi[k][0][j][ieta][iphi]);
	      	    }
	      	}
	    }
	}
    }

  TString legName[4] = {"MC input","Reconstructed in TPC","Match to MTD","Muon pid cuts"};
  const char *eff_type[4] = {"MC input","reconstructed in TPC","matched to MTD","fulfill muon pid cuts"};
  TString legName_charge[3] = {"All","Positive","Negative"};
  for(int k=0; k<nCentBins; k++)
    {
      // efficiency vs eta vs phi
      for(int i=1; i<4; i++)
	{
	  hMcTrkEtaPhiEff[k][i][0]->GetXaxis()->SetRangeUser(-1,1);
	  ScaleHistoTitle(hMcTrkEtaPhiEff[k][i][0],0.05,0.9,0.04,0.05,0.8,0.04,62);
	  char *title = Form("Efficiency of MC muon tracks %s (p_{T,mc}>1.5 GeV/c)",eff_type[i]);
	  if(nCentBins>1) title = Form("Efficiency of MC muon tracks %s (p_{T,mc}>1.5 GeV/c, %s%%)",eff_type[i],cent_Name[k]);
	  c = draw2D(hMcTrkEtaPhiEff[k][i][0],title,0.04,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaPhiEff_%s_cent%s.pdf",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaPhiEff_%s_cent%s.png",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	    }
	}

      // efficiency vs eta
      c = new TCanvas(Form("McTrkEffEta_Cent%d",k),Form("McTrkEffEta_Cent%d",k),800,600);
      double xmin = 0.15, xmax = 0.3, ymin = 0.7, ymax = 0.88;
      int counter = 0;
      for(int i=1; i<4; i++)
	{
	  leg = new TLegend(xmin+counter*0.2,ymin,xmax+counter*0.2,ymax);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.035);
	  if(i==1) leg->SetHeader("TPC tracking");
	  if(i==2) leg->SetHeader("MTD matched");
	  if(i==3) leg->SetHeader("Muon PID");
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkEtaEff[k][i][j]->SetMarkerStyle(19+i);
	      else     
		{
		  hMcTrkEtaEff[k][i][j]->SetMarkerStyle(23+i);
		  if(i==3) hMcTrkEtaEff[k][i][j]->SetMarkerStyle(32);
		}
	      hMcTrkEtaEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->SetLineColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->GetXaxis()->SetRangeUser(-0.6,0.6);
	      hMcTrkEtaEff[k][i][j]->GetYaxis()->SetRangeUser(0,1.4);
	      hMcTrkEtaEff[k][i][j]->SetTitle(";#eta;Efficiency");
	      if(j==0 && i==1) hMcTrkEtaEff[k][i][j]->Draw("P");
	      else hMcTrkEtaEff[k][i][j]->Draw("sames");
	      leg->AddEntry(hMcTrkEtaEff[k][i][j],legName_charge[j].Data(),"P");
	    }
	  leg->Draw();
	  counter++;
	}
      t1 = GetTitleText("Efficiency of MC muon track (p_{T,mc} > 1.5 GeV/c)");
      if(nCentBins>1)   t1 = GetTitleText(Form("Efficiency of MC muon track (p_{T,mc} > 1.5 GeV/c, %s%%)",cent_Name[k]));
      t1->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}
      
      // efficiency vs phi
      c = new TCanvas(Form("McTrkEffPhi_Cent%d",k),Form("McTrkEffPhi_Cent%d",k),800,600);
      leg = new TLegend(0.15,0.7,0.3,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
      int counter = 0;
      for(int i=1; i<4; i++)
	{
	  TH1F *htmp = (TH1F*)hMcTrkPhiEff[k][i][0]->Clone(Form("%s_clone",hMcTrkPhiEff[k][i][0]->GetName()));
	  htmp->SetLineColor(color[counter]);
	  htmp->GetYaxis()->SetRangeUser(0,1.4);
	  htmp->SetTitle(";#varphi_{mc};Efficiency");
	  if(i==1) htmp->DrawCopy("HIST");
	  else     htmp->DrawCopy("sames HIST");
	  leg->AddEntry(htmp,legName[counter+1].Data(),"L");
	  counter++;
	}
      leg->Draw();
      t1 = GetTitleText("Efficiency of MC muon track (p_{T,mc} > 1.5 GeV/c)");
      t1->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPhiEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPhiEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      for(int i=1; i<4; i++)
	{
	  c = new TCanvas(Form("McTrkEffPhi_%s_Cent%d",det_name[i],k),Form("McTrkEffPhi_%s_Cent%d",det_name[i],k),800,600);
	  leg = new TLegend(xmin,ymin,xmax,ymax);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.035);
	  if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkPhiEff[k][i][j]->SetMarkerStyle(19+i);
	      else     
		{
		  hMcTrkPhiEff[k][i][j]->SetMarkerStyle(23+i);
		  if(i==4) hMcTrkPhiEff[k][i][j]->SetMarkerStyle(32);
		}
	      hMcTrkPhiEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkPhiEff[k][i][j]->SetLineColor(color[j]);
	      hMcTrkPhiEff[k][i][j]->GetYaxis()->SetRangeUser(0,1.4);
	      hMcTrkPhiEff[k][i][j]->SetTitle("");
	      if(j==0 && i==1) hMcTrkPhiEff[k][i][j]->Draw("HIST");
	      else hMcTrkPhiEff[k][i][j]->Draw("sames HIST");
	      leg->AddEntry(hMcTrkPhiEff[k][i][j],legName_charge[j].Data(),"L");
	    }
	  leg->Draw();
	  t1 = GetTitleText(Form("Efficiency of MC muon track %s (p_{T,mc} > 1.5 GeV/c)",eff_type[i]));
	  t1->Draw();
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPhiEff_%s_cent%s.pdf",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPhiEff_%s_cent%s.png",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	    }
	}

      // efficiency vs pt
      for(int j=0; j<3; j++)
	{
	  list->Clear();
	  for(int i=0; i<4; i++)
	    {
	      list->Add(hMcTrkPt[k][i][j]);
	    }
	  char *centrality = Form(" (%s%%)",cent_Name[k]);
	  if(year==2013) centrality = "";
	  c = sysCompare(list,Form("McTrkEff_%s_Cent%d",charge_name[j],k),Form("p_{T} distribution of%smuon %s",charge_title[j],centrality),"Efficiency of single muon;p_{T}^{mc} (GeV/c);Efficiency",kTRUE,0,11,kFALSE,0.1,10,kTRUE,0,1,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.5,0.7,0.6,0.85,kTRUE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtEff%s_cent%s.pdf",run_type,run_cfg_name.Data(),charge_name[j],cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtEff%s_cent%s.png",run_type,run_cfg_name.Data(),charge_name[j],cent_Title[k]));
	    }
	}

      // efficiency vs pt in TPC sectors
      for(int ieta=0; ieta<2; ieta++)
	{
	  c = new TCanvas(Form("McTrkEffInPhi_eta%d_Cent%d",ieta,k),Form("McTrkEffInPhi_eta%d_Cent%d",ieta,k),1200,750);
	  c->Divide(4,3);
	  leg = new TLegend(0.4,0.25,0.8,0.5);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.06);
	  if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
	  for(int iphi=0; iphi<12; iphi++)
	    {
	      c->cd(iphi+1);
	      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	      for(int j=0; j<3; j++)
		{
		  hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->SetMarkerStyle(21);
		  hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->SetMarkerColor(color[j]);
		  hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->SetLineColor(color[j]);
		  hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->GetXaxis()->SetRangeUser(0,11);
		  hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->GetYaxis()->SetRangeUser(0,1.2);
		  hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->SetTitle(";p_{T,mc} [GeV/c];TPC Efficiency");
		  ScaleHistoTitle(hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi],0.06,1,0.05,0.06,1,0.05,62);
		  if(j==0) hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->Draw("P");
		  else     hMcTrkPtInEtaPhiEff[k][1][j][ieta][iphi]->Draw("sames");
		  int sector = -1;
		  if(ieta==0)
		    {
		      sector = 22 + iphi;
		      if(sector>24) sector -= 12;
		    }
		  else
		    {
		      sector = 2 - iphi;
		      if(sector<1) sector += 12;
		    }
		  TPaveText *t1 = GetTitleText(Form("TPC sector %d",sector),0.07);
		  t1->Draw();
		}
	      if(iphi==0)
		{
		  leg->AddEntry(hMcTrkPtInEtaPhiEff[k][1][0][ieta][iphi],"All","PL");
		  leg->AddEntry(hMcTrkPtInEtaPhiEff[k][1][1][ieta][iphi],"#mu^{+}","PL");
		  leg->AddEntry(hMcTrkPtInEtaPhiEff[k][1][2][ieta][iphi],"#mu^{-}","PL");
		}
	    }
	  c->cd(1);
	  leg->Draw();
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtEffInTpc_eta%d_cent%s.pdf",run_type,run_cfg_name.Data(),ieta,cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtEffInTpc_eta%d_cent%s.png",run_type,run_cfg_name.Data(),ieta,cent_Title[k]));
	    }
	}
    }

  list->Clear();
  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      TH1F *htmp = (TH1F*)hMcTrkPtEff[k][1][0]->Clone(Form("%s_clone",hMcTrkPtEff[k][1][0]->GetName()));
      list->Add(htmp);
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }
  c = drawHistos(list,Form("TrkEff"),Form("TPC tracking efficiency of muon tracks;p_{T,true} (GeV/c);Efficiency"),kTRUE,0,11,kTRUE,0,1.5,kFALSE,kTRUE,legName_cent,kTRUE,"",0.3,0.5,0.6,0.85,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtEff_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtEff_CentBins.png",run_type,run_cfg_name.Data()));
    }

  // track momentum resolution
  TH2F *hResVsTruePt[nCentBins];
  TH1F *hTrkPtRes[nCentBins];
  TH1F *hTrkPtShift[nCentBins];
  TF1 *funcRes[nCentBins];
  list->Clear();
  for(int k=0; k<nCentBins; k++)
    {
      hResVsTruePt[k] = (TH2F*) fin->Get(Form("PrimTrkRes_vs_TruePt_%s",cent_Title[k]));
      hResVsTruePt[k]->GetXaxis()->SetRangeUser(0,10);
      char *title = Form("Transverse momentum resolution of primary muon tracks (%s%%)",cent_Name[k]);
      c = draw2D(hResVsTruePt[k],title);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sDeltaPt_vs_Pt_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sDeltaPt_vs_Pt_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      hTrkPtRes[k] = (TH1F*)hResVsTruePt[k]->ProjectionX(Form("PrimTrkRes_%s",cent_Title[k]));
      hTrkPtRes[k]->Reset();
      hTrkPtShift[k] = (TH1F*)hResVsTruePt[k]->ProjectionX(Form("PrimTrkShift_%s",cent_Title[k]));
      hTrkPtShift[k]->Reset();

      TCanvas *c = new TCanvas(Form("FitTrkRes_cent%d",k),Form("FitTrkRes_cent%d",k),1200,800);
      c->Divide(10,10);
      for(int ibin=1; ibin<=hTrkPtRes[k]->GetNbinsX(); ibin++)
	{
	  TH1F *htmp = (TH1F*)hResVsTruePt[k]->ProjectionY(Form("TrkPt_bin%d_cent%d",ibin,k),ibin,ibin);
	  htmp->GetXaxis()->SetRangeUser(-0.2,0.2);
	  htmp->SetMarkerStyle(20);
	  TF1 *func = new TF1(Form("func_bin%d",ibin),"gaus",-0.15,0.15);
	  htmp->Fit(func,"IR0Q");
	  c->cd(ibin);
	  htmp->Draw("HIST");
	  func->SetLineColor(2);
	  func->Draw("same");
	  hTrkPtRes[k]->SetBinContent(ibin,func->GetParameter(2));
	  hTrkPtRes[k]->SetBinError(ibin,func->GetParError(2));
	  hTrkPtShift[k]->SetBinContent(ibin,func->GetParameter(1));
	  hTrkPtShift[k]->SetBinError(ibin,func->GetParError(1));
	}

      TH1F *htmp = (TH1F*)hTrkPtRes[k]->Clone(Form("%s_clone",hTrkPtRes[k]->GetName()));
      list->Add(htmp);

      funcRes[k] = new TF1(Form("FuncPrimTrkRes_vs_TruePt_%s",cent_Title[k]),"sqrt([0]^2*x^2+[1]^2)",0.8,10);
      funcRes[k]->SetParNames("a","b");
      funcRes[k]->SetParameter(0,0.005);
      hTrkPtRes[k]->Fit(funcRes[k],"IR0");
      hTrkPtRes[k]->SetMarkerStyle(21);
      hTrkPtRes[k]->GetYaxis()->SetRangeUser(0,0.08);
      c = draw1D(hTrkPtRes[k],Form("Transverse momentum resolution of primary tracks (%s%%);p_{T,true} (GeV/c);#sigma(p_{T})/p_{T}",cent_Name[k]));
      funcRes[k]->SetLineColor(2);
      funcRes[k]->Draw("same");
      leg = new TLegend(0.2,0.5,0.4,0.7);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hTrkPtRes[k],"Embedding data","P");
      leg->AddEntry(funcRes[k],"Fit: #sqrt{(a*p_{T})^{2}+b^{2}}","L");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sFitTrkPtRes_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sFitTrkPtRes_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}
    }
  bool drawLegend = true;
  if(year==2013) drawLegend = false;
  c = drawHistos(list,Form("TrkRes"),Form("p_{T} resolution of muon tracks;p_{T,true} (GeV/c);#sigma(p_{T})/p_{T}"),kTRUE,0,11,kFALSE,0.1,10,kFALSE,drawLegend,legName_cent,drawLegend,"",0.3,0.5,0.6,0.85,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtRes_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtRes_CentBins.png",run_type,run_cfg_name.Data()));
    }

  list->Clear();
  for(int k=0; k<nCentBins; k++)
    {
      TH1F *htmp = (TH1F*)hTrkPtShift[k]->Clone(Form("%s_clone",hTrkPtShift[k]->GetName()));
      htmp->Scale(100);
      list->Add(htmp);
    }
  c = drawHistos(list,Form("TrkPtShift"),Form("p_{T} shift of muon tracks;p_{T,true} (GeV/c);<#Deltap_{T}/p_{T}> (%%)"),kTRUE,0,11,kTRUE,-0.2,0.2,kFALSE,drawLegend,legName_cent,drawLegend,"",0.3,0.5,0.6,0.85,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtShift_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtShift_CentBins.png",run_type,run_cfg_name.Data()));
    }

  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  hResVsTruePt[k]->Write("",TObject::kOverwrite);
	  hTrkPtRes[k]->Write("",TObject::kOverwrite);
	  hTrkPtShift[k]->Write("",TObject::kOverwrite);
	  funcRes[k]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void makeHistos(TString outName, const bool save = 1)
{
  gStyle->SetOptStat(1);

  //const int nCentBins = 1;
  
  printf("+++ Single track resolution +++\n");
  // single track efficiency & resolution
  THnSparseF *hnTrkPtRes = (THnSparseF*)f->Get("mhpTrkPtRes_di_mu");
  TH2F *hResVsRecoPt[nCentBins];
  TH2F *hResVsTruePt[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      if(nCentBins>1) hnTrkPtRes->GetAxis(3)->SetRange(centBins_low[k],centBins_high[k]);
      hResVsRecoPt[k] = (TH2F*)hnTrkPtRes->Projection(0,1);
      hResVsRecoPt[k]->SetName(Form("PrimTrkRes_vs_RecoPt_%s",cent_Title[k]));

      hResVsTruePt[k] = (TH2F*)hnTrkPtRes->Projection(0,2);
      hResVsTruePt[k]->SetName(Form("PrimTrkRes_vs_TruePt_%s",cent_Title[k]));

      hnTrkPtRes->GetAxis(3)->SetRange(0,-1);
    }

  printf("+++ Single track efficiency +++\n");
  THnSparseF *hMcTrkInfo[nDet];
  TH1F *hMcTrkPt[nCentBins][nDet][3];
  TH1F *hMcTrkPtEff[nCentBins][nDet][3];
  TH3F *hMcTrkPtEtaPhi[nCentBins][nDet][3];
  TH3F *hMcTrkPtEtaPhiEff[nCentBins][nDet][3];
  const int nbinsX = 23, nbinsY = 8, nbinsZ = 360;
  const double xbinsX[nbinsX+1] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6,7,8,9,10,15,20};
  const double xbinsY[nbinsY+1] = {-0.7,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,0.7};
  double xbinsZ[nbinsZ+1];
  for(int i=0; i<=nbinsZ; i++)
    {
      xbinsZ[i] = pi/12 + 2*pi/360*i;
    }
  
  for(int i=0; i<nDet; i++)
    {
      hMcTrkInfo[i] = (THnSparseF*)f->Get(Form("mhMcTrkInfo%s_di_mu",det_name[i]));
      hMcTrkInfo[i]->SetTitle(Form("%s_jpsi",hMcTrkInfo[i]->GetName()));
      hMcTrkInfo[i]->Sumw2();

      for(int k=0; k<nCentBins; k++)
	{
	  if(nCentBins>1) hMcTrkInfo[i]->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkInfo[i]->GetAxis(3)->SetRange(1,3);
	      if(j==1) hMcTrkInfo[i]->GetAxis(3)->SetRange(3,3);
	      if(j==2) hMcTrkInfo[i]->GetAxis(3)->SetRange(1,1);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkInfo[i]->Projection(0);
	      hMcTrkPt[k][i][j]->Sumw2();
	      hMcTrkPt[k][i][j]->SetName(Form("McTrkPt%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);

	      TH3F *h3 = (TH3F*)hMcTrkInfo[i]->Projection(0,1,2);
	      h3->SetName(Form("McTrkPtEtaPhi%d%d%d",i,j,k));
	      //hMcTrkPtEtaPhi[k][i][j] = new TH3F(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),h3->GetTitle(),200,0,20,18,-0.9,0.9,360,pi/12,2*pi+pi/12);
	      hMcTrkPtEtaPhi[k][i][j] = new TH3F(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]),h3->GetTitle(),nbinsX,xbinsX,nbinsY,xbinsY,nbinsZ,xbinsZ);
	      hMcTrkPtEtaPhi[k][i][j]->Sumw2();
	      for(int binx=1; binx<=nbinsX; binx++)
		{
		  for(int biny=1; biny<=nbinsY; biny++)
		    {
		      for(int binz=1; binz<=nbinsZ; binz++)
			{
			  int new_binz = binz+15;
			  if(binz>360) new_binz = binz - 360;
			  int new_binx_low = h3->GetXaxis()->FindFixBin(xbinsX[binx-1]+1e-4);
			  int new_binx_hig = h3->GetXaxis()->FindFixBin(xbinsX[binx]-1e-4);
			  int new_biny_low = h3->GetYaxis()->FindFixBin(xbinsY[biny-1]+1e-4);
			  int new_biny_hig = h3->GetYaxis()->FindFixBin(xbinsY[biny]-1e-4);
			  double value = 0;
			  for(int new_binx=new_binx_low; new_binx<=new_binx_hig; new_binx++)
			    {
			      for(int new_biny=new_biny_low; new_biny<=new_biny_hig; new_biny++)
				{
				  value += h3->GetBinContent(new_binx, new_biny, new_binz);
				}
			    }
			  hMcTrkPtEtaPhi[k][i][j]->SetBinContent(binx,biny,binz,value);
			  hMcTrkPtEtaPhi[k][i][j]->SetBinError(binx,biny,binz,sqrt(value));
			}
		    }
		}
	      hMcTrkPtEtaPhi[k][i][j]->Sumw2();
	      hMcTrkPtEtaPhiEff[k][i][j] = (TH3F*)hMcTrkPtEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPtEtaPhi[k][i][j]->GetName()));
	      hMcTrkPtEtaPhiEff[k][i][j]->Divide(hMcTrkPtEtaPhi[k][0][j]);
	      hMcTrkInfo[i]->GetAxis(3)->SetRange(0,-1);
	    }
	  hMcTrkInfo[i]->GetAxis(4)->SetRange(0,-1);
	}
    }

  printf("+++ MTD efficiency +++\n");
  TH1F *hMtdEff[30][5];
  TH3F *hMtdInPtBlMod = (TH3F*)f->Get("mhTrkPtBlModMtd_di_mu");
  TH3F *hMtdOutPtBlMod = (TH3F*)f->Get("mhTrkPtBlModTrig_di_mu");
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  int bl = i+1; 
	  int mod = j+1;
	  hMtdEff[i][j] = (TH1F*)hMtdOutPtBlMod->ProjectionX(Form("RecoTrkPt_MtdEff_BL%d_Mod%d",bl,mod),bl,bl,mod,mod);
	  hMtdEff[i][j]->Rebin(5);
	  hMtdEff[i][j]->Sumw2();
	  TH1F *htmp = (TH1F*)hMtdInPtBlMod->ProjectionX(Form("Temp_BL%d_Mod%d",bl,mod),bl,bl,mod,mod);
	  htmp->Rebin(5);
	  hMtdEff[i][j]->Divide(htmp);
	}
    }
  
  if(save)
    {
      printf("+++ Save histograms +++\n");
      cout << outName.Data() << endl;
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outName.Data()),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  // single tracks
	  hResVsRecoPt[k]->Write();
	  hResVsTruePt[k]->Write();
	  for(int i=0; i<nDet; i++)
	    {
	      for(int j=0; j<3; j++)
		{
		  hMcTrkPt[k][i][j]->Write();
		  hMcTrkPtEff[k][i][j]->Write();
		  hMcTrkPtEtaPhi[k][i][j]->Write();
		  hMcTrkPtEtaPhiEff[k][i][j]->Write();
		}
	    }
	}

      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      hMtdEff[i][j]->Write();
	    }
	}
    }
}

//================================================
void TrkEff3D(TString inName, TString outPDFName, const bool savePlot = 0, const bool saveHisto = 0)
{
  if(year==2013) const int nCentBins = 1;

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);

  //----------------------------------------------------------------------------
  // Title page
  TPaveText *t1 = GetPaveText(0.28,0.7,0.5,0.7,0.07,62);
  t1->AddText("Efficiency of single muons by embedding");
  if(year==2013) t1->AddText("in pp 500 GeV from Run13");
  if(year==2014) t1->AddText("in Au+Au 200 GeV from Run14");
  t1->Draw();
  c1->Print(Form("%s(",outPDFName.Data()));

  //----------------------------------------------------------------------------
  // Get TPDF
  TPDF *pdf = 0;
  TSeqCollection *col = gROOT->GetListOfSpecials();
  for(Int_t i=0; i<col->GetEntries(); i++)
    {
      TObject *obj = (TObject*)col->At(i);
      if( obj->IsA()->InheritsFrom("TPDF") && strcmp(obj->GetName(),outPDFName.Data())==0 )
        pdf = dynamic_cast<TPDF*>obj;
    }
  if(!pdf) 
    {
      printf("No pointer to PDF file is available.\n");
      return;
    }
  pdf->Off();

  t1 = GetPaveText(0.28,0.7,0.5,0.7,0.05,62);
  t1->AddText("TPC tracking efficiency");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);


  // single tracks
  TFile *fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"read");
  TH3F *hMcTrkPtEtaPhi[nCentBins][nDet][3];
  TH1F *hMcTrkPt[nCentBins][nDet][3];
  TH1F *hMcTrkPtEff[nCentBins][nDet][3];
  TH2F *hMcTrkEtaPhi[nCentBins][nDet][3];
  TH2F *hMcTrkEtaPhiEff[nCentBins][nDet][3];
  TH1F *hMcTrkEta[nCentBins][nDet][3];
  TH1F *hMcTrkEtaEff[nCentBins][nDet][3];
  TH1F *hMcTrkPhi[nCentBins][nDet][3];
  TH1F *hMcTrkPhiEff[nCentBins][nDet][3];
  double eta_bounds[3] = {-0.5,0,0.5};
  int style[3] = {20,24,24};
  TString legName_charge[3] = {"#mu^{+}+#mu^{-}","#mu^{+}","#mu^{-}"};
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<nDet; i++)
	    {
	      hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      int low_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(-0.5+1e-4);
	      int high_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(0.5-1e-4);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_%d",det_name[i],charge_name[j],k),low_ybin,high_ybin,0,-1);
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);
	      hMcTrkPtEff[k][i][j]->SetMarkerStyle(style[j]);
	      hMcTrkPtEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkPtEff[k][i][j]->SetLineColor(color[j]);

	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRangeUser(1.5,20);
	      hMcTrkEtaPhi[k][i][j] = (TH2F*)hMcTrkPtEtaPhi[k][i][j]->Project3D("zy");
	      hMcTrkEtaPhi[k][i][j]->SetName(Form("McTrkEtaPhi%s%s_%d",det_name[i],charge_name[j],k));
	      hMcTrkEtaPhiEff[k][i][j] = (TH2F*)hMcTrkEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkEtaPhi[k][i][j]->GetName()));
	      hMcTrkEtaPhiEff[k][i][j]->Divide(hMcTrkEtaPhi[k][0][j]);
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRange(0,-1);

	      int low_xbin = hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->FindFixBin(1.5);
	      hMcTrkEta[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionY(Form("McTrkEta%s%s_%d",det_name[i],charge_name[j],k),low_xbin,-1,0,-1);
	      hMcTrkPhi[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionZ(Form("McTrkPhi%s%s_%d",det_name[i],charge_name[j],k),low_xbin,-1,low_ybin,high_ybin);
	      hMcTrkEtaEff[k][i][j] = (TH1F*)hMcTrkEta[k][i][j]->Clone(Form("%s_Eff",hMcTrkEta[k][i][j]->GetName()));
	      hMcTrkEtaEff[k][i][j]->Divide(hMcTrkEta[k][0][j]);
	      hMcTrkEtaEff[k][i][j]->SetMarkerStyle(style[j]);
	      hMcTrkEtaEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->SetLineColor(color[j]);

	      hMcTrkPhiEff[k][i][j] = (TH1F*)hMcTrkPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPhi[k][i][j]->GetName()));
	      hMcTrkPhiEff[k][i][j]->Divide(hMcTrkPhi[k][0][j]);
	      hMcTrkPhiEff[k][i][j]->SetMarkerStyle(style[j]);
	      hMcTrkPhiEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkPhiEff[k][i][j]->SetLineColor(color[j]);
	    }
	}
      c = new TCanvas(Form("McTrkPtEff_Cent%d",k),Form("McTrkPtEff_Cent%d",k),800,600);
      leg = new TLegend(0.3,0.65,0.5,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      if(nCentBins>1) leg->SetHeader(Form("%s%%",cent_Name[k]));
      for(int j=0; j<3; j++)
	{
	  hMcTrkPtEff[k][1][j]->GetXaxis()->SetRangeUser(0,10);
	  hMcTrkPtEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.4);
	  hMcTrkPtEff[k][1][j]->SetTitle(";p_{T,mc} (GeV/c);Efficiency");
	  if(j==0) hMcTrkPtEff[k][1][j]->Draw("P");
	  else     hMcTrkPtEff[k][1][j]->Draw("Psames");
	  leg->AddEntry(hMcTrkPtEff[k][1][j],legName_charge[j],"P");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (|#eta_{mc}| < 0.5)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPtTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      c = new TCanvas(Form("McTrkEtaEff_Cent%d",k),Form("McTrkEtaEff_Cent%d",k),800,600);
      for(int j=0; j<3; j++)
	{
	  hMcTrkEtaEff[k][1][j]->GetXaxis()->SetRangeUser(-1,1);
	  hMcTrkEtaEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.4);
	  hMcTrkEtaEff[k][1][j]->SetTitle(";#eta_{mc};Efficiency");
	  if(j==0) hMcTrkEtaEff[k][1][j]->Draw("P");
	  else     hMcTrkEtaEff[k][1][j]->Draw("PEsames");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (p_{T,mc} > 1.5 GeV/c)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      c = new TCanvas(Form("McTrkPhiEff_Cent%d",k),Form("McTrkPhiEff_Cent%d",k),800,600);
      for(int j=0; j<3; j++)
      	{
      	  hMcTrkPhiEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.4);
      	  hMcTrkPhiEff[k][1][j]->SetTitle(";#varphi_{mc};Efficiency");
      	  if(j==0) hMcTrkPhiEff[k][1][j]->Draw("HIST");
      	  else     hMcTrkPhiEff[k][1][j]->Draw("sames HIST");
      	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (|#eta_{mc}|<0.5)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPhiTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkPhiTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

     c = new TCanvas(Form("McTrkEtaPhiEff_Cent%d",k),Form("McTrkEtaPhiEff_Cent%d",k),1100,700);
     c->Divide(2,2);
     for(int j=0; j<3; j++)
      	{
	  c->cd(j+1);
     	  hMcTrkEtaPhiEff[k][1][j]->GetXaxis()->SetRangeUser(-1,1);
      	  hMcTrkEtaPhiEff[k][1][j]->GetZaxis()->SetRangeUser(0,1);
      	  hMcTrkEtaPhiEff[k][1][j]->SetTitle(";#eta_{mc};#varphi_{mc}");
	  ScaleHistoTitle(hMcTrkEtaPhiEff[k][1][j],0.06,0.7,0.05,0.06,0.7,0.05,62);
	  hMcTrkEtaPhiEff[k][1][j]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency for %s",legName_charge[j].Data()),0.06);
	  t1->Draw();
      	}
     c->cd(4);
     TPaveText *t1 = GetPaveText(0.3,0.5,0.4,0.6,0.07,62);
     t1->AddText(Form("%s%%",cent_Name[k]));
     t1->Draw();
     PrintCanvasToPDF(c,pdf);
     if(savePlot)
       {
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaPhiTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%sMcTrkEtaPhiTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
       }
    }

  TList *list = new TList;
  TString legName[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      list->Add(hMcTrkPtEff[k][3][0]);
      legName[k] = Form("%s%%",cent_Name[k]);
    }
  c = drawHistos(list,Form("McTrkEff"),Form("Efficiency of embedded tracks;p_{T} (GeV/c);Efficiency"),kTRUE,0,10,kTRUE,0,1.4,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.3,0.45,0.6,0.85,kTRUE);
  PaintCanvasToPDF(c,pdf);


  // in (phi,eta) bins
  const int index = 3; // 1 - tpc; 2 - mtd; 3 - final
  const int neta = 8;
  const int nphi = 180;
  const int nCanvas = nphi*neta/15;
  TCanvas *cTpcEff[nCanvas];
  for(int i=0; i<nCanvas; i++)
    {
      cTpcEff[i] = new TCanvas(Form("TpcEff_%d",i),Form("TpcEff_%d",i),1200,750);
      cTpcEff[i]->Divide(5,3);
    }
  TH3F *hMcTrkPtEtaPhiEff[nCentBins][3];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  hMcTrkPtEtaPhiEff[k][j] = (TH3F*)hMcTrkPtEtaPhi[k][index][j]->Clone(Form("%s_Eff",hMcTrkPtEtaPhi[k][index][j]->GetName()));
	  int nbinsy = hMcTrkPtEtaPhiEff[k][j]->GetNbinsY();
	  int nbinsz = hMcTrkPtEtaPhiEff[k][j]->GetNbinsZ();
	  hMcTrkPtEtaPhiEff[k][j]->RebinZ(nbinsz/nphi);
	  TH3F *h3tmp = (TH3F*)hMcTrkPtEtaPhi[k][0][j]->Clone(Form("%s_tmp",hMcTrkPtEtaPhi[k][0][j]->GetName()));
	  h3tmp->RebinZ(nbinsz/nphi);
	  hMcTrkPtEtaPhiEff[k][j]->Divide(h3tmp);
	}
    }

  for(int k=0; k<1; k++)
    {
      if(nCentBins>1)
	{
	  t1 = GetPaveText(0.28,0.7,0.5,0.7,0.05,62);
	  if(index==1) t1->AddText(Form("%s%%: TPC tracking",cent_Name[k]));
	  if(index==2) t1->AddText(Form("%s%%: TPC tracking + MTD matching",cent_Name[k]));
	  if(index==3) t1->AddText(Form("%s%%: total efficiency",cent_Name[k]));
	  c1->Clear();
	  t1->Draw();
	  PaintCanvasToPDF(c1,pdf);
	}

      leg = new TLegend(0.6,0.2,0.9,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      for(int ieta=1; ieta<=neta; ieta++)
	{
	  for(int iphi=1; iphi<=nphi; iphi++)
	    {
	      int ind = (ieta-1)*nphi + iphi - 1;
	      cTpcEff[ind/15]->cd(ind%15+1);
	      SetPadMargin(gPad,0.15,0.15,0.05,0.05);
	      for(int j=1; j<3; j++)
		{
		  TH1F *htmp = (TH1F*)hMcTrkPtEtaPhiEff[k][j]->ProjectionX(Form("%s_%d_%d",hMcTrkPtEtaPhiEff[k][j]->GetName(),ieta,iphi),ieta,ieta,iphi,iphi);
		  ScaleHistoTitle(htmp,0.06,1.1,0.05,0.06,1.1,0.05,62);
		  htmp->GetXaxis()->SetRangeUser(0,20);
		  htmp->GetYaxis()->SetRangeUser(0,1.2);
		  htmp->SetMarkerStyle(20+(j-1)*4);
		  if(index==1) htmp->SetTitle(";p_{T,mc} (GeV/c);Tpc Efficiency");
		  if(index==2) htmp->SetTitle(";p_{T,mc} (GeV/c);Tpc+Mtd Efficiency");
		  if(index==3) htmp->SetTitle(";p_{T,mc} (GeV/c);Total Efficiency");
		  if(j==1) htmp->Draw("PE");
		  else     htmp->Draw("PEsames");
		  TPaveText *t1 = GetPaveText(0.25,0.6,0.75,0.9,0.06);
		  t1->AddText(Form("%1.1f < #eta < %1.1f",hMcTrkPtEtaPhiEff[k][j]->GetYaxis()->GetBinLowEdge(ieta),hMcTrkPtEtaPhiEff[k][j]->GetYaxis()->GetBinUpEdge(ieta)));
		  t1->AddText(Form("%1.0f < #varphi < %1.0f",hMcTrkPtEtaPhiEff[k][j]->GetZaxis()->GetBinLowEdge(iphi)/pi*180,hMcTrkPtEtaPhiEff[k][j]->GetZaxis()->GetBinUpEdge(iphi)/pi*180));
		  t1->Draw();

		  if(ieta==1 && iphi==1)
		    {
		      if(j==1) leg->AddEntry(htmp,"#mu^{+}","P");
		      if(j==2) leg->AddEntry(htmp,"#mu^{-}","P");
		    }
		}
	      if(ind%15==0) leg->Draw();
	    }
	}
    }
  
  for(int i=0; i<nCanvas; i++)
    {
      PrintCanvasToPDF(cTpcEff[i],pdf);
    }

  // in (backleg,module)
  t1 = GetPaveText(0.28,0.7,0.5,0.7,0.05,62);
  t1->AddText(Form("%s%%: MTD efficiencies",cent_Name[0]));
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  TCanvas *cMtdEffBlMod[10];
  for(int i=0; i<10; i++)
    {
      cMtdEffBlMod[i] = new TCanvas(Form("MtdEff_%d",i),Form("MtdEff_%d",i),1200,750);
      cMtdEffBlMod[i]->Divide(5,3);
    }
  TH1F *hRcTrkMtdEff[30][5];
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  int bl = i+1; 
	  int mod = j+1;
	  hRcTrkMtdEff[i][j] = (TH1F*)fin->Get(Form("RecoTrkPt_MtdEff_BL%d_Mod%d",bl,mod));
	  cMtdEffBlMod[i/3]->cd((i%3)*5+j+1);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.05);
	  ScaleHistoTitle(hRcTrkMtdEff[i][j],0.06,1.1,0.05,0.06,1.1,0.05,62);
	  hRcTrkMtdEff[i][j]->GetXaxis()->SetRangeUser(0,20);
	  hRcTrkMtdEff[i][j]->GetYaxis()->SetRangeUser(0,1.2);
	  hRcTrkMtdEff[i][j]->SetMarkerStyle(20);
	  hRcTrkMtdEff[i][j]->SetTitle(";p_{T,rec} (GeV/c);Mtd Efficiency");
	  hRcTrkMtdEff[i][j]->Draw();
	  TPaveText *t1 = GetPaveText(0.25,0.6,0.75,0.9,0.06);
	  t1->AddText(Form("BL = %d, Mod = %d",bl,mod));
	  t1->Draw();
	}
    }
  for(int i=0; i<10; i++)
    {
      PrintCanvasToPDF(cMtdEffBlMod[i],pdf);
    }

  pdf->On();
  pdf->Close();


  TH2F *hMcTrkPtEtaMtd[nCentBins][3][2][3];
  TH2F *hMcTrkPtEtaMtdEff[nCentBins][3][2];
  for(int k=0; k<0; k++)
    {
      TCanvas *cMtdEff = new TCanvas(Form("MtdMthEff_%d",k),Form("MtdMthEff_%d",k),1200,700);
      cMtdEff->Divide(5,4);
      for(int j=0; j<3; j++)
	{
	  int nbin1 = hMcTrkPtEtaPhi[k][1][j]->GetZaxis()->FindFixBin(3.77);
	  int nbin2 = hMcTrkPtEtaPhi[k][1][j]->GetZaxis()->FindFixBin(5.65);
	  double bins[4] = {0,nbin1,nbin2,hMcTrkPtEtaPhi[k][1][j]->GetNbinsZ()};
	  for(int m=0; m<2; m++)
	    {
	      for(int l=0; l<3; l++)
		{
		  hMcTrkPtEtaPhi[k][m+1][j]->GetZaxis()->SetRange(bins[l]+1,bins[l+1]);
		  hMcTrkPtEtaMtd[k][j][m][l] = (TH2F*)hMcTrkPtEtaPhi[k][m+1][j]->Project3D("yx");
		  hMcTrkPtEtaMtd[k][j][m][l]->SetName(Form("%s_%d",hMcTrkPtEtaPhi[k][m+1][j]->GetName(),l));
		}
	    }

	  hMcTrkPtEtaMtdEff[k][j][1] = (TH2F*)hMcTrkPtEtaMtd[k][j][1][0]->Clone(Form("McTrkPtEtaMtdEff%s_5tray_%s",charge_name[j],cent_Title[k]));
	  hMcTrkPtEtaMtdEff[k][j][1]->Add(hMcTrkPtEtaMtd[k][j][1][2]);
	  hMcTrkPtEtaMtdEff[k][j][1]->Rebin2D(5,1);
	  h2tmp = (TH2F*)hMcTrkPtEtaMtd[k][j][0][0]->Clone(Form("McTrkPtEtaMtdEff%s_%s_tmp2",charge_name[j],cent_Title[k]));
	  h2tmp->Add(hMcTrkPtEtaMtd[k][j][0][2]);
	  h2tmp->Rebin2D(5,1);
	  hMcTrkPtEtaMtdEff[k][j][1]->Divide(h2tmp);
	  draw2D(hMcTrkPtEtaMtdEff[k][j][1]);

	  int begin = hMcTrkPtEtaMtdEff[k][j][1]->GetYaxis()->FindFixBin(-0.6);
	  int end = hMcTrkPtEtaMtdEff[k][j][1]->GetYaxis()->FindFixBin(0.6);
	  for(int ieta=begin; ieta<=end; ieta++)
	    {
	      cMtdEff->cd(ieta-begin+1);
	      TH1F *htmp = (TH1F*)hMcTrkPtEtaMtdEff[k][j][1]->ProjectionX(Form("%s_%d",hMcTrkPtEtaMtdEff[k][j][1]->GetName(),ieta),ieta,ieta);
	      htmp->SetMarkerStyle(style[j]+4);
	      htmp->SetMarkerColor(color[j]);
	      htmp->SetLineColor(color[j]);
	      if(j==0) htmp->Draw();
	      else     htmp->Draw("sames");
	    }
	  
	  hMcTrkPtEtaMtdEff[k][j][0] = (TH2F*)hMcTrkPtEtaMtd[k][j][1][1]->Clone(Form("McTrkPtEtaMtdEff%s_3tray_%s",charge_name[j],cent_Title[k]));
	  hMcTrkPtEtaMtdEff[k][j][0]->Rebin2D(5,1);
	  TH2F *h2tmp = (TH2F*)hMcTrkPtEtaMtd[k][j][0][1]->Clone(Form("McTrkPtEtaMtdEff%s_%s_tmp",charge_name[j],cent_Title[k]));
	  h2tmp->Rebin2D(5,1);
	  hMcTrkPtEtaMtdEff[k][j][0]->Divide(h2tmp);
	  draw2D(hMcTrkPtEtaMtdEff[k][j][0]);

	  for(int ieta=begin; ieta<=end; ieta++)
	    {
	      cMtdEff->cd(ieta-begin+1);
	      TH1F *htmp = (TH1F*)hMcTrkPtEtaMtdEff[k][j][0]->ProjectionX(Form("%s_%d",hMcTrkPtEtaMtdEff[k][j][0]->GetName(),ieta),ieta,ieta);
	      htmp->SetMarkerStyle(style[j]);
	      htmp->SetMarkerColor(color[j]);
	      htmp->SetLineColor(color[j]);
	      htmp->Draw("sames");
	    }
	}
      //PrintCanvasToPDF(cMtdEff,pdf);
    }

  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TH2F *hRes = (TH2F*)fin->Get("PrimTrkRes_vs_TruePt_0080");

      TFile *fout = TFile::Open(Form("Rootfiles/%s.ToyMC.root",run_type),"recreate");
      hRes->Write("PrimTrkRes_vs_TruePt",TObject::kOverwrite);
      for(int k=0; k<nCentBins; k++)
	{
	  for(int j=0; j<3; j++)
	    {
	      hMcTrkPtEtaPhiEff[k][j]->Write("",TObject::kOverwrite);
	      // for(int l=0; l<2; l++)
	      //   hMcTrkPtEtaMtdEff[k][j][l]->Write("",TObject::kOverwrite);
	    }
	}
      
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      int bl = i+1; 
	      int mod = j+1;
	      //hRcTrkMtdEff[i][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

