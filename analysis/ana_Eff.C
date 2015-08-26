const char *run_config = "";
const int year = 2013;
TString run_cfg_name;
const double pt1_cut = 1.5;
const double pt2_cut = 1.0;
const Double_t low_mass = 2.9;
const Double_t high_mass = 3.3;
const char *det_name[5] = {"","Tpc","Mtd","Muon","MtdFake"};
const char *charge_name[3] = {"","_pos","_neg"};
const char *charge_title[3] = {" "," positive "," negative "};

TFile *f;

//================================================
void ana_Eff()
{
  gStyle->SetOptStat(0);
  TString fileName;
  TString outName;
  TString outPDF;

  if(year==2013)
    {
      const int nCentBins = 1;
      run_type = "Run13_pp500";
      fileName = Form("Run13.pp500.jpsi.EmbedQA.MC.%sroot",run_config);
      outName = Form("Run13.pp500.jpsi.Eff.%sroot",run_config);
      outPDF = Form("Run13_pp500_Embed_Eff.pdf");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      fileName = Form("Run14.AuAu200.jpsi.EmbedQA.MC.%sroot",run_config);
      outName = Form("Run14.AuAu200.jpsi.Eff.%sroot",run_config);
      outPDF = Form("Run14_AuAu200_Embed_Eff.pdf");
    }
  run_cfg_name = run_config;

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("# of events: %4.4e\n",hStat->GetBinContent(3));

  //makeEffFromEmbed(outName);
  JpsiEffFromEmbed(outName, outPDF);
  //3DTrkEff(outName, outPDF);
}

//================================================
void JpsiEffFromEmbed(TString inName, TString outPDFName, const bool savePlot = 0)
{
 // single tracks
  TFile *fin = TFile::Open(Form("Rootfiles/%s",inName.Data()),"read");
  TH3F *hMcTrkPtEtaPhi[nCentBins][2][3];
  TH1F *hMcTrkPt[nCentBins][2][3];
  TH2F *hMcTrkEtaPhi[nCentBins][2][3];
  TH1F *hMcTrkEta[nCentBins][2][3];
  TH1F *hMcTrkPhi[nCentBins][2][3];
  double eta_bounds[3] = {-0.5,0,0.5};
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      if(nCentBins==1) hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s",det_name[i],charge_name[j]));
	      else             hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      int low_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(-0.5+1e-4);
	      int high_ybin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(0.5-1e-4);
	      hMcTrkPt[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_%d",det_name[i],charge_name[j],k),low_ybin,high_ybin,0,-1);

	      int low_xbin = hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->FindFixBin(1.5);
	      hMcTrkEta[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionY(Form("McTrkEta%s%s_%d",det_name[i],charge_name[j],k),low_xbin,-1,0,-1);
	      hMcTrkPhi[k][i][j] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionZ(Form("McTrkPhi%s%s_%d",det_name[i],charge_name[j],k),low_xbin,-1,low_ybin,high_ybin);
	      
	      // for(int ieta=0; ieta<2; ieta++)
	      // 	{
	      // 	  int low_etabin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(eta_bounds[ieta]+1e-4);
	      // 	  int high_etabin = hMcTrkPtEtaPhi[k][i][j]->GetYaxis()->FindFixBin(eta_bounds[ieta+1]-1e-4);
	      // 	  for(int iphi=0; iphi<12; iphi++)
	      // 	    {
	      // 	      int low_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12+pi/6*iphi + 1e-4);
	      // 	      int high_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12+pi/6*(iphi+1) - 1e-4);
	      // 	      hMcTrkPtInEtaPhi[k][i][j][ieta][iphi] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_eta%d_phi%d_cent%d",det_name[i],charge_name[j],ieta,iphi,k),low_etabin,high_etabin,low_phibin,high_phibin);
	      // 	      if(iphi==11)
	      // 		{
	      // 		  int low_phibin = 0;
	      // 		  int high_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/12 - 1e-4);
	      // 		  hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Add((TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_eta%d_phi%d_cent%d_2",det_name[i],charge_name[j],ieta,iphi,k),low_etabin,high_etabin,low_phibin,high_phibin));
	      // 		}
	      // 	      hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Rebin(2);
	      // 	    }
	      // 	}

	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRangeUser(1.5,20);
	      hMcTrkEtaPhi[k][i][j] = (TH2F*)hMcTrkPtEtaPhi[k][i][j]->Project3D("zy");
	      hMcTrkEtaPhi[k][i][j]->SetName(Form("McTrkEtaPhi%s%s_%d",det_name[i],charge_name[j],k));
	      hMcTrkPtEtaPhi[k][i][j]->GetXaxis()->SetRange(0,-1);

	      for(int iphi=0; iphi<30; iphi++)
		{
		  int low_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/15*iphi + 1e-4);
		  int high_phibin = hMcTrkPtEtaPhi[k][i][j]->GetZaxis()->FindFixBin(pi/15*(iphi+1) - 1e-4);
		  hMcTrkPtInBackleg[k][i][j][iphi] = (TH1F*)hMcTrkPtEtaPhi[k][i][j]->ProjectionX(Form("McTrkPt%s%s_phi%d_cent%d",det_name[i],charge_name[j],iphi,k),low_ybin,high_ybin,low_phibin,high_phibin);
		}
	    }
	}
    }

  TH1F *hMcTrkPtEff[nCentBins][5][3];
  TF1  *funcMcTrkPtEff[nCentBins][5][3];
  TH1F *hMcTrkPtInEtaPhiEff[nCentBins][5][3][2][12];
  TH2F *hMcTrkEtaPhiEff[nCentBins][5][3];
  TH1F *hMcTrkEtaEff[nCentBins][5][3];
  TH1F *hMcTrkPhiEff[nCentBins][5][3];
  TH1F *hMcTrkPtInBacklegEff[nCentBins][5][3][30];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<5; i++)
	    {
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);

	      hMcTrkEtaPhiEff[k][i][j] = (TH2F*)hMcTrkEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkEtaPhi[k][i][j]->GetName()));
	      hMcTrkEtaPhiEff[k][i][j]->Divide(hMcTrkEtaPhi[k][0][j]);

	      hMcTrkEtaEff[k][i][j] = (TH1F*)hMcTrkEta[k][i][j]->Clone(Form("%s_Eff",hMcTrkEta[k][i][j]->GetName()));
	      hMcTrkEtaEff[k][i][j]->Divide(hMcTrkEta[k][0][j]);

	      hMcTrkPhiEff[k][i][j] = (TH1F*)hMcTrkPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPhi[k][i][j]->GetName()));
	      hMcTrkPhiEff[k][i][j]->Divide(hMcTrkPhi[k][0][j]);

	      for(int ieta=0; ieta<2; ieta++)
		{
		  for(int iphi=0; iphi<12; iphi++)
		    {
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi] = (TH1F*)hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->Clone(Form("%s_Eff",hMcTrkPtInEtaPhi[k][i][j][ieta][iphi]->GetName()));
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->Divide(hMcTrkPtInEtaPhi[k][0][j][ieta][iphi]);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->SetMarkerStyle(21);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->SetMarkerColor(color[j]);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->SetLineColor(color[j]);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->GetXaxis()->SetRangeUser(0,11);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->GetYaxis()->SetRangeUser(0,1);
		      hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi]->SetTitle(";p_{T,mc} [GeV/c];TPC Efficiency");
		      ScaleHistoTitle(hMcTrkPtInEtaPhiEff[k][i][j][ieta][iphi],0.06,1,0.05,0.06,1,0.05,62);
		    }
		}

	      for(int iphi=0; iphi<30; iphi++)
		{
		  hMcTrkPtInBacklegEff[k][i][j][iphi] = (TH1F*)hMcTrkPtInBackleg[k][i][j][iphi]->Clone(Form("%s_Eff",hMcTrkPtInBackleg[k][i][j][iphi]->GetName()));
		  hMcTrkPtInBacklegEff[k][i][j][iphi]->Divide(hMcTrkPtInBackleg[k][0][j][iphi]);
		}
	    }
	  for(int i=0; i<5; i++)
	    {
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);
	    }

	  c = new TCanvas(Form("FitMcTrkEff_%s_Cent%d",charge_name[j],k),Form("FitMcTrkEff_%s_Cent%d",charge_name[j],k),1100,650);
	  c->Divide(2,2);
	  int counter = 1;
	  for(int i=1; i<5; i++)
	    {
	      if(i==3) continue;
	      funcMcTrkPtEff[k][i][j] = new TF1(Form("Fit_%s",hMcTrkPtEff[k][i][j]->GetName()),"pol3",1,10);
	      hMcTrkPtEff[k][i][j]->Fit(funcMcTrkPtEff[k][i][j],"IR0");
	      c->cd(counter);
	      counter++;
	      hMcTrkPtEff[k][i][j]->Draw();
	      funcMcTrkPtEff[k][i][j]->SetLineColor(2);
	      funcMcTrkPtEff[k][i][j]->Draw("sames");
	    }

	  
	}
    }

  TList *list = new TList;
  TString legName[4] = {"MC input","Reconstructed in TPC","Match to MTD","Muon pid cuts"};
  const char *eff_type[5] = {"MC input","reconstructed in TPC","matched to MTD","of fake MTD matches","fulfill muon pid cuts"};
  TString legName_charge[3] = {"All","Positive","Negative"};
  for(int k=0; k<nCentBins; k++)
    {
      // efficiency vs eta vs phi
      for(int i=1; i<5; i++)
	{
	  if(i==3) continue;
	  hMcTrkEtaPhiEff[k][i][0]->GetXaxis()->SetRangeUser(-1,1);
	  ScaleHistoTitle(hMcTrkEtaPhiEff[k][i][0],0.05,0.9,0.04,0.05,0.8,0.04,62);
	  c = draw2D(hMcTrkEtaPhiEff[k][i][0],Form("Efficiency of MC muon tracks %s (p_{T,mc}>1.5 GeV/c)",eff_type[i]),0.04,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaPhiEff_%s_cent%s.pdf",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaPhiEff_%s_cent%s.png",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	    }
	}

      // efficiency vs eta
      c = new TCanvas(Form("McTrkEffEta_Cent%d",k),Form("McTrkEffEta_Cent%d",k),800,600);
      double xmin = 0.15, xmax = 0.3, ymin = 0.7, ymax = 0.88;
      int counter = 0;
      for(int i=1; i<5; i++)
	{
	  if(i==3) continue;
	  leg = new TLegend(xmin+counter*0.2,ymin,xmax+counter*0.2,ymax);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.035);
	  if(i==1) leg->SetHeader("TPC tracking");
	  if(i==2) leg->SetHeader("MTD matched");
	  if(i==4) leg->SetHeader("Muon PID");
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hMcTrkEtaEff[k][i][j]->SetMarkerStyle(19+i);
	      else     
		{
		  hMcTrkEtaEff[k][i][j]->SetMarkerStyle(23+i);
		  if(i==4) hMcTrkEtaEff[k][i][j]->SetMarkerStyle(32);
		}
	      hMcTrkEtaEff[k][i][j]->SetMarkerColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->SetLineColor(color[j]);
	      hMcTrkEtaEff[k][i][j]->GetXaxis()->SetRangeUser(-0.6,0.6);
	      hMcTrkEtaEff[k][i][j]->GetYaxis()->SetRangeUser(0,1.1);
	      hMcTrkEtaEff[k][i][j]->SetTitle("");
	      if(j==0 && i==1) hMcTrkEtaEff[k][i][j]->Draw("P");
	      else hMcTrkEtaEff[k][i][j]->Draw("sames");
	      leg->AddEntry(hMcTrkEtaEff[k][i][j],legName_charge[j].Data(),"P");
	    }
	  leg->Draw();
	  counter++;
	}
      TPaveText *t1 = GetTitleText("Efficiency of single track (p_{T,mc} > 1.5 GeV/c)");
      t1->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}
      
      // efficiency vs phi
      c = new TCanvas(Form("McTrkEffPhi_Cent%d",k),Form("McTrkEffPhi_Cent%d",k),800,600);
      leg = new TLegend(0.15,0.7,0.3,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      int counter = 0;
      for(int i=1; i<5; i++)
	{
	  if(i==3) continue;
	  TH1F *htmp = (TH1F*)hMcTrkPhiEff[k][i][0]->Clone(Form("%s_clone",hMcTrkPhiEff[k][i][0]->GetName()));
	  htmp->SetLineColor(color[counter]);
	  htmp->GetYaxis()->SetRangeUser(0,1.1);
	  htmp->SetTitle(";#eta_{mc};Efficiency");
	  if(i==1) htmp->DrawCopy("HIST");
	  else     htmp->DrawCopy("sames HIST");
	  leg->AddEntry(htmp,legName[counter+1].Data(),"L");
	  counter++;
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText(Form("Efficiency of MC muon tracks (p_{T,mc} > 1.5 GeV/c)"));
      t1->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPhiEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPhiEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      for(int i=1; i<5; i++)
	{
	  if(i==3) continue;
	  c = new TCanvas(Form("McTrkEffPhi_%s_Cent%d",det_name[i],k),Form("McTrkEffPhi_%s_Cent%d",det_name[i],k),800,600);
	  leg = new TLegend(xmin,ymin,xmax,ymax);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.035);
	  if(i==1) leg->SetHeader("TPC tracking");
	  if(i==2) leg->SetHeader("MTD matched");
	  if(i==4) leg->SetHeader("Muon PID");
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
	      hMcTrkPhiEff[k][i][j]->GetYaxis()->SetRangeUser(0,1.1);
	      hMcTrkPhiEff[k][i][j]->SetTitle("");
	      if(j==0 && i==1) hMcTrkPhiEff[k][i][j]->Draw("HIST");
	      else hMcTrkPhiEff[k][i][j]->Draw("sames HIST");
	      leg->AddEntry(hMcTrkPhiEff[k][i][j],legName_charge[j].Data(),"L");
	    }
	  leg->Draw();
	  TPaveText *t1 = GetTitleText(Form("Efficiency of MC muon tracks %s (p_{T,mc} > 1.5 GeV/c)",eff_type[i]));
	  t1->Draw();
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPhiEff_%s_cent%s.pdf",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPhiEff_%s_cent%s.png",run_type,run_cfg_name.Data(),det_name[i],cent_Title[k]));
	    }
	}

      // efficiency vs pt
      for(int j=0; j<3; j++)
	{
	  list->Clear();
	  for(int i=0; i<5; i++)
	    {
	      if(i!=3) list->Add(hMcTrkPt[k][i][j]);
	    }
	  char *centrality = Form(" (%s%%)",cent_Name[k]);
	  if(year==2013) centrality = "";
	  c = sysCompare(list,Form("McTrkEff_%s_Cent%d",charge_name[j],k),Form("p_{T} distribution of%sMC tracks%s",charge_title[j],centrality),"Efficiency of single track;p_{T}^{mc} (GeV/c);Efficiency",kTRUE,0,11,kFALSE,0.1,10,kTRUE,0,1,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.5,0.7,0.6,0.85,kTRUE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtEff%s_cent%s.pdf",run_type,run_cfg_name.Data(),charge_name[j],cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtEff%s_cent%s.png",run_type,run_cfg_name.Data(),charge_name[j],cent_Title[k]));
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
	  for(int iphi=0; iphi<12; iphi++)
	    {
	      c->cd(iphi+1);
	      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	      for(int j=0; j<3; j++)
		{
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
		  leg->AddEntry(hMcTrkPtInEtaPhiEff[k][1][1][ieta][iphi],"Positive","PL");
		  leg->AddEntry(hMcTrkPtInEtaPhiEff[k][1][2][ieta][iphi],"Negative","PL");
		}
	    }
	  c->cd(1);
	  leg->Draw();
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtEffInTpc_eta%d_cent%s.pdf",run_type,run_cfg_name.Data(),ieta,cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtEffInTpc_eta%d_cent%s.png",run_type,run_cfg_name.Data(),ieta,cent_Title[k]));
	    }
	}

      // efficiency vs pt in backleg
      for(int iphi=0; iphi<30; iphi++)
	{
      	  if(iphi%15==0) 
	    {
	      c = new TCanvas(Form("McTrkEffInBackleg_Cent%d_%d",k,iphi/15),Form("McTrkEffInBackleg_Cent%d_%d",k,iphi/15),1250,700);
	      c->Divide(5,3);
	    }
	  c->cd(iphi%15+1);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.05);
	  for(int j=0; j<3; j++)
	    {
	      hMcTrkPtInBacklegEff[k][2][j][iphi]->SetMarkerStyle(21);
	      hMcTrkPtInBacklegEff[k][2][j][iphi]->SetMarkerColor(color[j]);
	      hMcTrkPtInBacklegEff[k][2][j][iphi]->SetLineColor(color[j]);
	      hMcTrkPtInBacklegEff[k][2][j][iphi]->GetXaxis()->SetRangeUser(0,11);
	      hMcTrkPtInBacklegEff[k][2][j][iphi]->GetYaxis()->SetRangeUser(0,0.7);
	      hMcTrkPtInBacklegEff[k][2][j][iphi]->SetTitle(";p_{T,mc} [GeV/c];MTD match");
	      ScaleHistoTitle(hMcTrkPtInBacklegEff[k][2][j][iphi],0.06,1,0.05,0.06,1,0.05,62);
	      if(j==0) hMcTrkPtInBacklegEff[k][2][j][iphi]->Draw("HIST");
	      else     hMcTrkPtInBacklegEff[k][2][j][iphi]->Draw("sames HIST");
	    }
	  int backleg = iphi + 24;
	  if(backleg>30) backleg -= 30;
	  TPaveText *t1 = GetPaveText(0.25,0.55,0.83,0.88,0.06,62);
	  t1->AddText(Form("%d < #varphi_{mc} < %d",12*iphi,12*(iphi+1)));
	  t1->Draw();
	  if(iphi%15==0)
	    {
	      leg = new TLegend(0.4,0.4,0.8,0.65);
	      leg->SetBorderSize(0);
	      leg->SetFillColor(0);
	      leg->SetTextSize(0.06);
	      leg->AddEntry(hMcTrkPtInBacklegEff[k][2][0][iphi],"All","L");
	      leg->AddEntry(hMcTrkPtInBacklegEff[k][2][1][iphi],"Positive","L");
	      leg->AddEntry(hMcTrkPtInBacklegEff[k][2][2][iphi],"Negative","L");
	      c->cd(1);
	      leg->Draw();
	    }
	  if(savePlot && iphi%15==14)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtEffInMtd_%d_cent%s.pdf",run_type,run_cfg_name.Data(),iphi/15,cent_Title[k]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtEffInMtd_%d_cent%s.png",run_type,run_cfg_name.Data(),iphi/15,cent_Title[k]));
	    }
	}
    }

  // single track momentum resolution

}



//================================================
double erf(double *x, double *par)
{
  return par[0] * TMath::Erf(x[0]);
}

//================================================
void makeEffFromEmbed(TString outName, const bool save = 1)
{
  gStyle->SetOptStat(1);
  THnSparseF *hJpsiUS_mc = (THnSparseF*)f->Get("hJpsiInfo_di_mu");
  if(year==2013) const int nCentBins = 1;
  printf("+++ Jpsi distribution +++\n");
  // j/psi distribution
  TH1F *hJpsi[nCentBins][3][4];
  TH2F *hJpsiMassVsPt[nCentBins][3];
  TString name[3] = {"MCinput","TPCreco","MTDreco"};
  TString variable[4] = {"mass","pT","rapidity","phi"};
  const char *mc_name[4] = {"invariant mass","p_{T}","y","#varphi"};
  for(int k=0; k<nCentBins; k++)
    {
      if(nCentBins>1) hJpsiUS_mc->GetAxis(9)->SetRange(centBins_low[k],centBins_high[k]);
      for(int i=0; i<3; i++)
	{
	  hJpsiUS_mc->GetAxis(7)->SetRange(i+4,i+4);
	  if(i>0)
	    {
	      hJpsiUS_mc->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
	      hJpsiUS_mc->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
	    }
	  hJpsiMassVsPt[k][i] = (TH2F*)hJpsiUS_mc->Projection(0,1);
	  if(nCentBins==1)
	    hJpsiMassVsPt[k][i]->SetName(Form("%s_Jpsi_MassVsPt",name[i].Data()));
	  else
	    hJpsiMassVsPt[k][i]->SetName(Form("%s_Jpsi_MassVsPt_%s",name[i].Data(),cent_Title[k]));
	  
	  for(int j=0; j<4; j++)
	    {
	      if(i>0 && j>0) hJpsiUS_mc->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
	      hJpsi[k][i][j] = (TH1F*)hJpsiUS_mc->Projection(j);
	      
	      if(nCentBins==1)
		{
		  hJpsi[k][i][j]->SetName(Form("%s_Jpsi_%s",name[i].Data(),variable[j].Data()));
		  hJpsi[k][i][j]->SetTitle(Form("%s distribution j/psi (%s, pt1 > %1.1f GeV/c)",mc_name[j],name[i].Data(),pt1_cut)); 
		}
	      else
		{
		  hJpsi[k][i][j]->SetName(Form("%s_Jpsi_%s_%s",name[i].Data(),variable[j].Data(),cent_Title[k]));
		  hJpsi[k][i][j]->SetTitle(Form("%s distribution j/psi (%s, pt1 > %1.1f GeV/c, %s%%)",mc_name[j],name[i].Data(),pt1_cut,cent_Name[k])); 
		}
	    }
	  hJpsiUS_mc->GetAxis(0)->SetRange(0,-1);
	  hJpsiUS_mc->GetAxis(4)->SetRange(0,-1);
	  hJpsiUS_mc->GetAxis(5)->SetRange(0,-1);
	}
    }
  hJpsiUS_mc->GetAxis(7)->SetRange(0,-1);

  // matched Jpsi for response
  THnSparseF *hnJpsiMatch = (THnSparseF*)f->Get("mhJpsiMatch_di_mu");
  hnJpsiMatch->GetAxis(2)->SetRangeUser(low_mass+0.001, high_mass-0.001);
  TH2F *hJpsiPtMatch[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      if(nCentBins>1) hnJpsiMatch->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
      hJpsiPtMatch[k] = (TH2F*)hnJpsiMatch->Projection(1,3);
      hJpsiPtMatch[k]->Sumw2();
      if(nCentBins==1)
	hJpsiPtMatch[k]->SetName(Form("JpsiPt_TrueVsReco"));
      else
	hJpsiPtMatch[k]->SetName(Form("JpsiPt_TrueVsReco_%s",cent_Title[k]));
    }

  printf("+++ Weight jpsi pt distribution +++\n");
  // weight input spectrum
  TFile *fWeight = TFile::Open("Rootfiles/HP2015/HTppSpectrum.root","read");
  TF1 *func = (TF1*)fWeight->Get("function");
  
  TH1F *hJpsiPt[nCentBins][3][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  hJpsiPt[k][i][0] = (TH1F*)hJpsi[k][i][1]->Clone(Form("%s_FlatPt",hJpsi[k][i][1]->GetName()));
	}
    }

  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  if(i<2)
	    {
	      hJpsiPt[k][i][1] = (TH1F*)hJpsi[k][i][1]->Clone(Form("%s_WeightPt",hJpsi[k][i][1]->GetName()));
	      for(int bin=1; bin<=hJpsiPt[k][i][1]->GetNbinsX(); bin++)
		{
		  double weight = func->Eval(hJpsiPt[k][i][1]->GetBinCenter(bin)/mean_pt);
		  hJpsiPt[k][i][1]->SetBinContent(bin,hJpsiPt[k][i][1]->GetBinContent(bin)*weight);
		  hJpsiPt[k][i][1]->SetBinError(bin,hJpsiPt[k][i][1]->GetBinError(bin)*weight);
		}
	    }
	  else
	    {
	      for(int biny=1; biny<=hJpsiPtMatch[k]->GetNbinsY(); biny++)
		{
		  double weight = func->Eval(hJpsiPtMatch[k]->GetYaxis()->GetBinCenter(biny)/mean_pt);
		  for(int binx=1; binx<=hJpsiPtMatch[k]->GetNbinsX(); binx++)
		    {
		      hJpsiPtMatch[k]->SetBinContent(binx,biny,weight*hJpsiPtMatch[k]->GetBinContent(binx,biny));
		      hJpsiPtMatch[k]->SetBinError(binx,biny,weight*hJpsiPtMatch[k]->GetBinError(binx,biny));
		    }
		}
	      hJpsiPt[k][i][1] = (TH1F*)hJpsiPtMatch[k]->ProjectionX(Form("%s_WeightPt",hJpsi[k][i][1]->GetName()));
	    }
	}
    }

  printf("+++ Rebin jpsi distribution +++\n");
  // rebin input spectrum
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  TH1F *hJpsiPtRebin[nCentBins][3][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hJpsiPtRebin[k][i][j] = (TH1F*)hJpsiPt[k][i][j]->Rebin(nbins,Form("%s_rebin",hJpsiPt[k][i][j]->GetName()),xbins);
	    }
	}
    }

  // j/psi efficiency
  TH1F *hJpsiEff[nCentBins][2][2], *hJpsiEffRebin[nCentBins][2][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hJpsiEff[k][i][j] = (TH1F*)hJpsiPt[k][i+1][j]->Clone(Form("%s_Eff",hJpsiPt[k][i+1][j]->GetName()));
	      hJpsiEff[k][i][j]->Divide(hJpsiPt[k][0][j]);
	      
	      hJpsiEffRebin[k][i][j] = (TH1F*)hJpsiPtRebin[k][i+1][j]->Clone(Form("%s_Eff_rebin",hJpsiPt[k][i+1][j]->GetName()));
	      hJpsiEffRebin[k][i][j]->Divide(hJpsiPtRebin[k][0][j]);
	    }
	}
    }

  printf("+++ Single track resolution +++\n");
  // single track efficiency & resolution
  THnSparseF *hnTrkPtRes = (THnSparseF*)f->Get("mhpTrkPtRes_di_mu");
  TH2F *hResVsRecoPt[nCentBins];
  TH2F *hResVsTruePt[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      if(nCentBins>1) hnTrkPtRes->GetAxis(3)->SetRange(centBins_low[k],centBins_high[k]);
      hResVsRecoPt[k] = (TH2F*)hnTrkPtRes->Projection(0,1);
      if(nCentBins==1) hResVsRecoPt[k]->SetName(Form("PrimTrkRes_vs_RecoPt"));
      else             hResVsRecoPt[k]->SetName(Form("PrimTrkRes_vs_RecoPt_%s",cent_Title[k]));

      hResVsTruePt[k] = (TH2F*)hnTrkPtRes->Projection(0,2);
      if(nCentBins==1) hResVsTruePt[k]->SetName(Form("PrimTrkRes_vs_TruePt"));
      else             hResVsTruePt[k]->SetName(Form("PrimTrkRes_vs_TruePt_%s",cent_Title[k]));

      hnTrkPtRes->GetAxis(3)->SetRange(0,-1);
    }

  printf("+++ Single track efficiency +++\n");
  THnSparseF *hMcTrkInfo[5], *hRcTrkInfo[4];
  TH1F *hMcTrkPt[nCentBins][5][3], *hRcTrkPt[nCentBins][4][3];
  TH1F *hMcTrkPtEff[nCentBins][5][3], *hRcTrkPtEff[nCentBins][4][3];
  TH3F *hMcTrkPtEtaPhi[nCentBins][5][3], *hRcTrkPtEtaPhi[nCentBins][4][3];
  TH3F *hMcTrkPtEtaPhiEff[nCentBins][5][3], *hRcTrkPtEtaPhiEff[nCentBins][4][3];
  for(int i=0; i<5; i++)
    {
      hMcTrkInfo[i] = (THnSparseF*)f->Get(Form("mhMcTrkInfo%s_di_mu",det_name[i]));
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
	      if(nCentBins==1) hMcTrkPt[k][i][j]->SetName(Form("McTrkPt%s%s",det_name[i],charge_name[j]));
	      else             hMcTrkPt[k][i][j]->SetName(Form("McTrkPt%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
	      hMcTrkPtEff[k][i][j] = (TH1F*)hMcTrkPt[k][i][j]->Clone(Form("%s_Eff",hMcTrkPt[k][i][j]->GetName()));
	      hMcTrkPtEff[k][i][j]->Divide(hMcTrkPt[k][0][j]);

	      TH3F *h3 = (TH3F*)hMcTrkInfo[i]->Projection(0,1,2);
	      h3->SetName(Form("McTrkPtEtaPhi%d%d%d",i,j,k));
	      TString name_tmp = Form("McTrkPtEtaPhi%s%s",det_name[i],charge_name[j]);
	      if(nCentBins>1) name_tmp = Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]);
	      hMcTrkPtEtaPhi[k][i][j] = new TH3F(name_tmp.Data(),h3->GetTitle(),200,0,20,40,-2,2,360,pi/12,2*pi+pi/12);
	      for(int binx=1; binx<=h3->GetNbinsX(); binx++)
		{
		  for(int biny=1; biny<=h3->GetNbinsY(); biny++)
		    {
		      for(int binz=1; binz<=h3->GetNbinsZ(); binz++)
			{
			  int new_binz = binz-15;
			  if(binz<=15) new_binz = binz + 360 - 15;
			  hMcTrkPtEtaPhi[k][i][j]->SetBinContent(binx,biny,new_binz,h3->GetBinContent(binx,biny,binz));
			}
		    }
		}
	      hMcTrkPtEtaPhi[k][i][j]->Sumw2();
	      hMcTrkPtEtaPhiEff[k][i][j] = (TH3F*)hMcTrkPtEtaPhi[k][i][j]->Clone(Form("%s_Eff",hMcTrkPtEtaPhi[k][i][j]->GetName()));
	      hMcTrkPtEtaPhiEff[k][i][j]->Divide(hMcTrkPtEtaPhi[k][0][j]);
	      hMcTrkInfo[i]->GetAxis(3)->SetRange(0,-1);
	    }
	}
    }

  for(int i=0; i<4; i++)
    {
      hRcTrkInfo[i] = (THnSparseF*)f->Get(Form("mhRcTrkInfo%s_di_mu",det_name[i+1]));
      for(int k=0; k<nCentBins; k++)
	{
	  if(nCentBins>1) hRcTrkInfo[i]->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
	  for(int j=0; j<3; j++)
	    {
	      if(j==0) hRcTrkInfo[i]->GetAxis(3)->SetRange(1,3);
	      if(j==1) hRcTrkInfo[i]->GetAxis(3)->SetRange(3,3);
	      if(j==2) hRcTrkInfo[i]->GetAxis(3)->SetRange(1,1);
	      hRcTrkPt[k][i][j] = (TH1F*)hRcTrkInfo[i]->Projection(0);
	      hRcTrkPt[k][i][j]->Sumw2();
	      if(nCentBins==1) hRcTrkPt[k][i][j]->SetName(Form("RcTrkPt%s%s",det_name[i+1],charge_name[j]));
	      else             hRcTrkPt[k][i][j]->SetName(Form("RcTrkPt%s%s_%s",det_name[i+1],charge_name[j],cent_Title[k]));
	      hRcTrkPtEff[k][i][j] = (TH1F*)hRcTrkPt[k][i][j]->Clone(Form("%s_Eff",hRcTrkPt[k][i][j]->GetName()));
	      hRcTrkPtEff[k][i][j]->Divide(hRcTrkPt[k][0][j]);

	      TH3F *h3 = (TH3F*)hRcTrkInfo[i]->Projection(0,1,2);
	      h3->SetName(Form("McTrkPtEtaPhi%d%d%d",i,j,k));
	      TString name_tmp = Form("RcTrkPtEtaPhi%s%s",det_name[i+1],charge_name[j]);
	      if(nCentBins>1) name_tmp = Form("RcTrkPtEtaPhi%s%s_%s",det_name[i+1],charge_name[j],cent_Title[k]);
	      hRcTrkPtEtaPhi[k][i][j] = new TH3F(name_tmp.Data(),h3->GetTitle(),200,0,20,40,-2,2,360,pi/12,2*pi+pi/12);
	      for(int binx=1; binx<=h3->GetNbinsX(); binx++)
		{
		  for(int biny=1; biny<=h3->GetNbinsY(); biny++)
		    {
		      for(int binz=1; binz<=h3->GetNbinsZ(); binz++)
			{
			  int new_binz = binz-15;
			  if(binz<=15) new_binz = binz + 360 - 15;
			  hRcTrkPtEtaPhi[k][i][j]->SetBinContent(binx,biny,new_binz,h3->GetBinContent(binx,biny,binz));
			}
		    }
		}
	      hRcTrkPtEtaPhi[k][i][j]->Sumw2();
	      hRcTrkPtEtaPhiEff[k][i][j] = (TH3F*)hRcTrkPtEtaPhi[k][i][j]->Clone(Form("%s_Eff",hRcTrkPtEtaPhi[k][i][j]->GetName()));
	      hRcTrkPtEtaPhiEff[k][i][j]->Divide(hRcTrkPtEtaPhi[k][0][j]);
	      hRcTrkInfo[i]->GetAxis(3)->SetRange(0,-1);
	    }
	}
    }
  
  if(save)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outName.Data()),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  // jpsi
	  hJpsiPtMatch[k]->Write();
	  for(int i=0; i<3; i++)
	    {
	      hJpsiMassVsPt[k][i]->Write();
	      for(int j=0; j<4; j++)
		{
		  hJpsi[k][i][j]->Write();
		}
	    }
	  
	  for(int i=0; i<3; i++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hJpsiPt[k][i][j]->Write();
		  if(i>0) hJpsiEff[k][i-1][j]->Write();
		  hJpsiPtRebin[k][i][j]->Write();
		  if(i>0) hJpsiEffRebin[k][i-1][j]->Write();
		}
	    }

	  // single tracks
	  hResVsRecoPt[k]->Write();
	  hResVsTruePt[k]->Write();
	  for(int i=0; i<5; i++)
	    {
	      for(int j=0; j<3; j++)
		{
		  hMcTrkPt[k][i][j]->Write();
		  hMcTrkPtEff[k][i][j]->Write();
		  hMcTrkPtEtaPhi[k][i][j]->Write();
		  hMcTrkPtEtaPhiEff[k][i][j]->Write();
		}
	    }

	  for(int i=0; i<4; i++)
	    {
	      for(int j=0; j<3; j++)
		{
		  hRcTrkPt[k][i][j]->Write();
		  hRcTrkPtEff[k][i][j]->Write();
		  hRcTrkPtEtaPhi[k][i][j]->Write();
		  hRcTrkPtEtaPhiEff[k][i][j]->Write();
		}
	    }
	}
    }
}

//-------------------------------------------------------
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter)
{
  float betax = parent.Px()/parent.E();
  float betay = parent.Py()/parent.E();
  float betaz = parent.Pz()/parent.E();
  daughter.Boost(betax,betay,betaz);
  return daughter;
}

//-------------------------------------------------------
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass) 
{
  Double_t e = parent.M()/2.;
  Double_t p = sqrt(e*e-dmass*dmass);
  Double_t costheta = myRandom->Uniform(-1.,1.);
  Double_t phi = myRandom->Uniform(0,TMath::Pi()*2);
  Double_t pz = p*costheta;
  Double_t px = p*sqrt(1.-costheta*costheta)*cos(phi);
  Double_t py = p*sqrt(1.-costheta*costheta)*sin(phi);
  TLorentzVector daughter(px,py,pz,e);
  return myBoost(parent,daughter);
}


//================================================
void 3DTrkEff(TString inName, TString outPDFName, const bool savePlot = 0)
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
  TH3F *hMcTrkPtEtaPhi[nCentBins][2][3];
  TH1F *hMcTrkPt[nCentBins][2][3];
  TH1F *hMcTrkPtEff[nCentBins][2][3];
  TH2F *hMcTrkEtaPhi[nCentBins][2][3];
  TH2F *hMcTrkEtaPhiEff[nCentBins][2][3];
  TH1F *hMcTrkEta[nCentBins][2][3];
  TH1F *hMcTrkEtaEff[nCentBins][2][3];
  TH1F *hMcTrkPhi[nCentBins][2][3];
  TH1F *hMcTrkPhiEff[nCentBins][2][3];
  double eta_bounds[3] = {-0.5,0,0.5};
  int style[3] = {20,24,24};
  TString legName_charge[3] = {"#mu^{+}+#mu^{-}","#mu^{+}","#mu^{-}"};
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<2; i++)
	    {
	      if(nCentBins==1) hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s",det_name[i],charge_name[j]));
	      else             hMcTrkPtEtaPhi[k][i][j] = (TH3F*)fin->Get(Form("McTrkPtEtaPhi%s%s_%s",det_name[i],charge_name[j],cent_Title[k]));
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
      for(int j=0; j<3; j++)
	{
	  hMcTrkPtEff[k][1][j]->GetXaxis()->SetRangeUser(0,10);
	  hMcTrkPtEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.2);
	  hMcTrkPtEff[k][1][j]->SetTitle(";p_{T,mc} (GeV/c);Efficiency");
	  if(j==0) hMcTrkPtEff[k][1][j]->Draw("");
	  else     hMcTrkPtEff[k][1][j]->Draw("sames");
	  leg->AddEntry(hMcTrkPtEff[k][1][j],legName_charge[j],"P");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (|#eta_{mc}| < 0.5)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPtTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      c = new TCanvas(Form("McTrkEtaEff_Cent%d",k),Form("McTrkEtaEff_Cent%d",k),800,600);
      for(int j=0; j<3; j++)
	{
	  hMcTrkEtaEff[k][1][j]->GetXaxis()->SetRangeUser(-1,1);
	  hMcTrkEtaEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.2);
	  hMcTrkEtaEff[k][1][j]->SetTitle(";#eta_{mc};Efficiency");
	  if(j==0) hMcTrkEtaEff[k][1][j]->Draw("");
	  else     hMcTrkEtaEff[k][1][j]->Draw("sames");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText("TPC tracking efficiency (p_{T,mc} > 1.5 GeV/c)");
      t1->Draw();
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}

      c = new TCanvas(Form("McTrkPhiEff_Cent%d",k),Form("McTrkPhiEff_Cent%d",k),800,600);
      for(int j=0; j<3; j++)
      	{
      	  hMcTrkPhiEff[k][1][j]->GetYaxis()->SetRangeUser(0,1.2);
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPhiTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkPhiTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
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
      PaintCanvasToPDF(c,pdf);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaPhiTpcEff_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[k]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sMcTrkEtaPhiTpcEff_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[k]));
	}
    }

  // in (phi,eta) bins
  const int neta = 10*2;
  const int nphi = 60;
  const int nCanvas = nphi*neta/2/15;
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
	  for(int i=0; i<2; i++)
	    {
	      int nbinsy = hMcTrkPtEtaPhi[k][i][j]->GetNbinsY();
	      int nbinsz = hMcTrkPtEtaPhi[k][i][j]->GetNbinsZ();
	      hMcTrkPtEtaPhi[k][i][j]->Rebin3D(10,nbinsy/neta,nbinsz/nphi);
	    }
	  hMcTrkPtEtaPhiEff[k][j] = (TH3F*)hMcTrkPtEtaPhi[k][1][j]->Clone(Form("%s_Eff",hMcTrkPtEtaPhi[k][1][j]->GetName()));
	  hMcTrkPtEtaPhiEff[k][j]->Divide(hMcTrkPtEtaPhi[k][0][j]);
	}

      for(int ieta=1+neta/4; ieta<=neta-neta/4; ieta++)
	{
	  for(int iphi=1; iphi<=nphi; iphi++)
	    {
	      int index = (ieta-1-neta/4)*nphi + iphi - 1;
	      cTpcEff[index/15]->cd(index%15+1);
	      for(int j=1; j<3; j++)
		{
		  TH1F *htmp = (TH1F*)hMcTrkPtEtaPhiEff[k][j]->ProjectionX(Form("%s_%d_%d",hMcTrkPtEtaPhiEff[k][j]->GetName(),ieta,iphi),ieta,ieta,iphi,iphi);
		  htmp->GetXaxis()->SetRangeUser(0,10);
		  htmp->SetMarkerStyle(20+(j-1)*4);
		  htmp->SetTitle(";p_{T,mc} (GeV/c);Tpc Efficiency");
		  if(j==1) htmp->Draw();
		  else     htmp->Draw("sames");
		}
	    }
	}
    }
  

  pdf->On();
  pdf->Close();
}

