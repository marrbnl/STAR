TFile *f;
const int year = YEAR;

//================================================
void ana_TrkEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Run14.AuAu200.Jpsi.Embed.root"),"read");
    }

  //checkEff();
  TpcAccLoss();
}

//================================================
void TpcAccLoss(const int savePlot = 0, const int saveHisto = 1)
{
  const int nHistos = 2;
  const char *name[nHistos] = {"McMuon","DataPion"};
  const char *title[nHistos] = {"recontructed MC track in embedding","pion candidates in data"};
  THnSparseF *hnTrk[nHistos];
  hnTrk[0] = (THnSparseF*)f->Get("hTrkEtaPhi_MCreco_di_mu");
  hnTrk[1] = (THnSparseF*)f->Get("hTrkEtaPhi_DataPion_di_mu");

  
  // === eta vs phi
  TH2F *hTrkEtaPhi[nHistos];
  for(int j=0; j<nHistos; j++)
    {
      hTrkEtaPhi[j] = (TH2F*)hnTrk[j]->Projection(2,1);
      hTrkEtaPhi[j]->RebinY(2);
      hTrkEtaPhi[j]->SetTitle(";#eta;#varphi");
      c = draw2D(hTrkEtaPhi[j], Form("#varphi vs #eta of %s",title[j]), 0.04, false);
      if(savePlot)  
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TrkEff/TrkPhivsEta_%s.pdf",run_type,name[j]));
    }

  const double eta_cuts[3] = {-1.0,0.2,1.0};
  const double pt_cuts[6] = {1.0,1.5,2.0,3.0,5.0,20.0};
  TH1F *hTrkPhi[2][nHistos][5];
  TCanvas *cPhi[2];
  for(int i=0; i<2; i++)
    {
      cPhi[i] = new TCanvas(Form("track_phi_%d",i),Form("track_phi_%d",i),1000,600);
      cPhi[i]->Divide(3,2);
      TLegend *leg = new TLegend(0.1,0.6,0.7,0.83);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(Form("%1.1f < #eta < %1.1f",eta_cuts[i], eta_cuts[i+1]));
      for(int j=0; j<nHistos; j++)
	{
	  hnTrk[j]->GetAxis(1)->SetRangeUser(eta_cuts[i]+0.01, eta_cuts[i+1]-0.01);
	  for(int k=0; k<5; k++)
	    {
	      hnTrk[j]->GetAxis(0)->SetRangeUser(pt_cuts[k]+0.01, pt_cuts[k+1]-0.01);
	      hTrkPhi[i][j][k] = (TH1F*)hnTrk[j]->Projection(2);
	      hTrkPhi[i][j][k]->SetName(Form("hTrkPhi_%s_%d_%d",name[j],i,k));
	      hTrkPhi[i][j][k]->SetMarkerStyle(21);
	      hTrkPhi[i][j][k]->SetMarkerColor(2-j);
	      hTrkPhi[i][j][k]->SetLineColor(2-j);
	      hTrkPhi[i][j][k]->Rebin(10);
	      hTrkPhi[i][j][k]->Scale(1./hTrkPhi[i][j][k]->Integral());
	      hTrkPhi[i][j][k]->GetYaxis()->SetRangeUser(0,0.04);
	      hTrkPhi[i][j][k]->SetTitle(";#varphi");
	      cPhi[i]->cd(k+1);
	      if(j==0) 
		{
		  hTrkPhi[i][j][k]->Draw("P");
		  t = GetTitleText(Form("#varphi distribution (%1.1f < p_{T} < %1.1f)",pt_cuts[k],pt_cuts[k+1]),0.05);
		  t->Draw();
		}
	      else     hTrkPhi[i][j][k]->Draw("samesHIST");
	      if(k==0)
		{
		  if(j==0) leg->AddEntry(hTrkPhi[i][j][k],"Reco MC #mu tracks","PE");
		  if(j==1) leg->AddEntry(hTrkPhi[i][j][k],"#pi candidates in data","L");
		}
	      hnTrk[j]->GetAxis(0)->SetRange(0,-1);
	    }
	  hnTrk[j]->GetAxis(1)->SetRange(0,-1);
	}
      cPhi[i]->cd(6);
      leg->Draw();
    }

  TH2F *hTpcCorr = new TH2F("hTpcCorr","Correction factor of TPC inefficiency for #eta < 0.2;p_{T} (GeV/c);#varphi",5,pt_cuts,36,0,2*pi);
  TH1F *hTrkPhiCorr[5];
  for(int k=0; k<5; k++)
    {
      hTrkPhiCorr[k] = (TH1F*)hTrkPhi[0][0][k]->Clone(Form("%s_corr",hTrkPhi[0][0][k]->GetName()));
      for(int bin=1; bin<=hTrkPhiCorr[k]->GetNbinsX(); bin++)
	{
	  if(bin<=31 || bin==36) hTpcCorr->SetBinContent(k+1,bin,1);
	  else
	    {
	      double eff = hTrkPhi[0][1][k]->GetBinContent(bin)/hTrkPhi[0][0][k]->GetBinContent(bin);
	      hTrkPhiCorr[k]->SetBinContent(bin,hTrkPhiCorr[k]->GetBinContent(bin)*eff);
	      hTpcCorr->SetBinContent(k+1,bin,eff);
	      if(k==4) hTpcCorr->SetBinContent(k+1,bin,hTpcCorr->GetBinContent(4,bin));
	    }
	}
      cPhi[0]->cd(k+1);
      hTrkPhiCorr[k]->SetMarkerStyle(25);
      hTrkPhiCorr[k]->SetMarkerColor(4);
      hTrkPhiCorr[k]->Draw("samesP");
    }
  c = draw2D(hTpcCorr,"",0.04,false);
  if(savePlot)  
    {
      cPhi[0]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TrkEff/CompTrkPhi_NegEta.pdf",run_type));
      cPhi[1]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TrkEff/CompTrkPhi_PosEta.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TrkEff/CorrFactor.pdf",run_type));
    }

  if(saveHisto)
    {
      TFile *fout =  TFile::Open("Rootfiles/Run14.AuAu200.Input.root","update");
      hTpcCorr->Write("",TObject::kOverwrite);
    }
}

//================================================
void checkEff(const int savePlot = 1)
{
  TList *list = new TList;

  // eta dependence
  const int nEtaCuts = 3;
  const double eta_cut[nEtaCuts] = {1.5,0.8,0.5};
  THnSparseF *hnTrkInfo[2];
  hnTrkInfo[0] = (THnSparseF*)f->Get("mhMcTrkInfo_di_mu");
  hnTrkInfo[1] = (THnSparseF*)f->Get("mhMcTrkInfoTpc_di_mu");

  TH2F *hMcTrkEtaVsPt = (TH2F*)hnTrkInfo[0]->Projection(1,0);
  hMcTrkEtaVsPt->GetXaxis()->SetRangeUser(0,20);
  c = draw2D(hMcTrkEtaVsPt,"#eta vs p_{T} for input MC tracks");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TrkEff/McTrkEtaVsPt.pdf",run_type));

  TH1F *hTrkPtEtaCut[2][nEtaCuts];
  TH1F *hTrkPtEffEtaCut[nEtaCuts];
  for(int i=0; i<2; i++)
    {
      hnTrkInfo[i]->Sumw2();
      hnTrkInfo[i]->GetAxis(4)->SetRange(13,16);
      for(int j=0; j<nEtaCuts; j++)
	{
	  hnTrkInfo[i]->GetAxis(1)->SetRangeUser(-1*eta_cut[j]+1e-4, eta_cut[j]-1e-4);
	  hTrkPtEtaCut[i][j] = (TH1F*)hnTrkInfo[i]->Projection(0);
	  hTrkPtEtaCut[i][j]->SetName(Form("hTrkPt_%d_etaCut%d",i,j));
	}
    }
  TString legName[nEtaCuts];
  for(int j=0; j<nEtaCuts; j++)
    {
      hTrkPtEffEtaCut[j] = (TH1F*)hTrkPtEtaCut[1][j]->Clone(Form("hTrkPtEff_etaCut%d",j));
      hTrkPtEffEtaCut[j]->Divide(hTrkPtEtaCut[0][j]);
      list->Add(hTrkPtEffEtaCut[j]);
      legName[j] = Form("|#eta_{MC}| < %1.1f",eta_cut[j]);
    }
  c = drawHistos(list,"TrkPtEff_EtaCut","TPC tracking efficiency for single muons;p_{T} (GeV/c);Efficiency",true,0,20,true,0,1,kFALSE,kTRUE,legName,kTRUE,"Au+Au 0-20%",0.5,0.7,0.25,0.5,kTRUE);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_TrkEff/McTrkEff_EtaCuts.pdf",run_type));


  
}

