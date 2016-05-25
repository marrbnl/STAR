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

  checkEff();
}

//================================================
void checkEff(const int savePlot = 1)
{
  TList *list = new TList;

  // eta dependence
  const int nEtaCuts = 3;
  const double eta_cut[nEtaCuts] = {1.0,0.8,0.5};
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

