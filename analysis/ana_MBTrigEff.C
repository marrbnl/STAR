const int year = YEAR;
TFile *f;

//================================================
void ana_TrigEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //trigBiasFromPythia();
}


//================================================
void trigBiasFromPythia(const int savePlot = 1, const int saveHisto = 1)
{
  TLegend *leg = new TLegend(0.4, 0.5, 0.6, 0.75);
  leg->SetFillColor(10);
  leg->SetFillStyle(10);
  leg->SetLineStyle(4000);
  leg->SetLineColor(10);
  leg->SetLineWidth(0.);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
	
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  //c1->DrawFrame(0,0.45,8,1);

  const int nPoints = nPtBins -1;
  double arrayOfLowEdges[nPoints+1];
  for(int i=0; i<nPoints; i++)
    arrayOfLowEdges[i] = ptBins_low[i+1];
  arrayOfLowEdges[nPoints] = ptBins_high[nPoints];

  TFile *fin = new TFile("Rootfiles/CharmedSample.root");
  TH2F *hPtEta_Mc = fin->Get("PtEta_Mc");
  hPtEta_Mc->ProjectionX("Mc",11,30,"");
  Mc->Rebin(nPoints,"hPtMc",arrayOfLowEdges);
	
  TH2F *hPtEta_Bbc = fin->Get("PtEta_Bbc");
  hPtEta_Bbc->ProjectionX("PtBbc",11,30,"");
  PtBbc->Rebin(nPoints,"hPtBbc",arrayOfLowEdges);
	
  TH2F *hPtEta_Vpd = fin->Get("PtEta_Vpd");
  hPtEta_Vpd->ProjectionX("PtVpd",11,30,"");
  PtVpd->Rebin(nPoints,"hPtVpd",arrayOfLowEdges);
    
  //************************ BBC and VDP efficiency ********************************
	
  TGraphAsymmErrors *gEffPtBbc = new TGraphAsymmErrors(hPtBbc,hPtMc,"w");
  TGraphAsymmErrors *gEffPtVpd = new TGraphAsymmErrors(hPtVpd,hPtMc,"w");
  gEffPtBbc->SetMarkerStyle(25);
  gEffPtBbc->SetMarkerColor(kGreen+2);
  gEffPtBbc->SetLineColor(kGreen+2);
  gEffPtVpd->SetMarkerStyle(25);
  gEffPtVpd->SetMarkerColor(2);
  gEffPtVpd->SetLineColor(2);
  gEffPtBbc->SetTitle(";p_{T,D} (GeV/c)");
  gEffPtBbc->GetYaxis()->SetRangeUser(0.45,1);
  gEffPtBbc->Draw("AP");
  gEffPtVpd->Draw("Psame");
    
  //**************************** BBC and VPD efficiency with vertexing **************************
	
  TH2F *hPtEta_Vtx = fin->Get("PtEta_Vtx");
  hPtEta_Vtx->ProjectionX("PtVtx",11,30,"");
  PtVtx->Rebin(nPoints,"hPtVtx",arrayOfLowEdges);
  TGraphAsymmErrors *gEffPtVtx = new TGraphAsymmErrors(hPtVtx,hPtMc,"w");
  gEffPtVtx->SetMarkerStyle(24);
  gEffPtVtx->SetMarkerColor(1);
  gEffPtVtx->SetLineColor(1);
  gEffPtVtx->Draw("Psame");
    
  TH2F *hPtEta_VpdVtx = fin->Get("PtEta_VpdVtx");
  hPtEta_VpdVtx->ProjectionX("PtVpdVtx",11,30,"");
  PtVpdVtx->Rebin(nPoints,"hPtVpdVtx",arrayOfLowEdges);
  TGraphAsymmErrors *gEffPtVpdVtx = new TGraphAsymmErrors(hPtVpdVtx,hPtMc,"w");

		
  gEffPtVpdVtx->SetMarkerStyle(20);
  gEffPtVpdVtx->SetMarkerColor(1);
  gEffPtVpdVtx->SetLineColor(1);
  gEffPtVpdVtx->Draw("Psame");

  double gx[100], gy[100], gy1[100], gy2[100];
  int n = 0;
  for (int i=0; i<gEffPtVpdVtx->GetN(); i++) {
    double x,y,y1,y2;
    gEffPtVpdVtx->GetPoint(i,x,y);
    y1 = gEffPtVpdVtx->GetErrorYhigh(i);
    y2 = gEffPtVpdVtx->GetErrorYlow(i);
    gx[i] = x;
    if(y>0){
      gx[i] = x;
      gy[i] = 0.388234/y;
      gy1[i] = gy[i]*y1/y;
      gy2[i] = gy[i]*y2/y;
      n++;
    }
  }
  TGraphAsymmErrors *gTrgBias = new TGraphAsymmErrors(n, gx, gy, 0, 0, gy2, gy1);
  gTrgBias->SetMarkerStyle(29);
  gTrgBias->SetMarkerColor(4);
  gTrgBias->SetLineColor(4);
  gTrgBias->SetMarkerSize(2);
  gTrgBias->Draw("Psame");

  TH1F *hTrgBias = new TH1F("VPDMB_TrigBias","VPDMB_TrigBias",nPoints,arrayOfLowEdges);
  double x,y;
  for(int i=0; i<nPoints; i++)
    {
      gTrgBias->GetPoint(i,x,y);
      hTrgBias->SetBinContent(i+1,y);
      hTrgBias->SetBinError(i+1,gy1[i]);
    }
    
  leg->AddEntry(gEffPtBbc,"BBC Coincidence","pl");
  leg->AddEntry(gEffPtVpd,"VPD Coincidence","pl");
  leg->AddEntry(gEffPtVtx,"Vertex Accepted","pl");
  leg->AddEntry(gEffPtVpdVtx,"VPD Coincidence & Vertex Accepted","pl");
  leg->AddEntry(gTrgBias,"Trigger Bias","pl");
  leg->Draw();	

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VPDMB_TrigBias.pdf",run_type));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VPDMB_TrigBias.png",run_type));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
      hTrgBias->Write("VPDMB_TrigBias",TObject::kOverwrite);
      fout->Close();
    }
}

