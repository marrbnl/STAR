TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "muon.FitDz.";
const Bool_t iPico = 1;
const int year = 2014;
const double low_mass = 2.95;
const double high_mass = 3.25;
TString run_cfg_name;
const char *outFileName = "Run14.AuAu200.jpsi.FitDz.root";

//================================================
void ana_DeltaZ()
{						
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.18); 

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  //compareFitDz();
  //MuonCandidates();
  //Compare2GausVs3Gaus();
  //JpsiMuon();
  HadronEmbed();
}

//================================================
void compareFitDz(const Int_t savePlot = 1)
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s",outFileName),"read");
  // embedding hadron
  const char *partName[4] = {"kaonP","kaonM","pionP","pionM"};
  TH1F *hEmbed[4][2];
  for(int j=0; j<4; j++)
    {
      for(int k=0; k<2; k++)
	{
	  hEmbed[j][k] = (TH1F*)fin->Get(Form("Embed_%s_dz_sigma%d",partName[j],k));
	  hEmbed[j][k]->SetMarkerStyle(20+j);
	  hEmbed[j][k]->SetMarkerColor(color[k+2]);
	  hEmbed[j][k]->SetLineColor(color[k+2]);
	}
    }

  // embedding muon
  TF1 *fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)",1,20);
  fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.444217);

  
  // data Jpsi muon
  TH1F *hDataJpsiMuon = (TH1F*)fin->Get("Data_JpsiMuon_sigma");
  hDataJpsiMuon->SetMarkerStyle(21);

  // data muon
  const int type = 1;
  const char *type_name[2] = {"2Gaus","3Gaus"};
  TH1F *hDataMuon[2+type];
  for(int i=0; i<2+type; i++)
    {
      hDataMuon[i] = (TH1F*)fin->Get(Form("Data_Muon_Sigma_%s_%d",type_name[type],i));
      hDataMuon[i]->SetMarkerStyle(25);
      hDataMuon[i]->SetMarkerColor(color[i]);
      hDataMuon[i]->SetLineColor(color[i]);
      hDataMuon[i]->GetYaxis()->SetRangeUser(0,60);
      hDataMuon[i]->GetXaxis()->SetRangeUser(0.5,20);
    }
  if(type==1)
    {
      for(int bin=hDataMuon[2]->FindFixBin(7.5); bin<=hDataMuon[2]->GetNbinsX(); bin++)
	{
	  hDataMuon[2]->SetBinContent(bin,-1);
	  hDataMuon[2]->SetBinError(bin,0);
	}
    }

  c = draw1D(hDataMuon[0],"Comparing #sigma of #Deltaz distributions");
  hDataMuon[1]->Draw("sames");
  if(type==1) hDataMuon[2]->Draw("sames");
  hDataJpsiMuon->Draw("sames");
  fResDzVsPt->Draw("same");
  for(int j=0; j<4; j++)
    {
      for(int k=0; k<2; k++)
	{
	  hEmbed[j][k]->Draw("same");
	}
    }

  double xmin = 0.3, xmax = 0.5, ymin = 0.65, ymax = 0.88;
  leg1 = new TLegend(xmin,ymin,xmax,ymax);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(hDataMuon[0],"Data muon","PL");
  leg1->AddEntry(hDataMuon[1],"Data muon","PL");
  if(type==1) leg1->AddEntry(hDataMuon[2],"Data muon","PL");
  leg1->AddEntry(hDataJpsiMuon,"Data J/psi muon","PL");
  leg1->Draw();

  leg2 = new TLegend(xmin+0.3,ymin,xmax+0.3,ymax);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(fResDzVsPt,"Embed J/psi muon","L");
  leg2->AddEntry(hEmbed[0][0],"Embed K^{+}","PL");
  leg2->AddEntry(hEmbed[1][0],"Embed K^{#font[122]{-}}","PL");
  leg2->AddEntry(hEmbed[2][0],"Embed #pi^{+}","PL");
  leg2->AddEntry(hEmbed[3][0],"Embed #pi^{#font[122]{-}}","PL");
  leg2->Draw();

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Compare_DzSigma_All.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Compare_DzSigma_All.pdf",run_type));
    }
}

//================================================
void HadronEmbed(const Int_t savePlot = 0, const Int_t saveHisto = 0)
{
  const int nPart = 6;
  const int geantId[nPart] = {8,9,11,12,14,15};
  const char *partName[nPart] = {"pionP","pionM","kaonP","kaonM","protonP","protonM"};
  const char *partTitle[nPart] = {"#pi^{+}","#pi^{#font[122]{-}}","K^{+}","K^{#font[122]{-}}","proton","anti-proton"};
  const char *partNameAll[nPart/2] = {"pion","kaon","pronton"};


  TFile *fin = TFile::Open("output/Run14.AuAu200.Hadron.embed.root","read");
  THnSparseF *hn = (THnSparseF*)fin->Get("mhHadronMatch");

  TList *list = new TList;

  // tracking efficiency
  TH2F *hMcTrkPtVsGeant[3];
  hMcTrkPtVsGeant[0] = (TH2F*)fin->Get("mhHadronMcPt");
  hMcTrkPtVsGeant[1] = (TH2F*)fin->Get("mhHadronMcPtTpc");
  hMcTrkPtVsGeant[2] = (TH2F*)fin->Get("mhHadronMcPtMtd");
  TH1F *hMcTrkPt[nPart][3];
  TString legName[3] = {"MC input","Reconstructed in TPC","Match to MTD"};
  for(int i=0; i<nPart; i++)
    {
      if(hMcTrkPtVsGeant[0]->GetBinContent(1,geantId[i]+1)<1) continue;
      list->Clear();
      for(int j=0; j<3; j++)
	{
	  hMcTrkPt[i][j] = (TH1F*)hMcTrkPtVsGeant[j]->ProjectionX(Form("hMcTrkPt_%s_%d",partName[i],j),geantId[i]+1,geantId[i]+1);
	  hMcTrkPt[i][j]->Rebin(5);
	  hMcTrkPt[i][j]->SetMaximum(1.2*hMcTrkPt[i][j]->GetMaximum());
	  hMcTrkPt[i][j]->SetMinimum(0.1);
	  hMcTrkPt[i][j]->Sumw2();

	  list->Add(hMcTrkPt[i][j]);
	}
      c = sysCompare(list,Form("McTrkEff_%s",partName[i]),Form("p_{T} distribution of embedded %s;p_{T}^{mc} (GeV/c);counts",partTitle[i]),Form("Efficiency of embedded %s;p_{T}^{mc} (GeV/c);Efficiency",partTitle[i]),kTRUE,0,13,kFALSE,0.1,10,kTRUE,0,1,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.5,0.7,0.25,0.5,kTRUE);
      c->cd(2);
      gPad->SetGridy();
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_McTrkPtEff.png",run_type,partName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_McTrkPtEff.pdf",run_type,partName[i]));
	}
    }

  // decay mode
  TH2F *hHadronDecay = (TH2F*)fin->Get("mhHadronDecay");
  for(int i=0; i<nPart; i++)
    {
      printf("=== %s\n",partName[i]);
      TH1F *h1 = (TH1F*)hHadronDecay->ProjectionX(Form("%s",partName[i]),geantId[i]+1,geantId[i]+1);
      for(int j=1; j<=h1->GetNbinsX(); j++)
	{
	  int entry = h1->GetBinContent(j);
	  if(entry<1) continue;
	  printf("Decay code %d with entries %1.0f (%1.2f%%)\n",j-1,entry,entry/h1->GetEntries()*100);
	}
    }

  TH2F *hDyVsPt[nPart];
  TH2F *hDzVsPt[nPart];
  TH2F *hDzVsDecay[nPart];
  TH1F *hAll[nPart];
  TH1F *hNoDecay[nPart];
  TH1F *hDecay[nPart];
  for(int i=0; i<nPart; i++)
    {
      hn->GetAxis(0)->SetRange(geantId[i]+1,geantId[i]+1);
      hDyVsPt[i] = (TH2F*)hn->Projection(3,1);
      hDyVsPt[i]->Sumw2();
      hDyVsPt[i]->Rebin2D(5,5);
      hDyVsPt[i]->SetName(Form("hDyVsPt_%s",partName[i]));
      hDyVsPt[i]->SetTitle(Form("#Deltay vs p_{T} for %s in embedding;p_{T} [GeV/c];#Deltay [cm]",partTitle[i]));
      hDzVsPt[i] = (TH2F*)hn->Projection(2,1);
      hDzVsPt[i]->Sumw2();
      hDzVsPt[i]->Rebin2D(5,5);
      hDzVsPt[i]->SetName(Form("hDzVsPt_%s",partName[i]));
      hDzVsPt[i]->SetTitle(Form("#Deltaz vs p_{T} for %s in embedding;p_{T} [GeV/c];#Deltaz [cm]",partTitle[i]));
      hDzVsDecay[i] = (TH2F*)hn->Projection(2,5);
      hDzVsDecay[i]->Sumw2();
      hDzVsDecay[i]->Rebin2D(5,5);
      hDzVsDecay[i]->SetName(Form("hDzVsDecay_%s",partName[i]));
      hDzVsDecay[i]->SetTitle(Form("#Deltaz vs decay radius for %s in embedding;r [cm];#Deltaz [cm]",partTitle[i]));
      hAll[i] = (TH1F*)hDzVsDecay[i]->ProjectionY(Form("hDz_All_%s",partName[i]),0,-1);
      hNoDecay[i] = (TH1F*)hDzVsDecay[i]->ProjectionY(Form("hDz_NoDecay_%s",partName[i]),0,0);
      hDecay[i] = (TH1F*)hDzVsDecay[i]->ProjectionY(Form("hDz_Decay_%s",partName[i]),1,hDzVsDecay[i]->GetXaxis()->FindFixBin(300));
      if(hDzVsPt[i]->GetEntries()<100) continue;
      hDzVsPt[i]->GetXaxis()->SetRangeUser(0,13);
      c = draw2D(hDzVsPt[i]);
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_DzVsPt.png",run_type,partName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_DzVsPt.pdf",run_type,partName[i]));
	}
      c = draw2D(hDzVsDecay[i]);
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_DzVsDecayR.png",run_type,partName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_DzVsDecayR.pdf",run_type,partName[i]));
	}
      hn->GetAxis(0)->SetRange(0,-1);
    }	

  // dz vs decay mode
  const int nDecay = 6;
  const int decay_mode[nDecay] = {0,5,13,16,20,40};
  const char *decay_process[nDecay] = {"No decay","Muonic decay","Elastic scattering","Nuclear Absorption","Inelastic scattering","Hadronic decay"};
  const char *var_title[3] = {"#Deltaz","#Deltay","#Deltatof"};
  const char *var_name[3] = {"dz","dy","dtof"};
  TH2F *hDVsMode[nPart][3];
  TString legName2[nDecay];
  printf("+++++ MTD match +++++\n");
  for(int i=0; i<nPart; i++)
    {
      hn->GetAxis(0)->SetRange(geantId[i]+1,geantId[i]+1);
      int all = 0;
      for(int j=0; j<3; j++)
	{
	  hDVsMode[i][j] = (TH2F*)hn->Projection(2+j,6);
	  hDVsMode[i][j]->SetName(Form("%s_%s_vs_mode",partName[i],var_name[j]));
	  int counter = 0;
	  list->Clear();
	  if(j==0)   all = hDVsMode[i][j]->GetEntries();
	  for(int k=0; k<nDecay; k++)
	    {
	      TH1F *htmp = (TH1F*)hDVsMode[i][j]->ProjectionY(Form("%s_%d",hDVsMode[i][j]->GetName(),decay_mode[k]),decay_mode[k]+1,decay_mode[k]+1);
	      if(htmp->GetEntries()<1) continue;
	      list->Add(htmp);
	      legName2[counter] = decay_process[k];
	      counter++;
	      if(j==0)
		{
		  printf("%s: decay code %d with entries %1.0f (%1.2f%%)\n",partName[i],decay_mode[k],htmp->GetEntries(),htmp->GetEntries()/all*100);
		}
	    }
	  if(j==0)       printf("%s: %1.0f\n",partName[i],all);
	  if(list->GetEntries()>0)
	    {
	      double min = -100, max = 100;
	      if(j==2) { min = -1; max = 5;}
	      c = drawHistos(list,Form("%s_%s",partName[i],var_name[j]),Form("%s distribution of %s",var_title[j],partTitle[i]),kTRUE,min,max,kFALSE,0.1,10,kFALSE,kTRUE,legName2,kTRUE,"",0.55,0.8,0.6,0.85,kFALSE);
	      TLine *line = GetLine(0,0,0,((TH1*)list->At(0))->GetMaximum(),1);
	      line->Draw();
	      if(savePlot) 
		{
		  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_%sDis.png",run_type,partName[i],var_name[j]));
		  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_%sDis.pdf",run_type,partName[i],var_name[j]));
		}
	    }

	  if(j>1) continue;
	  if(hDVsMode[i][0]->GetEntries()<1) continue;
	  TCanvas *c = new TCanvas(Form("Fit_%s_dz",partName[i]),Form("Fit_%s_dz",partName[i]),1100,500);
	  c->Divide(2,1);
	  for(int k=0; k<2; k++)
	    {
	      TH1F *htmp = (TH1F*)hDVsMode[i][j]->ProjectionY(Form("%s_%d",hDVsMode[i][j]->GetName(),decay_mode[k]),decay_mode[k]+1,decay_mode[k]+1);
	      htmp->Sumw2();
	      if(htmp->GetEntries()<1) continue;
	      double range = 60;
	      TF1 *func = new TF1(Form("func_%s_%d",partName[i],k),"gaus(0)+gaus(3)",-1*range,range);
	      func->SetParameter(1,0);
	      func->SetParameter(4,0);
	      func->SetParameter(2,10);
	      func->SetParameter(5,30);
	      htmp->Fit(func,"IR0");
	      htmp->GetYaxis()->SetNdivisions(505);
	      c->cd(k+1);
	      htmp->SetTitle(Form(";%s [cm]",var_title[j]));
	      htmp->GetXaxis()->SetRangeUser(-80,100);
	      htmp->SetMarkerStyle(20);
	      htmp->SetMarkerColor(1);
	      htmp->SetLineColor(1);
	      htmp->SetMaximum(1.5*htmp->GetMaximum());
	      htmp->Draw("P");
	      func->SetLineColor(2);
	      func->Draw("sames");
	      TF1 *func2 = new TF1(Form("func2_%s_%d",partName[i],k),"gaus",-1*range,range);
	      for(int ipar=0; ipar<3; ipar++) func2->SetParameter(ipar,func->GetParameter(ipar+3));
	      func2->SetLineColor(4);
	      func2->Draw("sames");
	      TPaveText *t = GetTitleText(Form("%s: %s",partTitle[i],decay_process[k]),0.06);
	      t->Draw();
	      TLine *line = GetLine(0,0,0,htmp->GetMaximum()*0.7,1);
	      line->Draw();
	    }
	  if(savePlot) 
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s.png",run_type,partName[i],var_name[j]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s.pdf",run_type,partName[i],var_name[j]));
	    }
	}
    }

  const int nTrkPtBin =6;
  const double trkPtBins[nTrkPtBin+1] = {1,1.5,2,3,5,8,10};
  double trkPtbins[nTrkPtBin+2];
  trkPtbins[0] = 0;
  for(int i=1; i<nTrkPtBin+2; i++) trkPtbins[i] = trkPtBins[i-1];
  TH1F *hMean[nPart][2][2];
  TH1F *hSigma[nPart][2][2];
  const int nCanvas = nTrkPtBin/6+1;
  TCanvas *canvas[nPart][2][nCanvas];
  TString legName1[2] = {"Narrow Gaussian","Wide Gaussian"};
  for(int i=0; i<nPart; i++)
    {
      TH2F *h2; 
      for(int j=0; j<2; j++)
	{
	  if(j==0) h2 = hDzVsPt[i];
	  else     h2 = hDyVsPt[i];
	  if(h2->GetEntries()<1) continue;
	  for(int k=0; k<2; k++)
	    {
	      hMean[i][j][k] = new TH1F(Form("Embed_%s_%s_mean%d",partName[i],var_name[j],k),Form("Mean of fitted %s distribution;p_{T} [GeV/c];mean [cm]",var_title[j]),nTrkPtBin+1,trkPtbins);
	      hSigma[i][j][k] = new TH1F(Form("Embed_%s_%s_sigma%d",partName[i],var_name[j],k),Form("Sigma of fitted %s distribution;p_{T} [GeV/c];#sigma [cm]",var_title[j]),nTrkPtBin+1,trkPtbins);
	    }
	  for(int k=0; k<nTrkPtBin; k++)
	    {
	      if(k%6==0) 
		{
		  canvas[i][j][k/6] = new TCanvas(Form("Fit_%s_%s_%d",var_name[j],partName[i],k/6),Form("Fit_%s_%s_%d",var_name[j],partName[i],k/6),1100,650);
		  canvas[i][j][k/6]->Divide(3,2);
		}

	      int low_bin  = h2->GetXaxis()->FindFixBin(trkPtBins[k]+1e-4);
	      int high_bin = h2->GetXaxis()->FindFixBin(trkPtBins[k+1]-1e-4);
      
	      TH1F *htmp = (TH1F*)h2->ProjectionY(Form("hTrk_%s_pt%1.1f-%1.1f_%s",var_name[j],trkPtBins[k],trkPtBins[k+1],partName[i]),low_bin,high_bin);
	      if(htmp->Integral(1,htmp->GetNbinsX())>20)
		{
		  htmp->SetTitle(Form("%s: %s of %s (%1.1f < p_{T} < %1.1f GeV/c);#Deltaz (cm)",partTitle[i],var_title[j],trkPtBins[k],trkPtBins[k+1]));
		  double range = 80;
		  if(i>2) range = 80;
		  TF1 *func = new TF1(Form("func_%s_pt%1.1f-%1.1f_%s",var_name[j],trkPtBins[k],trkPtBins[k+1],partName[i]),"gaus(0)+gaus(3)",-1*range,range);
		  func->SetParameter(1,0);
		  func->SetParameter(4,0);
		  func->SetParameter(2,5);
		  func->SetParameter(5,30);
		  htmp->Fit(func,"IR0");
		  htmp->GetYaxis()->SetNdivisions(505);
		  canvas[i][j][k/6]->cd(k%6+1);
		  htmp->SetTitle(Form(";%s [cm]",var_title[j]));
		  htmp->GetXaxis()->SetRangeUser(-80,100);
		  htmp->SetMarkerStyle(20);
		  htmp->SetMaximum(1.5*htmp->GetMaximum());
		  htmp->Draw("P");
		  func->SetLineColor(2);
		  func->Draw("sames");
		  TF1 *func2 = new TF1(Form("func2_%s_pt%1.1f-%1.1f_%s",var_name[j],trkPtBins[k],trkPtBins[k+1],partName[i]),"gaus",-1*range,range);
		  for(int ipar=0; ipar<3; ipar++) func2->SetParameter(ipar,func->GetParameter(ipar+3));
		  func2->SetLineColor(4);
		  func2->Draw("sames");
		  TPaveText *t = GetTitleText(Form("%s: %1.1f < p_{T} < %1.1f GeV/c",partTitle[i],trkPtBins[k],trkPtBins[k+1]),0.06);
		  t->Draw();
		  TLine *line = GetLine(0,0,0,htmp->GetMaximum()*0.7,1);
		  line->Draw();
		  for(int l=0; l<2; l++)
		    {
		      hMean[i][j][l]->SetBinContent(k+2,func->GetParameter(1+3*l));
		      hMean[i][j][l]->SetBinError(k+2, func->GetParError(1+3*l));
		      hSigma[i][j][l]->SetBinContent(k+2,fabs(func->GetParameter(2+3*l)));
		      hSigma[i][j][l]->SetBinError(k+2, func->GetParError(2+3*l));
		    }
		}
	      if(savePlot && (k%6==5 || k==nTrkPtBin-1))
		{
		  canvas[i][j][k/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s_Bin%d.png",run_type,partName[i],var_name[j],k/6));
		  canvas[i][j][k/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s_Bin%d.pdf",run_type,partName[i],var_name[j],k/6));
		}
	    }
	  list->Clear();
	  for(int k=0; k<2; k++)
	    {
	      list->Add(hMean[i][j][k]);
	    }
	  c = drawHistos(list,Form("FitMean_%s_%s",var_name[j],partTitle[i]),Form("%s: mean of %s distribution",partTitle[i],var_title[j]),true,0,10,kTRUE,-20,20,kFALSE,kTRUE,legName1,kTRUE,"",0.5,0.7,0.68,0.85,kTRUE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s_Mean.png",run_type,partName[i],var_name[j]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s_Mean.pdf",run_type,partName[i],var_name[j]));
	    }
	  
	  list->Clear();
	  for(int k=0; k<2; k++)
	    {
	      list->Add(hSigma[i][j][k]);
	    }
	  c = drawHistos(list,Form("FitSigma_%s_%s",var_name[j],partTitle[i]),Form("%s: #sigma of %s distribution",partTitle[i],var_title[j]),true,0,10,kTRUE,0,80,kFALSE,kTRUE,legName1,kTRUE,"",0.5,0.7,0.68,0.85,kTRUE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s_Sigma.png",run_type,partName[i],var_name[j]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/Embed_%s_Fit_%s_Sigma.pdf",run_type,partName[i],var_name[j]));
	    }
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outFileName),"update");
      for(int i=0; i<nPart; i++)
	{
	  TH2F *h2;
	  for(int j=0; j<2; j++)
	    {
	      if(j==0) h2 = hDzVsPt[i];
	      else     h2 = hDyVsPt[i];
	      if(h2->GetEntries()<1) continue;
	      for(int k=0; k<2; k++)
		{
		  hMean[i][j][k]->Write("",TObject::kOverwrite);
		  hSigma[i][j][k]->Write("",TObject::kOverwrite);
		}
	    }
	}
    }
}

//================================================
void JpsiMuon(const Int_t savePlot = 1, const Int_t saveHisto = 1)
{
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2); 
  const char *name = "JpsiMuon";

  TH2F *hMuonDzVsPt[2];
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhMuonDzDy_%s",trigName[kTrigType]));
  hn->GetAxis(0)->SetRangeUser(low_mass+0.001,high_mass-0.001);
  for(int i=0; i<2; i++)
    {
      hn->GetAxis(1)->SetRange(i+1,i+1);
      hMuonDzVsPt[i] = (TH2F*)hn->Projection(3,2);
      hMuonDzVsPt[i]->SetName(Form("hMuonDzVsPt_%d",i));
      hMuonDzVsPt[i]->Sumw2();
      cout << "# of di-muon pairs: " << hMuonDzVsPt[i]->Integral(0,-1)/2 << endl;
    }
  TH2F *hDzVsPt = (TH2F*)hMuonDzVsPt[0]->Clone("hDzVsPt_US-LS");
  hDzVsPt->Add(hMuonDzVsPt[1],-1);

  const int nTrkPtBin = 5;
  const double trkPtBins[nTrkPtBin+1] = {1,1.5,2.0,3.0,5.0,20.0};
  double trkPtbins[nTrkPtBin+2];
  trkPtbins[0] = 0;
  for(int i=1; i<nTrkPtBin+2; i++) trkPtbins[i] = trkPtBins[i-1];
  TH1F *hMean = new TH1F(Form("hMean_%s",name),Form("Mean of fitted #Deltaz distribution;p_{T} [GeV/c];mean [cm]"),nTrkPtBin+1,trkPtbins);
  TH1F *hSigma = new TH1F(Form("hSigma_%s",name),Form("Sigma of fitted #Deltaz distribution;p_{T} [GeV/c];#sigma [cm]"),nTrkPtBin+1,trkPtbins);
  const int nCanvas = nTrkPtBin/6+1;
  TCanvas *canvas[nCanvas];
  for(int i=0; i<nTrkPtBin; i++)
    {
      if(i%6==0) 
	{
	  canvas[i/6] = new TCanvas(Form("Fit_dz_%d",i/6),Form("Fit_dz_%d",i/6),1100,650);
	  canvas[i/6]->Divide(3,2);
	}
      int low_bin  = hDzVsPt->GetXaxis()->FindFixBin(trkPtBins[i]+1e-4);
      int high_bin = hDzVsPt->GetXaxis()->FindFixBin(trkPtBins[i+1]-1e-4);
      
      TH1F *htmp = (TH1F*)hDzVsPt->ProjectionY(Form("hTrkDz_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],trigName[kTrigType]),low_bin,high_bin);
      htmp->Rebin(5);
      htmp->SetTitle(Form("%s: #Deltaz of muons from J/#psi decay (%1.1f < p_{T} < %1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],trkPtBins[i],trkPtBins[i+1]));
      double range = 30;
      if(i>2) range = 20;
      TF1 *func = new TF1(Form("func_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],name),"gaus",-1*range,range);
      htmp->Fit(func,"IR0");
      htmp->GetYaxis()->SetNdivisions(505);
      canvas[i/6]->cd(i%6+1);
      htmp->SetTitle(";#Deltaz [cm]");
      htmp->GetXaxis()->SetRangeUser(-100,100);
      htmp->SetMarkerStyle(20);
      htmp->SetMaximum(1.5*htmp->GetMaximum());
      htmp->Draw("P");
      func->SetLineColor(2);
      func->Draw("sames");
      TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",trkPtBins[i],trkPtBins[i+1]),0.06);
      t->Draw();
      if(savePlot && (i%6==5 || i==nTrkPtBin-1))
	{
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%s.FitDz_Bin%d_%s.png",run_type,name,i/6,trigName[kTrigType]));
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%s.FitDz_Bin%d_%s.pdf",run_type,name,i/6,trigName[kTrigType]));
	}

      hMean->SetBinContent(i+2,func->GetParameter(1));
      hMean->SetBinError(i+2, func->GetParError(1));
      hSigma->SetBinContent(i+2,func->GetParameter(2));
      hSigma->SetBinError(i+2, func->GetParError(2));

    }

  hMean->SetMarkerStyle(21);
  c = draw1D(hMean);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%s.FitDz_Mean_%s.png",run_type,name,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%s.FitDz_Mean_%s.pdf",run_type,name,trigName[kTrigType]));
    }
  hSigma->SetMarkerStyle(21);
  c = draw1D(hSigma);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%s.FitDz_Sigma_%s.png",run_type,name,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%s.FitDz_Sigma_%s.pdf",run_type,name,trigName[kTrigType]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outFileName),"update");
      hMean->Write("Data_JpsiMuon_mean",TObject::kOverwrite);
      hSigma->Write("Data_JpsiMuon_sigma",TObject::kOverwrite);
      fout->Close();
    }
}

//================================================
void Compare2GausVs3Gaus(const Int_t savePlot = 1)
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s",outFileName),"read");
  const char *type_name[2] = {"2Gaus","3Gaus"};

  TList *list = new TList;
  TString legName[2] = {"Two-Gaussian fit","Three-Gaussian fit"};

  list->Clear();
  TH1F *hChi2[2];
  for(int i=0; i<2; i++)
    {
      hChi2[i] = (TH1F*)fin->Get(Form("Data_Muon_Chi2_%s",type_name[i]));
      list->Add(hChi2[i]);
    }
  c = drawHistos(list,"chi2","Chi2/NDF of #Deltaz fit;p_{T} [GeV/c];Chi2/ndf",kFALSE,2.0,3.8,kTRUE,0.1,300,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_FitDzChi2_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_FitDzChi2_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TH1F *hMean[2];
  for(int i=0; i<2; i++)
    {
      hMean[i] = (TH1F*)fin->Get(Form("Data_Muon_Mean_%s_0",type_name[i]));
      list->Add(hMean[i]);
    }
  c = drawHistos(list,"mean","Mean of fitted #Deltaz distribution",kFALSE,2.0,3.8,kFALSE,0.1,300,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.65,0.85,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_FitDzMean_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_FitDzMean_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  TH1F *hSigma[2];
  for(int i=0; i<2; i++)
    {
      hSigma[i] = (TH1F*)fin->Get(Form("Data_Muon_Sigma_%s_0",type_name[i]));
      hSigma[i]->SetMarkerStyle(20);
      hSigma[i]->SetMarkerColor(color[i]);
      hSigma[i]->SetLineColor(color[i]);
    }
  hSigma[0]->SetMaximum(16);
  c = draw1D(hSigma[0],"Sigma of fitted #Deltaz distribution");
  hSigma[1]->Draw("sames");
  TF1 *fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)",1,20);
  fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.444217);
  fResDzVsPt->SetLineColor(4);
  fResDzVsPt->Draw("sames");
  TLegend *leg = new TLegend(0.55,0.65,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hSigma[0],legName[0].Data(),"PL");
  leg->AddEntry(hSigma[1],legName[1].Data(),"PL");
  leg->AddEntry(fResDzVsPt,"Run13 embedding","L");
  leg->Draw();
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_FitDzSigma_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_FitDzSigma_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TH1F *hPurity[2];
  for(int i=0; i<2; i++)
    {
      hPurity[i] = (TH1F*)fin->Get(Form("Data_Muon_Purity_%s",type_name[i]));
      list->Add(hPurity[i]);
    }
  c = drawHistos(list,"hPurity","Purity of selected muon sample by #Deltaz cut",kFALSE,2.0,3.8,kFALSE,0.1,300,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.2,0.4,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_DzPurity_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_DzPurity_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TH1F *hEff[2];
  for(int i=0; i<2; i++)
    {
      hEff[i] = (TH1F*)fin->Get(Form("Data_Muon_Eff_%s",type_name[i]));
      list->Add(hEff[i]);
    }
  c = drawHistos(list,"hEff","Efficiency of #Deltaz cut",kFALSE,2.0,3.8,kFALSE,0.1,300,kFALSE,kTRUE,legName,kTRUE,"",0.2,0.4,0.2,0.4,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_DzCutEff_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sCompare_DzCutEff_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }
}

//================================================
void MuonCandidates(const Int_t savePlot = 1, const Int_t saveHisto = 1)
{
  TH2F *hDzVsPt = (TH2F*)f->Get(Form("mhMuonDzVsPt_%s",trigName[kTrigType]));
  c = draw2D(hDzVsPt,Form("%s: #Deltaz of muon candidates",trigName[kTrigType]));
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sDeltaZ_vs_pt_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sDeltaZ_vs_pt_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  // Fitting
  const int type = 1;
  const char *type_name[2] = {"2Gaus","3Gaus"};

  TH1F *hMthDz = (TH1F*)hDzVsPt->ProjectionY(Form("hTrkDz_%s_proj",trigName[kTrigType]));
  hMthDz->SetTitle(Form("%s: #Deltaz of matched track-hit pairs;#Deltaz (cm)",trigName[kTrigType]));
  Double_t range = 80;
  TF1 *func;
  if(type==0)
    {
      func = new TF1("func","gaus(0)+gaus(3)",-1*range,range);
      func->SetParameter(1,0);
      func->SetParameter(4,0);
      func->SetParameter(2,10);
      func->SetParameter(5,50);
    }
  else if(type==1)
    {
      func = new TF1("func","gaus(0)+gaus(3)+gaus(6)",-1*range,range);
      func->SetParameter(1,0);
      func->SetParameter(4,0);
      func->SetParameter(7,0);
      func->SetParameter(2,10);
      func->SetParameter(5,20);
      func->SetParameter(8,100);
    }
  c = FitDeltaZ(hMthDz,func,range,20.,type);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_%s_All_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_%s_All_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  //return;
  // pt dependence
  const int nTrkPtBin = 26;
  const double trkPtBins[nTrkPtBin+1] = {1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0};
  double trkPtbins[nTrkPtBin+2];
  trkPtbins[0] = 0;
  for(int i=1; i<nTrkPtBin+2; i++) trkPtbins[i] = trkPtBins[i-1];
  TH1F *hMean[3];
  TH1F *hSigma[3];
  for(int k=0; k<3; k++)
    {
      hMean[k] = new TH1F(Form("Data_Muon_Mean_%s_%d",type_name[type],k),Form("Mean of fitted #Deltaz distribution;p_{T} [GeV/c];mean [cm]"),nTrkPtBin,trkPtbins); 
      hSigma[k] = new TH1F(Form("Data_Muon_Sigma_%s_%d",type_name[type],k),Form("Sigma of fitted #Deltaz distribution;p_{T} [GeV/c];#sigma [cm]"),nTrkPtBin,trkPtbins);
    }
  TH1F *hPurity = new TH1F(Form("Data_Muon_Purity_%s",type_name[type]),Form("Purity of selected muon sample by #Deltaz cut;p_{T} [GeV/c];purity"),nTrkPtBin,trkPtbins);
  TH1F *hEff = new TH1F(Form("Data_Muon_Eff_%s",type_name[type]),Form("Efficiency of #Deltaz cut;p_{T} [GeV/c];Efficiency"),nTrkPtBin,trkPtbins);
  TH1F *hChi2 = new TH1F(Form("Data_Muon_Chi2_%s",type_name[type]),Form("Chi2/NDF of #Deltaz fit;p_{T} [GeV/c];Chi2/ndf"),nTrkPtBin,trkPtbins);
  TF1 *fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)");
  fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.444217);
  const int nCanvas = nTrkPtBin/6+1;
  TCanvas *canvas[nCanvas];
  for(int i=0; i<nTrkPtBin; i++)
    {
      if(i%6==0) 
	{
	  canvas[i/6] = new TCanvas(Form("Fit_dz_%d",i/6),Form("Fit_dz_%d",i/6),1100,650);
	  canvas[i/6]->Divide(3,2);
	}
      int low_bin  = hDzVsPt->GetXaxis()->FindFixBin(trkPtBins[i]+1e-4);
      int high_bin = hDzVsPt->GetXaxis()->FindFixBin(trkPtBins[i+1]-1e-4);
      
      TH1F *htmp = (TH1F*)hDzVsPt->ProjectionY(Form("hTrkDz_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],trigName[kTrigType]),low_bin,high_bin);
      htmp->SetTitle(Form("%s: #Deltaz of matched track-hit pairs (%1.1f < p_{T} < %1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],trkPtBins[i],trkPtBins[i+1]));

      TF1 *func;
      if(type==0)
	{
	  func = new TF1(Form("func_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],type_name[type]),"gaus(0)+gaus(3)",-1*range,range);
	  func->SetParameter(1,0);
	  func->SetParameter(4,0);
	  func->SetParameter(2,10);
	  func->SetParameter(5,50);
	}
      if(type==1)
	{
	  func = new TF1(Form("func_pt%1.1f-%1.1f_%s",trkPtBins[i],trkPtBins[i+1],type_name[type]),"gaus(0)+gaus(3)+gaus(6)",-1*range,range);
	  func->SetParameter(1,0);
	  func->SetParLimits(4,-3,3);
	  func->SetParameter(7,0);
	  func->SetParameter(2,8);
	  func->SetParameter(5,20);
	  func->SetParameter(8,40);
	  if(i>10) func->SetParameter(2,5);
	}
      htmp->Fit(func,"IR0");
      hChi2->SetBinContent(i+2, func->GetChisquare()/func->GetNDF());
      hChi2->SetBinError(i+2, 0);
      htmp->GetYaxis()->SetNdivisions(505);
      canvas[i/6]->cd(i%6+1);
      htmp->SetTitle(";#Deltaz [cm]");
      htmp->GetXaxis()->SetRangeUser(-100,100);
      htmp->SetMarkerStyle(20);
      htmp->Draw("HIST");
      func->SetLineColor(2);
      func->Draw("sames");

      TF1 *func1, *func2, *func3;
      func1 = new TF1("func1","gaus",-1*range,range);
      func1->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
      if(type==0)
	{
	  func2 = new TF1("func2","gaus",-1*range,range);
	  func2->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
	  func2->SetLineColor(4);
	  func2->Draw("sames");
	}
      else if(type==1)
	{
	  func2 = new TF1("func2","gaus(0)+gaus(3)",-1*range,range);
	  func2->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5),func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
	  func2->SetLineColor(kGreen+2);
	  func2->Draw("sames");

	  func3 = new TF1("func3","gaus",-1*range,range);
	  func3->SetParameters(func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
	  func3->SetLineColor(4);
	  func3->Draw("sames");
	}


      TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",trkPtBins[i],trkPtBins[i+1]),0.06);
      t->Draw();
      if(savePlot && (i%6==5 || i==nTrkPtBin-1))
	{
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_%s_Bin%d_%s.png",run_type,run_cfg_name.Data(),type_name[type],i/6,trigName[kTrigType]));
	  canvas[i/6]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_%s_Bin%d_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],i/6,trigName[kTrigType]));
	}

      for(int k=0; k<2+type; k++)
	{
	  hMean[k]->SetBinContent(i+2,func->GetParameter(1+k*3));
	  hMean[k]->SetBinError(i+2, func->GetParError(1+k*3));
	  hSigma[k]->SetBinContent(i+2,func->GetParameter(2+k*3));
	  hSigma[k]->SetBinError(i+2, func->GetParError(2+k*3));
	}

      // pt dependent dz cut
      double pt = hPurity->GetBinCenter(i+2);
      double dz_sigma = fResDzVsPt->Eval(pt);
      double cut = 999.;
      if(pt<3) cut = 2 * dz_sigma;
      else     cut = 2.5 * dz_sigma;
      
      double all = htmp->Integral(htmp->FindFixBin(-1*cut),htmp->FindFixBin(cut));
      double all_err = TMath::Sqrt(all);
      double bkg = (func->Integral(-1*cut,cut)-func1->Integral(-1*cut,cut)) * 1./htmp->GetBinWidth(1);
      double bkg_err = TMath::Sqrt(bkg);
      double signal = all - bkg;
      double signal_err = TMath::Sqrt(all_err*all_err+bkg_err*bkg_err);
      double purity = 1 - bkg/all;
      double purity_err = TMath::Sqrt(bkg_err*bkg_err/all/all+bkg*bkg/all/all/all/all*all_err*all_err);
      hPurity->SetBinContent(i+2,purity);
      hPurity->SetBinError(i+2,purity_err);
      
      double out = 
	func1->Integral(func1->GetParameter(1)-5*func1->GetParameter(2),-1*cut)* 1./htmp->GetBinWidth(1) +
	func1->Integral(cut, func1->GetParameter(1)+5*func1->GetParameter(2))* 1./htmp->GetBinWidth(1);
      double eff = signal/(signal+out);
      hEff->SetBinContent(i+2,eff);
      hEff->SetBinError(i+2,1e-10);

    }
  hChi2->SetMarkerStyle(21);
  hChi2->SetMaximum(1.2*hChi2->GetMaximum());
  c = draw1D(hChi2);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_Chi2_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_Chi2_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hMean[0]->SetMarkerStyle(21);
  hMean[0]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c = draw1D(hMean[0]);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_Mean_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_Mean_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hSigma[0]->SetMarkerStyle(21);
  hSigma[0]->GetYaxis()->SetRangeUser(0,15);
  c = draw1D(hSigma[0]);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_Sigma_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sFitDz_Sigma_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hPurity->SetMarkerStyle(21);
  hPurity->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hPurity);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sDzCut_Purity_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sDzCut_Purity_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }
  hEff->SetMarkerStyle(21);
  hEff->GetYaxis()->SetRangeUser(0,1.1);
  c = draw1D(hEff);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sDzCut_Efficiency_%s_%s.png",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DeltaZ/%sDzCut_Efficiency_%s_%s.pdf",run_type,run_cfg_name.Data(),type_name[type],trigName[kTrigType]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outFileName),"update");
      for(int k=0; k<2+type; k++)
	{
	  hMean[k]->Write("",TObject::kOverwrite);
	  hSigma[k]->Write("",TObject::kOverwrite);
	}
      hPurity->Write("",TObject::kOverwrite);
      hEff->Write("",TObject::kOverwrite);
      hChi2->Write("",TObject::kOverwrite);
      fout->Close();
    }
}

//================================================
TCanvas *FitDeltaZ(TH1 *histo, TF1 *func, const Double_t range1 = 50., Double_t range2 = 20., int type = 0)
{
  histo->Fit(func,"IR0");
  histo->GetYaxis()->SetNdivisions(505);
  c = draw1D(histo,"",kFALSE,kFALSE);
  func->SetLineColor(2);
  func->Draw("sames");

  TF1 *func1, *func2;
  if(type==0)
    {
      func1 = new TF1("func1","gaus",-1*range1,range1);
      func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
      func1->SetLineColor(4);
      func1->Draw("sames");
    }
  else if(type==1)
    {
      func1 = new TF1("func1","gaus(0)+gaus(3)",-1*range1,range1);
      func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5),func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
      func1->SetLineColor(kGreen+2);
      func1->Draw("sames");

      func2 = new TF1("func2","gaus",-1*range1,range1);
      func2->SetParameters(func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
      func2->SetLineColor(4);
      func2->Draw("sames");
    }
  TPaveText *t1 = GetPaveText(0.2,0.3,0.65,0.78,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f] cm",range1,range1));
  t1->AddText(Form("S/B ~ %1.3f:1",((func->Integral(-1*range1,range1))-(func1->Integral(-1*range1,range1)))/(func1->Integral(-1*range1,range1))));
  t1->Draw();
  t1 = GetPaveText(0.2,0.3,0.47,0.6,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f] cm",range2,range2));
  t1->AddText(Form("S/B ~ %1.3f:1",((func->Integral(-1*range2,range2))-(func1->Integral(-1*range2,range2)))/(func1->Integral(-1*range2,range2))));
  t1->Draw();
  return c;
}
