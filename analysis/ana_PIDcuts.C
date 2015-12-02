const int year = 2013;
TF1 *fResVsPt[3];

//================================================
void ana_PIDcuts()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      run_type = "Run13_pp500";
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
    }

  fResVsPt[0]= new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[0]->SetParameters(-32.6793, 32.6034, 0.444217);
  fResVsPt[1] = new TF1("fResDyVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[1]->SetParameters(-17.6867, 18.4528, 0.637142);
  fResVsPt[2] = new TF1("fResDtofVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[2]->SetParameters(0.0817528, 0.0169419, 4.34897);


  //makeData();
  anaEmbed();
  //anaData();
  //nSigmaPi();
}


//================================================
void nSigmaPi(const bool savePlot = 1, const bool saveHisto = 0)
{
  const char *name[2] = {"Data","Embed"};
  TH2F *hNsigmaPiVsPt[2];

  TFile *fMB = TFile::Open("Rootfiles/MB.V0.pion.root","read");
  hNsigmaPiVsPt[0] = (TH2F*)fMB->Get("hPinSigmaPivsPt");

  TFile *fEmbed = TFile::Open("output/Run14.AuAu200.Pion.Embed.root","read");
  hNsigmaPiVsPt[1] = (TH2F*)fEmbed->Get("mhNsigmaPiVsPt");

  const int nbins = 24;
  const double xbins[nbins+1] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0};
  TH1F *hEff[2];
  TList *list = new TList;
  TString legName[2] = {"Data","Embed"};
  for(int i=0; i<2; i++)
    {
      draw2D(hNsigmaPiVsPt[i]);
      hEff[i] = new TH1F(Form("_%s_pion_eff",name[i]),Form("%s: efficiency of n#sigma_{#pi} cut for #pi sample;p_{T} (GeV/c)",name[i]),nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  double pt1 = xbins[bin-1]; 
	  double pt2 = xbins[bin];
	  int start_bin = hNsigmaPiVsPt[i]->GetXaxis()->FindBin(pt1+1e-4);
	  int end_bin = hNsigmaPiVsPt[i]->GetXaxis()->FindBin(pt2-1e-4);
	  TH1F *htmp = (TH1F*)hNsigmaPiVsPt[i]->ProjectionY(Form("%s_%d",hNsigmaPiVsPt[i]->GetName(),bin),start_bin,end_bin);
	  double all = htmp->Integral(htmp->FindBin(-5),htmp->FindBin(5));
	  double acc = htmp->Integral(htmp->FindBin(-1.5),htmp->FindBin(2.5));
	  if(all>0 && acc>0)
	    {
	      hEff[i]->SetBinContent(bin,acc/all);
	      hEff[i]->SetBinError(bin,acc/all * sqrt(1./all + 1./acc));
	    }
	}
      list->Add(hEff[i]);
    }
  c = drawHistos(list,"nSigmaPi_eff","Efficiency of n#sigma_{#pi} cut for #pi sample;p_{T} (GeV/c)",kFALSE,0,10,kTRUE,0.5,1,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.3,0.5,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataVsEmbed_nSigmaPi_eff.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataVsEmbed_nSigmaPi_eff.png",run_type));
    }
  
}

//================================================
void anaData(const bool savePlot = 0, const bool saveHisto = 0)
{
  TFile *fMB = TFile::Open("Rootfiles/MB.V0.pion.root","read");
  TFile *fin = TFile::Open("Rootfiles/Run14.AuAu200.DataVsEmbed.root","update");
  const char *name[4] = {"nSigmaPi","dz","dy","dtof"};
  const char *title[4] = {"n#sigma_{#pi}","#Deltaz","#Deltay","#Deltatof"};
  const char *unit[4] = {""," (cm)"," (cm)", " (ns)"};
  TH2F *hMuonDisVsPt[4];
  hMuonDisVsPt[0] = (TH2F*)fMB->Get("hKsPinSigmaPivsPt");
  hMuonDisVsPt[1] = (TH2F*)fin->Get("Data_muon_dz_vs_pt");
  hMuonDisVsPt[2] = (TH2F*)fin->Get("Data_muon_dy_vs_pt");
  hMuonDisVsPt[3] = (TH2F*)fin->Get("Data_muon_dtof_vs_pt");
  const int nbins = 6;
  const double xbins[nbins+1] = {1.0,1.5,2.0,2.5,3.0,5.0,10.0};

  TFile *fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");
  TH2F *hEmbedDisVsPt[4];
  hEmbedDisVsPt[0] = (TH2F*)fEmbed->Get("hTrkNSigmaPi_MCreco_di_mu");
  hEmbedDisVsPt[1] = (TH2F*)fEmbed->Get("hTrkDzVsPt_MCreco_di_mu");
  hEmbedDisVsPt[2] = (TH2F*)fEmbed->Get("hTrkDyVsPt_MCreco_di_mu");
  hEmbedDisVsPt[3] = (TH2F*)fEmbed->Get("hMcDeltaTof_di_mu");

  TH1F *hEmbedEff[4];
  for(int i=0; i<4; i++)
    hEmbedEff[i] = (TH1F*)fin->Get(Form("Embed_muon_%s_Eff",name[i]));

  // check US and LS
  TH2F *hDisVsPt[2][3];
  TH1F *hDis[2][3][nbins], *hSub[3][nbins];
  TH1F *hEmbed[3][nbins];
  for(int i=0; i<3; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_distribution",name[i+1]),Form("%s_distribution",name[i+1]),1100,700);
      c->Divide(3,2);
      hDisVsPt[0][i] = (TH2F*)fin->Get(Form("US_%s_neg_muon",name[i+1]));
      hDisVsPt[0][i]->Add((TH2F*)fin->Get(Form("US_%s_pos_muon",name[i+1])));

      hDisVsPt[1][i] = (TH2F*)fin->Get(Form("LS_%s_neg_muon",name[i+1]));
      hDisVsPt[1][i]->Add((TH2F*)fin->Get(Form("LS_%s_pos_muon",name[i+1])));

      for(int bin=1; bin<=nbins; bin++)
	{
	  for(int j=0; j<2; j++)
	    {
	      int start_bin = hDisVsPt[j][i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	      int end_bin = hDisVsPt[j][i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	      hDis[j][i][bin-1] = (TH1F*)hDisVsPt[j][i]->ProjectionY(Form("%s_%d_bin%d",name[i+1],j,bin),start_bin,end_bin);
	      hDis[j][i][bin-1]->SetMarkerStyle(24+j);
	      hDis[j][i][bin-1]->SetMarkerColor(color[j]);
	      hDis[j][i][bin-1]->SetLineColor(color[j]);
	    }
	  hSub[i][bin-1] = (TH1F*)hDis[0][i][bin-1]->Clone(Form("%s_bin%d",name[i+1],bin));
	  hSub[i][bin-1]->Add(hDis[1][i][bin-1],-1);
	  hSub[i][bin-1]->Scale(1./hSub[i][bin-1]->Integral());
	  hSub[i][bin-1]->SetMarkerStyle(20);
	  hSub[i][bin-1]->SetMarkerColor(1);
	  hSub[i][bin-1]->SetLineColor(1);
	  hSub[i][bin-1]->SetTitle(Form(";%s%s",title[i+1],unit[i+1]));
	  if(i<2) hSub[i][bin-1]->Rebin(5);
	  else    hSub[i][bin-1]->Rebin(4);
	  if(i==0) hSub[i][bin-1]->GetXaxis()->SetRangeUser(-100,100);
	  else if(i==2) hSub[i][bin-1]->GetXaxis()->SetRangeUser(-3,4);
	  hSub[i][bin-1]->SetMinimum(3*hSub[i][bin-1]->GetMinimum());
	  c->cd(bin);
	  //hDis[0][i][bin-1]->Draw();
	  //hDis[1][i][bin-1]->Draw("sames");
	  hSub[i][bin-1]->Draw("");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();

	  double up = hSub[i][bin-1]->GetMaximum();
	  double low = hSub[i][bin-1]->GetMinimum();
	  double pt = (xbins[bin-1]+xbins[bin])/2;
	  double min = -999, max = -999;
	  if(i<2)
	    {
	      if(pt<3)
		{
		  min = -1 * fResVsPt[i]->Eval(pt) * 2;
		  max = -1 * min;
		}
	      else
		{
		  min = -1 * fResVsPt[i]->Eval(pt) * 2.5;
		  max = -1 * min;
		}
	      cout << min << "  " << pt << endl;
	      TLine *line = GetLine(min,low,min,up);
	      line->Draw();
	      TLine *line = GetLine(max,low,max,up);
	      line->Draw();
	    }
	  else
	    {
	      min = -100;
	      max = fResVsPt[i]->Eval(pt);
	      TLine *line = GetLine(max,low,max,up);
	      line->Draw();
	    }

	  // embedding data
	  int start_bin = hEmbedDisVsPt[i+1]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hEmbedDisVsPt[i+1]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hEmbed[i][bin-1] = (TH1F*)hEmbedDisVsPt[i+1]->ProjectionY(Form("Emebd_%s_bin%d",name[i+1],bin),start_bin,end_bin);
	  hEmbed[i][bin-1]->SetLineColor(4);
	  hEmbed[i][bin-1]->Scale(hSub[i][bin-1]->GetMaximum()/hEmbed[i][bin-1]->GetBinContent(hEmbed[i][bin-1]->FindBin(0)));
	  hEmbed[i][bin-1]->Draw("sames HIST");
	}

      c->cd(1);
      TLegend *leg = new TLegend(0.6,0.68,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(hSub[i][0],"Data","P");
      leg->AddEntry(hEmbed[i][0],"Embedding","L");
      leg->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataVsEmbed_%s_InPtBins.pdf",run_type,name[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataVsEmbed_%s_InPtBins.png",run_type,name[i+1]));
	}
    }

  TH1F *hDisAll[4], *hDisAcc[4], *hEff[4];
  for(int i=0; i<4; i++)
    {
      if(i==1 || i==2) hMuonDisVsPt[i]->RebinY(4);
      hDisAll[i] = new TH1F(Form("Data_%s_all",name[i]),Form("%s distribution of input muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      hDisAcc[i] = new TH1F(Form("Data_%s_Acc",name[i]),Form("%s distribution of accepted muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *htmp = (TH1F*)hMuonDisVsPt[i]->ProjectionY(Form("%s_%d",name[i],bin),start_bin,end_bin);
	  double all_err;
	  double all = htmp->IntegralAndError(0,-1,all_err);
	  hDisAll[i]->SetBinContent(bin,all);
	  hDisAll[i]->SetBinError(bin,all_err);

	  double pt = hDisAll[i]->GetBinCenter(bin);
	  double min = -999, max = -999;
	  if(i==0) { min = -1.5; max = 2.5; }
	  else if(i<3)
	    {
	      if(pt<3)
		{
		  min = -1 * fResVsPt[i-1]->Eval(pt) * 2;
		  max = -1 * min;
		}
	      else
		{
		  min = -1 * fResVsPt[i-1]->Eval(pt) * 2.5;
		  max = -1 * min;
		}
	    }
	  else
	    {
	      min = -100;
	      max = fResVsPt[i-1]->Eval(pt);
	    }
	  double acc_err;
	  double acc = htmp->IntegralAndError(htmp->FindBin(min),htmp->FindBin(max),acc_err);
	  if(i==1)
	    {
	      cout << all << "  " << all_err << endl;
	      cout << htmp->FindBin(min) << "  " << htmp->FindBin(max) << endl;
	      cout << acc << "  " << acc_err << endl;
	    }
	  hDisAcc[i]->SetBinContent(bin,acc);
	  hDisAcc[i]->SetBinError(bin,acc_err);
	}
      hEff[i] = (TH1F*)hDisAcc[i]->Clone(Form("Data_%s_Eff",name[i]));
      hEff[i]->Divide(hDisAll[i]);

      hMuonDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      if(i==2) hMuonDisVsPt[i]->GetYaxis()->SetRangeUser(-100,100);
      if(i==3) hMuonDisVsPt[i]->GetYaxis()->SetRangeUser(-1,2);
      if(i==0) c = draw2D(hMuonDisVsPt[i],Form("Data: %s distribution for K_{s}^{0} #pi;p_{T} (GeV/c);%s%s",title[i],title[i],unit[i]));
      else     c = draw2D(hMuonDisVsPt[i],Form("Data: %s distribution for J/#psi muon;p_{T} (GeV/c);%s%s",title[i],title[i],unit[i]));
      if(i==0)
	{
	  TLine *line = GetLine(0,-1.5,10,-1.5);
	  line->Draw();
	  line = GetLine(0,2.5,10,2.5);
	  line->Draw(); 
	}
      else if(i<3)
	{
	  TF1 *func11 = (TF1*)fResVsPt[i-1]->Clone(Form("%s_2sigma_upBound",fResVsPt[i-1]->GetName()));
	  func11->SetRange(0.5,3);
	  func11->SetParameter(0,2*func11->GetParameter(0));
	  func11->SetParameter(1,2*func11->GetParameter(1));
	  func11->SetLineColor(2);
	  func11->Draw("sames");
	  TF1 *func12 = (TF1*)fResVsPt[i-1]->Clone(Form("%s_2sigma_lowBound",fResVsPt[i-1]->GetName()));
	  func12->SetRange(0.5,3);
	  func12->SetParameter(0,-2*func12->GetParameter(0));
	  func12->SetParameter(1,-2*func12->GetParameter(1));
	  func12->SetLineColor(2);
	  func12->Draw("sames");
	  TF1 *func21 = (TF1*)fResVsPt[i-1]->Clone(Form("%s_2.5sigma_upBound",fResVsPt[i-1]->GetName()));
	  func21->SetRange(3,10);
	  func21->SetParameter(0,2.5*func21->GetParameter(0));
	  func21->SetParameter(1,2.5*func21->GetParameter(1));
	  func21->SetLineColor(2);
	  func21->Draw("sames");
	  TF1 *func22 = (TF1*)fResVsPt[i-1]->Clone(Form("%s_2.5sigma_lowBound",fResVsPt[i-1]->GetName()));
	  func22->SetRange(3,10);
	  func22->SetParameter(0,-2.5*func22->GetParameter(0));
	  func22->SetParameter(1,-2.5*func22->GetParameter(1));
	  func22->SetLineColor(2);
	  func22->Draw("sames");
	}
      else
	{
	  TF1 *func = (TF1*)fResVsPt[i-1]->Clone(Form("%s_upBound",fResVsPt[i-1]->GetName()));
	  func->SetRange(0.5,10);
	  func->Draw("sames");
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_muon_%s_vs_pt.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_muon_%s_vs_pt.png",run_type,name[i]));
	}

      hEff[i]->SetMarkerStyle(21);
      if(i==0) hEff[i]->GetYaxis()->SetRangeUser(0.6,1);
      else     hEff[i]->GetYaxis()->SetRangeUser(0.3,1.3);
      c = draw1D(hEff[i],Form("Embed: efficiency of %s cut;p_{T} (GeV/c);Efficiency",title[i]));
      gPad->SetGridy();
      hEmbedEff[i]->SetMarkerColor(2);
      hEmbedEff[i]->Draw("sames");
      TLegend *leg = new TLegend(0.2,0.15,0.4,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hEff[i],"Data","P");
      leg->AddEntry(hEmbedEff[i],"Embedding","P");
      leg->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_muon_%s_eff.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_muon_%s_eff.png",run_type,name[i]));
	}
    }
  return;

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14.AuAu200.DataVsEmbed.root","update");
      for(int i=0; i<4; i++)
	{
	  hMuonDisVsPt[i]->Write(Form("Embed_%s_vs_pt",name[i]),TObject::kOverwrite);
	  hDisAll[i]->Write("",TObject::kOverwrite);
	  hDisAcc[i]->Write("",TObject::kOverwrite);
	  hEff[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void anaEmbed(const bool savePlot = 0, const bool saveHisto = 0)
{
  TFile *fEmbed = 0x0;
  if(year==2013) fEmbed = TFile::Open("output/Run13.pp500.jpsi.Embed.root","read");
  if(year==2014) fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");
  const char *name[4] = {"Dz","Dy","Dtof","NsigmaPi"};
  const char *title[4] = {"#Deltaz","#Deltay","#Deltatof","n#sigma_{#pi}"};
  const char *unit[4] = {" (cm)"," (cm)", " (ns)",""};
  TH2F *hMuonDisVsPt[4];
  hMuonDisVsPt[0] = (TH2F*)fEmbed->Get("hTrkDzVsPt_MCreco_di_mu");
  hMuonDisVsPt[1] = (TH2F*)fEmbed->Get("hTrkDyVsPt_MCreco_di_mu");
  hMuonDisVsPt[2] = (TH2F*)fEmbed->Get("hMcDeltaTof_di_mu");
  hMuonDisVsPt[3] = (TH2F*)fEmbed->Get("hTrkNSigmaPi_MCreco_di_mu");
  const int nbins = 24;
  const double xbins[nbins+1] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0};

  TH1F *hDisAll[4], *hDisAcc[4], *hEff[4];
  TH1F *hFitSigma[2];
  for(int i=0; i<4; i++)
    {
      hDisAll[i] = new TH1F(Form("Embed_JpsiMuon_%s_all",name[i]),Form("%s distribution of input muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      hDisAcc[i] = new TH1F(Form("Embed_JpsiMuon_%s_Acc",name[i]),Form("%s distribution of accepted muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      TCanvas *cfit = 0x0;
      if(i<2)
	{
	  hFitSigma[i] = new TH1F(Form("Embed_JpsiMuon_Fit%s_Sigma",name[i]),Form("Embed: width of fitting %s distribution;p_{T} (GeV/c);#sigma (cm)",name[i]),nbins,xbins);
	  cfit = new TCanvas(Form("Fit_%s",name[i]),Form("Fit_%s",name[i]),1000,600);
	  cfit->Divide(6,4);
	}
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *htmp = (TH1F*)hMuonDisVsPt[i]->ProjectionY(Form("%s_%d",name[i],bin),start_bin,end_bin);
	  if(cfit) 
	    {
	      cfit->cd(bin);
	      htmp->SetTitle("");
	      htmp->GetXaxis()->SetRangeUser(-50,50);
	      htmp->Draw();
	      if(htmp->GetEntries()>0)
		{
		  TF1 *functmp = new TF1(Form("func_%s_%d",name[i],bin),"gaus",-30,30);
		  htmp->Fit(functmp,"IR0Q");
		  functmp->SetLineColor(2);
		  functmp->Draw("same");
		  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xbins[bin-1],xbins[bin]),0.07);
		  t1->Draw();
		  hFitSigma[i]->SetBinContent(bin,functmp->GetParameter(2));
		  hFitSigma[i]->SetBinError(bin,functmp->GetParError(2));
		}
	    }
	  double all = htmp->GetEntries();
	  hDisAll[i]->SetBinContent(bin,all);
	  hDisAll[i]->SetBinError(bin,sqrt(all));
	  double pt = hDisAll[i]->GetBinCenter(bin);
	  double min = -999, max = -999;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(year==2013) min = -2.5 * sigma;
	      else if (year==2014)
		{
		  if(pt<3) min = -2 * sigma;
		  else min = -2.5 * sigma;
		}
	      max = -1 * min;
	    }
	  else if (i==2)
	    {
	      min = -100;
	      max = fResVsPt[i]->Eval(pt);
	    }
	  else
	    { min = -1; max = 3; }
	  double acc = htmp->Integral(htmp->FindBin(min),htmp->FindBin(max));
	  //cout << htmp->FindBin(min) << "  " << htmp->FindBin(max) << endl;
	  hDisAcc[i]->SetBinContent(bin,acc);
	  hDisAcc[i]->SetBinError(bin,sqrt(acc));
	}
      if(i<2)
	{
	  hFitSigma[i]->SetMarkerStyle(21);
	  hFitSigma[i]->SetMarkerColor(2);
	  hFitSigma[i]->SetLineColor(2);
	  c = draw1D(hFitSigma[i],"");
	  fResVsPt[i]->Draw("sames");
	}
      hEff[i] = (TH1F*)hDisAcc[i]->Clone(Form("Embed_JpsiMuon_%s_Eff",name[i]));
      hEff[i]->Divide(hDisAll[i]);

      hMuonDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      if(i==1) hMuonDisVsPt[i]->GetYaxis()->SetRangeUser(-100,100);
      if(i==2) hMuonDisVsPt[i]->GetYaxis()->SetRangeUser(-1,2);
      c = draw2D(hMuonDisVsPt[i],Form("Embed: %s distribution for muon tracks;p_{T} (GeV/c);%s%s",title[i],title[i],unit[i]));
      if(i<2)
	{
	  double nsigma1 = 2.5, nsigma2 = 2.5;
	  if(year==2014) nsigma1 = 2;
	  TF1 *func11 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_upBound",fResVsPt[i]->GetName()));
	  func11->SetRange(0.5,3);
	  func11->SetParameter(0,nsigma1*func11->GetParameter(0));
	  func11->SetParameter(1,nsigma1*func11->GetParameter(1));
	  func11->SetLineColor(2);
	  func11->Draw("sames");
	  TF1 *func12 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_lowBound",fResVsPt[i]->GetName()));
	  func12->SetRange(0.5,3);
	  func12->SetParameter(0,-nsigma1*func12->GetParameter(0));
	  func12->SetParameter(1,-nsigma1*func12->GetParameter(1));
	  func12->SetLineColor(2);
	  func12->Draw("sames");
	  TF1 *func21 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_upBound",fResVsPt[i]->GetName()));
	  func21->SetRange(3,10);
	  func21->SetParameter(0,nsigma2*func21->GetParameter(0));
	  func21->SetParameter(1,nsigma2*func21->GetParameter(1));
	  func21->SetLineColor(2);
	  func21->Draw("sames");
	  TF1 *func22 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_lowBound",fResVsPt[i]->GetName()));
	  func22->SetRange(3,10);
	  func22->SetParameter(0,-nsigma2*func22->GetParameter(0));
	  func22->SetParameter(1,-nsigma2*func22->GetParameter(1));
	  func22->SetLineColor(2);
	  func22->Draw("sames");
	}
      else if(i==2)
	{
	  TF1 *func = (TF1*)fResVsPt[i]->Clone(Form("%s_upBound",fResVsPt[i]->GetName()));
	  func->SetRange(0.5,10);
	  func->Draw("sames");
	}
      else
	{
	  TLine *line = GetLine(0,-1,10,-1);
	  line->Draw();
	  line = GetLine(0,3,10,3);
	  line->Draw(); 
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_vs_pt.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_vs_pt.png",run_type,name[i]));
	}

      hEff[i]->SetMarkerStyle(21);
      hEff[i]->GetYaxis()->SetRangeUser(0.3,1.2);
      c = draw1D(hEff[i],Form("Embed: efficiency of %s cut;p_{T} (GeV/c);Efficiency",title[i]));
      gPad->SetGridy();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_eff.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_eff.png",run_type,name[i]));
	}
    }

  if(saveHisto)
    {
      TFile *fout = 0x0;
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.PIDcuts.root","update");
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.PIDcuts.root","update");
      for(int i=0; i<4; i++)
	{
	  hMuonDisVsPt[i]->Write(Form("Embed_JpsiMuon_%sVsPt",name[i]),TObject::kOverwrite);
	  hDisAll[i]->Write("",TObject::kOverwrite);
	  hDisAcc[i]->Write("",TObject::kOverwrite);
	  hEff[i]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void makeData()
{
  const int nHistos = 4;
  TH2F *hMuonDisVsPt[2][nHistos][2];
  const char *type[2] = {"US","LS"};
  const char *name[nHistos] = {"Dz","Dy","Dtof","NsigmaPi"};
  const char *charge[2] = {"neg","pos"};

  TFile *fin = 0x0;
  if(year==2013) fin = TFile::Open("output/Pico.Run13.pp500.jpsi.muon.root","read");

  THnSparseF *hnJpsiMuon[nHistos];
  for(int j=0; j<nHistos; j++)
    {
      hnJpsiMuon[j] = (THnSparseF*)fin->Get(Form("mhMuon%s_di_mu",name[j]));
      hnJpsiMuon[j]->GetAxis(0)->SetRangeUser(3.0,3.2);
      for(int i=0; i<2; i++)
	{
	  hnJpsiMuon[j]->GetAxis(1)->SetRange(i+1,i+1);
	  for(int k=0; k<2; k++)
	    {
	      for(int m=0; m<2; m++)
		{
		  hnJpsiMuon[j]->GetAxis(6+m)->SetRange(1+k*2, 1+k*2);
		  TH2F* h2 = (TH2F*)hnJpsiMuon[j]->Projection(4+m,2+m);
		  h2->SetName(Form("Data_JpsiMuon_%sVsPt_%s_%s_muon%d",name[j],charge[k],type[i],m+1));
		  if(m==0) hMuonDisVsPt[i][j][k] = (TH2F*)h2->Clone(Form("Data_JpsiMuon_%sVsPt_%s_%s",name[j],charge[k],type[i]));
		  else     hMuonDisVsPt[i][j][k]->Add(h2);
		  hnJpsiMuon[j]->GetAxis(6+m)->SetRange(0,-1);
		}
	    }
	  hnJpsiMuon[j]->GetAxis(1)->SetRange(0,-1);
	}   
    }

  TH2F *hDis[nHistos][3];
  for(int j=0; j<nHistos; j++)
    {
      for(int k=0; k<2; k++)
	{
	  hDis[j][k] = (TH2F*)hMuonDisVsPt[0][j][k]->Clone(Form("Data_JpsiMuon_%sVsPt_%s",name[j],charge[k]));
	  hDis[j][k]->Add(hMuonDisVsPt[1][j][k],-1);
	}
      hDis[j][2] = (TH2F*)hDis[j][0]->Clone(Form("Data_JpsiMuon_%sVsPt",name[j]));
      hDis[j][2]->Add((TH2F*)hDis[j][1]);
    }

  TFile *fout = 0x0;
  if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.PIDcuts.root","update");
  for(int j=0; j<nHistos; j++)
    {
      for(int k=2; k>-1; k--)
	{
	  hDis[j][k]->SetTitle("");
	  hDis[j][k]->Write("",TObject::kOverwrite);
	  if(k<2)
	    {
	      for(int i=0; i<2; i++)
		{
		  hMuonDisVsPt[i][j][k]->SetTitle("");
		  hMuonDisVsPt[i][j][k]->Write("",TObject::kOverwrite);
		}
	    }	    
	}
    }
}
