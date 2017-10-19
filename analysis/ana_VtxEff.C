TFile *f;
const int year = YEAR;
const char *setupName[3] = {"prod_low","prod_mid","prod_high"};

//================================================
void ana_VtxEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  f = TFile::Open(Form("./output/Run14.AuAu200.MB.VtxEff.root"),"read");

  anaDiffVzEff();
}

//================================================
void anaDiffVzEff(int savePlot = 1)
{
  gStyle->SetOptFit(0);
  THnSparseF *hn[2];
  const char* name[2] = {"mb","mtd"};
  TList *list = new TList;
  TString trigName[2] = {"VPD-ZDC-novtx-mon", "dimuon"};

  hn[0] = (THnSparseF*)f->Get("mhMbEvtEff");
  TFile *fmtd = TFile::Open("output/Run14_AuAu200.MTD.VtxEff.root", "read");
  THnSparseF *hn[1] = (THnSparseF*)fmtd->Get("mhEventCent_qa_di_mu");

  /*
  // dz vs vz in different centrality bins
  TH2F *hDzVsTpcVz[2][16];
  TCanvas *cDzVsTpcVz[2];
  for(int j=0; j<2; j++)
  {
  cDzVsTpcVz[j] = new TCanvas(Form("cDzVsTpcVz_%s",name[j]), Form("cDzVsTpcVz_%s",name[j]), 1100, 700);
  cDzVsTpcVz[j]->Divide(3,3);
  for(int i=0; i<9; i++)
  {
  if(j==0) 
  {
  hn[j]->GetAxis(5)->SetRange(i+1, i+1);
  hDzVsTpcVz[j][i] = (TH2F*)hn[j]->Projection(2, 1);
  hn[j]->GetAxis(5)->SetRange(0, -1);
  }
  else 
  {
  hn[j]->GetAxis(3)->SetRange(i+1, i+1);
  hDzVsTpcVz[j][i] = (TH2F*)hn[j]->Projection(0, 1);
  hn[j]->GetAxis(3)->SetRange(0, -1);
  }
  hDzVsTpcVz[j][i]->SetName(Form("hDzVsTpcVz_cent%d_%s",i,name[j]));
  if(i<9)
  {
  cDzVsTpcVz[j]->cd(i+1);
  gPad->SetLogz();
  hDzVsTpcVz[j][i]->SetTitle("");
  hDzVsTpcVz[j][i]->Draw("colz");
  if(j==0) t1 = GetTitleText(Form("VPD-ZDC-novtx-mon: %d-%d%%",75-i*5,80-i*5),0.08);
  if(j==1) t1 = GetTitleText(Form("dimuon: %d-%d%%",75-i*5,80-i*5),0.08);
  t1->Draw();
  }
  }
  if(savePlot)  cDzVsTpcVz[j]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/%s_VzDiffVsTpcVz.png",run_type,name[j]));
  }
  */

  // dz vs. luminosity
  TCanvas *c = new TCanvas("VzDiff_7580","VzDiff_7580",1100,600);
  c->Divide(2,1);
  hn[1]->GetAxis(3)->SetRange(1,1);
  TH2F *hVzDiffVsZdc = (TH2F*)hn[1]->Projection(0,4);
  TH2F *hVzDiffVsTofMult = (TH2F*)hn[1]->Projection(0,5);
  hVzDiffVsTofMult->SetName("hVzDiffVsTofMult");
  hn[1]->GetAxis(3)->SetRange(0,-1);
  TFile *flumi = TFile::Open("Rootfiles/Run14_AuAu200.Luminosity.root","read");
  TH1F *hDzDiff7580[2][3];
  for(int j=0; j<2; j++)
    {
      c->cd(j+1);
      TLegend *leg = new TLegend(0.15, 0.73, 0.35, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.035);
      leg->SetHeader(run_type);
      for(int i=0; i<3; i++)
	{
	  if(j==0) hDzDiff7580[j][i] = (TH1F*)hVzDiffVsZdc->ProjectionY(Form("hDzDiff7580_zdc%d",i+1),i+2,i+2);
	  else     hDzDiff7580[j][i] = (TH1F*)flumi->Get(Form("TpcVpdDz_cent7580_%s",setupName[i]));
	  hDzDiff7580[j][i]->Scale(1./hDzDiff7580[j][i]->GetBinContent(hDzDiff7580[j][i]->FindFixBin(0)));
	  hDzDiff7580[j][i]->GetXaxis()->SetRangeUser(-10,10);
	  hDzDiff7580[j][i]->SetMaximum(1.3);
	  hDzDiff7580[j][i]->SetLineColor(color[i]);
	  hDzDiff7580[j][i]->SetTitle("");
	  if(i==0) hDzDiff7580[j][i]->Draw();
	  else     hDzDiff7580[j][i]->Draw("samesHIST");
	  leg->AddEntry(hDzDiff7580[j][i],Form("%s",setupName[i]),"L");
	}
      t1 = GetTitleText(Form("%s: #Deltavz distribution",trigName[j].Data(),0.04));
      t1->Draw();
      leg->Draw();
    }
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/Mtd_VzDiffInZdc.pdf",run_type));

  c = draw2D(hVzDiffVsTofMult,"dimuon: #Deltavz vs. tofMult (75-80%);TofMult");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/Mtd_VzDiffVsTofMult.pdf",run_type));

  TH2F *hTofMultVsCent = (TH2F*)hn[1]->Projection(5, 3);
  hTofMultVsCent->SetName("hTofMultVsCent");
  for(int i=0; i<16; i++)
    {
      hTofMultVsCent->GetXaxis()->SetBinLabel(i+1, Form("%d-%d%%",75-i*5,80-i*5));
    }
  hTofMultVsCent->GetXaxis()->SetRangeUser(0,16);
  c = draw2D(hTofMultVsCent,"dimuon: tofMult vs. centrality;;tofMult");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/Mtd_TofMultVsCent.pdf",run_type));

  TH1F *hMtdVpdVz7580[3];
  TH1F *hMtdTpcVz7580[3];
  TH1F *hMtdVzDiff7580[3];
  hn[1]->GetAxis(3)->SetRange(1,1);
  const int tofMult_min[3] = {0, 0, 500};
  const int tofMult_max[3] = {3000, 500, 3000};
  for(int i=0; i<3; i++)
    {
      hn[1]->GetAxis(5)->SetRangeUser(tofMult_min[i], tofMult_max[i]);
      hMtdVzDiff7580[i] = (TH1F*)hn[1]->Projection(0);
      hMtdVzDiff7580[i]->SetName(Form("hMtdVzDiff7580_%d_%d",tofMult_min[i], tofMult_max[i]));

      hMtdTpcVz7580[i] = (TH1F*)hn[1]->Projection(1);
      hMtdTpcVz7580[i]->SetName(Form("hMtdTpcVz7580_%d_%d",tofMult_min[i], tofMult_max[i]));

      hMtdVpdVz7580[i] = (TH1F*)hn[1]->Projection(2);
      hMtdVpdVz7580[i]->SetName(Form("hMtdVpdVz7580_%d_%d",tofMult_min[i], tofMult_max[i]));
    }
  hn[1]->GetAxis(5)->SetRange(0,-1);
  hn[1]->GetAxis(3)->SetRange(0,-1);

  TString legName[3];
  for(int i=0; i<3; i++)
    {
      list->Add(hMtdTpcVz7580[i]);
      legName[i] = Form("%d < TofMult < %d",tofMult_min[i], tofMult_max[i]);
    }
  c = drawHistos(list,"hMtdTpcVz7580","dimuon: TPC vz distribution (75-80%)",false,0,120,false,0.2,1.1,false,true,legName,true,run_type,0.6,0.8,0.6,0.8,false,0.04,0.035);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/Mtd_TpcVzInTofMult.pdf",run_type));
  list->Clear();

  for(int i=0; i<3; i++)
    {
      list->Add(hMtdVpdVz7580[i]);
    }
  c = drawHistos(list,"hMtdVpdVz7580","dimuon: VPD vz distribution (75-80%)",false,0,120,false,0.2,1.1,false,true,legName,true,run_type,0.6,0.8,0.6,0.8,false,0.04,0.035);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/Mtd_VpdVzInTofMult.pdf",run_type));
  list->Clear();

  for(int i=0; i<3; i++)
    {
      list->Add(hMtdVzDiff7580[i]);
    }
  c = drawHistos(list,"hMtdVzDiff7580","dimuon: #Deltavz distribution (75-80%)",true,-10,10,false,0.2,1.1,false,true,legName,true,run_type,0.6,0.8,0.6,0.8,false,0.04,0.035);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/Mtd_VzDiffInTofMult.pdf",run_type));
  list->Clear();

  return;

  // dz distribution in different centrality bins
  TH2F *hVzDiffVsCent[2];
  hn[0]->GetAxis(1)->SetRangeUser(-100, 100);
  hn[0]->GetAxis(3)->SetRangeUser(0,2);
  hVzDiffVsCent[0] = (TH2F*)hn[0]->Projection(2,5);
  hn[0]->GetAxis(1)->SetRange(0,-1);
  hn[0]->GetAxis(3)->SetRange(0,-1);

  hn[1]->GetAxis(1)->SetRangeUser(-100, 100);
  hVzDiffVsCent[1] = (TH2F*)hn[1]->Projection(0,3);
  hn[1]->GetAxis(1)->SetRange(0,-1);
  TH1F *hVzDiff[2][16];
  TF1  *funcVzDiff[2][16];
  TH1F *hVzDiffEff[2];
  TH1F *hVzDiffSigma[2];
  TCanvas *cDz[2];
  for(int j=0; j<2; j++)
    {
      cDz[j] = new TCanvas(Form("VzDiff_%s",name[j]), Form("VzDiff_%s",name[j]), 1100, 700);
      cDz[j]->Divide(4,4);
      hVzDiffEff[j] = new TH1F(Form("hVzDiffEff_%s",name[j]),Form("hVzDiffEff_%s",name[j]),16,0,16);
      hVzDiffSigma[j] = new TH1F(Form("hVzDiffSigma_%s",name[j]),Form("hVzDiffSigma_%s",name[j]),16,0,16);
      for(int i=0; i<16; i++)
	{
	  hVzDiff[j][i] = (TH1F*)hVzDiffVsCent[j]->ProjectionY(Form("hVzDiff_cent%d_%s",i,name[j]),i+1, i+1);
	  hVzDiffEff[j]->GetXaxis()->SetBinLabel(i+1, Form("%d-%d%%",75-i*5,80-i*5));
	  hVzDiffEff[j]->SetBinContent(i+1, hVzDiff[j][i]->Integral(hVzDiff[j][i]->FindFixBin(-3),hVzDiff[j][i]->FindFixBin(3))/hVzDiff[j][i]->Integral(0,-1));
	  hVzDiffEff[j]->SetBinError(i+1, 1e-10);
	  printf("%s: %d-%d%%: %4.2f\n",name[j],75-i*5,80-i*5,hVzDiff[j][i]->Integral(hVzDiff[j][i]->FindFixBin(-3),hVzDiff[j][i]->FindFixBin(3))/hVzDiff[j][i]->Integral(0,-1));
	  cDz[j]->cd(i+1);
	  hVzDiff[j][i]->GetXaxis()->SetRangeUser(-10,10);
	  funcVzDiff[j][i] = new TF1(Form("funcVzDiff_cent%d_%s",i,name[j]), "gaus", -3, 3);
	  hVzDiff[j][i]->Fit(funcVzDiff[j][i], "IR0Q");
	  double mean = funcVzDiff[j][i]->GetParameter(1);
	  double sigma = funcVzDiff[j][i]->GetParameter(2);
	  double n = 2.5;
	  funcVzDiff[j][i]->SetRange(mean-n*sigma, mean+n*sigma);
	  hVzDiff[j][i]->Fit(funcVzDiff[j][i], "IR0Q");
	  hVzDiff[j][i]->SetTitle("");
	  hVzDiff[j][i]->Draw();
	  funcVzDiff[j][i]->SetLineStyle(2);
	  funcVzDiff[j][i]->SetLineColor(2);
	  funcVzDiff[j][i]->Draw("sames");
	  hVzDiffSigma[j]->GetXaxis()->SetBinLabel(i+1, Form("%d-%d%%",75-i*5,80-i*5));
	  hVzDiffSigma[j]->SetBinContent(i+1, funcVzDiff[j][i]->GetParameter(2));
	  hVzDiffSigma[j]->SetBinError(i+1, funcVzDiff[j][i]->GetParError(2));
	  if(j==0) t1 = GetTitleText(Form("VPD-ZDC-novtx-mon: %d-%d%%",75-i*5,80-i*5),0.08);
	  if(j==1) t1 = GetTitleText(Form("dimuon: %d-%d%%",75-i*5,80-i*5),0.08);
	  t1->Draw();
	  t1 = GetPaveText(0.15,0.35,0.6,0.7,0.08);
	  t1->AddText(Form("#sigma = %2.2fcm",funcVzDiff[j][i]->GetParameter(2)));
	  t1->Draw();
	}
      if(savePlot)  cDz[j]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/%s_FitVzDiffDis.pdf",run_type,name[j]));
    }

  for(int j=0; j<2; j++)
    {
      list->Add(hVzDiffEff[j]);
    }
  c = drawHistos(list,"VzDiffEff",Form("%s: efficiency of |#Deltavz| < 3 cm cut",run_type),false,0,120,true,0.2,1.1,false,true,trigName,true,"Bin-counting",0.5,0.7,0.5,0.65);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/MtdVsMb_VzDiffEff.pdf",run_type));
  list->Clear();

  for(int j=0; j<2; j++)
    {
      list->Add(hVzDiffSigma[j]);
    }
  c = drawHistos(list,"VzDiffSigma",Form("%s: width of #Deltavz distribution",run_type),false,0,120,true,0.3,1.7,false,true,trigName,true,"",0.5,0.7,0.5,0.65);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/MtdVsMb_VzDiffSigma.pdf",run_type));
  list->Clear();
}
