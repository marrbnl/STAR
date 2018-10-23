#include </Users/admin/Work/STAR/util/defs.h>
#include "/Users/admin/Work/STAR/util/drawHistos.C"

//******************* function definitions ************************
const double muMass = 0.1057;
const int year = 2014;
TString part_name = "Jpsi";
TRandom3 *myRandom;

double rotatePhi(double phi);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);
void getJpsiEff(const int nExpr, const int mass, TH1F *hJpsiTruth, TH2F *hTrkEtaVsPhiEffInPtBin[nCentBins_pt][gNTrgSetup], TH1F *hMcJpsi, TH1F *hJpsiCorrInPtBin[nCentBins_pt][gNTrgSetup]);

void anaTpcEff(const int isys = 0, const int savePlot = 0, const int saveHisto = 0);
void makeHistos(const int isys = 0, const int saveHisto = 0);

TH3F *hMuonEffPos;
TH3F *hMuonEffNeg;
TH1F *hInputJpsiPt;
TH1F *hInputJpsiEta;

//================================================
void toyMC_LumiDepEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  myRandom = new TRandom3();
  TDatime *clock = new TDatime();
  myRandom->SetSeed(clock->GetTime());

  TFile *feff = TFile::Open("Rootfiles/Run14_AuAu200.MbEmbTrkEff.root", "read");
  hMuonEffNeg = (TH3F*)feff->Get("McTrkPtEtaPhiFinal_neg_0080_Eff");
  hMuonEffPos = (TH3F*)feff->Get("McTrkPtEtaPhiFinal_pos_0080_Eff");

  for(int isys=0; isys<1; isys++)
    {
      //makeHistos(isys, 1);
      anaTpcEff(isys,1,1);
    }  

}

//================================================
void makeHistos(const int isys, const int saveHisto)
{
  // Get eta vs. phi distribution from data
  TFile *fdata = TFile::Open(Form("output/Run14_AuAu200.MB.P15ic.Study.TpcTracking.root"), "read");
  THnSparseF *hnDataTrkEtaPhi = (THnSparseF*)fdata->Get("mhnTrkEtaPhi_mb");
  if(isys==0) hnDataTrkEtaPhi->GetAxis(8)->SetRangeUser(-1*30, 30);
  hnDataTrkEtaPhi->GetAxis(7)->SetRange(1,1); // DCA < 1 cm
  hnDataTrkEtaPhi->GetAxis(0)->SetRangeUser(1.5,5); // pt > 1.5 GeV/c
  TH2F *hDataTrkEtaVsPhi[16][12];
  for(int j=0; j<16; j++)
    {
      hnDataTrkEtaPhi->GetAxis(3)->SetRange(j+1, j+1);
      TH3F* h3tmp = (TH3F*)hnDataTrkEtaPhi->Projection(2,1,6);
      h3tmp->SetName(Form("%s_%d",h3tmp->GetName(),j));
      for(int z=0; z<12; z++)
	{
	  h3tmp->GetZaxis()->SetRange(z+1, z+1);
	  hDataTrkEtaVsPhi[j][z] = (TH2F*)h3tmp->Project3D("yx");
	  hDataTrkEtaVsPhi[j][z]->SetName(Form("hDataTrkEtaVsPhi_Cent%d_Zdc%d_Sys%d",j,z,isys));
	  hDataTrkEtaVsPhi[j][z]->Sumw2();
	  if(hDataTrkEtaVsPhi[j][z]->Integral()>0)
	    hDataTrkEtaVsPhi[j][z]->Scale(1./hDataTrkEtaVsPhi[j][z]->Integral());
	}
    }

  hnDataTrkEtaPhi->GetAxis(3)->SetRange(15,16); // centrality
  hnDataTrkEtaPhi->GetAxis(6)->SetRangeUser(20,40); // Zdc
  TH2F *hRefTrkEtaVsPhi = (TH2F*)hnDataTrkEtaPhi->Projection(1,2);
  hRefTrkEtaVsPhi->SetName("hRefTrkEtaVsPhi_Cent0010_Zdc2040");
  hnDataTrkEtaPhi->GetAxis(3)->SetRange(0,-1);
  hnDataTrkEtaPhi->GetAxis(6)->SetRange(0,-1);

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Sys.LumiDepEff.root",run_type.Data()), "update");
      hRefTrkEtaVsPhi->SetTitle("");
      hRefTrkEtaVsPhi->Write("",TObject::kOverwrite);
      for(int j=0; j<16; j++)
	{
	  for(int z=0; z<12; z++)
	    {
	      hDataTrkEtaVsPhi[j][z]->SetTitle("");
	      hDataTrkEtaVsPhi[j][z]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void anaTpcEff(const int isys, const int savePlot, const int saveHisto)
{
  // this function deterimes luminsoity dependent 
  // additional TPC inefficiency 

 // Get weights from MB events
  TFile *fmb = TFile::Open("output/Run14_AuAu200.MB.P16id.VtxEff.root", "read");
  TH3F	*hEqMbEvtWeight = (TH3F*)fmb->Get("mhEqMbEvtWeight"); //Centrality vs. ZdcRate vs. TrgSetup
  hEqMbEvtWeight->Rebin3D(1,2,1);
  const double ncoll[16] = {10.54, 15.98, 23.96, 35.64, 51.62, 72.79, 101.51, 138.71, 186.69, 246.95, 320.78, 411.86, 524.31, 663.03, 838.41, 1048.11}; 
  const double raa[16] = {1, 1, 1, 0.87, 0.75, 0.68, 0.63, 0.60, 0.58, 0.56, 0.55, 0.51, 0.48, 0.44, 0.42, 0.41};

  TH2F *hEqMbEvt[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hEqMbEvtWeight->GetZaxis()->SetRange(k,k);
      hEqMbEvt[k] = (TH2F*)hEqMbEvtWeight->Project3D("yx");
      hEqMbEvt[k]->SetName(Form("hEqMbEvt%s",gTrgSetupTitle[k]));
      int nxbins = hEqMbEvt[k]->GetNbinsX();
      int nybins = hEqMbEvt[k]->GetNbinsY();
      for(int ibin=1; ibin<=nxbins; ibin++)
	{
	  double scale = ncoll[ibin-1]*raa[ibin-1];
	  for(int jbin=1; jbin<=nybins; jbin++)
	    {
	      hEqMbEvt[k]->SetBinContent(ibin, jbin, hEqMbEvt[k]->GetBinContent(ibin,jbin)*scale);
	      hEqMbEvt[k]->SetBinError(ibin, jbin, hEqMbEvt[k]->GetBinError(ibin, jbin)*scale);
	    }
	}
      printf("[i] %s has %1.0f counts\n",hEqMbEvt[k]->GetName(),hEqMbEvt[k]->Integral());
      hEqMbEvt[k]->Scale(1./hEqMbEvt[k]->Integral());
    }

  // Get eta vs. phi distribution from data
  TFile *fout = NULL;
  if(saveHisto) fout = TFile::Open(Form("Rootfiles/%s.Sys.LumiDepEff.root",run_type.Data()), "update");
  else          fout = TFile::Open(Form("Rootfiles/%s.Sys.LumiDepEff.root",run_type.Data()), "read");
  TH2F *hDataTrkEtaVsPhi[16][12];
  for(int j=0; j<16; j++)
    {
      for(int z=0; z<12; z++)
	{
	  hDataTrkEtaVsPhi[j][z] = (TH2F*)fout->Get(Form("hDataTrkEtaVsPhi_Cent%d_Zdc%d_Sys%d",j,z,isys));
	}
    }

  // Obtain track eta vs. phi in pt and npart bins
  TCanvas *c = new TCanvas("cEtaVsPhiInPtBin", "cEtaVsPhiInPtBin", 1200, 800);
  c->Divide(5,4);
  TH2F *hTrkEtaVsPhiInPtBin[nCentBins_pt][gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int icent=0; icent<nCentBins_pt; icent++)
	{
	  hTrkEtaVsPhiInPtBin[icent][k] = (TH2F*)hDataTrkEtaVsPhi[0][0]->Clone(Form("hTrkEtaVsPhi_cent%s%s",cent_Title_pt[icent],gTrgSetupTitle[k]));
	  hTrkEtaVsPhiInPtBin[icent][k]->Reset();
	  double total_w = 0;
	  for(int bin=centBins_low_pt[icent]; bin<=centBins_high_pt[icent]; bin++)
	    {
	      for(int z=0; z<12; z++)
		{
		  double weight = hEqMbEvt[k]->GetBinContent(bin, z+1);
		  hTrkEtaVsPhiInPtBin[icent][k]->Add(hDataTrkEtaVsPhi[bin-1][z], weight);
		  total_w += weight;
		}
	    }
	  hTrkEtaVsPhiInPtBin[icent][k]->Rebin2D(1,1);
	  hTrkEtaVsPhiInPtBin[icent][k]->Scale(1./total_w);
	  hTrkEtaVsPhiInPtBin[icent][k]->SetTitle("");

	  if(icent==0) continue;
	  c->cd((icent-1)*5+k+1);
	  hTrkEtaVsPhiInPtBin[icent][k]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("Data: %s%%, %s",cent_Name_pt[icent],gLegNameTrg[k].Data()),0.08);
	  t1->Draw();
	}  
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_DataTrkEtaVsPhi.pdf",run_type.Data()));

  TH2F *hTrkEtaVsPhiInNpartBin[nCentBins_npart[0]];

  // get reference track eta vs phi 
  TH2F *hRefTrkEtaVsPhi = (TH2F*)fout->Get("hRefTrkEtaVsPhi_Cent0010_Zdc2040");
  c = draw2D(hRefTrkEtaVsPhi);

  c = new TCanvas("cEtaVsPhiEffInPtBin", "cEtaVsPhiEffInPtBin", 1100, 700);
  c->Divide(5,4);
  TH2F *hTrkEtaVsPhiEffInPtBin[nCentBins_pt][gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int icent=0; icent<nCentBins_pt; icent++)
	{
	  hTrkEtaVsPhiEffInPtBin[icent][k] = (TH2F*)hTrkEtaVsPhiInPtBin[icent][k]->Clone(Form("hTrkEtaVsPhiEff_cent%s%s",cent_Title_pt[icent],gTrgSetupTitle[k]));
	  int nbinsx = hTrkEtaVsPhiEffInPtBin[icent][k]->GetNbinsX();
	  int nbinsy = hTrkEtaVsPhiEffInPtBin[icent][k]->GetNbinsY();
	  double scale = 0;
	  int nbins = 0;
	  for(int ibin=1; ibin<=nbinsx; ibin++)
	    {
	      for(int jbin=1; jbin<=nbinsy; jbin++)
		{
		  double data = hTrkEtaVsPhiEffInPtBin[icent][k]->GetBinContent(ibin, jbin);
		  double ref = hRefTrkEtaVsPhi->GetBinContent(ibin, jbin);
		  if(data<=0) continue;
		  hTrkEtaVsPhiEffInPtBin[icent][k]->SetBinContent(ibin, jbin, data/ref);
		  hTrkEtaVsPhiEffInPtBin[icent][k]->SetBinError(ibin, jbin, hTrkEtaVsPhiEffInPtBin[icent][k]->GetBinError(ibin,jbin)/ref);
		  if((ibin-1)%3!=0 && jbin>=15 && jbin<=18) 
		    {
		      scale += data/ref;
		      nbins++;
		    }
		}
	    }
	  hTrkEtaVsPhiEffInPtBin[icent][k]->Scale(1./(scale/nbins));
	  hTrkEtaVsPhiEffInPtBin[icent][k]->GetZaxis()->SetRangeUser(0,1.05);
	  if(icent==0) continue;
	  c->cd((icent-1)*5+k+1);
	  hTrkEtaVsPhiEffInPtBin[icent][k]->Draw("colz");
	  TPaveText *t1 = GetTitleText(Form("Data: %s%%, %s",cent_Name_pt[icent],gLegNameTrg[k].Data()),0.08);
	  t1->Draw();
	}
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_DataTrkEffEtaVsPhi.pdf",run_type.Data()));

  // ToyMC
  double jpsiMass = 0;
  int nPtBins = 0;
  double xPtBins[10] = {0};
  TString name = part_name;

  TH1F *hMcJpsiPt = NULL;
  if(name.Contains("Jpsi"))
    {
      jpsiMass = 3.096;
      if(year==2013)
	{
	  nPtBins = 5;
	  double xbins_tmp[6] = {0,1.5,3,5,7,9};
	  std::copy(std::begin(xbins_tmp), std::end(xbins_tmp), std::begin(xPtBins));
	}
      else if(year==2014)
	{
	  nPtBins = 9;
	  double xbins_tmp[10] = {0,1,2,3,4,5,6,8,10,15};
	  std::copy(std::begin(xbins_tmp), std::end(xbins_tmp), std::begin(xPtBins));
	}
      else if(year==2015)
	{
	  nPtBins = 8;
	  double xbins_tmp[9] = {0,1,2,3,4,5,6,8,10};
	  std::copy(std::begin(xbins_tmp), std::end(xbins_tmp), std::begin(xPtBins));
	}
      else if(year==2016)
	{
	  nPtBins = 6;
	  double xbins_tmp[7] = {0,1,2,3,4,6,10};
	  std::copy(std::begin(xbins_tmp), std::end(xbins_tmp), std::begin(xPtBins));
	}

      if(year==2013)
	{
	  TFile *fWeight = TFile::Open("Rootfiles/pp500GeVfit_new.root","read");
	  TF1 *funcJpsi = (TF1*)fWeight->Get("ffpt");
	  funcJpsi->SetNpx(1000);
	  hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
	  hMcJpsiPt->SetName(Form("GlobalFit_Jpsi_Yield_cent00100"));
	  for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
	    {
	      hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
	    }
	}
      else if(year==2014 || year==2016)
	{
          TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
          hMcJpsiPt = (TH1F*)fWeight->Get(Form("TBW_JpsiYield_AuAu200_cent0060"));
	}
      else if(year==2015)
	{
	  TFile *fWeight = TFile::Open("Rootfiles/JpsiSpectraShapepp200.root","read");
	  TF1 *funcJpsi = (TF1*)fWeight->Get("TsallisPowerLawFitJpsipp200");
	  funcJpsi->SetNpx(1000);
	  hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
	  hMcJpsiPt->SetName(Form("Jpsi_Yield_cent00100"));
	  for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
	    {
	      hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
	    }
	}
    }
  else if(name.Contains("Ups"))
    {
      jpsiMass = 9.46;
      nPtBins = 3;
      double xbins_tmp[4] = {0,2,4,10};
      std::copy(std::begin(xbins_tmp), std::end(xbins_tmp), std::begin(xPtBins));

      TF1 *fBol = new TF1("Boltzmann","x/(exp(x/[0]+1))",0,10);
      fBol->SetParameter(0,1.11);
      fBol->SetNpx(1000);
      hMcJpsiPt = (TH1F*)fBol->GetHistogram();

    }
  hMcJpsiPt->Scale(1./hMcJpsiPt->Integral()); // turn into PDF
  TH1F *hJpsiCorrInPtBin[nCentBins_pt][gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      for(int icent=0; icent<nCentBins_pt; icent++)
	{
	  hJpsiCorrInPtBin[icent][k] = new TH1F(Form("hJpsiCorrInPtBin_cent%s%s_Sys%d",cent_Title_pt[icent],gTrgSetupTitle[k],isys), "", nPtBins, xPtBins);
	  hJpsiCorrInPtBin[icent][k]->Sumw2();
	}
    }
  TH1F *hMcJpsi = new TH1F("hMcJpsi", "", nPtBins, xPtBins);

  const int nExpr = 5e6;
  getJpsiEff(nExpr, jpsiMass, hMcJpsiPt, hTrkEtaVsPhiEffInPtBin, hMcJpsi, hJpsiCorrInPtBin);

  // plot correction factor
  TH2F *hJpsiCorrInPt = new TH2F(Form("hJpsiCorrInPtBin_Sys%d",isys), "", nCentBins_pt, 0, nCentBins_pt, gNTrgSetup, 0, gNTrgSetup);
  for(int k=0; k<gNTrgSetup; k++)
    {
      hJpsiCorrInPt->GetYaxis()->SetBinLabel(gNTrgSetup-k, gLegNameTrg[k].Data());
      for(int icent=0; icent<nCentBins_pt; icent++)
	{
	  if(k==0) hJpsiCorrInPt->GetXaxis()->SetBinLabel(icent+1, Form("%s%%",cent_Name_pt[icent]));
	  hJpsiCorrInPt->SetBinContent(icent+1, gNTrgSetup-k, hJpsiCorrInPtBin[icent][k]->GetEntries()/hMcJpsi->GetEntries());
	  printf("[i] %s has %1.0f/%1.0f = %4.2f%%\n",hJpsiCorrInPtBin[icent][k]->GetName(),hJpsiCorrInPtBin[icent][k]->GetEntries(),hMcJpsi->GetEntries(),hJpsiCorrInPtBin[icent][k]->GetEntries()/hMcJpsi->GetEntries()*100);
	}
    }
  gStyle->SetPaintTextFormat("2.3f");
  hJpsiCorrInPt->SetMarkerSize(1.5);
  c = draw2D(hJpsiCorrInPt, Form("%s: additional TPC tracking efficiency correction",run_type.Data()), 0.04, false, "colzTEXT");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_JpsiEffInPtBin.pdf",run_type.Data()));

  if(saveHisto)
    {
      fout->cd();
      hJpsiCorrInPt->Write("",TObject::kOverwrite);
    }
}

//================================================
void getJpsiEff(const int nExpr, const int mass, TH1F *hJpsiTruth, TH2F *hTrkEtaVsPhiEffInPtBin[nCentBins_pt][gNTrgSetup], TH1F *hMcJpsi, TH1F *hJpsiCorrInPtBin[nCentBins_pt][gNTrgSetup])
{
  for(int i=0; i<nExpr; i++)
    {
      //printf("+++ Event %d +++\n",i+1);
      double mc_pt = myRandom->Uniform(0,10);
      double weight = hJpsiTruth->GetBinContent(hJpsiTruth->FindFixBin(mc_pt));
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);

      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = rotatePhi(daughter1.Phi());
      
      TLorentzVector daughter2 = parent - daughter1;
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = rotatePhi(daughter2.Phi());

      if(pt1<pt1_cut && pt2<pt1_cut) continue;
      if(pt1<pt2_cut || pt2<pt2_cut) continue;
      if(fabs(eta1)>0.6 || fabs(eta2)>0.6) continue;

      double prob1 = myRandom->Uniform(0,1);
      double prob2 = myRandom->Uniform(0,1);
      double eff1 = hMuonEffPos->GetBinContent(hMuonEffPos->GetXaxis()->FindBin(pt1),
					       hMuonEffPos->GetYaxis()->FindBin(eta1),
					       hMuonEffPos->GetZaxis()->FindBin(phi1));
      double eff2 = hMuonEffNeg->GetBinContent(hMuonEffNeg->GetXaxis()->FindBin(pt2),
					       hMuonEffNeg->GetYaxis()->FindBin(eta2),
					       hMuonEffNeg->GetZaxis()->FindBin(phi2));
      if(prob1>eff1 || prob2>eff2) continue;
      hMcJpsi->Fill(mc_pt);
      double prob11 = myRandom->Uniform(0,1);
      double prob21 = myRandom->Uniform(0,1);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  for(int icent=0; icent<nCentBins_pt; icent++)
	    {
	      double eff11 = hTrkEtaVsPhiEffInPtBin[icent][k]->GetBinContent(hTrkEtaVsPhiEffInPtBin[icent][k]->GetXaxis()->FindBin(phi1),
									     hTrkEtaVsPhiEffInPtBin[icent][k]->GetYaxis()->FindBin(eta1));
	      double eff21 = hTrkEtaVsPhiEffInPtBin[icent][k]->GetBinContent(hTrkEtaVsPhiEffInPtBin[icent][k]->GetXaxis()->FindBin(phi2),
									     hTrkEtaVsPhiEffInPtBin[icent][k]->GetYaxis()->FindBin(eta2));
	      if(prob11<eff11 && prob21<eff21)
		{
		  hJpsiCorrInPtBin[icent][k]->Fill(mc_pt);
		}
	    }
	}
    }
}

//_____________________________________________________________________________
double rotatePhi(double phi) 
{
  double outPhi = phi;
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
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
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass) 
{
  double e = parent.M()/2.;
  double p = sqrt(e*e-dmass*dmass);
  double costheta = myRandom->Uniform(-1.0,1.0);
  double phi = myRandom->Uniform(0,TMath::Pi()*2);
  double pz = p*costheta;
  double px = p*sqrt(1.-costheta*costheta)*cos(phi);
  double py = p*sqrt(1.-costheta*costheta)*sin(phi);
  TLorentzVector daughter(px,py,pz,e);
  return myBoost(parent,daughter);
}
