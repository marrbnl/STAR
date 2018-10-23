#include </Users/admin/Work/STAR/util/defs.h>
#include "/Users/admin/Work/STAR/util/drawHistos.C"

//******************* function definitions ************************
const double muMass = 0.1057;
const int year = 2014;
const int nPart = 4;
const char *part_name[nPart]  = {"Jpsi","Ups1S","Ups2S","Ups3S"};
const char *part_title[nPart] = {"Jpsi","Y(1S)","Y(2S)","Y(3S)"};
const double part_mass[nPart] = {3.097, 9.46, 10.023, 10.355};

TRandom3 *myRandom;
void anaSys(const char* type, const int saveHisto = 0);
void plotSys(const char* type, const int savePlot = 0, const int saveHisto = 0, const int saveAN = 0);
void mtdTrigEff(const int saveHisto = 0);
void mtdTrigEffLS(const int saveHisto = 0);
void mtdRespEff(const int saveHisto = 0);
void makeHisto(TString name, const double mass, const int nExpr = 1e4, const int nHistos = 3);
void toyMC(const double mass, const int nExpr, const int nHistos = 3, const int debug = 0);

TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);

TH2F *hTpcTrackRes;
TH1F *hTrkResBin[400];
TF1 *hMuonPtEff[5];
TH1F *hMcJpsiPt;
TH1F *hInJpsiPt;
TH1F *hInJpsiCent;
TH1F *hOutJpsiPt[5];
TH1F *hOutJpsiCent[5];

//================================================
void sys_ToyMC()
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

  const char* type = "MtdTrigEff";
  //const char* type = "Dtof";

  //anaSys(type, 1);
  //plotSys(type, 1, 1, 0);


  //const char* type = "MtdMthEff";
  //mtdTrigEff(1);
  mtdTrigEffLS(1);

  //const char* type = "MtdRespEff";
  //mtdRespEff(1);
}

//================================================
void plotSys(const char* type, const int savePlot, const int saveHisto, const int saveAN)
{
  char name_lumi[256], title_lumi[256], name_eff[256], title_eff[256];
  if(strcmp(type,"MtdTrigEff")==0)
    {
      sprintf(name_lumi, "");
      sprintf(title_lumi, "");
      sprintf(name_eff, "TacDiffEff");
      sprintf(title_eff, "MTD trigger efficiency");
    }

  else if(strcmp(type,"Dtof")==0)
    {
      sprintf(name_lumi, "");
      sprintf(title_lumi, "");
      if(year==2014) sprintf(name_eff, "Dtof0.75Eff");
      if(year==2015) sprintf(name_eff, "Dtof1.00Eff");
      sprintf(title_eff, "#Deltatof cut efficiency");
    }
  else if(strcmp(type,"MtdMthEff")==0)
    {
      sprintf(name_lumi, "");
      sprintf(title_lumi, "");
      sprintf(name_eff, "MtdMthEff");
      sprintf(title_eff, "MTD matching efficiency");
    }

  const char *sys_name[3] = {"", "_Sysup", "_Sysdown"};
  TH1F *hEffVsPt[nPart][3];
  TH1F *hEffVsCent[nPart][3];

  TH1F *hEffSysVsPt[nPart];
  TH1F *hEffSysVsCent[nPart];  

  TFile *fin = 0x0;
  if(strcmp(type,"MtdTrigEff")==0)
    {
      if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type.Data()),"update");
      else          fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type.Data()),"read");
    }
  if(strcmp(type,"Dtof")==0)
    {
      if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type.Data()),"update");
      else          fin = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type.Data()),"read");
    }
  if(strcmp(type,"MtdMthEff")==0)
    {
      if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type.Data()),"update");
      else          fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type.Data()),"read");
    }

  for(int i=0; i<nPart; i++)
    {
      for(int s=0; s<3; s++)
	{
	  if(strcmp(type,"MtdTrigEff")==0)
	    {
	      hEffVsPt[i][s] = (TH1F*)fin->Get(Form("%s_%sEffVsPt_%s%s",run_type.Data(),part_name[i],name_eff,sys_name[s]));
	    }
	  if(strcmp(type,"Dtof")==0)
	    {
	      char *tmp_name = Form("%s",sys_name[s]);
	      if(s==0) 
		{
		  if(year==2014) tmp_name = Form("_FitFunc");
		  if(year==2015) tmp_name = Form("_Syscenter");
		}
	      hEffVsPt[i][s] = (TH1F*)fin->Get(Form("TagAndProbe_%sEffVsPt_%s%s",part_name[i],name_eff,tmp_name));
	    }
	  if(strcmp(type,"MtdMthEff")==0)
	    {
	      cout << Form("%sEffVsPt_%s%s",part_name[i],name_eff,sys_name[s]) << endl;
	      hEffVsPt[i][s] = (TH1F*)fin->Get(Form("%sEffVsPt_%s%s",part_name[i],name_eff,sys_name[s]));
	    }
	}
      int nBins = hEffVsPt[i][0]->GetNbinsX();
      hEffSysVsPt[i] = (TH1F*)hEffVsPt[i][0]->Clone(Form("%s_%sEffVsPt_Sys_%s",run_type.Data(),part_name[i],name_eff));
      hEffSysVsPt[i]->Reset();
      for(int bin=1; bin<=nBins; bin++)
	{
	  double diff1 = fabs(hEffVsPt[i][1]->GetBinContent(bin)/hEffVsPt[i][0]->GetBinContent(bin) - 1);
	  double diff2 = fabs(hEffVsPt[i][2]->GetBinContent(bin)/hEffVsPt[i][0]->GetBinContent(bin) - 1);
	  double error = diff1 > diff2 ? diff1 : diff2;
	  hEffSysVsPt[i]->SetBinContent(bin, 1);
	  hEffSysVsPt[i]->SetBinError(bin, error);
	}

      for(int s=0; s<3; s++)
	{
	  if(strcmp(type,"MtdTrigEff")==0)
	    {
	      hEffVsCent[i][s] = (TH1F*)fin->Get(Form("%s_%sEffVsCent_%s%s",run_type.Data(),part_name[i],name_eff,sys_name[s]));
	    }
	  if(strcmp(type,"Dtof")==0)
	    {
	      char *tmp_name = Form("%s",sys_name[s]);
	      if(s==0) 
		{
		  if(year==2014) tmp_name = Form("_FitFunc");
		  if(year==2015) tmp_name = Form("_Syscenter");
		}
	      hEffVsCent[i][s] = (TH1F*)fin->Get(Form("TagAndProbe_%sEffVsCent_%s%s",part_name[i],name_eff,tmp_name));
	    }
	  if(strcmp(type,"MtdMthEff")==0)
	    {
	      hEffVsCent[i][s] = (TH1F*)fin->Get(Form("%sEffVsCent_%s%s",part_name[i],name_eff,sys_name[s]));
	    }
	}
      nBins = hEffVsCent[i][0]->GetNbinsX();
      hEffSysVsCent[i] = (TH1F*)hEffVsCent[i][0]->Clone(Form("%s_%sEffVsCent_Sys_%s",run_type.Data(),part_name[i],name_eff));
      hEffSysVsCent[i]->Reset();
      for(int bin=1; bin<=nBins; bin++)
	{
	  double diff1 = fabs(hEffVsCent[i][1]->GetBinContent(bin)/hEffVsCent[i][0]->GetBinContent(bin) - 1);
	  double diff2 = fabs(hEffVsCent[i][2]->GetBinContent(bin)/hEffVsCent[i][0]->GetBinContent(bin) - 1);
	  double error = diff1 > diff2 ? diff1 : diff2;
	  hEffSysVsCent[i]->SetBinContent(bin, 1);
	  hEffSysVsCent[i]->SetBinError(bin, error);
	}
    }

  TCanvas *c1[nPart]; 
  TCanvas *c2[nPart]; 
  char outName[256];
  for(int i=0; i<nPart; i++)
    {
      c1[i]= new TCanvas(Form("%s_%s_SysVsPt",part_name[i],name_eff),Form("%s_%s_SysVsPt",part_name[i],name_eff),800,600);
      SetPadMargin(gPad, 0.13, 0.13);
      ScaleHistoTitle(hEffSysVsPt[i],0.05,1,0.04,0.05,1.2,0.04,42);
      hEffSysVsPt[i]->SetTitle(";p_{T} (GeV/c);Sys. Uncert.");
      if(year==2014) 
	{
	  if(strcmp(type,"MtdMthEff")==0)
	    hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.8,1.2);
	  else 
	    hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.9,1.1);
	}
      if(year==2015) hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.95,1.05);
      if(year==2016) hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.9,1.1);
      hEffSysVsPt[i]->SetMarkerStyle(20);
      hEffSysVsPt[i]->Draw();
      TPaveText *t1 = GetTitleText(Form("%s: uncertainty for %s",part_title[i],title_eff),0.05,42);
      t1->Draw();
      if(savePlot) 
	{
	  sprintf(outName, "~/Work/STAR/analysis/Plots/%s/ana_%s/Sys.%s_%sVsPt.pdf",run_type.Data(),type,part_name[i],name_eff);
	  c1[i]->SaveAs(outName);
	}
      if(saveAN && i<2) 
	{
	  if(strcmp(type,"Dtof")==0) 
	    {
	      sprintf(outName, "~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch5_%sDtofEffSys.pdf",part_name[i]);
	    }
	  if(strcmp(type,"MtdTrigEff")==0) 
	    {
	      sprintf(outName, "~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch5_%sMtdTrigEffSys.pdf",part_name[i]);
	    }
	  if(strcmp(type,"MtdMthEff")==0) 
	    {
	      sprintf(outName, "~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch5_%sMtdMthEffSys.pdf",part_name[i]);
	    }
	  c1[i]->SaveAs(outName);
	}

      c2[i]= new TCanvas(Form("%s_%s_SysVsCent",part_name[i],name_eff),Form("%s_%s_SysVsCent",part_name[i],name_eff),800,600);
      SetPadMargin(gPad, 0.13, 0.13);
      ScaleHistoTitle(hEffSysVsCent[i],0.05,1,0.04,0.05,1.2,0.04,42);
      hEffSysVsCent[i]->GetXaxis()->SetBinLabel(1, "p_{T} > 0 GeV/c");
      hEffSysVsCent[i]->GetXaxis()->SetBinLabel(2, "p_{T} > 5 GeV/c");
      if(strcmp(type,"MtdMthEff")==0)
	hEffSysVsCent[i]->GetYaxis()->SetRangeUser(0.8,1.2);
      else
	hEffSysVsCent[i]->GetYaxis()->SetRangeUser(0.9,1.1);
      hEffSysVsCent[i]->GetXaxis()->SetLabelFont(42);
      hEffSysVsCent[i]->GetXaxis()->SetLabelOffset(0.01);
      hEffSysVsCent[i]->GetXaxis()->SetLabelSize(0.07);
      hEffSysVsCent[i]->SetMarkerStyle(20);
      hEffSysVsCent[i]->Draw();
      t1->Draw();
      if(savePlot) 
	{
	  sprintf(outName, "~/Work/STAR/analysis/Plots/%s/ana_%s/Sys.%s_%sVsCent.pdf",run_type.Data(),type,part_name[i],name_eff);
	  c2[i]->SaveAs(outName);
	}
    }

  if(saveHisto)
    {
      for(int i=0; i<nPart; i++)
	{
	  hEffSysVsPt[i]->Write("",TObject::kOverwrite);
	  hEffSysVsCent[i]->Write("",TObject::kOverwrite);
	}
    }

}

//================================================
void anaSys(const char* type, const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  // single muon efficiency
  // TFile *fMuonEff = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type.Data()),"read");
  // hMuonPtEff[0] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff46_BinCount_FitFunc");
  // hMuonPtEff[1] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff46_BinCount_FitFunc_Sysup");
  // hMuonPtEff[2] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff46_BinCount_FitFunc_Sysdown");

  TFile *fdata = 0x0;
  if(strcmp(type,"MtdTrigEff")==0)
    {
      if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type.Data()),"update");
      else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type.Data()),"read");
      hMuonPtEff[0] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff",run_type.Data()));
      hMuonPtEff[1] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff_Sysup",run_type.Data()));
      hMuonPtEff[2] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff_Sysdown",run_type.Data()));
    }
  else if(strcmp(type,"Dtof")==0)
    {
      double dtofCut = 0.75;
      if(year==2015)
	{
	  dtofCut = 1.0;
	}
      if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type.Data()),"update");
      else          fdata = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type.Data()),"read");
      if(year==2014) hMuonPtEff[0] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_FitFunc",dtofCut));
      if(year==2015) hMuonPtEff[0] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_Syscenter",dtofCut));
      hMuonPtEff[1] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_Sysup",dtofCut));
      hMuonPtEff[2] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_Sysdown",dtofCut));
    }
  else if(strcmp(type,"MtdMthEff")==0)
    {
      if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type.Data()),"update");
      else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type.Data()),"read");
      hMuonPtEff[0] = (TF1*)fdata->Get("funcTrkPtBtmBlEff_Type0");
      hMuonPtEff[0]->SetName("Muon_MtdMthEff");
      hMuonPtEff[1] = (TF1*)fdata->Get("funcTrkPtBtmBlEff_Type1");
      hMuonPtEff[1]->SetName("Muon_MtdMthEff_Sysup");
      hMuonPtEff[2] = (TF1*)fdata->Get("funcTrkPtBtmBlEff_Type1");
      hMuonPtEff[2]->SetName("Muon_MtdMthEff_Sysdown");
    }
  printf("[i] Process %s\n",fdata->GetName());

  for(int i=0; i<nPart; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e8);
      
      if(saveHisto)
	{
	  fdata->cd();
	  for(int e=0; e<3; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void mtdTrigEff(const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  TFile *fdata = 0x0;
  if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type.Data()),"update");
  else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type.Data()),"read");
  hMuonPtEff[0] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff",run_type.Data()));

  TFile *fEff = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type.Data()),"read");
  hMuonPtEff[1] = (TF1*)fEff->Get(Form("%s_gTacDiffEffFinalFit_BinCount_prod_Run15_pp200",run_type.Data()));
  hMuonPtEff[2] = (TF1*)fEff->Get(Form("%s_gTacDiffEffFinalFit_BinCount_prod_low_Run15_pp200",run_type.Data()));
  hMuonPtEff[1]->SetName(Form("%s_Muon_TacDiffEff_prod",run_type.Data()));
  hMuonPtEff[2]->SetName(Form("%s_Muon_TacDiffEff_prod_low",run_type.Data()));

  for(int e=0; e<3; e++)
    {
      hMuonPtEff[e]->SetName(Form("%s_TrigStudy",hMuonPtEff[e]->GetName()));
      cout <<  hMuonPtEff[e]->GetName() << endl;
    }

  for(int i=0; i<1; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e7);
      
      if(saveHisto)
	{
	  fdata->cd();
	  for(int e=0; e<3; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void mtdTrigEffLS(const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  for(int i=0; i<hTpcTrackRes->GetNbinsX(); i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  const int nHistos = 4;
  const int mode = 2;

  TFile *fEff = 0x0;
  if(saveHisto) fEff = TFile::Open(Form("Rootfiles/%s.StudyTrigEff.root",run_type.Data()),"update");
  else          fEff = TFile::Open(Form("Rootfiles/%s.StudyTrigEff.root",run_type.Data()),"read");

  const char* lumiName[nHistos] = {"prod", "prod_low", "prod_mid", "prod_high"};
  for(int e=0; e<4; e++)
    {
      hMuonPtEff[e] = (TF1*)fEff->Get(Form("%s_gTacDiffEff_LS_%s_func_M%d",run_type.Data(),lumiName[e],mode));
      hMuonPtEff[e]->SetName(Form("%s_MtdTrig_Muon_%s_M%d",run_type.Data(),lumiName[e],mode));
    }

  for(int i=0; i<1; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e7, nHistos);
      
      if(saveHisto)
	{
	  fEff->cd();
	  for(int e=0; e<nHistos; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void mtdRespEff(const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  TFile *fdata = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type.Data()),"update");
  hMuonPtEff[0] = (TF1*)fdata->Get("Cosmic_FitRespEff_Template");
  hMuonPtEff[1] = (TF1*)fdata->Get("Cosmic_RespEff_Btm_6300_fit");
  hMuonPtEff[2] = (TF1*)fdata->Get("Cosmic_RespEff_Btm_6400_fit");
 
  hMuonPtEff[0]->SetName("Cosmic_Muon_Btm_All");
  hMuonPtEff[1]->SetName("Cosmic_Muon_Btm_6300");
  hMuonPtEff[2]->SetName("Cosmic_Muon_Btm_6400");

  for(int e=0; e<3; e++)
    {
      hMuonPtEff[e]->SetName(Form("%s_TrigStudy",hMuonPtEff[e]->GetName()));
    }

  for(int i=0; i<1; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e7);
      
      if(saveHisto)
	{
	  fdata->cd();
	  for(int e=0; e<3; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void makeHisto(TString name, const double mass, const int nExpr, const int nHistos)
{
  int nbins_tmp = 0;
  double xbins_tmp[10] = {0};
  if(name.Contains("Jpsi"))
    {
      TFile *fin = 0x0;
      if(year==2013)
	{
	  fin = TFile::Open("Rootfiles/pp500GeVfit_new.root","read");
	  TF1 *funcJpsi = (TF1*)fin->Get("ffpt");
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
	  fin = TFile::Open("Rootfiles/models.root","read");
	  hMcJpsiPt = (TH1F*)fin->Get(Form("TBW_JpsiYield_AuAu200_cent0060"));
	}
      else if(year==2015)
        {
          fin = TFile::Open("Rootfiles/JpsiSpectraShapepp200.root","read");
          TF1 *funcJpsi = (TF1*)fin->Get("TsallisPowerLawFitJpsipp200");
          funcJpsi->SetNpx(1000);
          hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
          hMcJpsiPt->SetName(Form("Jpsi_Yield_cent00100"));
          for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
            {
              hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
            }
        }

      if(year==2013)
	{
	  nbins_tmp = 5;
	  double xbins_tmp_tmp[6] = {0,1.5,3,5,7,9};
	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
	}
      else if(year==2014)
	{
	  nbins_tmp = 9;
	  double xbins_tmp_tmp[10] = {0.15,1,2,3,4,5,6,8,10,15};
	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
	}
      else if(year==2015)
        {
          nbins_tmp = 8;
          double xbins_tmp_tmp[9] = {0,1,2,3,4,5,6,8,10};
    	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
        }
      else if(year==2016)
	{
	  nbins_tmp = 6;
	  double xbins_tmp_tmp[7] = {0,1,2,3,4,6,10};
	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
	}
    }
  else
    {
      TF1 *fBol = new TF1("Boltzmann","x/(exp(x/[0]+1))",0,10);
      fBol->SetParameter(0,1.11);
      fBol->SetNpx(1000);
      hMcJpsiPt = (TH1F*)fBol->GetHistogram();
      nbins_tmp = 3;
      double xbins_tmp_tmp[4] = {0,2,4,10};
      for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
    }
  hMcJpsiPt->Scale(1./hMcJpsiPt->Integral());

  
  // book histograms
  const int nbinsPt = nbins_tmp;
  double xbinsPt[nbinsPt+1];
  for(int i=0; i<nbinsPt+1; i++)
    {
      xbinsPt[i] = xbins_tmp[i];
    }
  TString hName;
  hInJpsiPt = new TH1F(Form("hInPartPt_%s",name.Data()), "", nbinsPt, xbinsPt);
  hInJpsiPt->Sumw2();
  hInJpsiCent = new TH1F(Form("hInPartCent_%s",name.Data()), "", 2, 0, 2);
  hInJpsiCent->Sumw2();
  for(int i=0; i<nHistos; i++)
    {
      hName =  hMuonPtEff[i]->GetName();
      hName.ReplaceAll("Muon",Form("%sEffVsPt",name.Data()));
      hOutJpsiPt[i] = new TH1F(hName.Data(), "", nbinsPt, xbinsPt);
      hOutJpsiPt[i]->Sumw2();

      hName.ReplaceAll("VsPt","VsCent");
      hOutJpsiCent[i] = new TH1F(hName.Data(), "", 2, 0, 2);
      hOutJpsiCent[i]->Sumw2();
    }
  toyMC(mass, nExpr, nHistos);

  TCanvas *c = 0x0;
  for(int e=0; e<nHistos; e++)
    {
      hOutJpsiPt[e]->Divide(hInJpsiPt);
      hOutJpsiPt[e]->SetMarkerStyle(21);
      hOutJpsiPt[e]->SetMarkerColor(TMath::Power(2, e));
      hOutJpsiPt[e]->SetLineColor(TMath::Power(2, e));
      hOutJpsiPt[e]->GetYaxis()->SetRangeUser(0.4,1);
    }
  c = draw1D(hOutJpsiPt[0],"");
  for(int e=1; e<nHistos; e++)
    hOutJpsiPt[e]->Draw("sames");

  for(int e=0; e<nHistos; e++)
    {
      hOutJpsiCent[e]->Divide(hInJpsiCent);
      hOutJpsiCent[e]->SetMarkerStyle(21);
      hOutJpsiCent[e]->SetMarkerColor(TMath::Power(2, e));
      hOutJpsiCent[e]->SetLineColor(TMath::Power(2, e));
      hOutJpsiCent[e]->GetYaxis()->SetRangeUser(0.4,1);
    }
  c = draw1D(hOutJpsiCent[0],"");
  for(int e=1; e<nHistos; e++)
    hOutJpsiCent[e]->Draw("sames");
}

//================================================
void toyMC(const double mass, const int nExpr, const int nHistos, const int debug)
{
  double pt1_cut = 1.5;
  double pt2_cut = 1.3;
  if(mass>8) //Upsilon
    {
      pt1_cut = 4.0;
      pt2_cut = 1.5;
    }
  int nHisto = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("+++ Event %d +++\n",i+1);
      double mc_pt  = myRandom->Uniform(0,20);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      double weight = hMcJpsiPt->GetBinContent(hMcJpsiPt->FindFixBin(mc_pt));

      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());

      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      
      TLorentzVector daughter2 = parent - daughter1;
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      if(debug) printf("daugther 2: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt2,eta2,phi2);

      // acceptance cut
      if(fabs(eta1)>0.5 || fabs(eta2)>0.5) continue;
      if(fabs(pt1)<0.5  || fabs(pt2)<0.5)  continue;

      // momentum resolution 
      int mom_index1 = hTpcTrackRes->GetXaxis()->FindBin(pt1)-1;
      if(mom_index1>=nHisto) mom_index1=nHisto-1;
      int mom_index2 = hTpcTrackRes->GetXaxis()->FindBin(pt2)-1;
      if(mom_index2>=nHisto) mom_index2=nHisto-1;
      double dpt1 = hTrkResBin[mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[mom_index2]->GetRandom();
      double rc_pt1 = (1-dpt1) * pt1;
      double rc_pt2 = (1-dpt2) * pt2;
      double leadpt = rc_pt1 > rc_pt2 ? rc_pt1 : rc_pt2;
      double subpt  = rc_pt1 < rc_pt2 ? rc_pt1 : rc_pt2;
      if(leadpt<pt1_cut || subpt<pt2_cut) continue;

      hInJpsiPt->Fill(mc_pt,weight);
      if(mc_pt>0) hInJpsiCent->Fill(0.5,weight);
      if(mc_pt>5) hInJpsiCent->Fill(1.5,weight);

      // single muon efficiency
      for(int e=0; e<nHistos; e++)
	{
	  double eff1 = hMuonPtEff[e]->Eval(rc_pt1);
	  double eff2 = hMuonPtEff[e]->Eval(rc_pt2);
	  if(myRandom->Uniform(0., 1.)<eff1 && myRandom->Uniform(0,1)<eff2)
	    {
	      hOutJpsiPt[e]->Fill(mc_pt,weight);
	      if(mc_pt>0) hOutJpsiCent[e]->Fill(0.5,weight);
	      if(mc_pt>5) hOutJpsiCent[e]->Fill(1.5,weight);
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

