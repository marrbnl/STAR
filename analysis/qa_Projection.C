TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const Double_t tpc_sector_low[12] = {0.78,  0.255, 6.01,  4.97, 4.45, 3.925, 3.4,   2.88, 2.355, 1.83,  1.305};
const Double_t tpc_sector_up[12]  = {1.305, 0.78,  0.255, 5.49, 4.97, 4.45,  3.925, 3.4,  2.88,  2.355, 1.83};

//const char *run_config = "NoEloss.Global.HLT";
const char *run_config = "projection.Global.HLT";

//================================================
void qa_Projection()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.Match.%s.root",),"read");
  //f = TFile::Open("/Users/admin/Work/STAR/analysis/output/Run13.pp500.jpsi.CutRanking.root","read");

  //makeHisto();
  //trackYZ();
  //trackProjection();
 EnergyLoss();
  //magneticField();
}


//================================================
void makeHisto()
{
  THnSparseF *hTrkYZ = (THnSparseF*)f->Get(Form("mhTrkProjYZ_qa_%s",trigName[kTrigType]));

  Int_t pt_cuts[4] = {1,2,5,20};

  const char *trkYZ_name[4] = {"OuterTPC","OuterEmc","OuterCoil","OuterMag"};
  TString charge_name[3] = {"all","negative","positive"};
 

  TH2F *hTrkPhiZ[4];
  TH1F *hTrkPhi[4][3][4];

  for(Int_t i=0; i<4; i++)
    {
      Int_t low_bin = i+1; 
      Int_t up_bin  = i+1;
      if(i==3)
	{
	  up_bin = 14;
	}

      hTrkYZ->GetAxis(4)->SetRange(low_bin,up_bin); // radius
      hTrkPhiZ[i] = (TH2F*)hTrkYZ->Projection(2,1);
      hTrkPhiZ[i]->SetName(Form("hTrkPhiZ_%s",trkYZ_name[i]));

      for(Int_t j=0; j<3; j++)
	{
	  if(j>0) hTrkYZ->GetAxis(3)->SetRange(2*j-1,2*j-1); // charge
	  for(Int_t k=0; k<4; k++)
	    {
	      if(k>0)  hTrkYZ->GetAxis(0)->SetRange(pt_cuts[k-1]+1,pt_cuts[k]);
	      hTrkPhi[i][j][k] = (TH1F*)hTrkYZ->Projection(2);
	      if(k==0) hTrkPhi[i][j][k]->SetName(Form("hTrkPhi_%s_%s",charge_name[j].Data(),trkYZ_name[i]));
	      else     hTrkPhi[i][j][k]->SetName(Form("hTrkPhi_%s_%s_pt_%d_%d",charge_name[j].Data(),trkYZ_name[i],pt_cuts[k-1],pt_cuts[k]));
	      hTrkYZ->GetAxis(0)->SetRange(0,-1);
	    }
	  hTrkYZ->GetAxis(3)->SetRange(0,-1);
	}
      hTrkYZ->GetAxis(4)->SetRange(0,-1);
    }

  // Backleg gap structure
  TH1F *hTrkPhiMag[11][3];
  for(Int_t i=0; i<11; i++)
    {
      hTrkYZ->GetAxis(4)->SetRange(4+i,4+i);
      for(Int_t j=0; j<3; j++)
	{
	  if(j>0) hTrkYZ->GetAxis(3)->SetRange(2*j-1,2*j-1);
	  hTrkPhiMag[i][j] = (TH1F*)hTrkYZ->Projection(2);
	  hTrkPhiMag[i][j]->SetName(Form("hTrkPhi_%s_OuterMag_zeroFieldStep%d",charge_name[j].Data(),i));
	  hTrkYZ->GetAxis(3)->SetRange(0,-1);
	}
      hTrkYZ->GetAxis(4)->SetRange(0,-1);
    }

  // pt dependence
  TH1F *hTrkPhiMag_pt[3];
  hTrkYZ->GetAxis(4)->SetRange(4,14);
  for(Int_t i=0; i<3; i++)
    {
      hTrkYZ->GetAxis(0)->SetRange(pt_cuts[i]+1,pt_cuts[i+1]);
      hTrkPhiMag_pt[i] = (TH1F*)hTrkYZ->Projection(2);
      hTrkPhiMag_pt[i]->SetName(Form("hTrkPhi_OuterMag_pt_%d_%d",pt_cuts[i],pt_cuts[i+1]));
      hTrkYZ->GetAxis(0)->SetRange(0,-1);
    } 
  hTrkYZ->GetAxis(4)->SetRange(0,-1);

  TFile *fout = TFile::Open("Rootfiles/Run13.pp500.TrkProj.NoELoss.root","update");   
  for(Int_t i=0; i<4; i++)
    {
      hTrkPhiZ[i]->Write("",TObject::kOverwrite);
      for(Int_t j=0; j<3; j++)
	{
	  for(Int_t k=0; k<4; k++)
	    {
	      hTrkPhi[i][j][k]->Write("",TObject::kOverwrite);
	    }
	}
    }
  for(Int_t i=0; i<11; i++)
    {
      for(Int_t j=0; j<3; j++)
	{
	  hTrkPhiMag[i][j]->Write("",TObject::kOverwrite);
	}
    }
  for(Int_t i=0; i<3; i++)
    {
      hTrkPhiMag_pt[i]->Write("",TObject::kOverwrite);
    }
  fout->Close();

}


//================================================
void trackProjection(const Int_t save = 0)
{
  f = TFile::Open("Rootfiles/Run13.pp500.TrkProj.NoELoss.root","read");   

  TList *list = new TList;

  THnSparseF *hTrkYZ = (THnSparseF*)f->Get(Form("mhTrkProjYZ_qa_%s",trigName[kTrigType]));

  const char *trkYZ_title[4] = {"outer TPC","outer BEMC","outer coil","outer magnet"};
  const char *trkYZ_name[4] = {"OuterTPC","OuterEmc","OuterCoil","OuterMag"};

  TH1F *hTrkPhiInCharge[4][2];
  TString legName_charge[2] = {"negative","positive"};
  
  for(Int_t i=0; i<4; i++)
    {
      TH2F *h2 = (TH2F*)f->Get(Form("hTrkPhiZ_%s",trkYZ_name[i]));
      if(save) 
	{
	  c = draw2D(h2,Form("%s: #varphi vs global z of %s tracks at %s",trigName[kTrigType],trk_name[trk_index],trkYZ_title[i]),0.04,kFALSE);
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/TrkPhiVsZ_%s_%s.png",run_type,trkYZ_name[i],trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/TrkPhiVsZ_%s_%s.pdf",run_type,trkYZ_name[i],trigName[kTrigType]));
	}

      for(Int_t j=0; j<2; j++)
	{
	  hTrkPhiInCharge[i][j] = (TH1F*)f->Get(Form("hTrkPhi_%s_%s",legName_charge[j].Data(),trkYZ_name[i]));
	} 
    }

  // charge dependence
  for(Int_t i=0; i<4; i++)
    {
      list->Clear();
      for(Int_t j=0; j<2; j++)
	{
	  list->Add(hTrkPhiInCharge[i][j]);
	}
      c = drawHistos(list,Form("hTrkPhi_%s",trkYZ_name[i]),Form("%s: #varphi of projected tracks at %s;#varphi",trigName[kTrigType],trkYZ_title[i]),kFALSE,0,100,kTRUE,1e-1,hTrkPhiInCharge[i][1]->GetMaximum()*1.5,kFALSE,kTRUE,legName_charge,kTRUE,"p_{T} > 1 GeV/c",0.6,0.8,0.65,0.85,kFALSE);

      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/TrkPhi_vs_charge_%s_%s.png",run_type,trkYZ_name[i],trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/TrkPhi_vs_charge_%s_%s.pdf",run_type,trkYZ_name[i],trigName[kTrigType]));
	}
    }

  // charge evolution
  Double_t low_boundary[4][2] = { {1.3, 1.32}, {1.26, 1.36}, {1.22, 1.40}, {1.21, 1.39}}; 
  TH1F *hTrkPhiInSector[4][2];
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hTrkPhiInSector[i][j] = new TH1F(Form("hTrkPhiInSector_%s_%s",legName_charge[j].Data(),trkYZ_name[i]),Form("%s: #varphi of projected tracks at %s;#varphi",trigName[kTrigType],trkYZ_title[i]), 600, pi/4, 0.75*pi);
	  Int_t period = hTrkPhiInCharge[i][j]->GetNbinsX()/12;
	  Int_t low_bin = hTrkPhiInCharge[i][j]->FindBin(low_boundary[i][j]) - period;
	  Int_t up_bin  = low_bin + period * 3;
	  for(Int_t ibin=low_bin; ibin<up_bin; ibin++)
	    {
	      Int_t jbin = hTrkPhiInSector[i][j]->FindBin(hTrkPhiInCharge[i][j]->GetBinCenter(ibin));
	      hTrkPhiInSector[i][j]->SetBinContent(jbin,hTrkPhiInCharge[i][j]->GetBinContent(ibin));
	      hTrkPhiInSector[i][j]->SetBinError(jbin,hTrkPhiInCharge[i][j]->GetBinError(ibin));
	    }
	}
    }

  for(Int_t j=0; j<2; j++)
    {
      TCanvas *c = new TCanvas(Form("TrkPhiInSector_%s",legName_charge[j].Data()),Form("TrkPhiInSector_%s",legName_charge[j].Data()),600,800);
      c->Divide(1,4);
      for(Int_t i=0; i<4; i++)
	{
	  c->cd(i+1);
	  SetPadMargin(gPad,0.1,0.1,0.1,0.01);
	  TH1F *htmp = (TH1F*)hTrkPhiInSector[i][j]->Clone(Form("%s_clone",hTrkPhiInSector[i][j]->GetName()));
	  htmp->SetTitle(";#varphi;");
	  htmp->SetLineColor(color[i]);
	  ScaleHistoTitle(htmp,0.08,0.5,0.06,0.045,1,0.035,62);
	  htmp->SetMaximum(1.2*htmp->GetMaximum());
	  htmp->Draw("HIST");
	  TPaveText *t1 = GetPaveText(0.2,0.25,0.8,0.9,0.1);
	  t1->AddText(trkYZ_title[i]);
	  t1->SetTextColor(color[i]);
	  t1->Draw();
	  TLine *line = GetLine(low_boundary[i][j],0,low_boundary[i][j],htmp->GetMaximum()*0.8,color[i]);
	  line->Draw();
	  line = GetLine(low_boundary[i][j]+hTrkPhiInCharge[i][j]->GetNbinsX()/12*hTrkPhiInCharge[i][j]->GetXaxis()->GetBinWidth(1),0,
			 low_boundary[i][j]+hTrkPhiInCharge[i][j]->GetNbinsX()/12*hTrkPhiInCharge[i][j]->GetXaxis()->GetBinWidth(1),htmp->GetMaximum()*0.8,color[i]);
	  line->Draw();
	}
      c->cd(1);
      t1 = GetPaveText(0.7,0.8,0.7,0.9,0.1);
      t1->AddText("p_{T} > 1 GeV/c");
      t1->AddText("No energy loss");
      t1->Draw();
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/TrkPhi_projection_%s_%s.png",run_type,legName_charge[j].Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/TrkPhi_projection_%s_%s.pdf",run_type,legName_charge[j].Data(),trigName[kTrigType]));
	}
    }

  return;

  // Backleg gap structure
  TH1F *hTrkPhiMag[3][11];
  for(Int_t i=0; i<3; i++)
    {
      if(i>0) hTrkYZ->GetAxis(3)->SetRange(2*i-1,2*i-1);
      for(Int_t j=0; j<11; j++)
	{
	  hTrkYZ->GetAxis(4)->SetRange(4+j,4+j);
	  hTrkPhiMag[i][j] = (TH1F*)hTrkYZ->Projection(2);
	  hTrkPhiMag[i][j]->SetName(Form("hTrkPhi_charge%d_zeroFiled%d",i,j));
	  hTrkYZ->GetAxis(4)->SetRange(0,-1);
	}
      hTrkYZ->GetAxis(3)->SetRange(0,-1);
    }

  list->Clear();
  list->Add(hTrkPhiMag[0][0]);
  list->Add(hTrkPhiMag[0][2]);
  list->Add(hTrkPhiMag[0][5]);
  list->Add(hTrkPhiMag[0][7]);
  list->Add(hTrkPhiMag[0][10]);
  TString legName_zeroMag[5] = {"nZeroFieldStep = 0","nZeroFieldStep = 2","nZeroFieldStep = 5","nZeroFieldStep = 7","nZeroFieldStep = 10"};
  c = drawHistos(list,"hTrkYPhi_ZeroFiled",Form("Au+Au %s: #varphi of projected tracks at outer magnet;#varphi",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e3,50*hTrkPhiMag[0][0]->GetMaximum(),kTRUE,kTRUE,legName_zeroMag,kTRUE,"",0.5,0.7,0.65,0.88,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_phi_OuterMag_ZeroField_%s.png",run_config,trigName[kTrigType]));

  c = drawHistos(list,"hTrkYPhi_ZeroFiled_zoomin",Form("Au+Au %s: #varphi of projected tracks at outer magnet;#varphi",trigName[kTrigType]),kTRUE,0,pi/4,kTRUE,1e3,hTrkPhiMag[0][0]->GetMaximum(),kTRUE,kTRUE,legName_zeroMag,kTRUE,"",0.5,0.7,0.65,0.88,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_phi_OuterMag_ZeroField_zoomin_%s.png",run_config,trigName[kTrigType]));
  hTrkPhiMag[0][0]->SetMaximum(1./50*hTrkPhiMag[0][0]->GetMaximum());

  TH1F *hTrkPhiMag_period[3][11];
  TH1F *hTrkPhiMag_period_all[3];
  for(Int_t i=0; i<3; i++)
    {
      for(Int_t j=0; j<11; j++)
	{
	  hTrkPhiMag_period[i][j] = new TH1F(Form("%s_period",hTrkPhiMag[i][j]->GetName()),hTrkPhiMag[i][j]->GetTitle(),hTrkPhiMag[i][j]->GetNbinsX()/2,0,pi);
	  hTrkPhiMag_period[i][j]->Sumw2();
	  for(Int_t ibin=1; ibin<=hTrkPhiMag[i][j]->GetNbinsX(); ibin++)
	    {
	      Int_t jbin = ibin-(ibin-1)/80*80;
	      hTrkPhiMag_period[i][j]->Fill(hTrkPhiMag_period[i][j]->GetBinCenter(jbin),hTrkPhiMag[i][j]->GetBinContent(ibin));
	    }

	  if(j==0) hTrkPhiMag_period_all[i] = (TH1F*)hTrkPhiMag_period[i][j]->Clone(Form("hTrkPhi_charge%d_period_all",i));
	  else     hTrkPhiMag_period_all[i]->Add(hTrkPhiMag_period[i][j]);
	}
    }
  list->Clear();
  list->Add(hTrkPhiMag_period_all[0]);
  list->Add(hTrkPhiMag_period_all[1]);
  list->Add(hTrkPhiMag_period_all[2]);
  TString legName_period[14];
  legName_period[0] = "Sum";
  legName_period[1] = "Negative";
  legName_period[2] = "Positive";
  for(Int_t j=0; j<11; j++)
    {
      list->Add(hTrkPhiMag_period[0][j]);
      legName_period[j+3] =  Form("nZeroFieldStep = %d", j);
    }
  c = drawHistos(list,"hTrkYPhi_ZeroFiled_period",Form("Au+Au %s: #varphi of projected tracks at outer magnet;#varphi",trigName[kTrigType]),kTRUE,0,0.35,kFALSE,1e3,50*hTrkPhiMag[0][0]->GetMaximum(),kFALSE,kTRUE,legName_period,kTRUE,"",0.7,0.85,0.3,0.88,kFALSE);
  gPad->SetRightMargin(0.01);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_phi_OuterMag_period_%s.png",run_config,trigName[kTrigType]));

  // pt dependence
  Int_t pt_cuts[4] = {1,2,5,20};
  TH1F *hTrkPhiMag_pt[3];
  hTrkYZ->GetAxis(4)->SetRange(4,14);
  list->Clear();
  TString legName_pt[3];
  for(Int_t i=0; i<3; i++)
    {
      hTrkYZ->GetAxis(0)->SetRange(pt_cuts[i]+1,pt_cuts[i+1]);
      TH1F *htmp = (TH1F*)hTrkYZ->Projection(2);
      htmp->SetName(Form("hTrkPhiMag_pt_%d_%d_tmp",pt_cuts[i],pt_cuts[i+1]));
      hTrkYZ->GetAxis(0)->SetRange(0,-1);

      hTrkPhiMag_pt[i] = new TH1F(Form("hTrkPhiMag_pt_%d_%d",pt_cuts[i],pt_cuts[i+1]),htmp->GetTitle(),htmp->GetNbinsX()/2,0,pi);
      for(Int_t ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
	{
	  Int_t jbin = ibin-(ibin-1)/80*80;
	  hTrkPhiMag_pt[i]->Fill(hTrkPhiMag_pt[i]->GetBinCenter(jbin),htmp->GetBinContent(ibin));
	}
      list->Add(hTrkPhiMag_pt[i]);;
      legName_pt[i] = Form("%d < p_{T} < %d GeV/c",pt_cuts[i],pt_cuts[i+1]);
    }  
  hTrkYZ->GetAxis(4)->SetRange(0,-1);
  c = drawHistos(list,"hTrkYPhi_ZeroFiled_pt",Form("Au+Au %s: #varphi of projected tracks at outer magnet;#varphi",trigName[kTrigType]),kTRUE,0,0.35,kTRUE,1e4,1e8,kTRUE,kTRUE,legName_pt,kTRUE,"",0.7,0.85,0.6,0.88,kFALSE);
  gPad->SetRightMargin(0.01);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_phi_OuterMag_period_pt_%s.png",run_config,trigName[kTrigType]));
  
}


//================================================
void trackYZ(const Int_t save = 1)
{
  TList *list = new TList;

  THnSparseF *hTrkProjYZ = (THnSparseF*)f->Get(Form("mhTrkYZatMtd_qa_%s",trigName[kTrigType]));

  const char *title[2] = {"global z","local y"};
  const char *name[2] = {"gz","ly"};
  TH2F *hBL[2];
  TH1F *hAll[2];

  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  for(Int_t i=0; i<2; i++)
    {
      hBL[i] = (TH2F*)hTrkProjYZ->Projection(2-i,3);
      hBL[i]->SetName(Form("Track_%s_vs_BL",name[i]));
      if(i==1) hBL[i]->GetYaxis()->SetRangeUser(-50,50);

      c = draw2D(hBL[i],Form("%s: %s distribution of projected tracks at MTD per module;backleg*5+module",trigName[kTrigType],title[i]));
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_%s_vs_Mod_%s.png",run_type,name[i],trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_%s_vs_Mod_%s.pdf",run_type,name[i],trigName[kTrigType]));
	}

      hAll[i] = (TH1F*)hBL[i]->ProjectionY(Form("Track_%s",name[i]));    

      c = draw1D(hAll[i],Form("%s: %s distribution of projected tracks at MTD",trigName[kTrigType],title[i]),kFALSE,kFALSE);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_%s_%s.png",run_type,name[i],trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_%s_%s.pdf",run_type,name[i],trigName[kTrigType]));
	}
    }

  TH2F *hTrkYZ = (TH2F*)hTrkProjYZ->Projection(1,2);
  hTrkYZ->SetName(Form("Track_y_vs_z"));
  hTrkYZ->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYZ,Form("%s: local y vs global z of projected tracks at MTD;global z (cm);local y (cm)",trigName[kTrigType]),0.04,kFALSE);
  gPad->SetRightMargin(0.13);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_gz_%s.png",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_gz_%s.pdf",run_type,trigName[kTrigType]));
    }
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);

  // non-overlapping
  hTrkProjYZ->GetAxis(5)->SetRange(1,1);
  TH2F *hTrkYvsZ = (TH2F*)hTrkProjYZ->Projection(1,2);
  hTrkYvsZ->SetName(Form("Track_y_vs_z_nonOverlapping"));
  hTrkYvsZ->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYvsZ,Form("%s:local y vs global z of projected tracks (non-overlapping);global z (cm);local y (cm)",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_gz_nonOverlapping_%s.png",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_gz_nonOverlapping_%s.pdf",run_type,trigName[kTrigType]));
    }

  TH1F *hTrkY = (TH1F*)hTrkYvsZ->ProjectionY(Form("Track_y_nonOverlapping"));    
  c = draw1D(hTrkY,Form("%s: %s of projected tracks at MTD (non-overlapping)",trigName[kTrigType],title[1]),kFALSE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_%s_nonOverlapping_%s.png",run_type,name[1],trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_%s_nonOverlapping_%s.pdf",run_type,name[1],trigName[kTrigType]));
    }
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);

  // pt dependence
  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  TH2F *hTrkYvsPt = (TH2F*)hTrkProjYZ->Projection(1,0);
  hTrkYvsPt->SetName(Form("Track_y_vs_pt"));
  hTrkYvsPt->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYvsPt,Form("%s:local y vs p_{T} of projected tracks at MTD radius;p_{T} (GeV/c);local y (cm)",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_pt_%s.png",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
    }

  Double_t pt_cut[5] = {1,1.5,2,5,20};
  TH1F *hTrkYinPt[4];
  list->Clear();
  TString legName_pt[4];
  for(Int_t i=0; i<4; i++)
    {
      Int_t low_bin = hTrkYvsPt->GetXaxis()->FindFixBin(pt_cut[i]+0.1);
      Int_t up_bin  = hTrkYvsPt->GetXaxis()->FindFixBin(pt_cut[i+1]-0.1);
      hTrkYinPt[i] = (TH1F*)hTrkYvsPt->ProjectionY(Form("hTrkY_pt_%1.1f_to_%1.1f",pt_cut[i],pt_cut[i+1]),low_bin,up_bin);
      //hTrkYinPt[i]->Scale(1./hTrkYinPt[i]->Integral());
      list->Add(hTrkYinPt[i]);
      legName_pt[i] = Form("%1.1f < p_{T} < %1.1f GeV/c",pt_cut[i],pt_cut[i+1]);
    }
  c = drawHistos(list,"hTrkY_in_pt",Form("%s: local y of projected tracks at MTD radius;y (cm);counts",trigName[kTrigType]),kFALSE,0,100,kFALSE,1e4,1e8,kFALSE,kTRUE,legName_pt,kTRUE,"",0.35,0.55,0.65,0.85,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_InPtBin_%s.png",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_InPtBin_%s.pdf",run_type,trigName[kTrigType]));
    }
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);
  list->Clear();

  // charge dependence
  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  TH2F *hTrkYvsCharge = (TH2F*)hTrkProjYZ->Projection(1,4);
  hTrkYvsCharge->SetName(Form("Track_y_vs_charge"));
  hTrkYvsCharge->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYvsCharge,Form("%s:local y vs charge of projected tracks at MTD radius;charge;local y (cm)",trigName[kTrigType]),0.04,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_charge_%s.png",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_vs_charge_%s.pdf",run_type,trigName[kTrigType]));
    }

  TH1F *hTrkYinCharge[2];
  TString legName_charge[2] = {"negative","positive"};
  for(Int_t i=0; i<2; i++)
    {
      hTrkYinCharge[i] = (TH1F*)hTrkYvsCharge->ProjectionY(Form("hTrkY_%s",legName_charge[i].Data()),1+i*2,1+i*2);
      list->Add(hTrkYinCharge[i]);
    }
  c = drawHistos(list,"hTrkY_in_charge",Form("%s: local y of projected tracks with different charges;y (cm);counts",trigName[kTrigType]),kFALSE,0,100,kFALSE,1e-1,6e6,kFALSE,kTRUE,legName_charge,kTRUE,"p_{T} > 1 GeV/c",0.4,0.6,0.65,0.85,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_InCharge_%s.png",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/Trk_ly_InCharge_%s.pdf",run_type,trigName[kTrigType]));
    }
  list->Clear();
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);
}


//================================================
void magneticField(const Int_t save = 0)
{
  TH1F *hCuts = (TH1F*)f->Get("hAnalysisCuts");
  Double_t scale = hCuts->GetBinContent(3)/1e4;

  TH2F *hMag = (TH2F*)f->Get("hMagneticMap");
  hMag->Scale(1./scale);
  hMag->SetZTitle("B (T)");
  hMag->SetTitleOffset(1.2,"Z");
  TCanvas *c = new TCanvas("MagneticMap","MagneticMap",880,800);
  TPaveText *t1 = GetTitleText(hMag->GetTitle());
  hMag->SetTitle("");
  hMag->GetYaxis()->SetTitleOffset(1.4);
  hMag->Draw("colz");
  t1->Draw();
  SetPadMargin(gPad,0.12,0.12,0.15,0.12);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/MagneticFieldMap.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Projection/MagneticFieldMap.pdf",run_type));
    }
}


//================================================
void EnergyLoss(const Int_t save = 0)
{
  TFile *f1 = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.Match.Global.HLT.root"),"read");

  const char *name[2] = {"wEloss","woEloss"};
  THnSparseF *hTrkProjYZ[2]; 
  TH1F *hTrkY[2];

  // With energy loss
  hTrkProjYZ[0] = (THnSparseF*)f1->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hTrkProjYZ[0]->SetName(Form("%s_%s",hTrkProjYZ[0]->GetName(),name[0]));
  hTrkY[0] = (TH1F*)hTrkProjYZ[0]->Projection(0);

  // Without energy loss
  hTrkProjYZ[1] = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hTrkProjYZ[1]->SetName(Form("%s_%s",hTrkProjYZ[1]->GetName(),name[1]));
  hTrkProjYZ[1]->GetAxis(5)->SetRange(2,2);
  hTrkY[1] = (TH1F*)hTrkProjYZ[1]->Projection(1);

  TList *list = new TList;
  TString legName[2] = {"With energy loss during projection","Without energy loss during projection"};
  for(Int_t i=0; i<2; i++)
    {
      hTrkY[i]->SetName(Form("hTrkY_%s",name[i]));
      hTrkY[i]->Scale(1./hTrkY[i]->Integral());
      list->Add(hTrkY[i]);
    }
  c = drawHistos(list,"hTrkY_Eloss",Form("Au+Au %s: local y of projected tracks at MTD radius;y (cm)",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-5,0.043,kFALSE,kTRUE,legName,kTRUE,"p_{T} > 1 GeV/c",0.15,0.6,0.7,0.88,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_Eloss_%s.png",run_config,trigName[kTrigType]));
}
