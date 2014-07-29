TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "RotateMTD.Global.HLT";

//================================================
void ana_RotateMTD()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  //randomMatch();
  MonteCarlo();
  //makeHisto();
}

//================================================
void randomMatch(const Int_t save = 0)
{
  TFile *f1 = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.GlobalTrk.NoDCA.HLT.root"),"read");
  TFile *f2 = TFile::Open("Rootfiles/AuAu200.RotateMTD.root","read");

  const char *title[2] = {"z","y"};
  const char *trkname[2] = {"projected tracks","MTD hits"};
  const char *name[2] = {"track","hit"};
  TString legName[2] = {"Standard","Rotated"};

  // hit map
  TH2F *hMtdHitMap = (TH2F*)f->Get(Form("hMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMtdHitMap,Form("Au+Au di-muon: channel vs backleg of good MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.MtdGoodHitMap_%s.png",run_config,trigName[kTrigType]));

  // matched hit map
  TH2F *hMthMtdHitMap = (TH2F*)f->Get(Form("hMthMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMthMtdHitMap,Form("Au+Au di-muon: channel vs backleg of matched MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.MthMtdHitMap_%s.png",run_config,trigName[kTrigType]));
  
  // y/z distribution of hits and projected tracks
  THnSparseF *hYZ[2];
  hYZ[0] = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hYZ[1] = (THnSparseF*)f->Get(Form("hHitYZ_qa_%s",trigName[kTrigType]));
  TH2F *hYBL[2], *hZBL[2];
  TH1F *hYAll[2], *hZAll[2];
  for(Int_t i=0; i<2; i++)
    {
      hYBL[i] = (TH2F*)hYZ[i]->Projection(0,2);
      hYBL[i]->SetName(Form("%s_YBL",trkname[i]));
      hYBL[i]->GetYaxis()->SetRangeUser(-50,50);

      hZBL[i] = (TH2F*)hYZ[i]->Projection(1,2);
      hZBL[i]->SetName(Form("%s_ZBL",trkname[i]));

      c1 = draw2D(hYBL[i],Form("Local y distribution of %s per module;backleg*5+module",trkname[i]));
      c2 = draw2D(hZBL[i],Form("Global z distribution of %s per module;backleg*5+module",trkname[i]));
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.%s_y_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.%s_z_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));

      hYAll[i] = (TH1F*)hYBL[i]->ProjectionY(Form("%s_Y_All",trkname[i]));
      hZAll[i] = (TH1F*)hZBL[i]->ProjectionY(Form("%s_Z_All",trkname[i]));     
      hYAll[i]->Scale(1./hYAll[i]->Integral());
      hZAll[i]->Scale(1./hZAll[i]->Integral());

      c1 = draw1D(hYAll[i],Form("Local y distribution of %s",trkname[i]),kFALSE,kFALSE);
      c2 = draw1D(hZAll[i],Form("Global z distribution of %s",trkname[i]),kFALSE,kFALSE);
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.%s_y_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.%s_z_%s.png",run_config,name[i],trigName[kTrigType]));
    }

  // local z distribution of MTD hits in each module
  TH1F *hHitZ[5];
  TList *list1 = new TList;
  TString legName1[5];
  for(Int_t i=0; i<5; i++)
    {
      for(Int_t j=0; j<30; j++)
	{
	  TH1F *htmp = (TH1F*)hZBL[1]->ProjectionY(Form("hHitZ_BL%d_Mod%d",j+1,i+1),j*5+i+1,j*5+i+1);
	  if(j==0) hHitZ[i] = (TH1F*)htmp->Clone(Form("hHitZ_Mod%d_AllBL",i+1));
	  else     hHitZ[i]->Add(htmp);
	}
      hHitZ[i]->Sumw2();
      hHitZ[i]->Scale(1./hHitZ[i]->Integral());
      list1->Add(hHitZ[i]);
      legName1[i] = Form("Module %d",i+1);
    }
  c = drawHistos(list1,list1->At(0)->GetName(),Form("Au+Au %s: global z distribution of MTD hits in each module;z (cm)",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-6,0.018,kFALSE,kTRUE,legName1,kTRUE,"",0.6,0.8,0.68,0.88,kFALSE);
  for(Int_t i=0; i<6; i++)
    {
      TLine *line = GetLine(-217.5+i*87, 0, -217.5+i*87, 0.0125, 2, 2, 2);
      line->Draw();
    }
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.hit_z_per_Module_%s.png",run_config,trigName[kTrigType]));

  // dz distribution vs Delta_module
  list1->Clear();
  THnSparseF *hn = (THnSparseF*)f->Get(Form("hTrkDzDy_%s",trigName[kTrigType]));
  TH1F *hDzDeltaMod[4];
  TString legName2[4];
  for(Int_t i=0; i<4; i++)
    {
      if(i==0) 
	{
	  hDzDeltaMod[i] = (TH1F*)hn->Projection(1);
	  legName2[i] = "Sum";
	}
      else
	{
	  hn->GetAxis(4)->SetRange(i,i);
	  hDzDeltaMod[i] = (TH1F*)hn->Projection(1);
	  hn->GetAxis(4)->SetRange(0,-1);
	  legName2[i] = Form("|Mod_{trk}-Mod_{hit}|=%d",i-1);
	}
      hDzDeltaMod[i]->SetName(Form("hDz_DeltaMod_%d",i-1));
      list1->Add(hDzDeltaMod[i]);
    }
  c = drawHistos(list1,list1->At(0)->GetName(),Form("Au+Au %s: #Deltaz distribution of matched track-hit pair;#Deltaz (cm)",trigName[kTrigType]),kFALSE,0,100,kFALSE,1e-6,0.018,kFALSE,kTRUE,legName2,kTRUE,"",0.6,0.8,0.6,0.88,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.dz_vs_Delta_Module_%s.png",run_config,trigName[kTrigType]));

  // dz/dy distribution in each backleg
  const Double_t pt_cut = 1;
  TCanvas *cRes[2];
  for(Int_t i=0; i<2; i++)
    {
      cRes[i] = new TCanvas(Form("D%s_Backleg_All",title[i]),Form("D%s_Backleg_All",title[i]),1550,800);
      cRes[i]->Divide(6,5);
      for(Int_t j=0; j<30; j++)
	{
	  cRes[i]->cd(j+1);
	  TH2F *h2 = (TH2F*)f2->Get(Form("hD%sVsPt_BL%d",title[i],j+1));
	  TH1F *h1 = (TH1F*)h2->ProjectionY(Form("%s_pro",h2->GetName()),h2->GetXaxis()->FindFixBin(pt_cut+0.1),-1);
	  h1->Draw();
	}
    }

  TH2F *hResi[2][2];
  hResi[0][0] = (TH2F*)f1->Get("hTrkDz_di-muon");
  hResi[0][1] = (TH2F*)f2->Get("hDzVsPt_All");
  hResi[1][0] = (TH2F*)f1->Get("hTrkDy_di-muon");
  hResi[1][1] = (TH2F*)f2->Get("hDyVsPt_All");
  TList *list[2];

  for(Int_t i=0; i<2; i++)
    {
      list[i] = new TList;
      for(Int_t j=0; j<2; j++)
	{
	  hResi[i][j]->SetName(Form("%s_%s",hResi[i][j]->GetName(),legName[j].Data()));
	  draw2D(hResi[i][j]);
	  TH1F *h = (TH1F*)hResi[i][j]->ProjectionY(Form("%s_pro",hResi[i][j]->GetName()),hResi[i][j]->GetXaxis()->FindFixBin(pt_cut+0.1),-1);
	  h->Sumw2();
	  h->Scale(1./h->Integral());
	  list[i]->Add(h);
	}
      c = drawHistos(list[i],list[i]->At(0)->GetName(),Form("Au+Au %s: #Delta%s distribution using %s tracks%s;#Delta%s (cm)",trigName[kTrigType],title[i],trk_name[trk_index],hlt_name[hlt_index],title[i]),kFALSE,0,100,kTRUE,0,0.02,kFALSE,kTRUE,legName,kTRUE,Form("p_{T} > %1.1f GeV/c",pt_cut),0.6,0.8,0.6,0.8,kFALSE);
      if(save)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/%s.d%s_Pt%1.0fGeV_%s.png",run_config,title[i],pt_cut,trigName[kTrigType]));
    }     
}

//================================================
void MonteCarlo(const Int_t save = 0)
{
  const char *title[2] = {"z","y"};
  const Int_t type = 1;

  TFile *fin = TFile::Open("Rootfiles/AuAu200.RotateMTD.root","read");

  TH2F *hHitYBL = (TH2F*)fin->Get(Form("h%sVsBL_hit",title[type]));
  TH2F *hTrkZBL = (TH2F*)fin->Get(Form("h%sVsBL_track",title[type]));
  if(type==1)
    {
      hHitYBL->GetYaxis()->SetRangeUser(-50,50);
      hTrkZBL->GetYaxis()->SetRangeUser(-50,50);
    }
  c = draw2D(hHitYBL,Form("%s distribution of MTD hits per module;backleg*5+module",title[type]));
  c = draw2D(hTrkZBL,Form("%s distribution of projected tracks per module;backleg*5+module",title[type]));

  hHitYBL->Sumw2();
  hTrkZBL->Sumw2();

  TH2F *h2 = (TH2F*)fin->Get(Form("hD%sVsPt_All",title[type]));
  draw2D(h2);
  TH1F *hDyMea = (TH1F*)h2->ProjectionY(Form("hD%s_Measured",title[type]),h2->GetXaxis()->FindFixBin(1+0.1),-1);
  hDyMea->Scale(1./hDyMea->Integral());

  TH1F *hHitY, *hTrkY;
  if(type==1)
    {
      hHitY = (TH1F*)fin->Get(Form("h%s_All_hit",title[type]));
      hTrkY = (TH1F*)fin->Get(Form("h%s_All_track",title[type]));
    }
  else
    {
      for(Int_t j=0; j<150; j++)
	{
	  TH1F *h11 = (TH1F*)hHitYBL->ProjectionY(Form("h%s_Mod%d_hit",title[type],j+1),j+1,j+1);
	  TH1F *h12 = (TH1F*)hTrkZBL->ProjectionY(Form("h%s_Mod%d_track",title[type],j+1),j+1,j+1);
	  Int_t bin_move = h11->FindFixBin(87.*TMath::Abs(j%5-2))-h11->FindFixBin(0);
	  
	  TH1F *h11tmp = (TH1F*)h11->Clone(Form("%s_clone",h11->GetName()));
	  TH1F *h12tmp = (TH1F*)h12->Clone(Form("%s_clone",h12->GetName()));
	  h11tmp->Reset();
	  h12tmp->Reset();
	  for(Int_t ibin=1; ibin<=h11->GetNbinsX(); ibin++)
	    {
	      Int_t jbin = ibin + bin_move * (j%5-2<0? 1 : -1);
	      if(jbin>=1 && jbin<=h11->GetNbinsX())
		{
		  h11tmp->SetBinContent(jbin,h11->GetBinContent(ibin));
		  h11tmp->SetBinError(jbin,h11->GetBinError(ibin));
		  
		  h12tmp->SetBinContent(jbin,h12->GetBinContent(ibin));
		  h12tmp->SetBinError(jbin,h12->GetBinError(ibin));
		}
	    }

	  if(j==0) 
	    {
	      hHitY = (TH1F*)h11tmp->Clone(Form("h%s_Sum_hit",title[type]));
	      hTrkY = (TH1F*)h12tmp->Clone(Form("h%s_Sum_track",title[type]));
	      // h11tmp->SetMaximum(1e5);
	      // draw1D(h11tmp,"",kTRUE,kFALSE);
	    }
	  else
	    {
	      hHitY->Add(h11tmp);
	      hTrkY->Add(h12tmp);
	      // h11tmp->Draw("sames HIST");
	    }
	}
    }

  hHitY->Scale(1./hHitY->Integral());
  hTrkY->Scale(1./hTrkY->Integral());
  if(type==1)
    {
      hHitY->GetXaxis()->SetRangeUser(-50,50);
      hTrkY->GetXaxis()->SetRangeUser(-50,50);
    }
  c1 = draw1D(hHitY,Form("%s distribution of MTD hits",title[type]),kFALSE,kFALSE);
  c2 = draw1D(hTrkY,Form("%s distribution of projected tracks",title[type]),kFALSE,kFALSE);

  TH2F *hDyHitY = new TH2F(Form("hD%s_hit%s_MC",title[type],title[type]),Form("Monte Carlo: hit %s vs #Delta%s distribution; #Delta%s (cm);%s_{hit} (cm)",title[type],title[type],title[type],title[type]),200,-100,100,500,-250,250);
  TH1F *hDy = new TH1F(Form("hD%s_MC",title[type]),Form("Monte Carlo: #Delta%s distribution; #Delta%s (cm)",title[type],title[type]),200,-100,100);
  //Double_t center[20] = {-24.2,-24.2,-24.2,-24.2,-19.8,-19.8,-15.4,-11,-6.6,-2.2,2.2,6.6,11,15.4,19.8,19.8,24.2,24.2,24.2,24.2};
  const Int_t nExpr = 1e5;

  gRandom->SetSeed(0);
  Double_t hit_y, trk_y;
  Double_t x1,y1,x2,y2;
  for(Int_t i=0; i<nExpr; i++)
    {
      if(type==1)
	{
	  hit_y = hHitY->GetRandom();
	  trk_y = hTrkY->GetRandom();
	}
      else
	{
	  hHitYBL->GetRandom2(x1,y1);
	  hTrkZBL->GetRandom2(x2,y2);
	  Int_t bl1 = (hHitYBL->GetXaxis()->FindFixBin(x1)-1)/5;
	  Int_t mo1 = (hHitYBL->GetXaxis()->FindFixBin(x1)-1)%5;
	  Int_t bl2 = (hTrkZBL->GetXaxis()->FindFixBin(x2)-1)/5;
	  Int_t mo2 = (hTrkZBL->GetXaxis()->FindFixBin(x2)-1)%5;
	  if(TMath::Abs(bl1-bl2)<=1 && TMath::Abs(mo1-mo2)<=1)
	    {
	      hit_y = y1;
	      trk_y = y2;
	    }
	  else
	    continue;
	}
      hDyHitY->Fill(trk_y-hit_y,hit_y);
      hDy->Fill(trk_y-hit_y);
    }
  if(type==1)
    {
      hDyHitY->GetYaxis()->SetRangeUser(-30,30);
    }
  c = draw2D(hDyHitY,"",0.04,kFALSE);
  gPad->SetGridx();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/MC_Rotate_d%s_vs_hit%s_%s.png",title[type],title[type],trigName[kTrigType]));

  TH1F *hHitYMC = (TH1F*)hDyHitY->ProjectionY(Form("h%s_hit_MC",title[type]));
  c1->cd();
  hHitYMC->Scale(1./hHitYMC->Integral());
  hHitYMC->SetLineColor(2);
  hHitYMC->Draw("HIST sames");

  c3 = draw1D(hDyMea,Form("#Delta%s distribution of matched track-hit pairs",title[type]),kFALSE,kFALSE);
  c3->cd();
  hDy->Scale(1./hDy->Integral());
  hDy->SetLineColor(2);
  hDy->Draw("HIST sames");
  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetHeader("Rotate MTD");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hDyMea,"Data","L");
  leg->AddEntry(hDy,"Monte Carlo","L");
  leg->Draw();
  if(save) c3->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_RotateMTD/MC_Rotate_compare_d%s_%s.png",title[type],trigName[kTrigType]));
}

//================================================
void makeHisto(const Int_t save = 0)
{
  const char *title[2] = {"z","y"};
  THnSparseF *hn = (THnSparseF*)f->Get("hTrkDzDy_di-muon");
  TH2F *hRes_All[2];
  TH2F *hRes_BL[2][30];
  for(Int_t i=0; i<2; i++)
    {
      hRes_All[i] = (TH2F*)hn->Projection(i+1,0);
      hRes_All[i]->SetName(Form("hD%sVsPt_All",title[i]));
      for(Int_t j=0; j<30; j++)
	{
	  hn->GetAxis(3)->SetRange(j*5+1, j*5+5);
	  hRes_BL[i][j] = (TH2F*)hn->Projection(i+1,0);
	  hRes_BL[i][j]->SetName(Form("hD%sVsPt_BL%d",title[i],j+1));
	  hn->GetAxis(3)->SetRange(0,-1);
	}
    }

  const char *trkname[2] = {"track","hit"};
  THnSparseF *hYZ[2];
  hYZ[0] = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hYZ[1] = (THnSparseF*)f->Get(Form("hHitYZ_qa_%s",trigName[kTrigType]));
  TH2F *hProjVsBL[2][2];
  TH1F *hProj_All[2][2];
  TH1F *hProj_BL[2][2][30];

  for(Int_t i=0; i<2; i++)
    {
      for(Int_t k=0; k<2; k++)
	{
	  hProjVsBL[i][k] = (TH2F*)hYZ[i]->Projection(1-k,2);
	  hProjVsBL[i][k]->SetName(Form("h%sVsBL_%s",title[k],trkname[i]));

	  hProj_All[i][k] = (TH1F*)hYZ[i]->Projection(1-k);
	  hProj_All[i][k]->SetName(Form("h%s_All_%s",title[k],trkname[i]));

	  for(Int_t j=0; j<30; j++)
	    {
	      hYZ[i]->GetAxis(2)->SetRange(j*5+1, j*5+5);
	      hProj_BL[i][k][j] = (TH1F*)hYZ[i]->Projection(1-k);
	      hProj_BL[i][k][j]->SetName(Form("h%s_BL%d_%s",title[k],j+1,trkname[i]));
	      hYZ[i]->GetAxis(2)->SetRange(0,-1);
	    }
	}
    }

  if(save)
    {
      TFile *fout = TFile::Open("Rootfiles/AuAu200.RotateMTD.root","recreate");
      for(Int_t i=0; i<2; i++)
	{
	  hRes_All[i]->Write();
	}

      for(Int_t i=0; i<2; i++)
	{
	  for(Int_t k=0; k<2; k++)
	    {
	      hProjVsBL[i][k]->Write();
	      hProj_All[i][k]->Write();
	    }
	}

      for(Int_t i=0; i<2; i++)
	{
	  for(Int_t j=0; j<30; j++)
	    {
	      hRes_BL[i][j]->Write();
	    }
	}

      for(Int_t i=0; i<2; i++)
	{
	  for(Int_t k=0; k<2; k++)
	    {
	      for(Int_t j=0; j<30; j++)
		{
		  hProj_BL[i][k][j]->Write();
		}
	    }
	}
    }
}



