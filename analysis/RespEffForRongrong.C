
const Int_t gMtdNBacklegs = 30;   // Total number of backlegs
const Int_t gMtdNModules  = 5;    // Total number of MTD modules per backleg

const Int_t gNPtbin = 15; // Pt bins for efficiency
const Double_t Ptbin[]={0.0,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,20.};

TH1D* hPtProBklMod[gMtdNBacklegs][gMtdNModules][3];
TH1D* hPtMatBklMod[gMtdNBacklegs][gMtdNModules][3];
TH1D* hPtEffBklMod[gMtdNBacklegs][gMtdNModules][3];
TH1D* hPtRefBklMod[gMtdNBacklegs][gMtdNModules][3];

TF1*  fitPtRefBklMod[gMtdNBacklegs][gMtdNModules][3];
TF1*  fitPtEffBklMod[gMtdNBacklegs][gMtdNModules][3];
TH1D* hFitChiSqBklMod[3];
TH1D* hFitSclFaBklMod[3];

void RespEffForRongrong(TString system="Good"){

	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	char name[256];
	bool printplot = 0;
	bool flag = -1;
	cout<< "Enter Flag : 0 : collision, 1 cosmic " << endl;
	cin >> flag;
	if(flag) sprintf(name,"Cosmic");
	else sprintf(name,"Collision");

	int mYear = -1;
	cout<< "Enter Year : 2013, 2014, or 2015 " << endl;
	cin >> mYear;
	if(mYear!=2013 && mYear!=2014 && mYear!=2015)
		return;

	TFile *fin;
	if(system.Contains("pp")){
		fin = TFile::Open(Form("./output/Run%d%sPPDay.root",mYear-2000,name));
	}
	else if(system.Contains("pAu")){
		fin = TFile::Open(Form("./output/Run%d%sPAuDay.root",mYear-2000,name));
	}
	else if(system.Contains("Good")) {
		fin = TFile::Open(Form("./output/Run%d%sGoodDay.root",mYear-2000,name));
	}
	else if(system.Contains("GoodWOCent")) {
		fin = TFile::Open(Form("./output/Run%d%sGoodDayWOCent.root",mYear-2000,name));
	}
	else {
		fin = TFile::Open(Form("./output/Run%d%s.root",mYear-2000,name));
	}

	for(int k=0;k<3;k++){
		for(int i=0;i<gMtdNBacklegs;i++){
			for(int j=0;j<gMtdNModules;j++){
				hPtProBklMod[i][j][k] = (TH1D*) fin->Get(Form("hPtMtdProBkl%d_Mod%d_%d",i,j,k));
				hPtMatBklMod[i][j][k] = (TH1D*) fin->Get(Form("hPtMtdMatBkl%d_Mod%d_%d",i,j,k));
				hPtProBklMod[i][j][k] = (TH1D*) hPtProBklMod[i][j][k]->Rebin(gNPtbin,Form("hPtMtdProBkl%d_Mod%d_%d",i,j,k),Ptbin);
				hPtMatBklMod[i][j][k] = (TH1D*) hPtMatBklMod[i][j][k]->Rebin(gNPtbin,Form("hPtMtdMatBkl%d_Mod%d_%d",i,j,k),Ptbin);
				hPtEffBklMod[i][j][k] = new TH1D(Form("hPtMtdEffBkl%d_Mod%d_%d",i,j,k),"",gNPtbin,Ptbin);
				hPtRefBklMod[i][j][k] = new TH1D(Form("hPtMtdRefBkl%d_Mod%d_%d",i,j,k),"",gNPtbin,Ptbin);
			}
		}
	}


	for(int k=0;k<3;k++){
		for(int i=0;i<gMtdNBacklegs;i++){
			for(int j=0;j<gMtdNModules;j++){
				hPtEffBklMod[i][j][k] ->Divide(hPtMatBklMod[i][j][k],hPtProBklMod[i][j][k],1.,1.,"B");
			}
		}
	}

	TH1D* hFitChiSqBklMod[3];
	TH1D* hFitSclFaBklMod[3];
	for(int k=0;k<3;k++){
		hFitChiSqBklMod[k] = new TH1D(Form("hFitChiSqBklMod%d",k),"Fit qualtiy vs Module;(BacklegID-1)x5+ModuleID;#chi^{2}/NDF",150,0.,150);
		hFitSclFaBklMod[k] = new TH1D(Form("hFitSclFaBklMod%d",k),"Cut Parameter Dependece of Efficiency ;(BacklegID-1)x5+ModuleID;Efficiency by Fit ((p_{T}>5 GeV/c))",150,0.,150);
		for(int i=0;i<gMtdNBacklegs;i++){
			for(int j=0;j<gMtdNModules;j++){
				fitPtRefBklMod[i][j][k] = new TF1(Form("fitPtMtdRefBkl%d_Mod%d_%d",i,j,k),"[0]",1.,20.);
				fitPtEffBklMod[i][j][k] = new TF1(Form("fitPtMtdRefBkl%d_Mod%d_%d",i,j,k),"[0]",5.,20.);

				bool isbklok=false;
				if(mYear==2013){
					if((i+1)<8||(i+1)==10||(i+1)==22||(i+1)>24) isbklok = true;
				}
				else if(mYear>=2014){
					if((i+1)<9||(i+1)==10||(i+1)==22||(i+1)>23) isbklok = true;
					else if(((i+1)>10&&(i+1)<22) && ((j+1)>1 && (j+1)<5) ) isbklok = true;
				}

				if(isbklok){
				  //hPtRefBklMod[i][j][k] ->Fit(fitPtRefBklMod[i][j][k],"R0QMI");
				  //double chindf = fitPtRefBklMod[i][j][k]->GetChisquare()/fitPtRefBklMod[i][j][k]->GetNDF();
				  //hFitChiSqBklMod[k] ->SetBinContent(i*5+j+1,chindf);

					hPtEffBklMod[i][j][k] ->Fit(fitPtEffBklMod[i][j][k],"R0QMI");
					hFitSclFaBklMod[k] ->SetBinContent(i*5+j+1,fitPtEffBklMod[i][j][k]->GetParameter(0));
					hFitSclFaBklMod[k] ->SetBinError  (i*5+j+1,fitPtEffBklMod[i][j][k]->GetParError (0));
				}
			}
		}
	}


	const int mcol[] = {1,2,4,3,6,7};
	const int msty[] = {26,25,24,3,6,7};
	const double effmax = 1.1;
	const double effmin = 0.0;
	const double refmax = 3.5;
	const double refmin = 0.0;
	for(int k=0;k<3;k++){
		for(int i=0;i<gMtdNBacklegs;i++){
			for(int j=0;j<gMtdNModules;j++){
				hPtEffBklMod[i][j][k]->SetMaximum(effmax);
				hPtEffBklMod[i][j][k]->SetMinimum(effmin);
				hPtEffBklMod[i][j][k]->SetLineColor(mcol[k]);
				hPtEffBklMod[i][j][k]->SetMarkerStyle(msty[k]);
				hPtEffBklMod[i][j][k]->SetMarkerColor(mcol[k]);
				hPtEffBklMod[i][j][k]->SetTitle(Form("Bkl=%d, Mod=%d;p_{T} GeV/c;Response Efficiency",i+1,j+1));
				hPtEffBklMod[i][j][k]->SetTitleSize(0.04,"xyz");

				hPtRefBklMod[i][j][k]->SetMaximum(refmax);
				hPtRefBklMod[i][j][k]->SetMinimum(refmin);
				hPtRefBklMod[i][j][k]->SetLineColor(mcol[k]);
				hPtRefBklMod[i][j][k]->SetMarkerStyle(msty[k]);
				hPtRefBklMod[i][j][k]->SetMarkerColor(mcol[k]);
				hPtRefBklMod[i][j][k]->SetTitle(Form("Bkl=%d, Mod=%d;p_{T} GeV/c;Relative Efficiency",i+1,j+1));
				hPtRefBklMod[i][j][k]->SetTitleSize(0.04,"xyz");
			}
		}
	}

	TLegend *leg = new TLegend(0.2,0.2,0.8,0.8);
	leg->SetFillStyle(0);
	leg->SetLineColor(0);
	leg->AddEntry(hPtEffBklMod[0][0][0],"Track Quality","p");
	leg->AddEntry(hPtEffBklMod[0][0][1],"+tof2tof","p");
	leg->AddEntry(hPtEffBklMod[0][0][2],"+track matching","p");

	TLegend *leg2 = new TLegend(0.4,0.4,0.6,0.6);
	leg2->SetFillStyle(0);
	leg2->SetLineColor(0);
	leg2->AddEntry(hPtEffBklMod[0][0][0],"Track Quality","p");
	leg2->AddEntry(hPtEffBklMod[0][0][1],"+tof2tof","p");
	leg2->AddEntry(hPtEffBklMod[0][0][2],"+track matching","p");

	TCanvas *canv[30];
	for(int i=0;i<gMtdNBacklegs;i++){

		bool isbklok=false;
		if(mYear==2013){
			if((i+1)<8||(i+1)==10||(i+1)==22||(i+1)>24) isbklok = true;
		}
		else if(mYear>=2014){
			if((i+1)<9|| ( (i+1)>=10 && (i+1)<=22) ||(i+1)>23) isbklok = true;
		}

		if(isbklok){
			canv[i] = new TCanvas(Form("canv%d",i+1),Form("canv%d",i+1),0,0,1200,600);
			canv[i] ->Divide(3,2);
			for(int j=0;j<gMtdNModules;j++){
				if((i+1)==7 && j==4) continue;
				canv[i] ->cd(j+1);
				for(int k=0;k<3;k++){
					hPtEffBklMod[i][j][k]->Draw("same"); 
				}
				canv[i] ->cd(6);
				leg->Draw("same");
			}
			if(printplot)canv[i] ->Print(Form("./fig/Run%dCosmicCutParDepBkl%d.pdf",mYear-2000,i+1));
		}
	}

	for(int k=0;k<3;k++){
		hFitSclFaBklMod[k]->SetMaximum(1.0);
		hFitSclFaBklMod[k]->SetMinimum(0.0);

		hFitSclFaBklMod[k]->SetLineColor  (mcol[k]);
		hFitSclFaBklMod[k]->SetMarkerColor(mcol[k]);
		hFitSclFaBklMod[k]->SetMarkerStyle(msty[k]);
	}

	TCanvas *canvscl = new TCanvas("canvscl","canvscl");
	canvscl->cd();
	hFitSclFaBklMod[0]->Draw("same");
	hFitSclFaBklMod[1]->Draw("same");
	hFitSclFaBklMod[2]->Draw("same");
	leg2->Draw("same");
	if(printplot)canvscl ->Print(Form("./fig/Run%dCosmicCutParDepModByMod.pdf",mYear-2000));


	//	TFile *out;
	//	if(system.Contains("pp")){
	//		out = TFile::Open(Form("./Run%dRawResponseEff%sPPDay.root",mYear-2000,name),"recreate");
	//	}
	//	else if(system.Contains("pAu")){
	//		out = TFile::Open(Form("./Run%dRawResponseEff%sPAuDay.root",mYear-2000,name),"recreate");
	//	}
	//	else if(system.Contains("Good")) {
	//		out = TFile::Open(Form("./Run%dRawResponseEff%sGoodDay.root",mYear-2000,name),"recreate");
	//	}
	//	else {
	//		out = TFile::Open(Form("./Run%dRawResponseEff%s.root",mYear-2000,name),"recreate");
	//	}
	//	out->cd();
	//
	//	for(int k=0;k<3;k++){
	//		hFitSclFaBklMod[k]->Write();
	//		for(int i=0;i<gMtdNBacklegs;i++){
	//			for(int j=0;j<gMtdNModules;j++){
	//				hPtProBklMod[i][j][k]->Write(); 
	//				hPtMatBklMod[i][j][k]->Write(); 
	//				hPtEffBklMod[i][j][k]->Write(); 
	//				hPtRefBklMod[i][j][k]->Write(); 
	//			}
	//		}
	//	}

}
