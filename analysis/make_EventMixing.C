const char *run_config = "EventMixing.PrimTrk.ClosePrimVtx";

//================================================
void make_EventMixing()
{
  ifstream file;
  file.open(Form("output/file.jpsi.AuAu200.Run14.%s.list",run_config));
  char filename[256];
  Double_t pt_cuts_1[6] = {1.0,1.2,1.5,2.0,2.5,3.0};
  Double_t pt_cuts_2[2] = {1.0,1.2};

  TH1F *hUS[6][2];
  TH1F *hLS[6][2];
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  TFile *fout = TFile::Open(Form("output/jpsi.AuAu200.Run14.%s.root",run_config),"recreate");

  Int_t counter = 0;
  char hname[256], hname2[256];
  while(!file.eof())
    {
      file.getline(filename,256);
      TString name = filename;
      if(name.Length()>0)
	{
	  TFile *f = TFile::Open(name.Data(),"read");
	  cout << "Open file " << counter << endl;
	  THnSparseF *hnInvMass[3];
	  for(Int_t j=0; j<3; j++)
	    {
	      sprintf(hname,"%s_%s",hName[j],trigName[kTrigType]);
	      hnInvMass[j] = (THnSparseF*)f->Get(hname);
	      for(Int_t i=0; i<6; i++)
		{
		  for(Int_t k=0; k<2; k++)
		    {
		      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt_cuts_1[i]+0.01,100);
		      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt_cuts_2[k]+0.01,100);
		      TH1F *h = (TH1F*)hnInvMass[j]->Projection(0);
		      sprintf(hname,"%s_%s_InvMass_pt1_%1.1f_pt2_%1.1f_tmp%d",hName[j],trigName[kTrigType],pt_cuts_1[i],pt_cuts_2[k],counter);
		      h->SetName(hname);
		      //cout << hname << "  =  " << h->GetName() << endl;
		      h->Sumw2();
		      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
		      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
		      if(counter==0)
		      	{
		      	  sprintf(hname,"US_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]);
		      	  sprintf(hname2,"LS_InvMass_pt1_%1.1f_pt2_%1.1f_%s",pt_cuts_1[i],pt_cuts_2[k],trigName[kTrigType]);
		      	  if(j==0) { hUS[i][k] = (TH1F*)h->Clone(hname); hUS[i][k]->SetDirectory(fout); }
		      	  else if(j==1) { hLS[i][k] = (TH1F*)h->Clone(hname2); hLS[i][k]->SetDirectory(fout); }
		      	  else hLS[i][k]->Add(h);
		      	}
		      else
		      	{
		      	  //cout << hUS[i][k]->GetName() << endl;
		      	  if(j==0) hUS[i][k]->Add(h);
		      	  else     hLS[i][k]->Add(h);
		      	}
		      //cout << hUS[i][k]->GetName() << endl;
		    }
		}
	    }
	  counter ++;
	  //if(counter==3) break;
	  f->Close();
	}
    }

  fout->cd();
  for(Int_t i=0; i<6; i++)
    {
      for(Int_t k=0; k<2; k++)
	{
	  hUS[i][k]->Write();
	  hLS[i][k]->Write();
	}
    }
  fout->Close();
}
