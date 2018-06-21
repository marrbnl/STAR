//================================================
void filter()
{
  mode1();
}

//================================================
void mode1()
{
  TString filename = "Rootfiles/Run14_AuAu200.StudyLumiDep.root";
  TFile *fin = TFile::Open(filename.Data());
  TFile *fout = TFile::Open(Form("%s_filter.root",filename.Data()),"recreate");
  TKey *key;
  TIter nextkey(fin->GetListOfKeys());
  while ((key = (TKey*)nextkey())) 
    {
      TObject *obj = key->ReadObj();   
      if ( obj->IsA()->InheritsFrom( "TH2" ) ) continue;
      obj->Write();
    }
  fout->Close();
  fin->Close();
}

