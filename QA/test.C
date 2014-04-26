void test()
{


  TTree *tree = new TTree("myTree","Testing tree");
  Double_t x,y,z;
  tree->Branch("x",&x,"x/D");
  tree->Branch("y",&y,"y/D");
  tree->Branch("z",&z,"z/D");
  for(Int_t i=0; i<1; i++)
    {
      x = gRandom->Rndm();
      y = gRandom->Rndm();
      z = gRandom->Rndm();
      tree->Fill();
    }

  TFile *fout = TFile::Open("tree.test.root","recreate");
  tree->Write();
  cout << fout->GetCompressionFactor() << endl;
  fout->Close();
}
