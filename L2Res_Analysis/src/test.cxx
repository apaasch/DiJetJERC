{

  TH2D* h2 = new TH2D("test", "test", 10,0,10,10,0,10);
  h2->Fill(1,3);
  h2->Fill(1,4);
  h2->Fill(1,7);
  h2->Fill(2,3);
  h2->Fill(2,3);
  h2->Fill(2,4);
  TProfile* prof = (TProfile*) h2->ProfileX();

  prof->Draw();

  TH1D* h1 = new TH1D();






}
