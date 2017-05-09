int epsilon2()
{
  gStyle->SetOptStat(0);
  TProfile *AuCu = new TProfile("AuCu", "AuCu", 29, 0, 29, 0, 1);  
  ifstream inputAuCu("./AuCu/epsilon_2.txt");
  double epsAuCu = 0.;
  for (int j = 0; j != 100; ++j)
    {
      for (int i = 0; i != 29; ++i)
        {
          inputAuCu >> epsAuCu;
          if (!inputAuCu.good()) break;
          AuCu->Fill(i, epsAuCu);
        }
    }
  inputAuCu.close();
  AuCu->SetMarkerStyle(31);
  AuCu->SetMarkerSize(2);
  AuCu->SetMarkerColor(1);
  AuCu->GetYaxis()->SetRangeUser(0, 1);
  AuCu->GetXaxis()->SetTitle("t [fm/c]");
  AuCu->GetXaxis()->SetTitleSize(0.035);
  AuCu->GetXaxis()->CenterTitle(true);
  AuCu->GetYaxis()->SetTitle("#epsilon_{2}");
  AuCu->GetYaxis()->SetTitleSize(0.035);
  AuCu->GetYaxis()->CenterTitle(true);
  AuCu->Draw(" P");

  TProfile *Au = new TProfile("Au", "Au", 29, 0, 29, 0, 1);
  ifstream inputAu("./Au/epsilon_2.txt");
  double epsAu = 0.;
  for (int j = 0; j != 100; ++j)
    {
      for (int i = 0; i != 29; ++i)
        {
          inputAu >> epsAu;
          if (!inputAu.good()) break;
          Au->Fill(i, epsAu);
        }
    }
  inputAu.close();
  Au->SetMarkerStyle(20);
  Au->SetMarkerSize(2);
  Au->SetMarkerColor(2);
  Au->Draw("P same");
  
  TProfile *Cu = new TProfile("Cu", "Cu", 29, 0, 29, 0, 1);  
  ifstream inputCu("./Cu/epsilon_2.txt");
  double epsCu = 0.;
  for (int j = 0; j != 100; ++j)
    {
      for (int i = 0; i != 29; ++i)
        {
          inputCu >> epsCu;
          if (!inputCu.good()) break;
          Cu->Fill(i, epsCu);
        }
    }
  inputCu.close();
  Cu->SetMarkerStyle(21);
  Cu->SetMarkerSize(2);
  Cu->SetMarkerColor(4);
  Cu->Draw("same P");

  TLegend *leg = new TLegend(0.67, 0.73, 0.86, 0.86);
  leg->AddEntry(AuCu, "AuCu", "lp");
  leg->AddEntry(Au, "Au-going", "lp");
  leg->AddEntry(Cu, "Cu-going", "lp");
  leg->Draw("same");

  AuCu->SetTitle("Au + Cu");
  return 0;
}
