#include <stdio.h>     
#include <stdlib.h>    

#include <TTree.h>

#include "TomiyaCorrection.h"
#include "pdmap.h"

void make_map2(TString inputname="1m/all1m.csv", TString outname="map.pdf", Double_t hv = 1600.){//2000.) {
 // hv is rougly in a range of 1600-2200. Taking 1600 V for the largest impact.

 gErrorIgnoreLevel = 5000;
  //gROOT->Reset();
  gStyle->SetOptStat(1);

  TTree *t = new TTree("t", "tree");
  t->ReadFile(inputname.Data(), "x/D:y/D:z/D:Bx/D:By/D:Bz/D:Bp/D:faceid/I",',');

  Double_t Bxyz[5][NPD]; // xyz, perpendicular, |B|
  Double_t Bpmt[3][NPD]; // xyz on PMT
  Double_t eff[NPD];
  Double_t effdir[3*2*1];
  Double_t xyz[NPD][3];
  Int_t faceid[NPD];

#ifdef takeminimum
  TH1D *heff = new TH1D("hEffLoss","Relative Detection Efficiency", 2020,0.6,1.01);
  heff->SetXTitle("Relative Detection Efficiency");
  TH1D *heffud[2];
  TH1D *hefftb[2];
  heffud[0] = new TH1D("hEffLossUp","Relative Detection Efficiency", 2020,0.6,1.01);
  heffud[1] = new TH1D("hEffLossDown","Relative Detection Efficiency", 2020,0.6,1.01);
  hefftb[0] = new TH1D("hEffLossTop","Relative Detection Efficiency", 2020,0.6,1.01);
  hefftb[1] = new TH1D("hEffLossBottom","Relative Detection Efficiency", 2020,0.6,1.01);
#else
  TH1F *heff = new TH1F("hEffLoss","Relative Detection Efficiency", 44,0.6,1.03);
  TH1F *hefftwoper = new TH1F("hEffLossdosper","Relative Detection Efficiency", 44,0.6,1.03);
  heff->SetXTitle("Relative Detection Efficiency");
  TH1F *heffud[2];
  TH1F *hefftb[2];
  heffud[0] = new TH1F("hEffLossUp","Relative Detection Efficiency", 44,0.6,1.03);
  heffud[1] = new TH1F("hEffLossDown","Relative Detection Efficiency", 44,0.6,1.03);
  hefftb[0] = new TH1F("hEffLossTop","Relative Detection Efficiency", 44,0.6,1.03);
  hefftb[1] = new TH1F("hEffLossBottom","Relative Detection Efficiency", 44,0.6,1.03);
#endif
  heff->SetLineColor(kBlack);
  heffud[0]->SetLineColor(kRed);
  heffud[1]->SetLineColor(kBlue);
  hefftb[0]->SetLineColor(kMagenta);
  hefftb[1]->SetLineColor(kCyan+2);

  TGraph * g[2];
  g[0] = new TGraph();
  g[1] = new TGraph();
  TGraph * gphi = new TGraph();

///   ifstream fin;
///   std::string str_buf;
///   std::string str_conma_buf;
///   fin.open (filename, ios::in);
///   if (!fin) {
///      std::cout << "Cannot open B-field file : " << filename << " " << std::endl;
///      exit(1);
///   }
///   Int_t pmtid = 0;
///   //while( ! fin.eof() ) {
///   while(  ) {
///     fin >> xyz[pmtid][0] >>  xyz[pmtid][1] >> xyz[pmtid][2] >> Bxyz[0][pmtid] >>  Bxyz[1][pmtid] >>  Bxyz[2][pmtid] >> Bxyz[3][pmtid] >> faceid[pmtid];
///     pmtid++;
///   }
///   fin.close();
///   fin.clear();

   TCanvas * c = new TCanvas("c", "c", 600, 600);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0);

   gStyle->SetPalette(1);
      c->SetTickx(1);
      c->SetTicky(1);
      c->SetRightMargin(0.16);
      //c->SetLeftMargin(0.056);

   c->Divide(2,2);

  TString outnamehead = outname;
  outnamehead.ReplaceAll(".pdf", "");
  //////outname = Form("%s.pdf", filename);
  c->Print(Form("%s[", outname.Data()));

   c->cd(1);
   t->Draw("faceid","","");
   for (int i = 0; i < t->GetEntries() ; i++ ){
      faceid[i] = t->GetV1()[i];
   }
   t->Draw("z:x:y:Bp","","colz");

   for (int i = 0; i < t->GetEntries() ; i++ ){
      xyz[i][2] = t->GetV1()[i]; // x
      xyz[i][0] = t->GetV2()[i]; // y
      xyz[i][1] = t->GetV3()[i]; // z
      //if ( xyz[2][i] < -1. * WallZ + 0.05 ) {
      //   faceid[i] = 3;
      //} else if ( xyz[2][i] > 1. * WallZ - 0.05 ) {
      //   faceid[i] = 1;
      //} else {
      //   faceid[i] = 2;
      //}
      Bxyz[3][i] = t->GetV4()[i];
   }

   c->cd(2);
   t->Draw("z:x:y:Bx","","colz");

   for (int i = 0; i < t->GetEntries() ; i++ ){
      Bxyz[0][i] = t->GetV4()[i];
   }

   c->cd(3);
   t->Draw("z:x:y:By","","colz");

   for (int i = 0; i < t->GetEntries() ; i++ ){
      Bxyz[1][i] = t->GetV4()[i];
   }

   c->cd(4);
   t->Draw("z:x:y:Bz","","colz");

   //Print B-map
   c->Print(Form("%s", outname.Data()));

   int updownid = -1;
   int topbottomid = -1;
   for (int i = 0; i < t->GetEntries() ; i++ ){
      Bxyz[2][i] = t->GetV4()[i];
      Bxyz[4][i] = TMath::Sqrt(Bxyz[0][i]*Bxyz[0][i]+Bxyz[1][i]*Bxyz[1][i]+Bxyz[2][i]*Bxyz[2][i]);
      //cout << "XYZ " << xyz[i][0] 
      //<< " " << xyz[i][1]
      //<< " " << xyz[i][2] << endl;
      switch(faceid[i]) {
          case 1: // top
             Bpmt[0][i] = Bxyz[0][i];
             Bpmt[1][i] = Bxyz[1][i];
             Bpmt[2][i] = Bxyz[2][i];
             updownid = 0;
             topbottomid = 0;
                 break;
          case 2: // barrel
             Bpmt[0][i] = TMath::Cos(TMath::ATan2(xyz[i][1],xyz[i][0]))*Bxyz[1][i]-TMath::Sin(TMath::ATan2(xyz[i][1],xyz[i][0]))*Bxyz[0][i];
             Bpmt[1][i] = Bxyz[2][i];
             Bpmt[2][i] = TMath::Cos(TMath::ATan2(xyz[i][1],xyz[i][0]))*Bxyz[0][i]+TMath::Sin(TMath::ATan2(xyz[i][1],xyz[i][0]))*Bxyz[1][i];
             topbottomid = -1;
             if (xyz[i][1] > 0) {
                updownid = 0;
             } else {
                updownid = 1;
             }
                 break;
          case 3: // bottom
             Bpmt[0][i] = Bxyz[0][i];
             Bpmt[1][i] = Bxyz[1][i];
             Bpmt[2][i] = -1.*Bxyz[2][i];
             updownid = 1;
             topbottomid = 1;
                 break;
      }
      //eff[i] = 1./ tomiyaCorrection(hv, Bpmt[0][i], Bpmt[1][i], Bpmt[2][i], false);
      effdir[0] = tomiyaCorrection(hv, Bpmt[0][i], Bpmt[1][i], Bpmt[2][i], false);
      effdir[1] = tomiyaCorrection(hv,-Bpmt[1][i], Bpmt[0][i], Bpmt[2][i], false);
      effdir[2] = tomiyaCorrection(hv,-Bpmt[0][i],-Bpmt[1][i], Bpmt[2][i], false);
      effdir[3] = tomiyaCorrection(hv, Bpmt[1][i],-Bpmt[0][i], Bpmt[2][i], false);
      //eff[i] = effdir[0];
      eff[i] = TMath::MinElement(4, effdir);
#ifdef takeminimum
      // Take minimum loss only in 4 directions on PMT
      heff->Fill(eff[i]);
      heffud[updownid]->Fill(eff[i]);
      if (topbottomid != -1) {
         hefftb[topbottomid]->Fill(eff[i]);
      }
#else
      // Fill 4 directions on PMT
      for (int idir = 0; idir < 4; idir++) {
         heff->Fill(effdir[idir]);
         heffud[updownid]->Fill(effdir[idir]);
         if (topbottomid != -1) {
            hefftb[topbottomid]->Fill(effdir[idir]);
         }
      }
#endif

      g[0]->SetPoint(g[0]->GetN(), Bxyz[3][i], eff[i]);
      g[1]->SetPoint(g[1]->GetN(), Bxyz[4][i], eff[i]);
      gphi->SetPoint(gphi->GetN(), TMath::ATan2(xyz[i][2],TMath::Sqrt(xyz[i][0]*xyz[i][0]+xyz[i][1]*xyz[i][1]))*TMath::RadToDeg(), eff[i]);
   }


   c->cd(1)->Clear();
   //t->Draw("TMath::Cos(TMath::ATan2(y,x))*Bx+TMath::Sin(TMath::ATan2(y,x))*By:Bp","faceid==2","");
   //t->Draw("TMath::Cos(TMath::ATan2(y,x))*By-TMath::Sin(TMath::ATan2(y,x))*Bx:Bp","faceid==2","");
   //t->Draw("TMath::Sqrt(TMath::Power(TMath::Cos(TMath::ATan2(y,x))*Bx+TMath::Sin(TMath::ATan2(y,x))*By, 2)+Bz*Bz):Bp","faceid==2","");

   //on vertical (PMT z) on barrel
   t->Draw("TMath::Cos(TMath::ATan2(y,x))*Bx+TMath::Sin(TMath::ATan2(y,x))*By:TMath::Sqrt(Bx*Bx+By*By+Bz*Bz-Bp*Bp)","faceid==2","");
   //on barrel face and Bp (|B| on PMT xy)
   //t->Draw("TMath::Sqrt(TMath::Power(TMath::Cos(TMath::ATan2(y,x))*By-TMath::Sin(TMath::ATan2(y,x))*Bx, 2)+Bz*Bz):Bp","faceid==2","");

   c->cd(2)->Clear();
   t->Draw("TMath::Sqrt(Bx*Bx+By*By):Bp","faceid==1","");

   c->cd(3)->Clear();
   t->Draw("TMath::Sqrt(Bx*Bx+By*By):Bp","faceid==3","");

   //c->Print(Form("%s", outname.Data()));
   c->cd(0)->Clear();
   Long64_t over100pmt = t->Draw("TMath::Sqrt(Bx*Bx+By*By+Bz*Bz)","TMath::Sqrt(Bx*Bx+By*By+Bz*Bz)>100");
   Long64_t totalpmt = t->Draw("TMath::Sqrt(Bx*Bx+By*By+Bz*Bz)>>hBtotal(250,0,500)", "", "");
   t->Draw("TMath::Sqrt(Bx*Bx+By*By+Bz*Bz)>>hBtotalTop(250,0,500)","faceid==1","same");
   t->Draw("TMath::Sqrt(Bx*Bx+By*By+Bz*Bz)>>hBtotalBarrel(250,0,500)","faceid==2","same");
   t->Draw("TMath::Sqrt(Bx*Bx+By*By+Bz*Bz)>>hBtotalBottom(250,0,500)","faceid==3","same");
   ((TH1*)(gDirectory->Get("hBtotalTop")))->SetLineColor(kRed);
   ((TH1*)(gDirectory->Get("hBtotalBarrel")))->SetLineColor(kGreen+2);
   ((TH1*)(gDirectory->Get("hBtotalBottom")))->SetLineColor(kBlue);
   c->Update();
   c->Print(Form("%s", outname.Data()));
   cout << "Num of PMTs exceeding 100 mG " << over100pmt << " / " << totalpmt << endl;

   t->Draw("Bp>>hBtotal", "", "");
   t->Draw("Bp>>hBtotalTop","faceid==1","same");
   t->Draw("Bp>>hBtotalBarrel","faceid==2","same");
   t->Draw("Bp>>hBtotalBottom","faceid==3","same");
   c->Print(Form("%s", outname.Data()));

   c->cd(0)->Clear();
   //DrawGeometry(10, eff, 0.97, 1.01, "", xyz);
   DrawGeometry(10, eff, 0.84, 1.01, "", xyz);
   c->Print(Form("%s", outname.Data()));

   c->cd(0)->Clear();
   //DrawGeometry(1, eff, 0.97, 1.01, "", xyz);
   DrawGeometry(1, eff, 0.84, 1.01, "", xyz);
   c->Print(Form("%s", outname.Data()));
   c->SetLogy();
   heff->Draw();
   heffud[0]->Draw("same");
   hefftb[0]->Draw("same");
   heffud[1]->Draw("same");
   hefftb[1]->Draw("same");
   c->Print(Form("%s", outname.Data()));
   c->SaveAs(Form("%s_eff.pdf", outnamehead.Data()));
   c->SetLogy(0);
   Double_t tbpeak[2], udpeak[2], tbmean[2], udmean[2];
   tbpeak[0] = hefftb[0]->GetBinCenter(hefftb[0]->GetMaximumBin());
   tbpeak[1] = hefftb[1]->GetBinCenter(hefftb[1]->GetMaximumBin());
   udpeak[0] = heffud[0]->GetBinCenter(heffud[0]->GetMaximumBin());
   udpeak[1] = heffud[1]->GetBinCenter(heffud[1]->GetMaximumBin());
   tbmean[0] = hefftb[0]->GetMean();
   tbmean[1] = hefftb[1]->GetMean();
   udmean[0] = heffud[0]->GetMean();
   udmean[1] = heffud[1]->GetMean();
   Double_t asymtbpeak, asymudpeak, asymtbmean, asymudmean;
   asymtbpeak = (hefftb[0]->GetBinCenter(hefftb[0]->GetMaximumBin()) - hefftb[1]->GetBinCenter(hefftb[1]->GetMaximumBin()))/(hefftb[0]->GetBinCenter(hefftb[0]->GetMaximumBin()) + hefftb[1]->GetBinCenter(hefftb[1]->GetMaximumBin()));
   asymudpeak = (heffud[0]->GetBinCenter(heffud[0]->GetMaximumBin()) - heffud[1]->GetBinCenter(heffud[1]->GetMaximumBin()))/(heffud[0]->GetBinCenter(heffud[0]->GetMaximumBin()) + heffud[1]->GetBinCenter(heffud[1]->GetMaximumBin()));
   asymtbmean = (hefftb[0]->GetMean() - hefftb[1]->GetMean())/(hefftb[0]->GetMean() + hefftb[1]->GetMean());
   asymudmean = (heffud[0]->GetMean() - heffud[1]->GetMean())/(heffud[0]->GetMean() + heffud[1]->GetMean());
   cout << endl;
   cout << inputname << endl;
   cout << "Assymmetry " << endl;
   cout << " Top-Bottom " << asymtbpeak << " by Peak,   " << asymtbmean << " by mean" << endl;
   cout << " Up-Down " << asymudpeak << " by Peak,   " << asymudmean << " by mean" << endl;
   cout << "Top/Bottom peaks, Up/Down peaks, Top/Bottom means, Up/Down means" << endl;
   cout <<
   tbpeak[0] << " " << 
   tbpeak[1] << " " << 
   udpeak[0] << " " << 
   udpeak[1] << " " << 
   tbmean[0] << " " << 
   tbmean[1] << " " << 
   udmean[0] << " " << 
   udmean[1] << " " <<  endl;
   cout << endl;

  cout << "Total Mean:    " << (1.-heff->GetMean()     )*100.  << "% RMS:" <<  100.*heff->GetRMS() << endl;
  cout << "  Up Mean:     " << (1.-heffud[0]->GetMean())*100.  << "% RMS:" <<  100.*heffud[0]->GetRMS() << endl;
  cout << "  Down Mean:   " << (1.-heffud[1]->GetMean())*100.  << "% RMS:" <<  100.*heffud[1]->GetRMS() << endl;
  cout << "  Top Mean:    " << (1.-hefftb[0]->GetMean())*100.  << "% RMS:" <<  100.*hefftb[0]->GetRMS() << endl;
  cout << "  Bottom Mean: " << (1.-hefftb[1]->GetMean())*100.  << "% RMS:" <<  100.*hefftb[1]->GetRMS() << endl;
  cout << heff->GetEntries() << endl;



   g[0]->SetMarkerStyle(1);//kCircle);
   g[0]->Draw("AP");
   c->Print(Form("%s", outname.Data()));

   g[1]->SetMarkerStyle(1);//kCircle);
   g[1]->Draw("AP");
   c->Print(Form("%s", outname.Data()));

   gphi->SetMarkerStyle(1);//kCircle);
   gphi->Draw("AP");
   c->Print(Form("%s", outname.Data()));

  c->Print(Form("%s]", outname.Data()));
  
  





return;

}
