// $Id$
// Version 1.0a

#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TVirtualHistPainter.h>
#include <TPad.h>
#include <TMarker.h>
#include <TBox.h>
#include <TArc.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TColor.h>
#include <TTimeStamp.h>
#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGaxis.h>
#include <TPave.h>
#include <TPaveLabel.h>
#include <TPDF.h>
#include <TPolyMarker3D.h>
#include <TRandom.h>

//#include "pdinfo.h"

const int NPMT   = 39844;
const int NPD    = 39844;

#define INNERF  0
#define OUTERF  1
#define UPF     2
#define DOWNF   3
#define TOPF    4
#define BOTTOMF 5

#define HKPDMARKER 25
//#define HKPDMARKER 21

#ifndef PDMAP_h
#define PDMAP_h

const Double_t WallR = 3240*2.; // in cm
const Double_t WallZ = 3290.00; // in cm

//______________________________________________________________________________
Double_t GetWallPos(Double_t pos[3], Double_t vec[3], Double_t wallpos[]) {
   const Double_t zwall = WallZ;
   const Double_t rwall = WallR*0.5;
   Double_t factor = 1;
   if (vec[0] != 0 && vec[1] !=0) {
      factor = TMath::Sqrt((pos[0]*vec[0]+pos[1]*vec[1])*(pos[0]*vec[0]+pos[1]*vec[1])
               - (vec[0]*vec[0]+vec[1]*vec[1]) * ((pos[0]*pos[0]+pos[1]*pos[1])-rwall*rwall))
               - (pos[0]*vec[0]+pos[1]*vec[1]);
   } else { return -1e6; }
   factor /= (vec[0]*vec[0]+vec[1]*vec[1]);
   if (TMath::Abs(factor * vec[2] + pos[2]) >= zwall) {
      factor = (zwall - TMath::Abs(pos[2])) / vec[2];
   }
   for (Int_t i = 0; i < 3; i++) {
      wallpos[i] = factor * vec[i] + pos[i];
   }
   //cout << "Wall " << rwall << " " << TMath::Sqrt(wallpos[0]*wallpos[0]+wallpos[1]*wallpos[1]) << endl;
   return factor;
}

//______________________________________________________________________________
Double_t GetXWall(Double_t x, Double_t y, Double_t z) {
   const Double_t zwall = WallZ;
   const Double_t rwall = WallR*0.5;
   if (z >= zwall - 0.05) {
      return y;
   } else if (z <= -1*zwall + 0.05) {
      return y;
   } else {
      return TMath::ATan2(y,-x)*rwall;
   }
}

//______________________________________________________________________________
Double_t GetYWall(Double_t x, Double_t y, Double_t z) {
   const Double_t zwall = WallZ;
   const Double_t rwall = WallR*0.5;
   if (z >= zwall - 0.05) {
      return x+zwall+rwall;
   } else if (z <= -1*zwall + 0.05) {
      return -x-zwall-rwall;
   } else {
      return z;
   }
}

/*
//______________________________________________________________________________
Double_t vecring(Double_t cang, Double_t phi, Double_t vecaxis[3], Double_t vec[])
{
   vec[0] = TMath::Sin(cang * TMath::DegToRad()) * TMath::Cos(phi*TMath::DegToRad());
   vec[1] = TMath::Sin(cang * TMath::DegToRad()) * TMath::Sin(phi*TMath::DegToRad());
   vec[2] = TMath::Cos(cang * TMath::DegToRad());
   TVector vecorg = TVector(vec[0], vec[1], vec[2]);
   TVector x1 = TVector(1, 0, 0);
   TVector x2 = TVector(0, 1, 0);
   TVector x3 = TVector(0, 0, 1);
   TRotation r;
   r.MakeBasis();
}
*/

//______________________________________________________________________________
Int_t GetValueColor(Double_t zc = 0, Double_t wmin = 0, Double_t wmax = 1)
{
   // for old ROOT and to avoid bug in current ROOT revision
   // Returns the color index of the given z value
   //
   // This function should be used after an histogram has been plotted with the
   // option COL or COLZ like in the following example:
   //
   //   h2->Draw("COLZ");
   //   gPad->Update();
   //   TPaletteAxis *palette =
   //      (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
   //   Int_t ci = palette->GetValueColor(30.);
   //
   // Then it is possible to retrieve the RGB components in the following way:
   //
   //   TColor *c = gROOT->GetColor(ci);
   //   float x,y,z;
   //   c->GetRGB(x,y,z);

   TColor::SetPalette(1, 0);
   TH2F      *fPadHist = new TH2F();
   //if (fPadHistMinimum) wmin = fPadHistMinimum;
   //if (fPadHistMaximum) wmax = fPadHistMaximum;
   Double_t wlmin = wmin;
   Double_t wlmax = wmax;
   //Int_t fContourLevel   = 40;
   Int_t fContourLevel   = gStyle->GetNumberContours();

   if (gPad->GetLogz()) {
      if (wmin <= 0 && wmax > 0) wmin = TMath::Min((Double_t)1,
                                                   (Double_t)0.001*wmax);
      wlmin = TMath::Log10(wmin);
      wlmax = TMath::Log10(wmax);
   }

   Int_t ncolors = gStyle->GetNumberOfColors();
   //Int_t ndivz   = TMath::Abs(fPadHist->GetContour());
   Int_t ndivz   = TMath::Abs(fContourLevel);
   Int_t theColor = 0;
   Double_t scale = wmax == wlmin ? 1 : ndivz/(wlmax - wlmin);

   if (fPadHist->TestBit(TH1::kUserContour) && gPad->GetLogz()) zc = TMath::Log10(zc);
   if (zc < wlmin) zc = wlmin;
   if (zc > wlmax) zc = wlmax;

   Int_t color = Int_t(-0.01+(zc-wlmin)*scale);
   delete fPadHist;

   if (ndivz)
      theColor = Int_t((color+0.99)*Double_t(ncolors)/Double_t(ndivz));
   return gStyle->GetColorPalette(theColor);
}

//______________________________________________________________________________
void GetGeometry(Double_t xyz_[][3])
{
   //Double_t xyz[NPD][3];
   const char*  filename = "geofile_Cyl_74x60_fix.txt";
   ifstream fin;
   fin.open (filename, ios::in);
   if (!fin) {
      std::cout << "Cannot open Geometry file : " << filename << ". in GetGeometry()" << std::endl;
      exit(1);
   }
   Int_t pmtid;
   while( ! fin.eof() ) {
     fin >> pmtid;
     fin >> xyz_[pmtid-1][0] >> xyz_[pmtid-1][1] >> xyz_[pmtid-1][2];
   }
   fin.close();
   fin.clear();
}

//______________________________________________________________________________
void DrawGeometry(Int_t type = 0, Double_t val[] = 0, Double_t min = 0, Double_t max = 1, const char *title = "HYPERK PMT map", Double_t xyz[][3] = 0)
{
/*
  Draw a PD map in a current campas.
  type 1  : Draw 2D map
  type 1x : Draw 3D map
*/

   Double_t pmtid = 1, theta, q, x, y, z, dir[3];
   //TString  filename[1][2]; // [sample][x,y,z]
   //   //filename[0][0] = "PMTData2.txt";
   //   //filename[0][0] = "PMTData.txt";
   //   filename[0][0] = "geofile_Cyl_74x60_fix.txt";

   //TArc * wallcirc = new TArc(0,0,WallR*0.5,0,360);
   //TArc * circrange[2];
   //wallcirc->SetFillStyle(0);
   TBox *wall_bg = new TBox(-1. * WallR * TMath::Pi() * 0.5 - 0.5 * 71., -1. * WallZ, WallR * TMath::Pi() * 0.5 + 0.5 * 71., WallZ);
   TArc *wall_top    = new TArc(0, GetYWall(0, 0, WallZ),  WallR*0.5);
   TArc *wall_bottom = new TArc(0, GetYWall(0, 0, -WallZ), WallR*0.5);
      wall_bg->SetFillStyle(1001);
      wall_bg->SetLineWidth(0);
      wall_bg->SetFillColor(kBlack);
      wall_top->SetLineWidth(0);
      wall_top->SetFillColor(kBlack);
      wall_bottom->SetLineWidth(0);
      wall_bottom->SetFillColor(kBlack);

   //TFrame * fFrame0 = new TFrame(-1.75 * WallR, 1.75 * WallR, -3.4 * WallZ, 3.4 * WallZ);
   if (gDirectory->Get("HYPERKPM")) { gDirectory->Get("HYPERKPM")->Delete(); }
   TH2D * fFrame = new TH2D(TString::Format("HYPERKPM"),  title, 4, -1.7 * WallR, 1.7 * WallR, 4, -4. * WallZ, 4. * WallZ);
   fFrame->SetFillStyle(0);
   fFrame->SetMinimum(min);
   fFrame->SetMaximum(max);
   //TH2D * fFrameCirc = new TH2D(TString::Format("HYPERKXY"),  "HYPERK PMT Map", 4, -1*WallR*0.5*1.1, WallR*0.5*1.1, 4, -1*WallR*0.5*1.1, WallR*0.5*1.1);
   fFrame->SetEntries(NPMT);
   fFrame->GetXaxis()->SetLabelSize(0.02);
   fFrame->GetYaxis()->SetLabelSize(0.02);

   //const Int_t npmt = 11146; // 11053
   const Int_t npmt = NPMT; // 

   TPolyMarker3D **polpmt[1][1];
   polpmt[0][0] = new TPolyMarker3D*[npmt];
   TMarker **pmt[1][1];
   pmt[0][0] = new TMarker*[npmt];
   for (Int_t k = 0; k < 1; k++) {
   for (Int_t i = 0; i < 1; i++) {
   for (Int_t ipmt = 0; ipmt < npmt; ipmt++) {
      polpmt[k][i][ipmt] = new TPolyMarker3D();
      polpmt[k][i][ipmt]->SetMarkerStyle(20);
      pmt[k][i][ipmt] = new TMarker();
      pmt[k][i][ipmt]->SetMarkerStyle(20);
   }
   }
   }
   //ifstream fin;
   //ofstream fout;
   //fout.open("PMTData3.txt");
   Int_t k = 0;
   for (Int_t j = 0; j < 1; j++) {
   for (Int_t i = 0; i < 1; i++) {
      k = 0;
      //fin.open (filename[j][i].Data(), ios::in);
      //cout << "OPEN " << filename[j][i].Data() << endl;
      while( (int)pmtid-1 < NPD ) {
        Int_t ipmt = (Int_t)(pmtid-1);
        x = xyz[ipmt][0] * 100.;
        y = xyz[ipmt][1] * 100.;
        z = xyz[ipmt][2] * 100.;
        q = val[ipmt];

        // 3D
        polpmt[0][i][(Int_t)(pmtid)-1]->SetNextPoint(x,y,z);
        polpmt[0][i][(Int_t)(pmtid)-1]->SetMarkerColor(GetValueColor(q,min,max));
        ////polpmt[0][i][(Int_t)(pmtid)-1]->SetMarkerColor(GetValueColor(q,0,100.0));
        ////polpmt[0][i][(Int_t)(pmtid)-1]->SetMarkerSize(q*0.003);
        polpmt[0][i][(Int_t)(pmtid)-1]->SetMarkerSize(0.15);

        // 2D
        pmt[0][i][(Int_t)(pmtid)-1]->SetX(GetXWall(x,y,z));
        pmt[0][i][(Int_t)(pmtid)-1]->SetY(GetYWall(x,y,z));
        pmt[0][i][(Int_t)(pmtid)-1]->SetMarkerColor(GetValueColor(q,min,max));
        //pmt[0][i][(Int_t)(pmtid)-1]->SetMarkerSize(q*0.05);
        pmt[0][i][(Int_t)(pmtid)-1]->SetMarkerSize(0.15);
        k++;

        //if(fin.eof() ) { break;}
        //cout << pmtid << " " << q <<" "<< x <<" "<< y << endl;
        //cout << pmtid << " " << pmt[0][i][(int)pmtid-1]->GetX() <<" "<< x <<" "<< y << endl;
        pmtid++;
      }
      //fin.close();
      //fin.clear();
   }
   }

   //fout.close();
   //cout << "DRAW" << endl;



   for (Int_t i = 0; i < 1 ; i++) {
      // 3D view
      if (type >= 10 && type < 100) {
         for (Int_t ipmt = 0; ipmt < npmt; ipmt++) {
           polpmt[0][i][ipmt]->Draw();
         }
      }

      // 2D view
      if (type >= 1 && type < 10) {
         fFrame->Draw("colz");
         gPad->Update();
         wall_bg->Draw();
         wall_top->Draw("");
         wall_bottom->Draw();
         ////TPaletteAxis *palette = (TPaletteAxis*)fFrame->GetListOfFunctions()->FindObject("palette");
         //TPaletteAxis *palette = 0;
         //if (fFrame->FindObject("palette")) {
         //   palette = ((TPaletteAxis*)(fFrame->FindObject("palette"))->Clone("PD_palette"));
         //   palette->Draw("same");
         //} else {
         //   cout << "Could not find pallete " << endl;
         //}
         //fFrame->Reset();
         for (Int_t ipmt = 0; ipmt < npmt; ipmt++) {
           pmt[0][i][ipmt]->Draw();
         }
      }
   }



}


//______________________________________________________________________________
TCanvas* pdmap_sample(Double_t *q = 0, const char* file = "HYPERKPMTMap.pdf")
{
   //Sample function to draw PD geometry using DrawGeometry()

   TRandom * gr = new TRandom();
   Double_t qq[NPMT];
   Double_t qqmin, qqmax;
   for (Int_t ipmt = 0; ipmt < NPMT; ipmt++) {
     if (q==0) { qq[ipmt] = gr->Gaus(); } 
     else {
        qq[ipmt] = q[ipmt];
     }
     if (ipmt == 0) { qqmin = qqmax = qq[ipmt]; }
     else if (qqmin>qq[ipmt]) { qqmin = qq[ipmt]; }
     else if (qqmax<qq[ipmt]) { qqmax = qq[ipmt]; }
   }

   TCanvas * c = new TCanvas("c", "c", 600, 600);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0);


   gStyle->SetPalette(1);
   //c->SetGridx(1);
   //c->SetGridy(1);
   if (kTRUE) { // Draw 2d or 3d
      c->SetTickx(1);
      c->SetTicky(1);
      DrawGeometry(1, qq, qqmin, qqmax, "Test");
      c->SaveAs(file);
   } else {
      DrawGeometry(10, qq, qqmin, qqmax, "Test");
      c->SaveAs(file);
   }
#ifndef __CINT__
   delete gr;
#endif
   return c;
}


//______________________________________________________________________________
Double_t *array_from_file(const char* file = "")
{
   //Input file should have one value in one line by PMT ID ordering. 

   Double_t *qq = new Double_t(NPMT);

   ifstream inputfile( file ) ; 
   int idx = 0;

   while (inputfile >> qq[idx]){
     idx++;
     if (idx > NPMT) { std::cout << "ERROR : Num. of lines in input file exceeds the limit " << NPMT << std::endl; }
   }
   if (idx != NPMT) { std::cout << "ERROR : Num. of lines in input file does not agree with the number of photodetectors = " << idx << std::endl; }
   inputfile.close();

   return qq;
}

//______________________________________________________________________________
void pdmap_from_file(const char* file = "", double min = 0, double max = 1, const char* outfile = "HYPERK_PDmap.png",  Bool_t Is3D = kFALSE)
{
   //Input file should have one value in one line by PMT ID ordering. 

   Double_t qq[NPMT];

   ifstream inputfile( file ) ; 
   int idx = 0;

   while (inputfile >> qq[idx]){
     idx++;
     if (idx > NPMT) { std::cout << "ERROR : Num. of lines in input file exceeds the limit " << NPMT << std::endl; }
   }
   if (idx != NPMT) { std::cout << "ERROR : Num. of lines in input file does not agree with the number of photodetectors = " << idx << std::endl; }
   inputfile.close();

   TCanvas * c = new TCanvas("c", "c", 600, 600);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0);


   const UInt_t Number = 3; 
   Int_t nb = 3; 
   Int_t MyPalette[Number];
   MyPalette[0] = kGreen;
   MyPalette[1] = kBlue;
   MyPalette[2] = kRed;
   gStyle->SetPalette(nb, MyPalette);
   //gStyle->SetPalette(1);
   //c->SetGridx(1);
   //c->SetGridy(1);
   if (!Is3D) { // Draw 2d or 3d
      c->SetTickx(1);
      c->SetTicky(1);
      DrawGeometry(1, qq, min, max, file);
      c->SaveAs(outfile);
   } else {
      DrawGeometry(10, qq, min, max, file);
      c->SaveAs(outfile);
   }
#ifndef __CINT__
   delete c;
#endif
}


//______________________________________________________________________________
void usage()
{
   std::cout << "Usage : " << std::endl 
             <<  "       pdmap_from_hv(const char* file = \"HVTable.dat\", const char* outfile = \"HYPERK_HVmap.png\",  Bool_t Is3D = kFALSE)" << std::endl 
             <<  "       pdmap_from_file(const char* file = \"\", double min = 0, double max = 1, const char* outfile = \"HYPERK_PDmap.png\",  Bool_t Is3D = kFALSE)" << std::endl 
             <<  "       pdmap_sample(Double_t *q = 0)" << std::endl 
             <<  "       DrawGeometry(Int_t type = 0, Double_t val[] = 0, Double_t min = 0, Double_t max = 1, const char *title = \"HYPERK PMT map\")" << std::endl <<
                              std::endl ;
}


//______________________________________________________________________________
void test()
{
   gROOT->SetStyle("Plain");
   //pdmap_from_hv();
   pdmap_sample();
   Double_t q[NPMT];
   for (int i=0; i<NPMT; i++) {q[i] = (Double_t)i;} 
   pdmap_sample(q, "test2.pdf");
   
/*
   Int_t MyPalette[100];
   Double_t r[]    = {0., 0.0, 1.0, 1.0, 1.0};
   Double_t g[]    = {0., 0.0, 0.0, 1.0, 1.0};
   Double_t b[]    = {0., 1.0, 0.0, 0.0, 1.0};
   Double_t stop[] = {0., .25, .50, .75, 1.0};
   Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 100);
   for (int i=0;i<100;i++) MyPalette[i] = FI+i;
   gStyle->SetPalette(100, MyPalette);
*/
   const UInt_t Number = 3; 
   Int_t MyPalette[Number];
   Double_t Red[Number] = { 1.00, 0.00, 0.00}; 
   Double_t Green[Number] = { 0.00, 1.00, 0.00}; 
   Double_t Blue[Number] = { 1.00, 0.00, 1.00}; 
   Double_t Stops[Number] = { 0.00, 0.50, 1.00 }; 
   Int_t nb = 3; 
   Int_t FI = TColor::CreateGradientColorTable(Number+1,Stops,Red,Blue,Green,nb); 
   for (int i=0;i<nb;i++) MyPalette[i] = FI+i;
//   MyPalette[0] = 0;
//   pdmap_from_file("PMTtype.dat", 0, 3, "PMTType.root");
}


//______________________________________________________________________________
void pdmap()
{
}

#endif
