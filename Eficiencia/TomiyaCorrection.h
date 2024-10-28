#include "TMath.h"
#include <vector>

using namespace std;
vector<string> split(string str, string s_cut){
	vector<string> v_str;
	str = str + s_cut;
	long unsigned int l = str.length();
	long sl = s_cut.length();
	string::size_type pos = 0, prev = 0;
	for (;pos < l && (pos = str.find(s_cut, pos)) != string::npos; prev = (pos += sl)) {
		string item = str.substr(prev, pos - prev);
		v_str.push_back(item);
	}
	return v_str;
}

void rotateBfield(double &Bx, double &By, char Dir){
  if( Dir<'A' || Dir>'H'){
    std::cout<<"PMT Direction should be in the ragne from A to H"<<std::endl;
    return;
  }

  static bool isFirstTime = true;
  if(isFirstTime){
    std::cout<<"##############################"<<std::endl
	     <<"Warning in rotateBfield."<<std::endl
      	     <<"This function is tentative. You should edit the relation of A-H and the dynode direction correctly."<<std::endl
	     <<"##############################"<<std::endl;
    isFirstTime = false;
  }
    
  int irot = Dir - 'A';
  double const deg = 1./TMath::RadToDeg();
  double const offset = 0 * deg; // tentative
  double rot = 45*deg * irot - offset; //tentative
  double newx = Bx*cos(rot) - By*sin(rot); //tentative
  double newy = Bx*sin(rot) + By*cos(rot); //tentative

  Bx = newx;
  By = newy;
}

double tomiyaCorrection(double ebb/*HV in Volt*/, double Bx, double By, double Bz, bool HVdependence = false){ // HV dependence : if true, use interpolation or extrapolation. Otherwise not.
  int const nDir=3;
  int const nPar=7;
  int const nHV = 5;
  double const HV[nHV]={1540,1640,1740,1840,1940};

  /*  //// based on tables in X,Y,Z tags of https://hep.phys.keio.ac.jp/~nishimura/share/PMTRateCorrectionFactor.xlsx
  double const par[nDir][nHV][nPar]={
    { // Bx
      { 0.999812, -0.0000841147, -0.00000271294, 0.00000000104102, 0.000000000000352248, -2.87749E-15, 0 }, // 1540V
      {1, -0.000131902, -0.0000018079, 0.00000000155544,-0.00000000000221921, -3.80754E-15, 0}, // 1640V
      {1, -0.00010322, -0.00000166296, 0.00000000138139, -0.00000000000214074, -3.40495E-15, 0}, // 1740V
      {1, -0.0000169, -0.00000144742, 0.000000000247498, -0.00000000000253004, -1.91244E-16, 0}, // 1840V
      {0.999993, -0.0000997472, -0.00000107409, 0.00000000135765, -0.00000000000310736, -2.78903E-15, 0} // 1940V
    },{ // By
      {0.999998, 0.000156965, -0.00000142785, -0.00000000153435, -0.00000000000240921, 2.98969E-15, 5.67711E-18}, // 1540V
      {0.997922, 0.0000694082, -0.00000076869, -0.000000000876493, -0.00000000000533161, 1.42977E-15, 1.04802E-17}, // 1640V
      {0.989872, 0.00000440899, -0.000000423596, -0.000000000714794, -0.00000000000592093, 1.36499E-15, 1.02494E-17}, // 1740V
      {1, 0.000062676, -0.000000188694, -0.0000000012973, -0.00000000000686332, 2.51542E-15, 1.20204E-17}, // 1840V
      {0.997301, 0.0000360363, -0.0000000642535, -0.000000000675914, -0.00000000000721265, 6.76981E-16, 1.30657E-17} // 1940V
    },{ // Bz
      {0.984078, 0.0000374851, 0.000000105514, -0.00000000103782, -0.00000000000632785, 2.02829E-15, 8.24998E-18}, // 1540V
      {0.99234, 0.0000444225, 0.00000015511, -0.0000000010849, -0.00000000000533636, 2.23955E-15, 6.66422E-18}, // 1640V
      {0.999993, 0.0000349757, -0.000000326224, -0.000000000287445, -0.000000000000247255, -6.01343E-16, -3.19635E-18}, // 1740V
      {0.999238, 0.00000326132, -0.0000000948637, -0.000000000257745, -0.000000000000891122, -2.64933E-16, -2.04591E-18}, // 1840V
      {1, 0.0000170068, 0.00000017373, -0.000000000564378, -0.0000000000026247, 9.33524E-16, 1.66666E-18} // 1940V
    }
  };
  ///// */

  //    /*  //// based on cells I~M,3~6 in "Correction" tag of https://hep.phys.keio.ac.jp/~nishimura/share/PMTRateCorrectionFactor.xlsx. Almost similar to functions in X,Y,Z tags, but par[0] is fixed to 0.
    double const par[nDir][nHV][nPar]={
    { // Bx
      { 1, -0.0000841147, -0.00000271294, 0.00000000104102, 0.000000000000352248, -2.87749E-15, 0 }, // 1540V
      {1, -0.000131902, -0.0000018079, 0.00000000155544,-0.00000000000221921, -3.80754E-15, 0}, // 1640V
      {1, -0.00010322, -0.00000166296, 0.00000000138139, -0.00000000000214074, -3.40495E-15, 0}, // 1740V
      {1, -0.0000169, -0.00000144742, 0.000000000247498, -0.00000000000253004, -1.91244E-16, 0}, // 1840V
      {1, -0.0000997472, -0.00000107409, 0.00000000135765, -0.00000000000310736, -2.78903E-15, 0} // 1940V
    },{ // By
      {1, 0.000156965, -0.00000142785, -0.00000000153435, -0.00000000000240921, 2.98969E-15, 5.67711E-18}, // 1540V
      {1, 0.0000694082, -0.00000076869, -0.000000000876493, -0.00000000000533161, 1.42977E-15, 1.04802E-17}, // 1640V
      {1, 0.00000440899, -0.000000423596, -0.000000000714794, -0.00000000000592093, 1.36499E-15, 1.02494E-17}, // 1740V
      {1, 0.000062676, -0.000000188694, -0.0000000012973, -0.00000000000686332, 2.51542E-15, 1.20204E-17}, // 1840V
      {1, 0.0000360363, -0.0000000642535, -0.000000000675914, -0.00000000000721265, 6.76981E-16, 1.30657E-17} // 1940V
    },{ // Bz
      {1, 0.0000374851, 0.000000105514, -0.00000000103782, -0.00000000000632785, 2.02829E-15, 8.24998E-18}, // 1540V
      {1, 0.0000444225, 0.00000015511, -0.0000000010849, -0.00000000000533636, 2.23955E-15, 6.66422E-18}, // 1640V
      {1, 0.0000349757, -0.000000326224, -0.000000000287445, -0.000000000000247255, -6.01343E-16, -3.19635E-18}, // 1740V
      {1, 0.00000326132, -0.0000000948637, -0.000000000257745, -0.000000000000891122, -2.64933E-16, -2.04591E-18}, // 1840V
      {1, 0.0000170068, 0.00000017373, -0.000000000564378, -0.0000000000026247, 9.33524E-16, 1.66666E-18} // 1940V
    }
  };
  ///// */

  double B[3]={Bx,By,Bz};
  double const fitRangeMin[3] = {-500.,-450.,-450.};
  double const fitRangeMax[3] = {500.,550.,550.};
  for(int ix=0;ix<3;ix++){
    if( B[ix]<fitRangeMin[ix] || B[ix]>fitRangeMax[ix]){
      static int nWarning = 0;
      if(nWarning<50){
	std::cout<<"warning in TomiyaCorrection."<<std::endl
		 <<"B-field is too large : B["<<ix<<"] = "<<B[ix]<<std::endl
		 <<"It is out of fitting range("<<fitRangeMin[ix]<<"--"<<fitRangeMax[ix]<<")"<<std::endl;
	
	nWarning++;
      }
    }
  }
  

  double cor[nDir]={0.};
  if(HVdependence){ // use interpolation/extrapolation for HV value
    
    int iHVlow = 0;
    if(ebb<HV[0]) iHVlow = 0;
    else if(ebb>=HV[nHV-1]) iHVlow = nHV-2;
    else{
      for(int i=1;i<nHV;i++){
	if(ebb<HV[i]){
	  iHVlow = i-1;
	  break;
	}
      }
    }

    double w = (ebb-HV[iHVlow])/(HV[iHVlow+1]-HV[iHVlow]);
    for(int dir = 0; dir<3; dir++){
      double f[2]={0};
      for(int i=0;i<2;i++){
	int iHV = iHVlow+i;
	for(int ipar = 0;ipar<nPar;ipar++){
	  f[i]+=par[dir][iHV][ipar]*pow(B[dir],ipar);
	}
      }
      cor[dir] = (1-w)*f[0] + w*f[1];
      if(cor[dir]<0) cor[dir] = 0;
    }
    
  }else{ 
    
    double const boarder[nHV]={1600,1700,1800,1900,2650};
    int iHV = -1;
    for(int i=0;i<nHV;i++){
      if(ebb<boarder[i]){
	iHV = i;
	break;
      }
    }
    if(iHV==-1){
      static int nWarning = 0;
      if(nWarning<50){
	std::cout<<"warning in TomiyaCorrection."<<std::endl
		 <<"ebb(="<<ebb<<") is too high."<<std::endl;
	return 1.;
	nWarning++;
      }
    }
    for(int dir=0;dir<3;dir++){
      cor[dir]=0;
      for(int ipar = 0;ipar<nPar;ipar++){
	cor[dir] += par[dir][iHV][ipar]*pow(B[dir],ipar);
      }
    }
  }
  return cor[0]*cor[1]*cor[2];
}
