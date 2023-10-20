#include "TMath.h"
#include <assert.h>
#include <iostream>
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

int getnJetnBJet(int nJ, int nBJ) {
  int bin = -1;

  if      (nJ == 0)             bin = 0;
  else if (nJ == 1 && nBJ == 0) bin = 1;
  else if (nJ == 1 && nBJ == 1) bin = 2;
  else if (nJ == 2 && nBJ == 0) bin = 3;
  else if (nJ == 2 && nBJ == 1) bin = 4;
  else if (nJ == 2 && nBJ == 2) bin = 5;
  else if (nJ == 3 && nBJ == 0) bin = 6;
  else if (nJ == 3 && nBJ == 1) bin = 7;
  else if (nJ == 3 && nBJ == 2) bin = 8;
  else if (nJ == 3 && nBJ == 3) bin = 9;
  else if (nJ >= 4 && nBJ >= 0) bin = 10;
  return bin;
}


int getnLooseJets(int nlj, int nj) {
  int bin = -1;
  bin = nlj - nj;

  return bin;
}


Float_t getBDtW(Double_t BDT){
  if      (BDT < -0.346575   ) return 1;
  else if (BDT < -0.264520) return 2;
  else if (BDT < -0.167959) return 3;
  else if (BDT < -0.073513) return 4;
  else if (BDT < 0.024968) return 5;
  else if (BDT < 0.119442) return 6;
  else if (BDT < 0.201342) return 7;
  else if (BDT < 0.258906) return 8;
  else if (BDT < 0.317417) return 9;
  else                       return 10;
}

Float_t getBDtWOther(Double_t BDT){
  if      (BDT < -0.238182   ) return 1;
  else if (BDT < -0.093824) return 2;
  else if (BDT < 0.002727) return 3;
  else if (BDT < 0.073376) return 4;
  else if (BDT < 0.144134) return 5;
  else                       return 6;
}

// Random Forests
Float_t getRFtW(Double_t BDT){
  if      (BDT < 0.362929   ) return 1;
  else if (BDT < 0.393709) return 2;
  else if (BDT < 0.411112) return 3;
  else if (BDT < 0.424242) return 4;
  else if (BDT < 0.452956) return 5;
  else if (BDT < 0.485076) return 6;
  else if (BDT < 0.518692) return 7;
  else if (BDT < 0.534609) return 8;
  else if (BDT < 0.546438) return 9;
  else                       return 10;
}

Float_t getRFtWOther(Double_t BDT){
  if      (BDT < 0.335388   ) return 1;
  else if (BDT < 0.361837) return 2;
  else if (BDT < 0.387806) return 3;
  else if (BDT < 0.440645) return 4;
  else if (BDT < 0.496462) return 5;
  else                       return 6;
}

Float_t getRFtW1j0b(Double_t BDT){
  if      (BDT < 0.395016   ) return 1;
  else if (BDT < 0.415115) return 2;
  else if (BDT < 0.430495) return 3;
  else if (BDT < 0.445367) return 4;
  else if (BDT < 0.460616) return 5;
  else if (BDT < 0.482564) return 6;
  else if (BDT < 0.507189) return 7;
  else if (BDT < 0.529880) return 8;
  else if (BDT < 0.547096) return 9;
  else                       return 10;
} 

Float_t getRFtW_mm(Double_t BDT){
  if      (BDT < 0.273437   ) return 1;
  else if (BDT < 0.326876) return 2;
  else if (BDT < 0.344809) return 3;
  else if (BDT < 0.360357) return 4;
  else if (BDT < 0.376982) return 5;
  else if (BDT < 0.411057) return 6;
  else if (BDT < 0.429613) return 7;
  else if (BDT < 0.446195) return 8;
  else if (BDT < 0.459559) return 9;
  else                       return 10;
}

Float_t getRFtW_ee(Double_t BDT){
  if      (BDT < 0.272420   ) return 1;
  else if (BDT < 0.328963) return 2;
  else if (BDT < 0.347061) return 3;
  else if (BDT < 0.359887) return 4;
  else if (BDT < 0.380686) return 5;
  else if (BDT < 0.414465) return 6;
  else if (BDT < 0.432662) return 7;
  else if (BDT < 0.445169) return 8;
  else if (BDT < 0.456785) return 9;
  else                       return 10;
}

// Function to obtain the electron in the emu channel from the LepGood collection
Float_t getElePt(int l1_pdg, int l2_pdg, Float_t l1_pt, Float_t l2_pt){
  if (abs(l1_pdg) == 11) return l1_pt;
  else if (abs(l2_pdg) == 11) return l2_pt;
  else return -999;
}

// For muon
Float_t getMuPt(int l1_pdg, int l2_pdg, Float_t l1_pt, Float_t l2_pt){
  if (abs(l1_pdg) == 13) return l1_pt;
  else if (abs(l2_pdg) == 13) return l2_pt;
  else return -999;
}

///////////


    
Float_t getRF1j1b_opt_10bins(Double_t BDT){
  if      (BDT < 0.392348   ) return 1;
  else if (BDT < 0.414540) return 2;
  else if (BDT < 0.431385) return 3;
  else if (BDT < 0.448032) return 4;
  else if (BDT < 0.464931) return 5;
  else if (BDT < 0.487771) return 6;
  else if (BDT < 0.510953) return 7;
  else if (BDT < 0.532430) return 8;
  else if (BDT < 0.548111) return 9;
  else                       return 10;
} 

    
Float_t getRF1j1b_opt_20bins(Double_t BDT){
  if      (BDT < 0.349516   ) return 1;
  else if (BDT < 0.362929) return 2;
  else if (BDT < 0.378330) return 3;
  else if (BDT < 0.393709) return 4;
  else if (BDT < 0.404067) return 5;
  else if (BDT < 0.411112) return 6;
  else if (BDT < 0.417000) return 7;
  else if (BDT < 0.424242) return 8;
  else if (BDT < 0.440702) return 9;
  else if (BDT < 0.452956) return 10;
  else if (BDT < 0.467729) return 11;
  else if (BDT < 0.485076) return 12;
  else if (BDT < 0.503881) return 13;
  else if (BDT < 0.518692) return 14;
  else if (BDT < 0.527935) return 15;
  else if (BDT < 0.534609) return 16;
  else if (BDT < 0.540526) return 17;
  else if (BDT < 0.546438) return 18;
  else if (BDT < 0.553782) return 19;
  else                       return 20;
} 
Float_t getRF2j1b_opt_12bins(Double_t BDT){
  if      (BDT < 0.319423   ) return 1;
  else if (BDT < 0.335388) return 2;
  else if (BDT < 0.347965) return 3;
  else if (BDT < 0.361837) return 4;
  else if (BDT < 0.373409) return 5;
  else if (BDT < 0.387806) return 6;
  else if (BDT < 0.409925) return 7;
  else if (BDT < 0.440645) return 8;
  else if (BDT < 0.461806) return 9;
  else if (BDT < 0.496462) return 10;
  else if (BDT < 0.516066) return 11;
  else                       return 12;
}
    
Float_t getRF1j1b_opt_30bins(Double_t BDT){
  if      (BDT < 0.359345   ) return 1;
  else if (BDT < 0.379355) return 2;
  else if (BDT < 0.392348) return 3;
  else if (BDT < 0.401583) return 4;
  else if (BDT < 0.408705) return 5;
  else if (BDT < 0.414540) return 6;
  else if (BDT < 0.420179) return 7;
  else if (BDT < 0.425511) return 8;
  else if (BDT < 0.431385) return 9;
  else if (BDT < 0.437229) return 10;
  else if (BDT < 0.442798) return 11;
  else if (BDT < 0.448032) return 12;
  else if (BDT < 0.453191) return 13;
  else if (BDT < 0.459034) return 14;
  else if (BDT < 0.464931) return 15;
  else if (BDT < 0.471469) return 16;
  else if (BDT < 0.479569) return 17;
  else if (BDT < 0.487771) return 18;
  else if (BDT < 0.495838) return 19;
  else if (BDT < 0.503584) return 20;
  else if (BDT < 0.510953) return 21;
  else if (BDT < 0.518260) return 22;
  else if (BDT < 0.525770) return 23;
  else if (BDT < 0.532430) return 24;
  else if (BDT < 0.538578) return 25;
  else if (BDT < 0.544030) return 26;
  else if (BDT < 0.548111) return 27;
  else if (BDT < 0.552943) return 28;
  else if (BDT < 0.558362) return 29;
  else                       return 30;
} 

    
Float_t getRF1j1b_opt_40bins(Double_t BDT){
  if      (BDT < 0.352456   ) return 1;
  else if (BDT < 0.370521) return 2;
  else if (BDT < 0.383012) return 3;
  else if (BDT < 0.392348) return 4;
  else if (BDT < 0.399568) return 5;
  else if (BDT < 0.405315) return 6;
  else if (BDT < 0.410230) return 7;
  else if (BDT < 0.414540) return 8;
  else if (BDT < 0.418893) return 9;
  else if (BDT < 0.422710) return 10;
  else if (BDT < 0.426950) return 11;
  else if (BDT < 0.431385) return 12;
  else if (BDT < 0.435789) return 13;
  else if (BDT < 0.440062) return 14;
  else if (BDT < 0.444143) return 15;
  else if (BDT < 0.448032) return 16;
  else if (BDT < 0.451842) return 17;
  else if (BDT < 0.456059) return 18;
  else if (BDT < 0.460491) return 19;
  else if (BDT < 0.464931) return 20;
  else if (BDT < 0.469612) return 21;
  else if (BDT < 0.475433) return 22;
  else if (BDT < 0.481651) return 23;
  else if (BDT < 0.487771) return 24;
  else if (BDT < 0.493889) return 25;
  else if (BDT < 0.499823) return 26;
  else if (BDT < 0.505413) return 27;
  else if (BDT < 0.510953) return 28;
  else if (BDT < 0.516391) return 29;
  else if (BDT < 0.522045) return 30;
  else if (BDT < 0.527552) return 31;
  else if (BDT < 0.532430) return 32;
  else if (BDT < 0.537136) return 33;
  else if (BDT < 0.541378) return 34;
  else if (BDT < 0.545154) return 35;
  else if (BDT < 0.548111) return 36;
  else if (BDT < 0.551517) return 37;
  else if (BDT < 0.555574) return 38;
  else if (BDT < 0.559614) return 39;
  else                       return 40;
} 

    
Float_t getRF1j1b_opt_50bins(Double_t BDT){
  if      (BDT < 0.347603   ) return 1;
  else if (BDT < 0.364172) return 2;
  else if (BDT < 0.376086) return 3;
  else if (BDT < 0.385071) return 4;
  else if (BDT < 0.392348) return 5;
  else if (BDT < 0.398295) return 6;
  else if (BDT < 0.403167) return 7;
  else if (BDT < 0.407411) return 8;
  else if (BDT < 0.411105) return 9;
  else if (BDT < 0.414540) return 10;
  else if (BDT < 0.418006) return 11;
  else if (BDT < 0.421140) return 12;
  else if (BDT < 0.424364) return 13;
  else if (BDT < 0.427832) return 14;
  else if (BDT < 0.431385) return 15;
  else if (BDT < 0.434929) return 16;
  else if (BDT < 0.438377) return 17;
  else if (BDT < 0.441715) return 18;
  else if (BDT < 0.444944) return 19;
  else if (BDT < 0.448032) return 20;
  else if (BDT < 0.451061) return 21;
  else if (BDT < 0.454310) return 22;
  else if (BDT < 0.457817) return 23;
  else if (BDT < 0.461379) return 24;
  else if (BDT < 0.464931) return 25;
  else if (BDT < 0.468562) return 26;
  else if (BDT < 0.473033) return 27;
  else if (BDT < 0.477900) return 28;
  else if (BDT < 0.482905) return 29;
  else if (BDT < 0.487771) return 30;
  else if (BDT < 0.492700) return 31;
  else if (BDT < 0.497427) return 32;
  else if (BDT < 0.502123) return 33;
  else if (BDT < 0.506516) return 34;
  else if (BDT < 0.510953) return 35;
  else if (BDT < 0.515308) return 36;
  else if (BDT < 0.519754) return 37;
  else if (BDT < 0.524295) return 38;
  else if (BDT < 0.528573) return 39;
  else if (BDT < 0.532430) return 40;
  else if (BDT < 0.536230) return 41;
  else if (BDT < 0.539724) return 42;
  else if (BDT < 0.543014) return 43;
  else if (BDT < 0.545777) return 44;
  else if (BDT < 0.548111) return 45;
  else if (BDT < 0.550708) return 46;
  else if (BDT < 0.554022) return 47;
  else if (BDT < 0.557306) return 48;
  else if (BDT < 0.560457) return 49;
  else                       return 50;
} 

    
Float_t getRF1j1b_opt_60bins(Double_t BDT){
  if      (BDT < 0.344595   ) return 1;
  else if (BDT < 0.359345) return 2;
  else if (BDT < 0.370521) return 3;
  else if (BDT < 0.379355) return 4;
  else if (BDT < 0.386373) return 5;
  else if (BDT < 0.392348) return 6;
  else if (BDT < 0.397416) return 7;
  else if (BDT < 0.401583) return 8;
  else if (BDT < 0.405315) return 9;
  else if (BDT < 0.408705) return 10;
  else if (BDT < 0.411677) return 11;
  else if (BDT < 0.414540) return 12;
  else if (BDT < 0.417417) return 13;
  else if (BDT < 0.420179) return 14;
  else if (BDT < 0.422710) return 15;
  else if (BDT < 0.425511) return 16;
  else if (BDT < 0.428423) return 17;
  else if (BDT < 0.431385) return 18;
  else if (BDT < 0.434344) return 19;
  else if (BDT < 0.437229) return 20;
  else if (BDT < 0.440062) return 21;
  else if (BDT < 0.442798) return 22;
  else if (BDT < 0.445469) return 23;
  else if (BDT < 0.448032) return 24;
  else if (BDT < 0.450561) return 25;
  else if (BDT < 0.453191) return 26;
  else if (BDT < 0.456059) return 27;
  else if (BDT < 0.459034) return 28;
  else if (BDT < 0.461976) return 29;
  else if (BDT < 0.464931) return 30;
  else if (BDT < 0.467884) return 31;
  else if (BDT < 0.471469) return 32;
  else if (BDT < 0.475433) return 33;
  else if (BDT < 0.479569) return 34;
  else if (BDT < 0.483721) return 35;
  else if (BDT < 0.487771) return 36;
  else if (BDT < 0.491891) return 37;
  else if (BDT < 0.495838) return 38;
  else if (BDT < 0.499823) return 39;
  else if (BDT < 0.503584) return 40;
  else if (BDT < 0.507267) return 41;
  else if (BDT < 0.510953) return 42;
  else if (BDT < 0.514589) return 43;
  else if (BDT < 0.518260) return 44;
  else if (BDT < 0.522045) return 45;
  else if (BDT < 0.525770) return 46;
  else if (BDT < 0.529233) return 47;
  else if (BDT < 0.532430) return 48;
  else if (BDT < 0.535623) return 49;
  else if (BDT < 0.538578) return 50;
  else if (BDT < 0.541378) return 51;
  else if (BDT < 0.544030) return 52;
  else if (BDT < 0.546183) return 53;
  else if (BDT < 0.548111) return 54;
  else if (BDT < 0.550199) return 55;
  else if (BDT < 0.552943) return 56;
  else if (BDT < 0.555574) return 57;
  else if (BDT < 0.558362) return 58;
  else if (BDT < 0.561096) return 59;
  else                       return 60;
} 

    
Float_t getRF1j1b_opt_70bins(Double_t BDT){
  if      (BDT < 0.342410   ) return 1;
  else if (BDT < 0.355623) return 2;
  else if (BDT < 0.366159) return 3;
  else if (BDT < 0.374556) return 4;
  else if (BDT < 0.381509) return 5;
  else if (BDT < 0.387278) return 6;
  else if (BDT < 0.392348) return 7;
  else if (BDT < 0.396768) return 8;
  else if (BDT < 0.400445) return 9;
  else if (BDT < 0.403798) return 10;
  else if (BDT < 0.406831) return 11;
  else if (BDT < 0.409587) return 12;
  else if (BDT < 0.412087) return 13;
  else if (BDT < 0.414540) return 14;
  else if (BDT < 0.416996) return 15;
  else if (BDT < 0.419458) return 16;
  else if (BDT < 0.421572) return 17;
  else if (BDT < 0.423883) return 18;
  else if (BDT < 0.426338) return 19;
  else if (BDT < 0.428846) return 20;
  else if (BDT < 0.431385) return 21;
  else if (BDT < 0.433926) return 22;
  else if (BDT < 0.436405) return 23;
  else if (BDT < 0.438865) return 24;
  else if (BDT < 0.441244) return 25;
  else if (BDT < 0.443568) return 26;
  else if (BDT < 0.445844) return 27;
  else if (BDT < 0.448032) return 28;
  else if (BDT < 0.450206) return 29;
  else if (BDT < 0.452415) return 30;
  else if (BDT < 0.454807) return 31;
  else if (BDT < 0.457308) return 32;
  else if (BDT < 0.459862) return 33;
  else if (BDT < 0.462397) return 34;
  else if (BDT < 0.464931) return 35;
  else if (BDT < 0.467433) return 36;
  else if (BDT < 0.470387) return 37;
  else if (BDT < 0.473717) return 38;
  else if (BDT < 0.477192) return 39;
  else if (BDT < 0.480747) return 40;
  else if (BDT < 0.484300) return 41;
  else if (BDT < 0.487771) return 42;
  else if (BDT < 0.491317) return 43;
  else if (BDT < 0.494728) return 44;
  else if (BDT < 0.498108) return 45;
  else if (BDT < 0.501476) return 46;
  else if (BDT < 0.504627) return 47;
  else if (BDT < 0.507815) return 48;
  else if (BDT < 0.510953) return 49;
  else if (BDT < 0.514070) return 50;
  else if (BDT < 0.517191) return 51;
  else if (BDT < 0.520404) return 52;
  else if (BDT < 0.523657) return 53;
  else if (BDT < 0.526804) return 54;
  else if (BDT < 0.529687) return 55;
  else if (BDT < 0.532430) return 56;
  else if (BDT < 0.535190) return 57;
  else if (BDT < 0.537757) return 58;
  else if (BDT < 0.540204) return 59;
  else if (BDT < 0.542552) return 60;
  else if (BDT < 0.544681) return 61;
  else if (BDT < 0.546466) return 62;
  else if (BDT < 0.548111) return 63;
  else if (BDT < 0.549861) return 64;
  else if (BDT < 0.552120) return 65;
  else if (BDT < 0.554455) return 66;
  else if (BDT < 0.556816) return 67;
  else if (BDT < 0.559059) return 68;
  else if (BDT < 0.561773) return 69;
  else                       return 70;
} 

    
Float_t getRF1j1b_opt_80bins(Double_t BDT){
  if      (BDT < 0.340732   ) return 1;
  else if (BDT < 0.352456) return 2;
  else if (BDT < 0.362347) return 3;
  else if (BDT < 0.370521) return 4;
  else if (BDT < 0.377345) return 5;
  else if (BDT < 0.383012) return 6;
  else if (BDT < 0.387977) return 7;
  else if (BDT < 0.392348) return 8;
  else if (BDT < 0.396264) return 9;
  else if (BDT < 0.399568) return 10;
  else if (BDT < 0.402581) return 11;
  else if (BDT < 0.405315) return 12;
  else if (BDT < 0.407910) return 13;
  else if (BDT < 0.410230) return 14;
  else if (BDT < 0.412402) return 15;
  else if (BDT < 0.414540) return 16;
  else if (BDT < 0.416682) return 17;
  else if (BDT < 0.418893) return 18;
  else if (BDT < 0.420775) return 19;
  else if (BDT < 0.422710) return 20;
  else if (BDT < 0.424788) return 21;
  else if (BDT < 0.426950) return 22;
  else if (BDT < 0.429165) return 23;
  else if (BDT < 0.431385) return 24;
  else if (BDT < 0.433612) return 25;
  else if (BDT < 0.435789) return 26;
  else if (BDT < 0.437949) return 27;
  else if (BDT < 0.440062) return 28;
  else if (BDT < 0.442124) return 29;
  else if (BDT < 0.444143) return 30;
  else if (BDT < 0.446122) return 31;
  else if (BDT < 0.448032) return 32;
  else if (BDT < 0.449939) return 33;
  else if (BDT < 0.451842) return 34;
  else if (BDT < 0.453881) return 35;
  else if (BDT < 0.456059) return 36;
  else if (BDT < 0.458270) return 37;
  else if (BDT < 0.460491) return 38;
  else if (BDT < 0.462712) return 39;
  else if (BDT < 0.464931) return 40;
  else if (BDT < 0.467095) return 41;
  else if (BDT < 0.469612) return 42;
  else if (BDT < 0.472447) return 43;
  else if (BDT < 0.475433) return 44;
  else if (BDT < 0.478532) return 45;
  else if (BDT < 0.481651) return 46;
  else if (BDT < 0.484733) return 47;
  else if (BDT < 0.487771) return 48;
  else if (BDT < 0.490883) return 49;
  else if (BDT < 0.493889) return 50;
  else if (BDT < 0.496831) return 51;
  else if (BDT < 0.499823) return 52;
  else if (BDT < 0.502676) return 53;
  else if (BDT < 0.505413) return 54;
  else if (BDT < 0.508220) return 55;
  else if (BDT < 0.510953) return 56;
  else if (BDT < 0.513679) return 57;
  else if (BDT < 0.516391) return 58;
  else if (BDT < 0.519193) return 59;
  else if (BDT < 0.522045) return 60;
  else if (BDT < 0.524844) return 61;
  else if (BDT < 0.527552) return 62;
  else if (BDT < 0.530028) return 63;
  else if (BDT < 0.532430) return 64;
  else if (BDT < 0.534856) return 65;
  else if (BDT < 0.537136) return 66;
  else if (BDT < 0.539294) return 67;
  else if (BDT < 0.541378) return 68;
  else if (BDT < 0.543396) return 69;
  else if (BDT < 0.545154) return 70;
  else if (BDT < 0.546678) return 71;
  else if (BDT < 0.548111) return 72;
  else if (BDT < 0.549631) return 73;
  else if (BDT < 0.551517) return 74;
  else if (BDT < 0.553619) return 75;
  else if (BDT < 0.555574) return 76;
  else if (BDT < 0.557721) return 77;
  else if (BDT < 0.559614) return 78;
  else if (BDT < 0.562413) return 79;
  else                       return 80;
} 

    
Float_t getRF1j1b_opt_90bins(Double_t BDT){
  if      (BDT < 0.339265   ) return 1;
  else if (BDT < 0.349748) return 2;
  else if (BDT < 0.359345) return 3;
  else if (BDT < 0.367198) return 4;
  else if (BDT < 0.373656) return 5;
  else if (BDT < 0.379355) return 6;
  else if (BDT < 0.384154) return 7;
  else if (BDT < 0.388485) return 8;
  else if (BDT < 0.392348) return 9;
  else if (BDT < 0.395858) return 10;
  else if (BDT < 0.398868) return 11;
  else if (BDT < 0.401583) return 12;
  else if (BDT < 0.404141) return 13;
  else if (BDT < 0.406498) return 14;
  else if (BDT < 0.408705) return 15;
  else if (BDT < 0.410717) return 16;
  else if (BDT < 0.412647) return 17;
  else if (BDT < 0.414540) return 18;
  else if (BDT < 0.416437) return 19;
  else if (BDT < 0.418401) return 20;
  else if (BDT < 0.420179) return 21;
  else if (BDT < 0.421813) return 22;
  else if (BDT < 0.423623) return 23;
  else if (BDT < 0.425511) return 24;
  else if (BDT < 0.427439) return 25;
  else if (BDT < 0.429413) return 26;
  else if (BDT < 0.431385) return 27;
  else if (BDT < 0.433367) return 28;
  else if (BDT < 0.435312) return 29;
  else if (BDT < 0.437229) return 30;
  else if (BDT < 0.439133) return 31;
  else if (BDT < 0.440983) return 32;
  else if (BDT < 0.442798) return 33;
  else if (BDT < 0.444588) return 34;
  else if (BDT < 0.446335) return 35;
  else if (BDT < 0.448032) return 36;
  else if (BDT < 0.449729) return 37;
  else if (BDT < 0.451408) return 38;
  else if (BDT < 0.453191) return 39;
  else if (BDT < 0.455084) return 40;
  else if (BDT < 0.457026) return 41;
  else if (BDT < 0.459034) return 42;
  else if (BDT < 0.460982) return 43;
  else if (BDT < 0.462957) return 44;
  else if (BDT < 0.464931) return 45;
  else if (BDT < 0.466846) return 46;
  else if (BDT < 0.469022) return 47;
  else if (BDT < 0.471469) return 48;
  else if (BDT < 0.474097) return 49;
  else if (BDT < 0.476797) return 50;
  else if (BDT < 0.479569) return 51;
  else if (BDT < 0.482351) return 52;
  else if (BDT < 0.485071) return 53;
  else if (BDT < 0.487771) return 54;
  else if (BDT < 0.490539) return 55;
  else if (BDT < 0.493233) return 56;
  else if (BDT < 0.495838) return 57;
  else if (BDT < 0.498487) return 58;
  else if (BDT < 0.501114) return 59;
  else if (BDT < 0.503584) return 60;
  else if (BDT < 0.506022) return 61;
  else if (BDT < 0.508531) return 62;
  else if (BDT < 0.510953) return 63;
  else if (BDT < 0.513375) return 64;
  else if (BDT < 0.515781) return 65;
  else if (BDT < 0.518260) return 66;
  else if (BDT < 0.520768) return 67;
  else if (BDT < 0.523297) return 68;
  else if (BDT < 0.525770) return 69;
  else if (BDT < 0.528125) return 70;
  else if (BDT < 0.530294) return 71;
  else if (BDT < 0.532430) return 72;
  else if (BDT < 0.534586) return 73;
  else if (BDT < 0.536635) return 74;
  else if (BDT < 0.538578) return 75;
  else if (BDT < 0.540463) return 76;
  else if (BDT < 0.542294) return 77;
  else if (BDT < 0.544030) return 78;
  else if (BDT < 0.545500) return 79;
  else if (BDT < 0.546843) return 80;
  else if (BDT < 0.548111) return 81;
  else if (BDT < 0.549452) return 82;
  else if (BDT < 0.551051) return 83;
  else if (BDT < 0.552943) return 84;
  else if (BDT < 0.554695) return 85;
  else if (BDT < 0.556529) return 86;
  else if (BDT < 0.558362) return 87;
  else if (BDT < 0.560054) return 88;
  else if (BDT < 0.562994) return 89;
  else                       return 90;
} 

    
Float_t getRF1j1b_opt_100bins(Double_t BDT){
  if      (BDT < 0.338093   ) return 1;
  else if (BDT < 0.347603) return 2;
  else if (BDT < 0.356761) return 3;
  else if (BDT < 0.364172) return 4;
  else if (BDT < 0.370521) return 5;
  else if (BDT < 0.376086) return 6;
  else if (BDT < 0.380898) return 7;
  else if (BDT < 0.385071) return 8;
  else if (BDT < 0.388891) return 9;
  else if (BDT < 0.392348) return 10;
  else if (BDT < 0.395510) return 11;
  else if (BDT < 0.398295) return 12;
  else if (BDT < 0.400789) return 13;
  else if (BDT < 0.403167) return 14;
  else if (BDT < 0.405315) return 15;
  else if (BDT < 0.407411) return 16;
  else if (BDT < 0.409325) return 17;
  else if (BDT < 0.411105) return 18;
  else if (BDT < 0.412843) return 19;
  else if (BDT < 0.414540) return 20;
  else if (BDT < 0.416241) return 21;
  else if (BDT < 0.418006) return 22;
  else if (BDT < 0.419678) return 23;
  else if (BDT < 0.421140) return 24;
  else if (BDT < 0.422710) return 25;
  else if (BDT < 0.424364) return 26;
  else if (BDT < 0.426093) return 27;
  else if (BDT < 0.427832) return 28;
  else if (BDT < 0.429612) return 29;
  else if (BDT < 0.431385) return 30;
  else if (BDT < 0.433172) return 31;
  else if (BDT < 0.434929) return 32;
  else if (BDT < 0.436653) return 33;
  else if (BDT < 0.438377) return 34;
  else if (BDT < 0.440062) return 35;
  else if (BDT < 0.441715) return 36;
  else if (BDT < 0.443337) return 37;
  else if (BDT < 0.444944) return 38;
  else if (BDT < 0.446506) return 39;
  else if (BDT < 0.448032) return 40;
  else if (BDT < 0.449561) return 41;
  else if (BDT < 0.451061) return 42;
  else if (BDT < 0.452646) return 43;
  else if (BDT < 0.454310) return 44;
  else if (BDT < 0.456059) return 45;
  else if (BDT < 0.457817) return 46;
  else if (BDT < 0.459613) return 47;
  else if (BDT < 0.461379) return 48;
  else if (BDT < 0.463155) return 49;
  else if (BDT < 0.464931) return 50;
  else if (BDT < 0.466654) return 51;
  else if (BDT < 0.468562) return 52;
  else if (BDT < 0.470703) return 53;
  else if (BDT < 0.473033) return 54;
  else if (BDT < 0.475433) return 55;
  else if (BDT < 0.477900) return 56;
  else if (BDT < 0.480391) return 57;
  else if (BDT < 0.482905) return 58;
  else if (BDT < 0.485343) return 59;
  else if (BDT < 0.487771) return 60;
  else if (BDT < 0.490264) return 61;
  else if (BDT < 0.492700) return 62;
  else if (BDT < 0.495063) return 63;
  else if (BDT < 0.497427) return 64;
  else if (BDT < 0.499823) return 65;
  else if (BDT < 0.502123) return 66;
  else if (BDT < 0.504311) return 67;
  else if (BDT < 0.506516) return 68;
  else if (BDT < 0.508779) return 69;
  else if (BDT < 0.510953) return 70;
  else if (BDT < 0.513131) return 71;
  else if (BDT < 0.515308) return 72;
  else if (BDT < 0.517512) return 73;
  else if (BDT < 0.519754) return 74;
  else if (BDT < 0.522045) return 75;
  else if (BDT < 0.524295) return 76;
  else if (BDT < 0.526496) return 77;
  else if (BDT < 0.528573) return 78;
  else if (BDT < 0.530507) return 79;
  else if (BDT < 0.532430) return 80;
  else if (BDT < 0.534371) return 81;
  else if (BDT < 0.536230) return 82;
  else if (BDT < 0.538005) return 83;
  else if (BDT < 0.539724) return 84;
  else if (BDT < 0.541378) return 85;
  else if (BDT < 0.543014) return 86;
  else if (BDT < 0.544486) return 87;
  else if (BDT < 0.545777) return 88;
  else if (BDT < 0.546975) return 89;
  else if (BDT < 0.548111) return 90;
  else if (BDT < 0.549309) return 91;
  else if (BDT < 0.550708) return 92;
  else if (BDT < 0.552367) return 93;
  else if (BDT < 0.554022) return 94;
  else if (BDT < 0.555574) return 95;
  else if (BDT < 0.557306) return 96;
  else if (BDT < 0.558847) return 97;
  else if (BDT < 0.560457) return 98;
  else if (BDT < 0.563487) return 99;
  else                       return 100;
} 
