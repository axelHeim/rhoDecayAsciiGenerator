#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include <TFile.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Generator.cc"

using namespace std;

#define E_BEAM 190.0
#define PROT_MASS 0.938
#define PI_MASS 0.13957
#define PI0_MASS 0.13497
#define RHO_MASS 0.7755

#define PI_MASS_SQ (PI_MASS*PI_MASS)
#define RHO_MASS_SQ (RHO_MASS*RHO_MASS)
#define PROT_MASS_SQ (PROT_MASS*PROT_MASS)

double invariantMass(std::vector<double>, std::vector<double>);
void approx_HitPosECAL2(vector<double>, double hitpos[3]);
double approx_distance_between2Gamma_ECAL2(vector<double>, vector<double>);

int main() {
  TFile file1("Rho_gammas.root", "RECREATE");
  TH2D *histo_xyHitsECAL2;
  histo_xyHitsECAL2 = new TH2D("xyHitsECAL2", "xyHitsECAL2", 65, -1225.6,
           1263.9, 48, -919.2, 919.2);
  TH1D *histo_2gamma_dist_ECAL2;
  histo_2gamma_dist_ECAL2 = new TH1D("2gamma_dist_ECAL2 [mm]","2gamma_dist_ECAL2 [mm]",1000,0,3500);




  std::random_device rd;
  std::mt19937 en(rd());
  std::uniform_real_distribution<> dist_cosT(-1.0, 1.0);
  std::uniform_real_distribution<> dist_phi(-M_PI, M_PI);

  double s0 = PI_MASS_SQ+PROT_MASS_SQ + 2*E_BEAM*PROT_MASS_SQ;
  double t = -0.4;  // GeV^2
  double cosTh_from_t = (2*s0*(t-PI_MASS_SQ - RHO_MASS_SQ) +
                         (s0 + PI_MASS_SQ - PROT_MASS_SQ)*
                         (s0+ RHO_MASS_SQ - PROT_MASS_SQ) ) /
                        sqrt( LAMBDA(s0,  PI_MASS_SQ, PROT_MASS_SQ)*
                              LAMBDA(s0, RHO_MASS_SQ, PROT_MASS_SQ));

for(int o = 1; o <= 1000 ; o++)
{
  fstream file;
  file.open(TString::Format("ascii_events%d",o), ios::out);

  const int Nevents = 10;
  for (int i = 0; i < Nevents; i++) {
          std::vector<double> piProt {0.0, 0.0, sqrt(E_BEAM*E_BEAM-PI_MASS_SQ), E_BEAM+PROT_MASS};
          auto rhoprot = Generator::decay_p(piProt,RHO_MASS,PROT_MASS, cosTh_from_t, dist_phi(en));
          auto rho = rhoprot[0];
          auto pions = Generator::decay_p(rho,PI_MASS,PI0_MASS, dist_cosT(en), dist_phi(en));
          auto pi0 = pions[1];
          auto gammas = Generator::decay_p(pi0,0.0,0.0, dist_cosT(en), dist_phi(en));
          // print photons
          file << "22" << "   "
            << gammas[0][0] << "  "
            << gammas[0][1] << "  "
            << gammas[0][2] << "  "
            << endl
            << "22" << "   "
            << gammas[1][0] << "  "
            << gammas[1][1] << "  "
            << gammas[1][2] << "  "
            << endl;
          file << "--------------------------------------------------" << endl;

          




          // double hitpos[3];
          // approx_HitPosECAL2(gammas[0], hitpos);
          // histo_xyHitsECAL2->Fill(hitpos[0], hitpos[1]);
          // approx_HitPosECAL2(gammas[1], hitpos);
          // histo_xyHitsECAL2->Fill(hitpos[0], hitpos[1]);

          //cout <<  approx_distance_between2Gamma_ECAL2(gammas[0], gammas[1]) << '\n';
          //histo_2gamma_dist_ECAL2->Fill(approx_distance_between2Gamma_ECAL2(gammas[0], gammas[1]));


  }

  file.close();
}


        histo_2gamma_dist_ECAL2->Write();
        histo_xyHitsECAL2->Write();
        return 0.0;
}

double approx_distance_between2Gamma_ECAL2(vector<double> fourVec1, vector<double> fourVec2)
{
  double ECAL2_zPos = 33252; //ECAL2 z position in mm
  double hitpos1[3];
   hitpos1[0] = 0;
   hitpos1[1] = 0;
   hitpos1[2] = -30.0;  //approx. the primary vertex pos

  double mom_absValue = sqrt(pow(fourVec1[0],2)+pow(fourVec1[1],2)+pow(fourVec1[2],2));

  for(int j = 0; ; j++)
  {
    for(int i = 0; i<3; i++)
    {
       hitpos1[i] += 0.01 *  fourVec1[i]/mom_absValue;
    }
    if(fabs(hitpos1[2] - ECAL2_zPos) < 0.5)
    {
      //cout << "hier" << '\n';
      break;}
  }

  double  hitpos2[3];
    hitpos2[0] = 0;
    hitpos2[1] = 0;
    hitpos2[2] = -30.0;  //approx. the primary vertex pos

  mom_absValue = sqrt(pow(fourVec2[0],2)+pow(fourVec2[1],2)+pow(fourVec2[2],2));

  for(int j = 0; ; j++)
  {
    for(int i = 0; i<3; i++)
    {
        hitpos2[i] += 0.01 * fourVec2[i]/mom_absValue;
    }
    if(fabs(hitpos2[2] - ECAL2_zPos) < 0.5)
    {
      //cout << "hier" << '\n';
      break;}
  }

  double distance = sqrt(pow(hitpos1[0] - hitpos2[0],2) +
                        pow(hitpos1[1] - hitpos2[1],2));

  return distance;
}



double invariantMass(std::vector<double> fourVec1, std::vector<double> fourVec2)
{
  double invariant_mass = sqrt(pow(fourVec1[0] + fourVec2[0],2) +
                              pow(fourVec1[1] + fourVec2[1],2)  +
                              pow(fourVec1[2] + fourVec2[2],2)  +
                              pow(fourVec1[3] + fourVec2[3],2));

  return invariant_mass;
}

void approx_HitPosECAL2(vector<double> fourVec, double hitpos[3])
{
  double ECAL2_zPos = 33252; //ECAL2 z position in mm
  hitpos[0] = 0;
  hitpos[1] = 0;
  hitpos[2] = -30.0;  //approx. the primary vertex pos

  double mom_absValue = sqrt(pow(fourVec[0],2)+pow(fourVec[1],2)+pow(fourVec[2],2));

  for(int j = 0; ; j++)
  {
    for(int i = 0; i<3; i++)
    {
      hitpos[i] += 0.01 * fourVec[i]/mom_absValue;
    }
    if(fabs(hitpos[2] - ECAL2_zPos) < 0.5)
    {
      //cout << "hier" << '\n';
      break;}
  }

}
