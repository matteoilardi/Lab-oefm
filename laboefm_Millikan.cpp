#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "TH1F.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"

using namespace std;

// ===========================================================
// Read data from file and store them in a vector
// ===========================================================

vector<double> ParseFile(string filename){
  
  ifstream fin(filename.c_str());

  vector<double> v;
  double val;
  
  if (!fin){
    cout << "Cannot open file " << filename << endl;
    exit(11);
  };
  
  while(fin >> val) v.push_back(val);
  
  fin.close();
  return v;
}

// ===========================================================                                                                                                                                              
// Compute S(q)                                                                                                                                                             
// ===========================================================     

double fun ( double q , vector<double> params ) {
  double sum = 0;
  for (int k = 0; k < params.size(); k++) sum+= pow(q - params[k]/(round(params[k]/q)),  2);
  return sum;
}

// ===========================================================                                                                                                                                                
// Compute q_min                                                                                                                                                                                        
// ===========================================================           

double deriv(double qmin, vector<double> params){
  double sum = 0;
  for (int k = 0; k < params.size(); k++) sum+= (params[k]/round(params[k]/qmin));  
  return sum/params.size();
}

// ===========================================================                                                                                                                                                      
// This code estimates the best value of qe from a set of 
// measurements (drop charges)                                                                                                                                                                       
// ===========================================================      

int main() {

  TApplication app(0,0,0);

  // read charges from file
  ofstream out("Millikan.output");
  vector<double> charges = ParseFile("data_millikan.dat");
  
  // show charges distribution

  TCanvas can1 ;
  can1.cd();
  TH1F histo("cariche","Charges distribution", 100, 0 , 20E-19);
  for(auto i = 0; i < charges.size(); i++ ) histo.Fill(charges[i]);
  histo.Draw();  
  histo.GetXaxis()->SetTitle("Charge [C]");
  
  TGraph g ;
  int counter = 0;
  double qmin = 0;
  double sqmin = DBL_MAX;

  for(double value = 1.4e-19; value < 1.8E-19; value+=0.0001E-19){
    g.SetPoint(counter, value,fun( value , charges ) );
    if (fun(value, charges) < sqmin){sqmin = fun(value, charges ); qmin = value;};
    counter++;
  }

  cout << "Found approximate minimum at q = " << qmin << endl;
  out << "Found approximate minimum at q = " << qmin << endl;

  TCanvas can2;
  can2.cd();
  g.Draw("ALP");
  g.SetMarkerStyle(20);
  g.SetMarkerSize(0.5);
  g.SetTitle("Best charge value");
  g.GetXaxis()->SetTitle("Charge (C)");
  g.GetYaxis()->SetTitle("S(q) (C^{2})");

  double mycharge = deriv(qmin, charges);
  double uncer = sqrt(fun(mycharge, charges)/(charges.size() * (charges.size()-1))); 
  cout << "Measured charge = " << mycharge << " +/- " << uncer << "(stat only)" << endl;
  out << "Measured charge = " << mycharge << " +/- " << uncer << "(stat only)" << endl;

  out << "Charges (ascending order) :" << endl;
  sort(charges.begin(), charges.end());
  for(int i=0; i<charges.size(); i++){
    out << charges[i] << endl;
  }

  out.close();
  app.Run();

}