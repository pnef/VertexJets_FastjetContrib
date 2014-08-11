// VertexJets package
// Quesstions / Comments: pascal.nef@slac.stanford.edu, sch@slac.stanford.edu
//
//
// Copyright (c) 2014
// Pascal Nef, Ariel Schwartzman
// 
// Compile this code as: 
//   make example 
// And run it with:
//   ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//
// $Id$
//
// Copyright (c) -, 
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <sstream>
#include "VertexJets.hh" // In external code, this should be fastjet/contrib/VertexJets.hh
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

// example user info class 
class MyUserInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
    MyUserInfo(const bool &is_pileup, const int charge):  _is_pileup(is_pileup), _charge(charge){}
    bool is_pileup()  const { return _is_pileup;}
    int  charge   ()  const { return _charge;}
 protected:
    bool _is_pileup; // true if pileup particle
    int  _charge;    // charge: -1, 0, 1
};


// Selector worker for charged pileup particles
class SelectorWorkerChargedPileup : public fastjet::SelectorWorker {
public:

  virtual bool pass(const fastjet::PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero
    return (particle.has_user_info<MyUserInfo>()  
            &&      particle.user_info<MyUserInfo>().is_pileup() == true 
            &&  abs(particle.user_info<MyUserInfo>().charge())   == 1     );
  }
  
  virtual string description() const {return "charged particle from pileup interaction";}
};

// Selector worker for charged hard-scatter particles
class SelectorWorkerChargedHardScatter : public fastjet::SelectorWorker {
public:

  virtual bool pass(const fastjet::PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero
    return (particle.has_user_info<MyUserInfo>()  
            &&      particle.user_info<MyUserInfo>().is_pileup() == false 
            &&  abs(particle.user_info<MyUserInfo>().charge())   == 1     );

  }
  
  virtual string description() const {return "charged particle from hard-scatter";}
};

// now the two selector functions
fastjet::Selector SelectorChargedPileup() {  
  return new SelectorWorkerChargedPileup();
}
fastjet::Selector SelectorChargedHardScatter() {  
  return new SelectorWorkerChargedHardScatter();
}


//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  // calculate pileup pT density
  JetMedianBackgroundEstimator bge(SelectorAbsRapMax(4.), JetDefinition(kt_algorithm,0.4), AreaDefinition(active_area_explicit_ghosts));
  bge.set_particles(event);
  Subtractor theSubtractor(&bge); 

  //----------------------------------------------------------
  // area definition
  AreaDefinition      area_def(active_area_explicit_ghosts);

  // cluster particles into small-R jets (anti-kt 0.4)
  JetDefinition         jet_def_small(antikt_algorithm, 0.4);
  ClusterSequenceArea   cs_small(event, jet_def_small, area_def);
  Selector sel_jets_small      = SelectorNHardest(4) * SelectorAbsRapMax(2.5);
  vector<PseudoJet> jets_small = sel_jets_small(theSubtractor(cs_small.inclusive_jets()));

  // cluster particles into large-R jets (anti-kt 1.0)
  JetDefinition       jet_def_large(antikt_algorithm, 1.0);
  ClusterSequenceArea cs_large(event, jet_def_large, area_def);
  Selector sel_jets_large       = SelectorNHardest(2) * SelectorAbsRapMax(1.5);
  vector<PseudoJet> jets_large  = sel_jets_large(cs_large.inclusive_jets());


  //----------------------------------------------------------
  // now use VertexJets to calculate corrJVF for small-R jets

  // get total number of charged pileup particles within the tracker coverage
  Selector sel_charged_pu_particles = SelectorChargedPileup();
  int n_pileup_tracks = sel_charged_pu_particles(event).size(); 

  // here comes VertexJets for small R jets 
  contrib::VertexJets    vertexjet(SelectorChargedPileup(), SelectorChargedHardScatter());
  vertexjet.set_tot_n_pu_tracks      (n_pileup_tracks);
  vertexjet.set_corrJVF_scale_factor (0.01);
  vertexjet.set_corrJVF_cut          (0.6);
  cout << vertexjet.description() << endl;
  cout << setprecision(3);
  
  // loop over small R jets
  cout << ">> anti-kt 0.4 jets: " << endl;
  for(unsigned int ij=0; ij < jets_small.size(); ++ij){
      PseudoJet resultjet = vertexjet(jets_small[ij]);
      cout << "      jet " << ij << ": pt "    << resultjet.pt() 
                                 << " eta "    << resultjet.eta() 
                                 << " phi "    << resultjet.phi()
                                 << " m "      << resultjet.m() 
                                 << " corrJVF "<< resultjet.structure_of<contrib::VertexJets>().corrJVF() << endl;
  }



  //----------------------------------------------------------
  // and for large Rjets in grooming mode
  contrib::VertexJets    vertexjet_largeR(SelectorChargedPileup(), SelectorChargedHardScatter(), JetDefinition(kt_algorithm, 0.3));
  vertexjet_largeR.set_tot_n_pu_tracks      (n_pileup_tracks);
  vertexjet_largeR.set_corrJVF_scale_factor (0.01);
  vertexjet_largeR.set_corrJVF_cut          (0.6);
  vertexjet_largeR.set_trimming_fcut        (0.05);
  vertexjet_largeR.set_subtractor           (&theSubtractor);
  cout << vertexjet_largeR.description() << endl;

  // loop over large R jets 
  cout << ">> anti-kt 1.0 jets: " << endl;
  for(unsigned int ij=0; ij< jets_large.size(); ++ij){
      PseudoJet resultjet = vertexjet_largeR(jets_large[ij]);
      cout << "      jet " << ij << ": pt "  << resultjet.pt() 
                                 << " eta "  << resultjet.eta()
                                 << " phi "  << resultjet.phi()
                                 << " mass " << resultjet.m() << endl;

      vector<PseudoJet> subjets = resultjet.pieces();
      for (unsigned int is=0; is<subjets.size(); ++is){
        cout << "         > subjet " << is << ": pt " << subjets[is].pt() << " corrJVF " << subjets[is].structure_of<contrib::VertexJets>().corrJVF() << endl;
      }
  }





  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
    string line;
    int  nsub  = 0;    // counter to keep track of which sub-event we're reading
    while (getline(cin, line)) {
        istringstream linestream(line);
        // characters (e.g. line-feed).
        if (line.substr(0,4) == "#END")      {break;}
        if (line.substr(0,9) == "#SUBSTART") {nsub += 1;}
        if (line.substr(0,1) == "#")         {continue;}
        double px,py,pz,E;
        int    pdg_id, charge;
        linestream >> px >> py >> pz >> E >> pdg_id >> charge;
        PseudoJet particle(px,py,pz,E);
        particle.reset_PtYPhiM(particle.pt(), particle.rapidity(), particle.phi(), 0.); // set particles massless
        

        // set charge to 0 if outside tracker
        if(fabs(particle.eta())> 2.4) charge = 0;                            
        
        // set user info
        particle.set_user_info(new MyUserInfo(nsub==1?false:true, charge));  // is_pileup, charge

        // push event into event vector
        event.push_back(particle);
    }
    
    // if there was nothing in the event
    if (nsub == 0) {
    cerr << "Error: read empty event\n";
    exit(-1);
    }
    
    cout << "# " << nsub-1 << " pileup events on top of the hard event" << endl;
}
