// VertexJets package
// Quesstions / Comments: pascal.nef@slac.stanford.edu, sch@slac.stanford.edu
//
//
// Copyright (c) 2014
// Pascal Nef, Ariel Schwartzman
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

#include "VertexJets.hh"
#include <algorithm> 
#include <numeric>
#include <sstream>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

    // ---------------------------
    // result method: actual worker 
    PseudoJet VertexJets::result(const PseudoJet &jet) const {

        // first perform some checks
        if (! jet.has_constituents()){
            throw Error("result(): jet does not have constituents!");
        }
        if(_tot_n_pu_tracks <0){
            throw Error("result(): tot_n_pu_tracks must be >= 0, use setter function set_tot_n_pu_tracks() to set the variable for each event");
        }
        if(_corrJVF_cut <-1 || _corrJVF_cut > 1){
            throw Error("result(): _corrJVF_cut must be between -1 and 1, set to -1 if no cut should be applied");
        }
        if(_corrJVF_scale_factor <= 0){
            throw Error("result(): _corrJVF_scale_factor must be > 0"); 
        }



        // either grooming_mode or tagging_mode
        if(_run_mode == grooming_mode){

            // use filter with SelectorIdentity() to reconstruct all area-corrected subjets
            PseudoJet tmp_jet = _filter(jet);

            // loop over subjets and compute corrJVF
            std::vector<PseudoJet> subjets = tmp_jet.pieces();
            for(unsigned int isub=0; isub < subjets.size(); ++isub){
                subjets[isub] = CalculateCorrJVF(subjets[isub]);
            }
            
            // setup trimming selector, defined as MUTABLE selector
            _selector_pTfracMin.set_reference(jet); 
            
            // setup combined selector: corrJVF and trimming
            Selector selector_combined = _selector_corrJVF && _selector_pTfracMin;

            // get subjets that pass both selectors
            std::vector<PseudoJet> selected_subjets = selector_combined(subjets); // subjets that pass corrJVF cut
            PseudoJet joined = join(selected_subjets, *_def_subjet.recombiner()); // recombine subjets 

            return joined;

        }else{
            PseudoJet resultjet = jet;
            resultjet = CalculateCorrJVF(resultjet);

            if(  _selector_corrJVF.pass(resultjet) ) return resultjet;           // jet passes corrJVF cut, return jet
            else                                     return resultjet*0;         // jet fails corrJVF cut,  return 0*jet
        }
        
    }

    // ---------------------------
    // Description
    std::string VertexJets::description() const {
         std::ostringstream ostr;
         ostr << std::endl;
         ostr << "VertexJets >> class for corrJVF based tagging or grooming of jets." << std::endl;
         if(_run_mode == tagging_mode){
         ostr << "           >> running in tagging mode : corrJVF calculated for full jet" << std::endl;
         }else{
         ostr << "           >> running in grooming mode: corrJVF calculated for subjets" << std::endl;;
         }
         ostr << "           >> _selector_pu_ghosttracks: " << _selector_pu_ghosttracks.description() << std::endl; 
         ostr << "           >> _selector_hs_ghosttracks: " << _selector_hs_ghosttracks.description() << std::endl; 
         ostr << "           >> _selector_pTfracMin     : " << _selector_pTfracMin.description()      << std::endl; 
         ostr << "           >> _selector_corrJVF       : " << _selector_corrJVF.description()        << std::endl; 
         ostr << "           >> _def_subjet             : " << _def_subjet.description()              << std::endl;
         if(_subtractor!=0){
         ostr << "           >> _subtractor             : " << _subtractor->description()             << std::endl;
         }else{
         ostr << "           >> _subtractor             : " << "NULL"                                 << std::endl;
         }

         ostr << "           >> _ghost_scale_factor     : " << _ghost_scale_factor                    << std::endl;
         ostr << "           >> _tot_n_pu_tracks        : " << _tot_n_pu_tracks                       << std::endl;
         ostr << "           >> _corrJVF_scale_factor   : " << _corrJVF_scale_factor                  << std::endl;
         ostr << "           >> _corrJVF_cut            : " << _corrJVF_cut  
                                                            << "   (-1 -> no cut. Note: jets failing the cut get their momenta rescaled by 0)"    
                                                            << std::endl;
         return ostr.str();
    }
        
    //-------------------------------
    // CalculateCorrJVF method
    PseudoJet VertexJets::CalculateCorrJVF(PseudoJet &jet) const{
        // define new vector<PseudoJet> to contain PU tracks, HS tracks and other(discarded)
        std::vector<PseudoJet> pu_ghost_tracks, hs_ghost_tracks, other;

        // get pu tracks
        _selector_pu_ghosttracks.sift(jet.constituents(), pu_ghost_tracks, other);
        // get hs tracks
        _selector_hs_ghosttracks.sift(jet.constituents(), hs_ghost_tracks, other);

    
        // sum of pT of pu tracks: correct momenta for ghost scale factor (if applicable)
        float pu_sum_pt = 0, hs_sum_pt = 0;
        for(unsigned int it=0;it<pu_ghost_tracks.size();++it){
            pu_sum_pt += 1./_ghost_scale_factor * pu_ghost_tracks[it].pt();
        }
        for(unsigned int it=0;it<hs_ghost_tracks.size();++it){
            hs_sum_pt += 1./_ghost_scale_factor * hs_ghost_tracks[it].pt();
        }
        
        // calculate corrJVF 
        float corrJVF = -1;
        if(pu_sum_pt + hs_sum_pt >0){
            corrJVF = hs_sum_pt / (hs_sum_pt + (pu_sum_pt / (std::max(1, _tot_n_pu_tracks) * _corrJVF_scale_factor) ) );
        }

        // create jet with VertexJetStructure structure
        PseudoJet resultjet    = join<VertexJetStructure>(jet);
        StructureType *s = (StructureType *) resultjet.structure_non_const_ptr();
        s->_corrJVF = corrJVF;
        
        // and return the jet. 
        return resultjet;
    }

    //--------------------------
    // Selector based on corrJVF cut
    Selector VertexJets::SelectorCorrJVF(float corrJVF_cut) {
        return new SelectorWorkerCorrJVF(corrJVF_cut);
    }

} // namespace contrib

FASTJET_END_NAMESPACE
