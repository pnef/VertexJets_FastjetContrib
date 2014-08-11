// VertexJets package
// Quesstions / Comments: pascal.nef@slac.stanford.edu, sch@slac.stanford.edu
//
//
// Copyright (c) 2014
// Pascal Nef, Ariel Schwartzman
//
// $Id$
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

#ifndef __FASTJET_CONTRIB_VERTEXJETS_HH__
#define __FASTJET_CONTRIB_VERTEXJETS_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Error.hh"
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/tools/Transformer.hh>     // to derive VertexJets from Transformer
#include <fastjet/CompositeJetStructure.hh>
#include <fastjet/tools/Subtractor.hh>
#include <algorithm>                        // std::max
#include <sstream>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

class VertexJets;
class VertexJetStructure;
class SelectorWorkerCorrJVF;



//------------------------------------------------------------------------
/// \class VertexJets
/// Class for corrJVF based tagging or grooming of jets.
/// 
/// corrJVF is a modified version of the jet-vertex fraction JVF,
/// aiming at correcting the variable in average for its pileup dependence
/// (see below).
///
/// VertexJets has two run_modes:
/// 1: Calculate corrJVF for the input PseudoJet:
///    If a corrJVF cut value is specified and the jet fails the cut, 
///    the original jet * 0 is returned.
/// 2: corrJVF based grooming:
///    The constituents of the jet are reclustered into subjets.
///    A corrJVF cut and an optional f_cut (trimming) is applied 
///    to the subjets. The ones passing the selection are recombined and
///    returned as a PseudoJet. 
///    
/// In the grooming mode, the subjets that were kept are accessible through:
/// \code
///    vector<PseudoJet> kept_subjets = jet.pieces();
/// \endcode
/// 
/// The corrJVF value of the jet (tagging mode) or subjets (grooming mode) 
/// is accessible through:
/// \code
///    float corrJVF = jet.structure_of<VertexJets>().corrJVF();
/// \endcode
/// 
/// ----------- corrJVF definition ----------------------------------------
/// corrJVF is defined jet-by-jet as:
///
/// corrJVF = HS_pt / (HS_pt + PU_pt/n_PU_part*1/scale_fact),
/// where 
/// HS_pt       = scalar pt sum of charged particles from the hard-scatter
/// PU_pt       = scalar pt sum of charged particles from pileup interactions
/// n_PU_part   = total number of charged PU particles in the *event*
/// scale_fact  = optional scaling parameter. default value is 0.01
/// 
/// For details, please refer to ATLAS-CONF-2014-018 or send us an email. 
/// 
/// ----------- Usage -----------------------------------------------------
/// The charged hard-scatter and pileup constituents can either be:
/// 1: true constituents (as in CMS-style Pflow jets)
/// 2: ghost associated tracks (as in ATLAS-style calo jets)
///    in this case, set_ghost_scale_factor() must be used to 
///    specify the scale factor that was used to scale the 4 momentum 
///    of the ghost particles. 
///
/// Constructor: 
/// The constructor takes as input two fastjet Selectors. These are 
/// implemented by the user and specify which jet constituents are
/// charged pileup particles and which are charged particles from the
/// hard-scatter. See example.cc for for an example implementation. 
/// 
/// The alternative constructor has an additional argument: the JetDefinition
/// to reconstruct subjets. If this constructor is used, the grooming mode
/// is enabed. 
///
/// Setter functions:
///
/// set_tot_n_pu_tracks():      To set the total number of charged pileup 
///                             particles in the event
///
/// set_corrJVF_scale_factor(): Optional corrJVF scale factor
///
/// set_ghost_scale_factor():   Only used if tracks are ghost associated. 
///                             specify scale factor that was used to 
///                             rescale the 4 momentum of the ghosts.
///                             The particles original momentum should be
///                             p * 1./scale_factor
///
/// set_trimming_fcut():        Optional trimming cut. Only functional in 
///                             grooming mode. 
///
/// set_subtractor():           Optionally, a Subtractor can be provided 
///                             to perform a background subtraction on the 
///                             subjets in grooming mode
/// 
//------------------------------------------------------------------------
class VertexJets : public Transformer{
public:

  /// default ctor
  VertexJets(const Selector selector_pu_ghosttracks, const Selector selector_hs_ghosttracks )
      :
        _run_mode(tagging_mode),
        _selector_pu_ghosttracks(selector_pu_ghosttracks),
        _selector_hs_ghosttracks(selector_hs_ghosttracks),
        _selector_corrJVF(SelectorCorrJVF(-1)),
        _selector_pTfracMin(SelectorPtFractionMin(0)),
        _def_subjet(JetDefinition ()),
        _ghost_scale_factor(1.),
        _tot_n_pu_tracks(0),
        _corrJVF_scale_factor(1.),
        _corrJVF_cut(-1.),
        _subtractor(0),
        _filter(_def_subjet, SelectorIdentity()){};

  VertexJets(const Selector selector_pu_ghosttracks, const Selector selector_hs_ghosttracks, const JetDefinition def_subjet)
      :
        _run_mode(grooming_mode),
        _selector_pu_ghosttracks(selector_pu_ghosttracks),
        _selector_hs_ghosttracks(selector_hs_ghosttracks),
        _selector_corrJVF(SelectorCorrJVF(-1)),
        _selector_pTfracMin(SelectorPtFractionMin(0)),
        _def_subjet(def_subjet),
        _ghost_scale_factor(1.),
        _tot_n_pu_tracks(0),
        _corrJVF_scale_factor(1.),
        _corrJVF_cut(-1),
        _subtractor(0),
        _filter(_def_subjet, SelectorIdentity()){};


  /// default dtor
  virtual ~VertexJets(){};
  
  // setter fuction: tot_n_pu_tracks
  inline void set_tot_n_pu_tracks(int tot_n_pu_tracks){
      _tot_n_pu_tracks = tot_n_pu_tracks;
  }

  // setter function: corrJVF_scale_factor
  inline void set_corrJVF_scale_factor(float corrJVF_scale_factor){
      _corrJVF_scale_factor = corrJVF_scale_factor;
  }
  
  // setter function: _corrJVF_cut
  inline void set_corrJVF_cut(float corrJVF_cut){
      _corrJVF_cut      = corrJVF_cut;
      _selector_corrJVF = SelectorCorrJVF(corrJVF_cut);
  }
  
  // setter function: ghost_scale_factor
  inline void set_ghost_scale_factor(float ghost_scale_factor){
      _ghost_scale_factor = ghost_scale_factor;
  }

  // optional trimming for subjet mode
  inline void set_trimming_fcut(float fcut){
      _selector_pTfracMin = SelectorPtFractionMin(fcut);
  }

  // optionally specify subtractor
  inline void set_subtractor(Subtractor *subtractor){
      _subtractor = subtractor;
      _filter.set_subtractor( _subtractor);
  }

  /// the result of the Transformer acting on the PseudoJet.
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// Description
  virtual std::string description() const;

  // Structure of resulting jet
  typedef VertexJetStructure StructureType;
  
protected:
  
  // Selector based on _corrJVF_cut
  Selector SelectorCorrJVF(float corrJVF_cut);
  
  /// calculate corrJVF
  virtual PseudoJet CalculateCorrJVF(PseudoJet & jet) const;

  enum run_mode{
    grooming_mode,
    tagging_mode
  }  _run_mode;

  Selector           _selector_pu_ghosttracks;
  Selector           _selector_hs_ghosttracks;
  Selector           _selector_corrJVF;
  mutable Selector   _selector_pTfracMin; // needs to be mutable since reference jet only set in result method
  JetDefinition      _def_subjet;
  float              _ghost_scale_factor;
  int                _tot_n_pu_tracks;
  float              _corrJVF_scale_factor;
  float              _corrJVF_cut;
  const Transformer  *_subtractor;
  mutable Filter     _filter;
  

};


//------------------------------------------------------------------------
/// \class VertexJetStructure
///
class VertexJetStructure : public CompositeJetStructure{
public:
    /// ctor with initialisation
    VertexJetStructure(const std::vector<PseudoJet> & pieces, const JetDefinition::Recombiner *recombiner = 0)
        : CompositeJetStructure(pieces, recombiner), _corrJVF(-1){}

    /// return corrJVF of the jet
    inline double corrJVF() const {return _corrJVF;}


protected:
    double _corrJVF;               /// corrJVF

    // allow VertexJets to access these
    friend class VertexJets;
};




// -------------------------------------------------------------------------
/// \class SelectorWorkerCorrJVF
class SelectorWorkerCorrJVF : public fastjet::SelectorWorker {
public:

    SelectorWorkerCorrJVF(float corrJVFcut) : _corrJVFcut(corrJVFcut) {}

    virtual bool pass(const fastjet::PseudoJet & jet) const {
        return (jet.has_structure_of<VertexJets>() && jet.structure_of<VertexJets>().corrJVF() >= _corrJVFcut );
    }
  
    virtual std::string description() const {
        std::ostringstream ostr;
        ostr << "selecting jets passing corrJVF >= " << _corrJVFcut;
        return ostr.str();
    }
private:
    float _corrJVFcut;
};


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_VERTEXJETS_HH__
