#ifndef __sbnanalysis_core_Event__
#define __sbnanalysis_core_Event__

/**
 * \file Event.hh
 *
 * The standard minimum output tree.
 *
 * This event structure is written out by every Processor subclass.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include <map>
#include <string>
#include <vector>
#include <TTree.h>
#include <TVector3.h>

/**
 * \class Event
 * \brief The standard event data definition.
 */
class Event {
public:
  /**
   * \class Event::Metadata
   * \brief Event-level information
   */
  class Metadata {
  public:
    int run;  //!< Run ID
    int subrun;  //!< Subrun ID
    int eventID;  //!< Event ID
  };

  /**
   * \class Event::Neutrino
   * \brief Neutrino interaction information
   */
  class Neutrino {
  public:
    bool ccnc;  //!< CC (true) or NC (false)
    int pdg;  //!< PDG code of neutrino
    int targetPDG;  //!< PDG code of struck target
    int intcode;  //!< Interaction code (as for LArSoft MCNeutrino)
    double bjorkenX;  //!< Bjorken x
    double inelasticityY;  //!< Inelasticity y
    double q2;  //!< Q squared
    double w;  //!< Hadronic invariant mass W
    double t;  //!< Kinematic t
    double energy;  //!< Neutrino energy
    TVector3 momentum;  //!< Neutrino three-momentum
  };

  /**
   * \class Event::FinalStateParticle
   * \brief Final state particle information
   */
  class FinalStateParticle {
  public:
    int pdg;  //!< PDG Code
    double energy;  //!< Energy
    TVector3 momentum;  //!< Three-momentum
  };

  /**
   * \class Event::Interaction
   * \brief All truth information associated with one neutrino interaction
   */
  class Interaction {
  public:
    Neutrino neutrino;  //!< The neutrino
    FinalStateParticle lepton;  //!< The primary final state lepton
    std::vector<FinalStateParticle> finalstate;  //!< The other final state particles
  };

  Metadata metadata;  //!< Event metadata
  std::vector<Interaction> interactions;  //!< All interactions

  /**
   * Event weights.
   *
   * This is a map from the weight calculator name to the list of weights
   * for all the sampled universes.
   */
  std::map<std::string, std::vector<double> > weights;
};

#endif  // __sbnanalysis_core_Event__

