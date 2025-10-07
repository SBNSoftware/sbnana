#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/Pi0StudyVars.h"
#include "sbnana/SBNAna/Vars/NuGraph2Vars.h"
#include "sbnana/SBNAna/Cuts/Pi0StudyCuts.h"

#include "TVector3.h"

namespace ana {

  const Var kMuonCandidate_SemCat([](const caf::SRSliceProxy* slc) -> int {
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return -9999.f;
    if ( std::isnan(slc->reco.pfp[Idx].ngscore.sem_cat) || std::isinf(slc->reco.pfp[Idx].ngscore.sem_cat) ) return -9999.f;
    return slc->reco.pfp[Idx].ngscore.sem_cat;
  });

  const Var kPi0LeadingPhoton_SemCat([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    //int sem_cat = slc->reco.pfp[idx].ngscore.sem_cat;
    if ( std::isnan(slc->reco.pfp[idx].ngscore.sem_cat) || std::isinf(slc->reco.pfp[idx].ngscore.sem_cat) ) return -5.f;
    return slc->reco.pfp[idx].ngscore.sem_cat;
  });

  const Var kPi0SubLeadingPhoton_SemCat([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    //int sem_cat = slc->reco.pfp[idx].ngscore.sem_cat;
    if ( std::isnan(slc->reco.pfp[idx].ngscore.sem_cat) || std::isinf(slc->reco.pfp[idx].ngscore.sem_cat) ) return -5.f;
    return slc->reco.pfp[idx].ngscore.sem_cat;
  });
}





