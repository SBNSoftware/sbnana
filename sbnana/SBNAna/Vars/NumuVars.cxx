#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{


  const Var kPrimMuonIdx([](const caf::SRSliceProxy *slc) -> double
                        {       //Find the most muon-like track
			  if( (int)slc->reco.npfp == 0 ) return -5.0;

                          double best_idx   = 0;
                          //double best_score = -5.0;
                          double best_len   = -5.0;
			  for( unsigned int trkIdx = 0; trkIdx < slc->reco.npfp; trkIdx++ ){
			    auto &trk = slc->reco.pfp[trkIdx].trk;
			    if(trk.chi2pid[2].pid_ndof < 0 || slc->reco.pfp[trkIdx].trackScore < 0.5) return -5.0;

			    //Find longest trk w/Chi2 for muon < Chi2 for pion
			    bool isMuonLike = trk.chi2pid[2].chi2_pion > trk.chi2pid[2].chi2_muon;
			    if( isMuonLike && trk.len > best_len ){ //Chi2 muon < Chi2 pion
			      // best_score = score;
			      best_len = trk.len;
			      best_idx = trkIdx;
			    }
			  } // end loop over trks

                          return best_idx;
                        });


  const Var kPrimTrkLen(
			[](const caf::SRSliceProxy *slc) -> double
			{
			  double len = -5.0;
			  if ( slc->reco.npfp > 0 ){
			    int muIdx = (int)kPrimMuonIdx(slc);
                            if (muIdx >= 0) {
  			      len = slc->reco.pfp[muIdx].trk.len;
                            }
			  }
			  return len;
			});

}
