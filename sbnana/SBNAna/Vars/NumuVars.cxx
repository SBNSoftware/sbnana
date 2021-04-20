#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{


  const Var kPrimMuonIdx([](const caf::SRSliceProxy *slc) -> double
                        {       //Find the most muon-like track
			  if( (int)slc->reco.ntrk == 0 ) return -5.0;

                          double best_idx   = 0;
                          //double best_score = -5.0;
                          double best_len   = -5.0;
			  for( unsigned int trkIdx = 0; trkIdx < slc->reco.ntrk; trkIdx++ ){
			    auto &trk = slc->reco.trk[trkIdx];
			    if(trk.chi2pid2.pid_ndof < 0 ) return -5.0;

			    //Find longest trk w/Chi2 for muon < Chi2 for pion
			    bool isMuonLike = trk.chi2pid2.chi2_pion > trk.chi2pid2.chi2_muon;
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
			  if ( slc->reco.ntrk > 0 ){
			    int muIdx = (int)kPrimMuonIdx(slc);
                            if (muIdx >= 0) {
  			      len = slc->reco.trk[muIdx].len;
                            }
			  }
			  return len;
			});

}
