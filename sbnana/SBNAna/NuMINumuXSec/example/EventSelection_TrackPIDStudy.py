from SelectionInfo import *

EventSelectionInfos = []

## For muon study

EventSelectionInfos.append(
  SelectionInfo(
    Name = "ContainedMuonPresel",
    Cut = "kNoCut && cutRFiducial && cutNotClearCosmic && cutCRLongestTrackDirYHard && ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackContained && !ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackIsochronous && ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLengthCut",
    SpillCut = "kNoSpillCut",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "ExitingMuonPresel",
    Cut = "kNoCut && cutRFiducial && cutNotClearCosmic && cutCRLongestTrackDirYHard && !ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackContained && !ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackIsochronous && ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLengthCut",
    SpillCut = "kNoSpillCut",
  )
)

## For proton study

EventSelectionInfos.append(
  SelectionInfo(
    Name = "ProtonPresel",
    Cut = "kNoCut && cutRFiducial && cutNotClearCosmic && cutCRLongestTrackDirYHard && ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackContained && !ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackIsochronous",
    SpillCut = "kNoSpillCut",
  )
)


