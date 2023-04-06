from SelectionInfo import *

SampleSelectionInfos = []

trackTypes = [
"ContainedNuMuon", "ExitingNuMuon", "CosmicMuon", "StoppingProton", "InelProton", "OtherProton", "StoppingChargedPion", "ExitingChargedPion", "Other",
]
for trackType in trackTypes:
  SampleSelectionInfos.append(
    SelectionInfo(
      Name = "RelaxedMuonTrackTruth%s"%(trackType),
      Cut = "ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruth%s"%(trackType),
    )
  )
  SampleSelectionInfos.append(
    SelectionInfo(
      Name = "RelaxedProtonTrackTruth%s"%(trackType),
      Cut = "ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruth%s"%(trackType),
    )
  )


SampleSelectionInfos.append(
  SelectionInfo(
    Name = "RelaxedMuonTrackTruthStoppingProtonFromNC",
    Cut = "ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthStoppingProton && cutIsNC",
  )
)
SampleSelectionInfos.append(
  SelectionInfo(
    Name = "RelaxedMuonTrackTruthStoppingProtonFromCC",
    Cut = "ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthStoppingProton && cutIsCC",
  )
)
