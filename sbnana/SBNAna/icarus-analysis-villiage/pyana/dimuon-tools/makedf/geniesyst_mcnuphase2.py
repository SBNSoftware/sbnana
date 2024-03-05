from . import getsyst

# MCNuPhase2 systematic variations
mcnuphase2_systematics = [
'GENIEReWeight_ICARUS_v1_multisigma_MaNCEL',
 'GENIEReWeight_ICARUS_v1_multisigma_EtaNCEL',
 'GENIEReWeight_ICARUS_v1_multisigma_MaCCRES',
 'GENIEReWeight_ICARUS_v1_multisigma_MvCCRES',
 'GENIEReWeight_ICARUS_v1_multisigma_MaNCRES',
 'GENIEReWeight_ICARUS_v1_multisigma_MvNCRES',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvpCC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvpCC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvpNC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvpNC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvnCC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvnCC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvnNC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvnNC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarpCC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarpCC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarpNC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarpNC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarnCC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarnCC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarnNC1pi',
 'GENIEReWeight_ICARUS_v1_multisigma_NonRESBGvbarnNC2pi',
 'GENIEReWeight_ICARUS_v1_multisigma_RDecBR1gamma',
 'GENIEReWeight_ICARUS_v1_multisigma_RDecBR1eta',
 'GENIEReWeight_ICARUS_v1_multisigma_Theta_Delta2Npi',
 'GENIEReWeight_ICARUS_v1_multisigma_AhtBY',
 'GENIEReWeight_ICARUS_v1_multisigma_BhtBY',
 'GENIEReWeight_ICARUS_v1_multisigma_CV1uBY',
 'GENIEReWeight_ICARUS_v1_multisigma_CV2uBY',
 'GENIEReWeight_ICARUS_v1_multisigma_MFP_pi',
 'GENIEReWeight_ICARUS_v1_multisigma_FrCEx_pi',
 'GENIEReWeight_ICARUS_v1_multisigma_FrInel_pi',
 'GENIEReWeight_ICARUS_v1_multisigma_FrAbs_pi',
 'GENIEReWeight_ICARUS_v1_multisigma_FrPiProd_pi',
 'GENIEReWeight_ICARUS_v1_multisigma_MFP_N',
 'GENIEReWeight_ICARUS_v1_multisigma_FrCEx_N',
 'GENIEReWeight_ICARUS_v1_multisigma_FrInel_N',
 'GENIEReWeight_ICARUS_v1_multisigma_FrAbs_N',
 'GENIEReWeight_ICARUS_v1_multisigma_FrPiProd_N',
]

def geniesyst(f, nuind):
    return getsyst.getsyst(f, mcnuphase2_systematics, nuind)

