from . import getsyst

# Regen systematic variations
regen_systematics = [
    # CCQE
    "GENIEReWeight_ICARUS_v2_multisigma_VecFFCCQEshape",
    "GENIEReWeight_ICARUS_v2_multisigma_RPA_CCQE",
    "GENIEReWeight_ICARUS_v2_multisigma_CoulombCCQE",
    "GENIEReWeight_ICARUS_v2_multisim_ZExpAVariationResponse",
    
    # MEC
    "GENIEReWeight_ICARUS_v2_multisigma_NormCCMEC",
    "GENIEReWeight_ICARUS_v2_multisigma_NormNCMEC",
    "GENIEReWeight_ICARUS_v2_multisigma_DecayAngMEC",

    # RES
    "GENIEReWeight_ICARUS_v2_multisim_CCRESVariationResponse",
    "GENIEReWeight_ICARUS_v2_multisim_NCRESVariationResponse",
    "GENIEReWeight_ICARUS_v2_multisigma_RDecBR1gamma",
    "GENIEReWeight_ICARUS_v2_multisigma_RDecBR1eta",
    "GENIEReWeight_ICARUS_v2_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_ICARUS_v2_multisigma_ThetaDelta2NRad",

    # Non-Res
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvpCC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvpCC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvpNC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvpNC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvnCC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvnCC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvnNC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvnNC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarpCC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarpCC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarpNC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarpNC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarnCC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarnCC2pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarnNC1pi",
    "GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarnNC2pi",

    # DIS
    "GENIEReWeight_ICARUS_v2_multisim_DISBYVariationResponse",

    # COH
    "GENIEReWeight_ICARUS_v2_multisigma_NormCCCOH",
    "GENIEReWeight_ICARUS_v2_multisigma_NormNCCOH",

    # FSI
    "GENIEReWeight_ICARUS_v2_multisim_FSI_pi_VariationResponse",
    "GENIEReWeight_ICARUS_v2_multisim_FSI_N_VariationResponse",

    # NCEL
    "GENIEReWeight_ICARUS_v2_multisim_NCELVariationResponse",
]

def geniesyst(f, nuind):
    return getsyst.getsyst(f, regen_systematics, nuind)

