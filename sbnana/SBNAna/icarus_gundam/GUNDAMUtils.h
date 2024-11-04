namespace ICARUSGUNDAM{

inline std::vector<std::string> GetGENIEMultisigmaKnobNames(){

  return {
"ZExpA1CCQE",
"ZExpA2CCQE",
"ZExpA3CCQE",
"ZExpA4CCQE",
"RPA_CCQE",
"CoulombCCQE",
"NormCCMEC",
"NormNCMEC",
"MaNCEL",
"EtaNCEL",
"MaCCRES",
"MvCCRES",
"MaNCRES",
"MvNCRES",
"NonRESBGvpCC1pi",
"NonRESBGvpCC2pi",
"NonRESBGvpNC1pi",
"NonRESBGvpNC2pi",
"NonRESBGvnCC1pi",
"NonRESBGvnCC2pi",
"NonRESBGvnNC1pi",
"NonRESBGvnNC2pi",
"NonRESBGvbarpCC1pi",
"NonRESBGvbarpCC2pi",
"NonRESBGvbarpNC1pi",
"NonRESBGvbarpNC2pi",
"NonRESBGvbarnCC1pi",
"NonRESBGvbarnCC2pi",
"NonRESBGvbarnNC1pi",
"NonRESBGvbarnNC2pi",
"RDecBR1gamma",
"RDecBR1eta",
"NormCCCOH",
"NormNCCOH",
"AhtBY",
"BhtBY",
"CV1uBY",
"CV2uBY",
"MFP_pi",
"FrCEx_pi",
"FrInel_pi",
"FrAbs_pi",
"FrPiProd_pi",
"MFP_N",
"FrCEx_N",
"FrInel_N",
"FrAbs_N",
"FrPiProd_N",
  };

}

inline std::vector<std::string> GetGENIEDependentKnobNames(){

  return {
"ZExpAVariationResponse",
"NCELVariationResponse",
"CCRESVariationResponse",
"NCRESVariationResponse",
"DISBYVariationResponse",
"FSI_pi_VariationResponse",
"FSI_N_VariationResponse",
  };

}

inline std::vector<std::string> GetGENIEMorphKnobNames(){
  return {
"VecFFCCQEshape",
"DecayAngMEC",
"Theta_Delta2Npi",
"ThetaDelta2NRad",
  };
}

} // END name space
