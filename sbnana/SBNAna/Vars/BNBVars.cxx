//! ////////////////////////////////////////////////////////////////////////////
//!  @file: BNBVars.cxx
//!  @author: Jacob Smith (smithja)
//!  @email: jacob.a.smith@stonybrook.edu
//!
//!  @brief Explicit definitions of variables used in BNB Quality Cuts for the
//!  ICARUS experiment. 
//!  
//!  @details Most SpillVars here are declared with Bruce Howard's
//!  SIMPLESPILLVAR() function allowing SpillVars to be easily integrated with
//!  creating Spectra (like what is done with regular Vars). Other composite
//!  SpillVars are created and some SpillVars have aliases.
//! ////////////////////////////////////////////////////////////////////////////

#include "sbnana/SBNAna/Vars/BNBVars.h"
#include "sbnana/SBNAna/Vars/getBNBFoM.h"

#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <iostream>

namespace ana {
    /** Multi-wire Readout Fit Parameters:
     * @note We determine beam sigma primarily by fitting gaussian curves to
     * multi-wire device data; we can also extract beam position from fits.
     * Beam position should be taken from beam position monitor variables
     * above, but can be taken from below if needed.
     *
     * @note Sigmas and positions provided for both horizontal and vertical
     * directions perpendicular to the beam direction.
     * 
     * units: mm
    */

    /** @var kSpillMW875HorSigma */
    const SpillVar kSpillMW875HorSigma = SIMPLESPILLVAR( hdr.spillbnbinfo.M875HS);

    /** @var kSpillMW875VerSigma */
    const SpillVar kSpillMW875VerSigma = SIMPLESPILLVAR( hdr.spillbnbinfo.M875VS);

    /** @var kSpillMW876HorSigma */
    const SpillVar kSpillMW876HorSigma = SIMPLESPILLVAR( hdr.spillbnbinfo.M876HS);

    /** @var kSpillMW876VerSigma */
    const SpillVar kSpillMW876VerSigma = SIMPLESPILLVAR( hdr.spillbnbinfo.M876VS);

    //! ////////////////////////////////////////////////////////////////////////
    //!   ^--- Multi-wire Fit Parameters ; Other BNB Variables ---v
    //! ////////////////////////////////////////////////////////////////////////

    //! Readout time for a given spill.
    //! units: s (i.e. seconds)

    /** @var kSpillTimeSec */
    const SpillVar kSpillTimeSec = SIMPLESPILLVAR( hdr.spillbnbinfo.spill_time_sec);


    //! Toroid Devices: monitor protons on target (POT)
    //! units: POT

    /** @var kSpillTOR860 */
    const SpillVar kSpillTOR860 = SIMPLESPILLVAR( hdr.spillbnbinfo.TOR860);
    /** @var kSpillTOR875 */
    const SpillVar kSpillTOR875 = SIMPLESPILLVAR( hdr.spillbnbinfo.TOR875);


    //! Use closest toroid device to target to get POT reading.

    /** @var kSpillTOR */
    const SpillVar kSpillTOR([](const caf::SRSpillProxy *sr) {
        if (kSpillTOR875(sr) > 0) { //! use TOR875 if it has data
            return kSpillTOR875(sr);
        }
        else { //! otherwise, use TOR860
            return kSpillTOR860(sr);
        }
    });


    //! Loss Monitors: measure backscattering of beam just upstream of target
    //! 
    //! @note higher reading is better, implies beam is hitting more of target
    //!
    //! units: rad / s (i.e. rads per second)

    /** @var kSpillLM875A */
    const SpillVar kSpillLM875A = SIMPLESPILLVAR( hdr.spillbnbinfo.LM875A);
    /** @var kSpillLM875B */
    const SpillVar kSpillLM875B = SIMPLESPILLVAR( hdr.spillbnbinfo.LM875B);
    /** @var kSpillLM875C */
    const SpillVar kSpillLM875C = SIMPLESPILLVAR( hdr.spillbnbinfo.LM875C);


    //! Use closest loss monitor to the target to get back-scattering radiation
    //! reading.

    /** @var kSpillLM875 */
    const SpillVar kSpillLM875([](const caf::SRSpillProxy *sr) {
        if (kSpillLM875C(sr) > 0) {       //! use LM875C if it has data
            return kSpillLM875C(sr);
        }
        else if (kSpillLM875B(sr) > 0) {  //! otherwise, use LM875B if it has data
            return kSpillLM875B(sr);
        }
        else {
            return kSpillLM875A(sr);      //! other-otherwise, default to LM875A
        }
    });


    //! Beam Position Monitors: measure beam position at different points 
    //! upstream of target
    //!
    //! @note not always centered at (0, 0); @see SBNAna/Cuts/BNBQualityCuts.cxx
    //! for information on nominal beam position for different runs
    //! 
    //! @note VPTG1 and VPTG2 unreliable for ICARUS (and presumably SBND) runs;
    //! VP873 identified as a suitable replacement
    //! 
    //! units: mm

    /** @var kSpillVP873 */
    const SpillVar kSpillVP873 = SIMPLESPILLVAR( hdr.spillbnbinfo.VP873);

    /** @var kSpillHP875 */
    const SpillVar kSpillHP875 = SIMPLESPILLVAR( hdr.spillbnbinfo.HP875);
    /** @var kSpillVP875 */
    const SpillVar kSpillVP875 = SIMPLESPILLVAR( hdr.spillbnbinfo.VP875);

    /** @var kSpillHPTG1 */
    const SpillVar kSpillHPTG1 = SIMPLESPILLVAR( hdr.spillbnbinfo.HPTG1);
    /** @var kSpillVPTG1 */
//!    const SpillVar kSpillVPTG1 = SIMPLESPILLVAR( hdr.spillbnbinfo.VPTG1);

    /** @var kSpillHPTG2 */
    const SpillVar kSpillHPTG2 = SIMPLESPILLVAR( hdr.spillbnbinfo.HPTG2);
    /** @var kSpillVPTG2 */
//!    const SpillVar kSpillVPTG2 = SIMPLESPILLVAR( hdr.spillbnbinfo.VPTG2);


    //! @deprecated
    //! Temperature reading near target.
    //!
    //! @note this beam monitoring device provided diminishing returns when 
    //! combined with all of the other devices listed here. additionally, 
    //! temperature readings can be correlated between spills if they are close 
    //! enough in time
    //!
    //! units: degrees celsius

    /** @var kSpillBTJT2 */
//!    const SpillVar kSpillBTJT2 = SIMPLESPILLVAR( hdr.spillbnbinfo.BTJT2);


    //! Focusing Horn Current
    //! units: kA

    /** @var kSpillTHCURR */
    const SpillVar kSpillTHCURR = SIMPLESPILLVAR( hdr.spillbnbinfo.THCURR);





    //!/////////////////////////////////////////////////////////////////////////
    //! Figure(s) of Merit: @see getBNBFoM.cxx for details
    //!/////////////////////////////////////////////////////////////////////////

    //! @note This Figure of Merit DOES NOT take BNB width into account.
    /** @var kSpillFoM_noMultiWire */
    const SpillVar kSpillFoM_noMultiWire([](const caf::SRSpillProxy *sr) {
        //! Extract variables used in FoM calculations and cast as double
        //! for compatibility with getBNBFoM_noMultiWire().
        double kSpillTimeSecVal = kSpillTimeSec(sr);
        double kSpillTOR860Val  = kSpillTOR860(sr);
        double kSpillTOR875Val  = kSpillTOR875(sr);
        double kSpillHP875Val   = kSpillHP875(sr);
        double kSpillHPTG1Val   = kSpillHPTG1(sr);
        double kSpillHPTG2Val   = kSpillHPTG2(sr);
        double kSpillVP873Val   = kSpillVP873(sr);
        double kSpillVP875Val   = kSpillVP875(sr);

        double fom = getBNBFoM_noMultiWire( kSpillTimeSecVal, 
            kSpillTOR860Val, kSpillTOR875Val,
            kSpillHP875Val, kSpillHPTG1Val, kSpillHPTG2Val,
            kSpillVP873Val, kSpillVP875Val);
        return fom;
    });

    //! @note This Figure of Merit DOES take BNB width into account.
    /** @var kSpillFoM */
    const SpillVar kSpillFoM([](const caf::SRSpillProxy *sr) {
        //! Extract variables used in FoM calculations and cast as double
        //! for compatibility with getBNBFoM().
        double kSpillTimeSecVal       = kSpillTimeSec(sr);
        double kSpillTOR860Val        = kSpillTOR860(sr);
        double kSpillTOR875Val        = kSpillTOR875(sr);
        double kSpillHP875Val         = kSpillHP875(sr);
        double kSpillHPTG1Val         = kSpillHPTG1(sr);
        double kSpillHPTG2Val         = kSpillHPTG2(sr);
        double kSpillVP873Val         = kSpillVP873(sr);
        double kSpillVP875Val         = kSpillVP875(sr);
        double kSpillMW875HorSigmaVal = kSpillMW875HorSigma(sr);
        double kSpillMW875VerSigmaVal = kSpillMW875VerSigma(sr);
        double kSpillMW876HorSigmaVal = kSpillMW876HorSigma(sr);
        double kSpillMW876VerSigmaVal = kSpillMW876VerSigma(sr);

        double fom = getBNBFoM( kSpillTimeSecVal, 
            kSpillTOR860Val, kSpillTOR875Val,
            kSpillHP875Val, kSpillHPTG1Val, kSpillHPTG2Val,
            kSpillVP873Val, kSpillVP875Val,
            kSpillMW875HorSigmaVal, kSpillMW875VerSigmaVal,
            kSpillMW876HorSigmaVal, kSpillMW876VerSigmaVal
        );
        return fom;
    });
}