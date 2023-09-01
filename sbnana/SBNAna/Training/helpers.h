#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"

#include "TMath.h"

using namespace ana;


const double ClosestDistanceBetweenLines(TVector3 a0, TVector3 a1, TVector3 b0, TVector3 b1, bool clampAll):

    bool campA0 = false;
    bool campA1 = false;
    bool campB0 = false;
    bool campB1 = false;

    if (clampAll) {
        clampA0=True;
        clampA1=True;
        clampB0=True;
        clampB1=True;
    }

    // Calculate denomitator
    TVector3 A = a1 - a0;
    TVector3 B = b1 - b0;
    double magA = A.Mag();
    double magB = B.Mag();
    
    TVectro3 _A = A / magA;
    TVectro3 _B = B / magB;
    
    TVector3 cross = (_A.Cross(_B);
    double denom = cross.Mag()*cross.Mag();
    
    
    // If lines are parallel (denom=0) test if lines overlap.
    // If they don't overlap then there is a closest point solution.
    // If they do overlap, there are infinite closest positions, but there is a closest distance
if (denom == 0)
        {
            TVector3 d0 = _A.Dot((b0 - a0));
 
            // Overlap only possible with clamping
 
            TVector3 d1 = _A.Dot((b1 - a0));
 
            // Is segment B before A?
            if (d0 <= 0 && 0 >= d1)
            {
                if (std::abs(d0) < std::abs(d1))
                {
                    result.Line1Closest = a0;
                    result.Line2Closest = b0;
                    result.Distance = (a0 - b0).magnitude;
 
                    return result;
                }
                result.Line1Closest = a0;
                result.Line2Closest = b1;
                result.Distance = (a0 - b1).magnitude;
 
                return result;
            }
            // Is segment B after A?
            else if (d0 >= magA && magA <= d1)
            {
                if (std::abs(d0) < std::abs(d1))
    
            {
                    result.Line1Closest = a1;
                    result.Line2Closest = b0;
                    result.Distance = (a1 - b0).magnitude;
 
                    return result;
                }
                result.Line1Closest = a1;
                result.Line2Closest = b1;
                result.Distance = (a1 - b1).magnitude;
 
                return result;
            }
 
            // Segments overlap, return distance between parallel segments
            result.Line1Closest = Vector3.zero;
            result.Line2Closest = Vector3.zero;
            result.Distance = (((d0 * _A) + a0) - b0).magnitude;
            return result;
        }
 
 
        // Lines criss-cross: Calculate the projected closest points
        var t = (b0 - a0);
        var detA = Determinant(t, _B, cross);
        var detB = Determinant(t, _A, cross);
 
        var t0 = detA / denom;
        var t1 = detB / denom;
 
        var pA = a0 + (_A * t0); // Projected closest point on segment A
        var pB = b0 + (_B * t1); // Projected closest point on segment B
 
 
        // Clamp projections
        if (t0 < 0)
            pA = a0;
        else if (t0 > magA)
            pA = a1;
 
        if (t1 < 0)
            pB = b0;
        else if (t1 > magB)
            pB = b1;
 
        float dot;
        // Clamp projection A
        if (t0 < 0 || t0 > magA)
        {
            dot = Vector3.Dot(_B, (pA - b0));
            if (dot < 0)
                dot = 0;
            else if (dot > magB)
                dot = magB;
            pB = b0 + (_B * dot);
        }
        // Clamp projection B
        if (t1 < 0 || t1 > magB)
        {
            dot = Vector3.Dot(_A, (pB - a0));
            if (dot < 0)
                dot = 0;
            else if (dot > magA)
                dot = magA;
            pA = a0 + (_A * dot);
        }
 
        result.Line1Closest = pA;
        result.Line2Closest = pB;
        result.Distance = (pA - pB).magnitude;
        return result;
    }
