#pragma once

#include "CAFAnaCore/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// \brief Representation of a variable to be retrieved from a \ref
  /// caf::StandardRecord object
  ///
  /// A Var consists of a function, taking a StandardRecord and returning the
  /// value of the variable (which may be some complicated function).
  typedef _Var<caf::SRSliceProxy> Var;

  /// \brief Equivalent of \ref Var acting on \ref caf::SRSpill
  typedef _Var<caf::SRSpillProxy> SpillVar;

  /// \brief For Vars where literally all you need is a single CAF variable
  ///
  /// eg Var myVar = SIMPLEVAR(my.var.str);
  /// NB lack of quotes around my.var.str
#define SIMPLEVAR(CAFNAME) Var([](const caf::SRSliceProxy* sr) -> double {return sr->CAFNAME;})

#define SIMPLESPILLVAR(CAFNAME) SpillVar([](const caf::SRSpillProxy* sr) -> double {return sr->CAFNAME;})


  typedef _Var<caf::SRTrackProxy> TrackVar;
  typedef _Var<caf::SRShowerProxy> ShowerVar;
  typedef _Var<caf::SRStubProxy> StubVar;

#define SIMPLETRACKVAR(CAFNAME) TrackVar([](const caf::SRTrackProxy* trk) -> double {return trk->CAFNAME;})
#define SIMPLESHOWERVAR(CAFNAME) ShowerVar([](const caf::SRShowerProxy* shw) -> double {return shw->CAFNAME;})
#define SIMPLESTUBVAR(CAFNAME) StubVar([](const caf::SRStubProxy* stub) -> double {return stub->CAFNAME;})
} // namespace
