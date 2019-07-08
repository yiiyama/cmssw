#include "CLHEP/Random/RandomEngine.h"

CLHEP::HepRandomEngine* _BeamHalo_randomEngine;

extern "C" {
  // C++ function to be called from FORTRAN - at some point we should just port all from FORTRAN to C++
  float bhgpyr_(int* idummy)
  {
    return (float)_BeamHalo_randomEngine->flat();
  }
}

