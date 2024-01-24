#ifndef MICROPHYSIcSVARIABLES_H
#define MICROPHYSIcSVARIABLES_H

namespace MicVar {
   enum {
      // independent variables
      qt = 0,
      qp,
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      rho,   // density
      pres,  // pressure
      // derived variables
      qv,   // water vapor
      qn,   // cloud condensate (liquid+ice), initial to zero
      qci,  // cloud ice
      qcl,  // cloud water
      qpl,  // precip rain
      qpi,  // precip ice
      qg,   // graupel
      // temporary variable
      omega,
      NumVars
  };
}

namespace MicVar_Kess {
   enum {
      // independent variables
      qt = 0,
      qp,
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      rho,   // density
      pres,  // pressure
      // derived variables
      qr,   // rain water
      qv,   // water vapor
      qn,   // cloud condensate (liquid+ice), initial to zero
      qci,  // cloud ice
      qcl,  // cloud water
      qpl,  // precip rain
      qpi,  // precip ice
      qg,   // graupel
      // temporary variable
      omega,
      NumVars
  };
}

namespace MicVar_FE {
   enum {
      // independent variables
      qv = 0,
      qc,
      qt,
      rho,   // density
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      pres,  // pressure
      NumVars
  };
}

#endif
