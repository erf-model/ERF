    ParmParse pp_prob("prob");

    Real advection_heating_rate = zero;
    pp_prob.query("advection_heating_rate", advection_heating_rate);

    src->setVal(advection_heating_rate);
