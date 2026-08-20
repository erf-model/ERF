import numpy as np

def compute_max_reflectivity_dbz(
    rho_air,
    tmk,
    qra,
    qsn,
    qgr,
    in0r=1,
    in0s=1,
    in0g=1,
    iliqskin=1,
):
    """
    Compute radar reflectivity (dBZ).

    Parameters
    ----------
    rho_air : ndarray
        Air density [kg/m^3]

    tmk : ndarray
        Temperature [K]

    qra : ndarray
        Rain mixing ratio [kg/kg]

    qsn : ndarray
        Snow mixing ratio [kg/kg]

    qgr : ndarray
        Graupel mixing ratio [kg/kg]

    Returns
    -------
    dbz : ndarray
        Reflectivity [dBZ]
    """

    # ------------------------------------------------------------------
    # Constants
    # ------------------------------------------------------------------

    pi = np.pi

    rhowat = 1000.0
    rho_r = 1000.0
    rho_s = 100.0
    rho_g = 400.0

    celkel = 273.15

    alpha = 0.224
    gamma_seven = 720.0

    rn0_r = 8.0e6
    rn0_s = 2.0e7
    rn0_g = 4.0e6

    r1 = 1.0e-15

    ron2 = 1.0e10
    gon = 5.0e7

    ron_min = 8.0e6
    ron_qr0 = 1.0e-4
    ron_delqr0 = 0.25 * ron_qr0

    ron_const1r = 0.5 * (ron2 - ron_min)
    ron_const2r = 0.5 * (ron2 + ron_min)

    # ------------------------------------------------------------------
    # Precompute factors
    # ------------------------------------------------------------------

    factor_r = (
        gamma_seven
        * 1.0e18
        * (1.0 / (pi * rho_r)) ** 1.75
    )

    factor_s = (
        gamma_seven
        * 1.0e18
        * (1.0 / (pi * rho_s)) ** 1.75
        * (rho_s / rhowat) ** 2
        * alpha
    )

    factor_g = (
        gamma_seven
        * 1.0e18
        * (1.0 / (pi * rho_g)) ** 1.75
        * (rho_g / rhowat) ** 2
        * alpha
    )

    # ------------------------------------------------------------------
    # Bright-band correction
    # ------------------------------------------------------------------

    factorb_s = np.full_like(tmk, factor_s)
    factorb_g = np.full_like(tmk, factor_g)

    if iliqskin == 1:
        warm = tmk > celkel
        factorb_s[warm] = factor_s / alpha
        factorb_g[warm] = factor_g / alpha

    # ------------------------------------------------------------------
    # Snow intercept
    # ------------------------------------------------------------------

    sonv = np.full_like(tmk, rn0_s)

    if in0s == 1:
        temp_c = np.minimum(-0.001, tmk - celkel)

        sonv = np.minimum(
            2.0e8,
            2.0e6 * np.exp(-0.12 * temp_c)
        )

    # ------------------------------------------------------------------
    # Graupel intercept
    # ------------------------------------------------------------------

    gonv = np.full_like(tmk, rn0_g)

    if in0g == 1:
        mask = qgr > r1

        gonv[mask] = 2.38 * (
            pi * rho_g /
            (rho_air[mask] * qgr[mask])
        ) ** 0.92

        gonv[mask] = np.clip(gonv[mask], 1.0e4, gon)

    # ------------------------------------------------------------------
    # Rain intercept
    # ------------------------------------------------------------------

    ronv = np.full_like(tmk, rn0_r)

    if in0r == 1:
        ronv[:] = ron2

        mask = qra > r1

        ronv_calc = (
            ron_const1r
            * np.tanh((ron_qr0 - qra) / ron_delqr0)
            + ron_const2r
        )

        ronv[mask] = ronv_calc[mask]

    # ------------------------------------------------------------------
    # Equivalent reflectivity
    # ------------------------------------------------------------------

    z_e = (
        factor_r
        * (rho_air * qra) ** 1.75
        / ronv ** 0.75
        +
        factorb_s
        * (rho_air * qsn) ** 1.75
        / sonv ** 0.75
        +
        factorb_g
        * (rho_air * qgr) ** 1.75
        / gonv ** 0.75
    )

    z_e = np.maximum(z_e, 0.01)

    dbz = 10.0 * np.log10(z_e)

    return dbz


def compute_Thunderstorm_Squall_HeavyRain (nz, rho_air, tmk, qr, qs, qg):
    dbz_3d = compute_max_reflectivity_dbz(rho_air, tmk, qr, qs, qg)
    max_dbz_2d = np.max(dbz_3d, axis=2)
    max_dbz_3d = np.repeat(max_dbz_2d[:, :, np.newaxis], nz, axis=2)

    return max_dbz_2d, max_dbz_3d
