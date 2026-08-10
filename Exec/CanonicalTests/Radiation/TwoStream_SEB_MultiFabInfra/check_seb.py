#!/usr/bin/env python3
import os, sys, numpy as np

def read_diag(path):
    if not os.path.exists(path):
        return None, [f"missing diagnostics file: {path}"]
    try:
        data = np.genfromtxt(path, delimiter=",", skip_header=1)
    except Exception as exc:
        return None, [f"failed reading {path}: {exc}"]
    if data.size == 0:
        return None, [f"empty diagnostics file: {path}"]
    if not np.all(np.isfinite(data[:, [1,3,4,5,6,7]])):
        return None, [f"non-finite diagnostic values in {path}"]
    return data, []

def main():
    errors = []
    disabled, err = read_diag('radiation_seb_disabled_diag.dat'); errors += err
    enabled, err = read_diag('radiation_seb_enabled_diag.dat'); errors += err
    if disabled is not None and enabled is not None:
        if disabled.shape != enabled.shape:
            errors.append('enabled/disabled diagnostic shapes differ')
        else:
            if not np.allclose(disabled[:,3:8], enabled[:,3:8], rtol=0.0, atol=0.0):
                errors.append('seb_enable changed legacy radiation diagnostics; expected bitwise-identical outputs')
    if errors:
        print('SEB validation FAILED')
        for e in errors: print(' -', e)
        return 1
    print('SEB validation PASSED')
    return 0

if __name__ == '__main__':
    sys.exit(main())
