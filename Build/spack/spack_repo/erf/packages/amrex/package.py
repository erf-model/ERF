# Spack's builtin amrex recipe, plus the one thing ERF needs that it lacks: a
# particle-data precision of its own (AMReX_PARTICLES_PRECISION), so that an
# ERF built with particles_precision different from precision, as ERF's own
# Windows CI does, has an AMReX to link against. The same change is proposed
# for spack/spack-packages (branch larenspear/spack-packages:amrex-particles-
# precision); once it is merged there, delete this file and the erf recipe's
# `amrex particles_precision=` dependencies keep working unchanged.
#
# Registering this repository therefore shadows builtin.amrex; the override
# changes nothing unless particles_precision is set to double or single.

from spack_repo.builtin.packages.amrex.package import Amrex as BuiltinAmrex

from spack.package import *


class Amrex(BuiltinAmrex):
    """AMReX, Spack's builtin recipe plus a particle-precision variant (see above)."""

    # AMReX_PARTICLES_PRECISION; AMReX's CMake defaults it to the mesh precision.
    variant(
        "particles_precision",
        default="same",
        values=("same", "double", "single"),
        multi=False,
        description="Precision of the particle data ('same' follows precision)",
        when="+particles",
    )

    def cmake_args(self):
        args = super().cmake_args()
        if self.spec.satisfies("+particles") and not self.spec.satisfies(
            "particles_precision=same"
        ):
            # The builtin recipe configures AMReX in xSDK mode, where AMReX's
            # CMake assigns the particle precision from XSDK_PRECISION with a
            # plain set() that shadows any -DAMReX_PARTICLES_PRECISION. Leave
            # xSDK mode for this build and state both precisions directly;
            # everything else xSDK mode would default (BUILD_SHARED_LIBS) the
            # builtin recipe already passes explicitly.
            args = [a for a in args if not a.startswith("-DUSE_XSDK_DEFAULTS")]
            args.extend(
                [
                    self.define("USE_XSDK_DEFAULTS", False),
                    self.define("AMReX_PRECISION", self.spec.variants["precision"].value.upper()),
                    self.define(
                        "AMReX_PARTICLES_PRECISION",
                        self.spec.variants["particles_precision"].value.upper(),
                    ),
                ]
            )
        return args
