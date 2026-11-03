# ERF's Spack package; Build/spack/README.md says how to register it. It uses Spack's
# builtin repository (spack/spack-packages), which must be registered too.

import itertools
import os

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage

from spack.package import *


def _amrex_submodule(package):
    """Only the AMReX submodule (a plain list would mean every submodule)."""
    return ["Submodules/AMReX"]


class Erf(CMakePackage, CudaPackage, ROCmPackage):
    """ERF solves the compressible Navier-Stokes equations on an Arakawa C-grid
    for large-scale weather modeling.
    """

    homepage = "https://erf.readthedocs.io/en/latest/index.html"
    url = "https://github.com/erf-model/ERF/archive/refs/tags/26.09.tar.gz"
    git = "https://github.com/erf-model/ERF.git"

    maintainers("larenspear", "jmsexton03")

    license("BSD-3-Clause", checked_by="larenspear")

    # @develop builds the AMReX commit its submodule points at, as ERF's own build does: the
    # commit variant on amrex@develop (a git-commit version in a depends_on is not a
    # concretizer candidate; the variant is). check_pins.py compares the pin with the gitlink
    # and rewrites it with --fix. A release builds an AMReX release, so that the amrex
    # recipe's version-keyed logic (patches, the tracer-header fix for 24.10-26.07) applies:
    # the AMReX release of the same month, which contains the development commits ERF's
    # submodule pointed at, unless the build matrix showed ERF needs the next one (24.10
    # calls amrex::decompose(), added after amrex@24.10).
    _amrex_develop = "1d9bd90c35837dd54b20e5e8f94651ebfdcf5f96"
    _amrex_newer = {"24.10": "24.11"}
    version("develop", branch="development")
    depends_on(f"amrex@develop commit={_amrex_develop}", when="@develop", type=("build", "link"))
    for _v, _sha in (
        ("26.09", "d86372f6b9d85996f1b30beb7353c01b33f1bf0e454c0e63c913b255e17dadc0"),
        ("26.06", "379ffdebb30242cacddde95357511d5de1da29332c6624d5470088927743db78"),
        ("26.05", "198f8093b7afb2e405f509beca3b6a6da94c299af564e889e3c80daecc376e65"),
        ("26.04", "93a7ef8e6b55c4f830fac85129d6cf36c0b0e931eb060d7f8a1579afea8f909d"),
        ("26.03", "90399f4a7904e2b0d3c8e785484bebb89f2a5c46826fa4d751f939c0678d2d47"),
        ("26.02", "22b6570b9f022e58e0952801a85d4634b3076b40953422e0bd4df3de1d22e44d"),
        ("26.01", "1b40ec52b225e6bed1455f074cec04db1337121683ae1300c3adb3d048596d4d"),
        ("25.12", "0d11389c59561dea51f6b12253e1af7ab157599e0657bde603db7c102713e158"),
        ("25.11", "bd4160686eb52a49173ba07fe68ef86c2cd20b17283eaacf01b6c2d61a82982f"),
        ("25.10", "92575d25a8e87266f1b3a9a0c16a5755b84c45d489a060abeec6ddea8b8d8fe0"),
        ("25.08", "c8723384e00fb4bbb694385bdbc5187f38ee1509a0a18a23efc1094212fd21f4"),
        ("25.07", "bee8c136872c4e5400fdf2f6a432c8dc17f9b8384715158d9780a2c561a51642"),
        ("25.06", "dd1d627faf67477e3fa1f58b006f7e2b00316c66e8c99baf8b6cfbedce37b5f2"),
        ("25.05", "eb56a44d915d6af3ad30438dff6b96d9c489ae7c610fe4cdd0c17f2f26cc4b28"),
        ("25.04", "fc6837b252c1cebd0f2d38cdb3cea8ad55c837e4163c76e5d3508d973b282a42"),
        ("25.03", "3068c0c5f66538bcc1a12ee8c3eb142e5d02b73b71ab1cfbb244a053e8d5cf4a"),
        ("25.01", "11e48aadfc420c3f9adb05a0c5778776d6cb0915476532fd7a095d9f167ef584"),
        ("24.11", "19e7ce4829a46c98baf55b1dc4b2cefc508f2db41182e47c0ca788fe9280783f"),
        ("24.10", "7c9ba35374e71103fab8c5ba3ea107b4574e3270038bc3e8cbdcbe30eba26e7a"),
    ):
        version(_v, sha256=_sha)
        depends_on(f"amrex@{_amrex_newer.get(_v, _v)}", when=f"@{_v}", type=("build", "link"))

    # Releases before 24.10 have no ERF_USE_INTERNAL_AMREX: they build the AMReX submodule
    # they vendor, so they are git checkouts with that submodule and depend on no amrex.
    for _v, _commit in (
        ("24.09", "16e75b3e4de96ba9c4ca38b945d2a452d028244a"),
        ("24.08", "e3ff1bf8b76da36adac9bbfe83fc06eafa4cbb16"),
        ("24.06", "6f794f2fa7ac618a2dff77a573f973946809bab9"),
        ("24.05", "22931d4f402ade07b60440bcc7f678f585e8adec"),
        ("24.04", "5eb01f4a611d5eb80ab323468657db8ec53be9e5"),
        ("24.03", "5621a11955c27f57b9955a11d611b41c70cb42d3"),
        ("24.02", "a7762e8530d86b1cc654c0c1d5e45f07196fc350"),
        ("24.01", "850aa4079b6b71eeb376be0558e1927d3b999775"),
        ("23.12", "c0a1c9f212511c200a337c1b55f65e157a92df48"),
        ("23.11", "885a8c432fd8c81dc32863bfc51c6ee1542531bc"),
        ("23.10", "de410c8d03f61a008cafb2326c44a9a739785359"),
    ):
        version(_v, tag=_v, commit=_commit, submodules=_amrex_submodule)

    variant("mpi", default=True, description="Enable MPI support")
    variant("openmp", default=False, description="Enable OpenMP support")
    variant("sycl", default=False, description="Enable the SYCL GPU backend (needs %oneapi)")
    variant("netcdf", default=False, description="Enable NetCDF I/O (parallel, needs MPI)")
    # ERF's CMake defaults ERF_ENABLE_HDF5 to the NetCDF setting.
    variant("hdf5", default=True, description="Enable AMReX HDF5 plotfile I/O", when="+netcdf")
    variant("hdf5", default=False, description="Enable AMReX HDF5 plotfile I/O", when="~netcdf")
    variant("particles", default=False, description="Enable Lagrangian particles", when="@24.02:")
    variant("multiblock", default=False, description="Enable multiblock support")
    variant("fft", default=False, description="Enable the AMReX FFT-based solvers", when="@25.01:")
    variant("tests", default=False, description="Build the regression and unit tests")
    variant("fcompare", default=False, description="Build AMReX's fcompare alongside ERF")
    # Single precision exists from 26.06; earlier releases hard-code double.
    variant(
        "precision",
        default="double",
        values=("double",),
        multi=False,
        description="Floating point precision of the mesh and particle data",
        when="@:26.05",
    )
    variant(
        "precision",
        default="double",
        values=("double", "single"),
        multi=False,
        description="Floating point precision of the mesh and particle data",
        when="@26.06:",
    )
    variant("tools", default=False, description="Build the additional tools", when="@25.07:")
    variant(
        "rrtmgp",
        default=False,
        description="Enable RTE-RRTMGP radiation (builds EKAT and Kokkos)",
        when="@25.10:",
    )
    variant(
        "noahmp",
        default=False,
        description="Enable the Noah-MP land surface model",
        when="@25.12:",
    )
    variant(
        "windfarm",
        default=False,
        description="Enable the wind farm parameterizations",
        when="@26.02:",
    )
    variant(
        "ml_uphys_diagnostics",
        default=False,
        description="Enable the ML microphysics diagnostics in the SDM particle microphysics",
        when="@26.02:",
    )
    variant(
        "wsm6_fort",
        default=False,
        description="Enable the WSM6 Fortran microphysics bridge",
        when="@26.05:",
    )
    variant(
        "wdm6_fort",
        default=False,
        description="Enable the WDM6 Fortran microphysics bridge",
        when="@26.09:",
    )
    variant(
        "eamxx_shoc",
        default=False,
        description="Build the EAMxx SHOC turbulence parameterization (EKAT, Kokkos and "
        "a sparse checkout of E3SM's EAMxx sources at the commit ERF pins)",
        when="@26.09:",
    )
    # ERF_ENABLE_FIRE is on development only; 26.10 is the first release that can carry it.
    variant(
        "fire",
        default=False,
        description="Enable the fire coupling (atmosphere side; fire engines are external)",
        when="@26.10:",
    )
    variant(
        "unit_tests",
        default=True,
        description="Build the GoogleTest unit tests alongside the regression tests",
        when="@26.06:+tests",
    )
    variant(
        "particles_precision",
        default="same",
        values=("same", "double", "single"),
        multi=False,
        description="Floating point precision of the particle data ('same' follows precision)",
        when="@26.06:",
    )

    with default_args(type="build"):
        depends_on("cmake@3.20:")
        depends_on("c")
        depends_on("cxx")
        depends_on("fortran", when="+noahmp")
        depends_on("fortran", when="+wsm6_fort")
        depends_on("fortran", when="+wdm6_fort")
        depends_on("fortran", when="+rrtmgp")
        depends_on("fortran", when="+eamxx_shoc")
        # The Noah-MP ERF driver generates part of its interface with Python (26.09+).
        depends_on("python@3:", when="+noahmp")
        # ERF's FindNetCDF tries pkg-config before nc-config.
        depends_on("pkgconfig", when="+netcdf")
        # Replaces the bundled googletest submodule, which is not fetched (see patch()).
        depends_on("googletest", when="@26.06:+tests+unit_tests")

    with default_args(type=("build", "link")):
        # From 24.10 ERF links an installed AMReX (ERF_USE_INTERNAL_AMREX=OFF).
        with when("@24.10:"):
            for _v in ("mpi", "openmp", "cuda", "rocm", "sycl"):
                depends_on(f"amrex+{_v}", when=f"+{_v}")
                depends_on(f"amrex~{_v}", when=f"~{_v}")
            for _sm in CudaPackage.cuda_arch_values:
                depends_on(f"amrex+cuda cuda_arch={_sm}", when=f"+cuda cuda_arch={_sm}")
            for _tgt in ROCmPackage.amdgpu_targets:
                depends_on(f"amrex+rocm amdgpu_target={_tgt}", when=f"+rocm amdgpu_target={_tgt}")

            # ERF uses only AMReX's HDF5 and FFT headers; the libraries are AMReX's dependencies.
            depends_on("amrex+hdf5", when="+hdf5")
            depends_on("amrex+fft", when="+fft")
            # ERF's own build turns AMReX_EB on from 25.03, when ERF.H starts including
            # AMReX_EBFabFactory.H; before that its EB-only code does not build, and AMReX's
            # EB-only MLNodeLaplacian overload makes 25.01's call ambiguous.
            depends_on("amrex~eb", when="@24.10:25.02")
            depends_on("amrex+eb", when="@25.03:")

            # ERF's find_package(AMReX) always asks for the PIC, PARTICLES and LSOLVERS
            # components, even ~particles, and for 3D only.
            depends_on("amrex dimensions=3 +pic+particles+linear_solvers")

            depends_on("amrex precision=double", when="@:26.05")
            for _p in ("double", "single"):
                depends_on(f"amrex precision={_p}", when=f"@26.06:precision={_p}")
            # Spack's builtin amrex has no particles_precision variant yet; packages/amrex in
            # this repository adds it until spack-packages does.
            for _p in ("double", "single"):
                depends_on(
                    f"amrex particles_precision={_p}", when=f"@26.06:particles_precision={_p}"
                )
            # The tests and +fcompare need <amrex prefix>/bin/amrex_fcompare, which only
            # +plotfile_tools installs.
            depends_on("amrex+plotfile_tools", when="+tests")
            depends_on("amrex+plotfile_tools", when="+fcompare")
        # The releases that build their vendored AMReX need its HDF5 themselves.
        depends_on("hdf5@1.10.4:+mpi", when="@:24.09+hdf5+mpi")
        depends_on("hdf5@1.10.4:~mpi", when="@:24.09+hdf5~mpi")

        depends_on("mpi", when="+mpi")
        # ERF opens NetCDF files with nc_create_par in NC_NETCDF4 mode: HDF5-backed parallel I/O.
        depends_on("netcdf-c+mpi", when="+netcdf")
        depends_on("netcdf-fortran", when="+noahmp")
        # EKAT's bundled Kokkos imports rocThrust (HIP) and oneDPL (SYCL), as Spack's kokkos does.
        for _ekat in ("+rrtmgp", "+eamxx_shoc"):
            depends_on("rocthrust", when=f"{_ekat}+rocm")
            depends_on("intel-oneapi-dpl", when=f"{_ekat}+sycl")

    # Spack cuda_arch/amdgpu_target -> Kokkos_ARCH_<name> for the Kokkos EKAT bundles, one
    # map per EKAT generation (the ranges of _rrtmgp_stacks).
    _kokkos_cuda_arch = {
        "50": "MAXWELL50",
        "52": "MAXWELL52",
        "53": "MAXWELL53",
        "60": "PASCAL60",
        "61": "PASCAL61",
        "70": "VOLTA70",
        "72": "VOLTA72",
        "75": "TURING75",
        "80": "AMPERE80",
        "86": "AMPERE86",
        "89": "ADA89",
        "90": "HOPPER90",
    }
    _kokkos_amdgpu_arch = {
        "gfx906": "VEGA906",
        "gfx908": "VEGA908",
        "gfx90a": "VEGA90A",
        "gfx940": "AMD_GFX940",
        "gfx942": "AMD_GFX942",
        "gfx1030": "NAVI1030",
        "gfx1100": "AMD_GFX1100",
        "gfx1103": "AMD_GFX1103",
    }
    _kokkos_arch = {
        "@25.10:26.09": (_kokkos_cuda_arch, _kokkos_amdgpu_arch),
        "@26.10:": (
            {
                **_kokkos_cuda_arch,
                "87": "AMPERE87",
                "100": "BLACKWELL100",
                "103": "BLACKWELL103",
                "120": "BLACKWELL120",
                "121": "BLACKWELL121",
            },
            _kokkos_amdgpu_arch,
        ),
    }

    conflicts("+openmp", when="+cuda", msg="ERF supports one of OpenMP, CUDA, HIP and SYCL")
    conflicts("+openmp", when="+rocm", msg="ERF supports one of OpenMP, CUDA, HIP and SYCL")
    conflicts("+openmp", when="+sycl", msg="ERF supports one of OpenMP, CUDA, HIP and SYCL")
    conflicts("+cuda", when="+rocm", msg="ERF supports one of OpenMP, CUDA, HIP and SYCL")
    conflicts("+cuda", when="+sycl", msg="ERF supports one of OpenMP, CUDA, HIP and SYCL")
    conflicts("+rocm", when="+sycl", msg="ERF supports one of OpenMP, CUDA, HIP and SYCL")
    # CudaPackage allows cuda_arch=none; cmake_args needs a real one.
    conflicts("cuda_arch=none", when="+cuda", msg="ERF needs an explicit cuda_arch")
    # Apple Clang ships no OpenMP runtime, so AMReX's CMake cannot find OpenMP with it.
    conflicts(
        "+openmp",
        when="platform=darwin %c,cxx=apple-clang",
        msg="OpenMP on macOS needs gcc for C and C++ everywhere (%%c,cxx=gcc); Apple Clang has none",
    )
    requires("%c,cxx=oneapi", when="+sycl", msg="the SYCL backend needs the oneAPI compilers")
    conflicts("+netcdf", when="~mpi", msg="ERF's NetCDF I/O is parallel and needs MPI")
    conflicts("+noahmp", when="~netcdf", msg="Noah-MP requires NetCDF")
    conflicts("+rrtmgp", when="~netcdf", msg="RRTMGP requires NetCDF")
    conflicts("+rrtmgp", when="~mpi", msg="RRTMGP (EKAT) requires MPI")
    conflicts("+eamxx_shoc", when="~mpi", msg="SHOC (EKAT) requires MPI")
    conflicts("+eamxx_shoc", when="precision=single", msg="EAMxx SHOC is double precision only")
    conflicts("+tools", when="~netcdf", msg="The tools require NetCDF")
    # Backport the CUDA source-language fix to the current release. Older
    # releases still compile the tools as C++, which CUDA AMReX headers reject.
    patch("cuda-tools.patch", when="@26.09+cuda+tools")
    conflicts("+tools", when="@:26.06+cuda", msg="CUDA plotfile tools require ERF 26.09 or newer")
    conflicts(
        "+ml_uphys_diagnostics",
        when="~particles",
        msg="the ML microphysics diagnostics live in the particle (SDM) microphysics",
    )
    patch("ml-diagnostics-single.patch", when="@26.09+ml_uphys_diagnostics precision=single")
    conflicts(
        "+wsm6_fort",
        when="precision=single",
        msg="the WSM6 Fortran bridge is double precision only",
    )
    conflicts(
        "+wdm6_fort",
        when="precision=single",
        msg="the WDM6 Fortran bridge is double precision only",
    )
    # EKAT's Kokkos takes exactly one GPU architecture, as Kokkos_ARCH_<name> (see cmake_args);
    # it ignores CMAKE_CUDA_ARCHITECTURES and otherwise autodetects, which fails on a node
    # without the GPU.
    for _range, (_cuda_map, _amd_map) in _kokkos_arch.items():
        for _ekat in ("+rrtmgp", "+eamxx_shoc"):
            for _sm in CudaPackage.cuda_arch_values:
                if _sm not in _cuda_map:
                    conflicts(
                        f"cuda_arch={_sm}",
                        when=f"{_range}{_ekat}+cuda",
                        msg="EKAT's bundled Kokkos has no architecture for this cuda_arch",
                    )
            for _tgt in ROCmPackage.amdgpu_targets:
                if _tgt not in _amd_map:
                    conflicts(
                        f"amdgpu_target={_tgt}",
                        when=f"{_range}{_ekat}+rocm",
                        msg="EKAT's bundled Kokkos has no architecture for this amdgpu_target",
                    )
            for _a, _b in itertools.combinations(sorted(_cuda_map), 2):
                conflicts(
                    f"cuda_arch={_a},{_b}",
                    when=f"{_range}{_ekat}+cuda",
                    msg="Kokkos takes one cuda_arch",
                )
            for _a, _b in itertools.combinations(sorted(_amd_map), 2):
                conflicts(
                    f"amdgpu_target={_a},{_b}",
                    when=f"{_range}{_ekat}+rocm",
                    msg="Kokkos takes one amdgpu_target",
                )

    # RRTMGP, EKAT and EKAT's third-party sources at the commits ERF and EKAT pin, one set
    # per EKAT generation; check_pins.py compares them with the gitlinks and FetchContent
    # tags. Spack's kokkos, spdlog and yaml-cpp packages are not substitutes: ERF
    # add_subdirectory()s all of these, forces EKAT's spdlog, passes in-tree Kokkos options,
    # and E3SM's Kokkos carries E3SM patches. They are resources rather than submodules
    # because EKAT lists its submodules over SSH, which an unattended fetch cannot follow;
    # patch() moves them onto the submodule paths.
    _rrtmgp_stacks = {
        "@25.10:26.09": (
            (
                "rrtmgp",
                "rte-rrtmgp",
                "71693e7490db0908a625eadb6ac3f23967de2c7d",
                "Submodules/RRTMGP",
            ),
            ("ekat", "EKAT", "10ab9d34ceba4eead542702a29676adf77ad91af", "Submodules/ekat"),
            (
                "kokkos",
                "kokkos",
                "6a19eb9454da490dbf785443d8f021028baeb161",
                "Submodules/ekat/extern/kokkos",
            ),
            (
                "yaml-cpp",
                "yaml-cpp",
                "95088a0a2b6f2dec0b3e6e59020cdcc0d4f3c658",
                "Submodules/ekat/extern/yaml-cpp",
            ),
            (
                "spdlog",
                "spdlog",
                "bdd1dff3788ebfe520f48f9ad216c60da6dd8f00",
                "Submodules/ekat/extern/spdlog",
            ),
            (
                "catch2",
                "Catch2",
                "958944d27a2d2fb82aa008377bf4f8752f6b848e",
                "Submodules/ekat/extern/Catch2",
            ),
        ),
        "@26.10:": (
            (
                "rrtmgp",
                "rte-rrtmgp",
                "71693e7490db0908a625eadb6ac3f23967de2c7d",
                "Submodules/RRTMGP",
            ),
            ("ekat", "EKAT", "ae5e7f8227943a25c66b1a1822973ec665d224e4", "Submodules/ekat"),
            (
                "kokkos",
                "kokkos",
                "3ec81abe1816109f6f62ac48cef41921f91a4d00",
                "Submodules/ekat/extern/kokkos",
            ),
            (
                "yaml-cpp",
                "yaml-cpp",
                "95088a0a2b6f2dec0b3e6e59020cdcc0d4f3c658",
                "Submodules/ekat/extern/yaml-cpp",
            ),
            (
                "spdlog",
                "spdlog",
                "79524ddd08a4ec981b7fea76afd08ee05f83755d",
                "Submodules/ekat/extern/spdlog",
            ),
            (
                "catch2",
                "Catch2",
                "0d70154ddf178cf99bca140669e0618267e51682",
                "Submodules/ekat/extern/Catch2",
            ),
        ),
    }
    # SHOC needs EKAT too; the disjoint conditions fetch each once for +rrtmgp+eamxx_shoc.
    for _range, _stack in _rrtmgp_stacks.items():
        for _name, _repo, _commit, _ in _stack:
            for _when in ("+rrtmgp", "~rrtmgp+eamxx_shoc"):
                resource(
                    name=_name,
                    git=f"https://github.com/E3SM-Project/{_repo}.git",
                    commit=_commit,
                    destination="spack-resources",
                    placement=_name,
                    when=f"{_range}{_when}",
                )

    # The EAMxx sources SHOC compiles, at the E3SM commit Build/GNU_Ekat/eamxx_clone.sh
    # checks out (check_pins.py compares them); Spack 1.1+ clones them sparsely.
    _e3sm_commit = "2eb52d9b3d452e593f3825f9408507cae3519bf6"
    resource(
        name="e3sm-eamxx",
        git="https://github.com/E3SM-Project/E3SM.git",
        commit=_e3sm_commit,
        git_sparse_paths=["components/eamxx/src/physics", "components/eamxx/src/share"],
        destination="spack-resources",
        placement="e3sm",
        when="+eamxx_shoc",
    )

    # Noah-MP's ERF driver (drivers/erf) is versioned with ERF, past the last Noah-MP tag,
    # and compiled in-tree with ERF's precision macro, so ERF builds its pinned commit rather
    # than Spack's noahmp package. Each pin is that release's Submodules/Noah-MP gitlink;
    # check_pins.py --fix keeps the develop one.
    _noahmp_pins = (
        ("25.12:26.04", "6dc16027de56f2abc9fe430016acd086625e9ae2"),
        ("26.05", "c08dd406761330ac823a3e268bd2e98eb6889b40"),
        ("26.06", "41587582bf9b73786ce693660b0a56ff412d1eff"),
        ("26.09", "e03243d2c5eac19b1c6a02a2c3ca12f03f92cceb"),
        ("develop", "fbaa38f0753a0f199dcebd508e2ef64a5dc74bfd"),
    )
    for _versions, _commit in _noahmp_pins:
        resource(
            name="noahmp",
            git="https://github.com/erf-model/noahmp.git",
            commit=_commit,
            destination="spack-resources",
            placement="noahmp",
            when=f"@{_versions}+noahmp",
        )

    @staticmethod
    def _move_into_place(src, dst):
        """Move a staged resource onto the submodule path ERF's CMake expects."""
        # An uninitialised submodule leaves an empty directory behind.
        if os.path.isdir(dst) and not os.listdir(dst):
            os.rmdir(dst)
        if os.path.exists(dst):
            raise InstallError(f"{dst} is already populated; refusing to replace it")
        mkdirp(os.path.dirname(dst))
        os.rename(src, dst)

    def patch(self):
        if self.spec.satisfies("+rrtmgp") or self.spec.satisfies("+eamxx_shoc"):
            stack = next(s for r, s in self._rrtmgp_stacks.items() if self.spec.satisfies(r))
            for name, _, _, target in stack:
                self._move_into_place(join_path("spack-resources", name), target)
        if self.spec.satisfies("+eamxx_shoc"):
            self._move_into_place(
                join_path("spack-resources", "e3sm"), join_path("external", "E3SM")
            )
        if self.spec.satisfies("+noahmp"):
            self._move_into_place(join_path("spack-resources", "noahmp"), "Submodules/Noah-MP")

        # The bundled googletest submodule is not fetched; point the unit tests at Spack's.
        # They already use the GTest::gtest target and include(GoogleTest).
        if self.spec.satisfies("@26.06:+tests+unit_tests"):
            unit_cmakelists = join_path("Tests", "Unit", "CMakeLists.txt")
            marker = (
                "add_subdirectory(${PROJECT_SOURCE_DIR}/Submodules/googletest\n"
                "                 ${PROJECT_BINARY_DIR}/Submodules/googletest\n"
                "                 EXCLUDE_FROM_ALL)"
            )
            with open(unit_cmakelists, encoding="utf-8") as f:
                contents = f.read()
            if marker not in contents:
                raise InstallError("ERF unit tests no longer match the googletest patch marker")
            with open(unit_cmakelists, "w", encoding="utf-8") as f:
                f.write(contents.replace(marker, "find_package(GTest REQUIRED)"))

    def setup_build_environment(self, env):
        if self.spec.satisfies("+rocm"):
            # Spack-built rocprim (through rocthrust) puts llvm-amdgpu among the link
            # dependencies, so the compiler wrapper adds its include directory, where gfortran
            # then opens flang's iso_c_binding.mod instead of its own. hipcc finds it itself.
            llvm_include = self.spec["llvm-amdgpu"].prefix.include
            for var in ("SPACK_INCLUDE_DIRS", "SPACK_STORE_INCLUDE_DIRS"):
                env.remove_path(var, llvm_include)

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define("ERF_USE_INTERNAL_AMREX", spec.satisfies("@:24.09")),
            self.define_from_variant("ERF_ENABLE_MPI", "mpi"),
            self.define_from_variant("ERF_ENABLE_OPENMP", "openmp"),
            self.define_from_variant("ERF_ENABLE_CUDA", "cuda"),
            self.define_from_variant("ERF_ENABLE_HIP", "rocm"),
            self.define_from_variant("ERF_ENABLE_SYCL", "sycl"),
            self.define_from_variant("ERF_ENABLE_NETCDF", "netcdf"),
            self.define_from_variant("ERF_ENABLE_HDF5", "hdf5"),
            self.define_from_variant("ERF_ENABLE_PARTICLES", "particles"),
            self.define_from_variant("ERF_ENABLE_MULTIBLOCK", "multiblock"),
            self.define_from_variant("ERF_ENABLE_FFT", "fft"),
            self.define_from_variant("ERF_ENABLE_TESTS", "tests"),
            # Before 24.10 the tests' amrex_fcompare comes from the vendored AMReX, built
            # only when ERF_ENABLE_FCOMPARE is on.
            self.define(
                "ERF_ENABLE_FCOMPARE",
                spec.satisfies("+fcompare") or spec.satisfies("@:24.09+tests"),
            ),
            # Conditional variants absent from this version define nothing.
            self.define_from_variant("ERF_ENABLE_TOOLS", "tools"),
            self.define_from_variant("ERF_ENABLE_RRTMGP", "rrtmgp"),
            self.define_from_variant("ERF_ENABLE_NOAHMP", "noahmp"),
            self.define_from_variant("ERF_ENABLE_WINDFARM", "windfarm"),
            self.define_from_variant("ERF_ENABLE_ML_UPHYS_DIAGNOSTICS", "ml_uphys_diagnostics"),
            self.define_from_variant("ERF_ENABLE_WSM6_FORT", "wsm6_fort"),
            self.define_from_variant("ERF_ENABLE_WDM6_FORT", "wdm6_fort"),
            self.define_from_variant("ERF_ENABLE_EAMXX_SHOC", "eamxx_shoc"),
            self.define_from_variant("ERF_ENABLE_FIRE", "fire"),
            self.define("ERF_PRECISION", spec.variants["precision"].value.upper()),
        ]

        if spec.satisfies("@26.06:"):
            args.append(self.define("ERF_ENABLE_UNIT_TESTS", spec.satisfies("+tests+unit_tests")))
            particles_precision = spec.variants["particles_precision"].value
            if particles_precision == "same":
                particles_precision = spec.variants["precision"].value
            args.append(self.define("ERF_PARTICLES_PRECISION", particles_precision.upper()))

        if spec.satisfies("@25.11:+cuda%nvhpc"):
            args.append(self.define("ERF_ENABLE_NVHPC", True))

        # --oversubscribe is Open MPI's spelling; only ERF's test harness reads MPIEXEC_PREFLAGS.
        if spec.satisfies("+tests+mpi ^[virtuals=mpi] openmpi"):
            args.append(self.define("MPIEXEC_PREFLAGS", "--oversubscribe"))

        # 25.12 and later derive ERF_ENABLE_EKAT from ERF_ENABLE_RRTMGP.
        if spec.satisfies("@25.10:25.11+rrtmgp"):
            args.append(self.define("ERF_ENABLE_EKAT", True))

        if spec.satisfies("+rrtmgp") or spec.satisfies("+eamxx_shoc"):
            cuda_map, amd_map = next(
                maps for rng, maps in self._kokkos_arch.items() if spec.satisfies(rng)
            )
            if spec.satisfies("+cuda"):
                (sm,) = spec.variants["cuda_arch"].value
                args.append(self.define("Kokkos_ARCH_" + cuda_map[sm], True))
            if spec.satisfies("+rocm"):
                (tgt,) = spec.variants["amdgpu_target"].value
                args.append(self.define("Kokkos_ARCH_" + amd_map[tgt], True))
            if spec.satisfies("@26.10:"):
                # This EKAT fetches yaml-cpp, spdlog and Catch2 with FetchContent and trusts a
                # pre-populated directory only with a .git in it; the override makes it
                # add_subdirectory() the resource as it is.
                extern = join_path(self.stage.source_path, "Submodules", "ekat", "extern")
                for name, var in (
                    ("yaml-cpp", "YAML_CPP"),
                    ("spdlog", "SPDLOG"),
                    ("Catch2", "CATCH2"),
                ):
                    args.append(
                        self.define(f"FETCHCONTENT_SOURCE_DIR_{var}", join_path(extern, name))
                    )
            if spec.satisfies("+sycl"):
                # Kokkos's FindTPLONEDPL falls back to find_package(oneDPL) when the compiler
                # does not see the headers; the config lives in the component.
                dpl = spec["intel-oneapi-dpl"].package.component_prefix
                args.append(self.define("oneDPL_DIR", dpl.lib.cmake.oneDPL))

        if spec.satisfies("+netcdf"):
            args.append(self.define("NETCDF_DIR", spec["netcdf-c"].prefix))

        if spec.satisfies("+noahmp"):
            # The Noah-MP driver's FindNetCDF looks for the Fortran interface under
            # netcdf-c's prefix only.
            netcdf_fortran = spec["netcdf-fortran"]
            args.extend(
                [
                    self.define("NETCDF_INCLUDES_F90", netcdf_fortran.prefix.include),
                    self.define("NETCDF_LIBRARIES_F90", netcdf_fortran.libs.joined(";")),
                ]
            )

        if spec.satisfies("+cuda"):
            # SASS for each architecture, PTX for the newest only (AMReX's "all" policy):
            # SASS+PTX for three architectures pushed erf_exec past the 2 GB small code model
            # at link time. Compressed fatbins cost nothing at load time.
            archs = sorted(spec.variants["cuda_arch"].value, key=self._cuda_arch_number)
            spelled = [a if a == archs[-1] else f"{a}-real" for a in archs]
            args.extend(
                [
                    self.define("CMAKE_CUDA_ARCHITECTURES", ";".join(spelled)),
                    self.define("CMAKE_CUDA_FLAGS", "-Xfatbin=-compress-all"),
                    self.define("CUDAToolkit_ROOT", spec["cuda"].prefix),
                ]
            )

        if spec.satisfies("+rocm"):
            args.append(self.define("CMAKE_CXX_COMPILER", spec["hip"].hipcc))
            if spec.satisfies("@:24.09"):
                targets = spec.variants["amdgpu_target"].value
                args.append(self.define("AMReX_AMD_ARCH", ";".join(targets)))

        return args

    def check(self):
        """Run CTest directly so Make's job server cannot override test concurrency."""
        if not self.spec.satisfies("+tests"):
            return
        ctest = Executable(join_path(self.spec["cmake"].prefix.bin, "ctest"))
        # Spack's default CMake check invokes make test and replaces
        # CTEST_PARALLEL_LEVEL with the build worker count. CTest 3.31.6
        # crashed while using that job server. Keep tests serial by default,
        # honor an explicit test worker count, and preserve failure output.
        with working_dir(self.build_directory):
            ctest(
                "--output-on-failure",
                "--no-tests=error",
                "--parallel",
                os.environ.get("CTEST_PARALLEL_LEVEL") or "1",
                extra_env={"MAKEFLAGS": "", "MFLAGS": ""},
            )

    @staticmethod
    def _cuda_arch_number(arch):
        """Numeric part of a compute capability ("90a" -> 90), for ordering."""
        return int("".join(c for c in arch if c.isdigit()))

    @run_after("install")
    def install_examples(self):
        """Ship two self-contained cases so an installed ERF can be run without a
        checkout: prefix/share/erf/examples/<case>/inputs_*."""
        examples = join_path(self.prefix.share, "erf", "examples")
        for case in ("ScalarAdvDiff", "IsentropicVortex"):
            # The cases sat under DryRegTests from 24.11 to 26.03.
            for group in ("RegTests", "DryRegTests"):
                src = join_path(self.stage.source_path, "Exec", group, case)
                if os.path.isdir(src):
                    install_tree(src, join_path(examples, case))
                    break

    def _erf_exe(self):
        """The executable that runs the examples: one per case through 26.02, one for the
        dry regression tests at 26.03, erf_exec from 26.04."""
        if self.spec.satisfies("@26.04:"):
            name = "erf_exec"
        elif self.spec.satisfies("@26.03"):
            name = "erf_dryregtests"
        else:
            name = "erf_scalar_advdiff"
        return Executable(join_path(self.prefix.bin, name))

    def test_help(self):
        """erf_exec prints its usage"""
        if not self.spec.satisfies("@26.04:"):
            raise SkipTest("before 26.04 the executables exit with an error without an input file")
        erf_exec = self._erf_exe()
        out = erf_exec("--help", output=str.split, error=str.split)
        assert "usage" in out.lower()

    def test_scalar_advection(self):
        """run two steps of the uniform-velocity scalar advection case"""
        inputs = join_path(
            self.prefix.share, "erf", "examples", "ScalarAdvDiff", "inputs_advect_uniformU"
        )
        if not os.path.isfile(inputs):
            raise SkipTest("no example inputs installed for this version")
        if self.spec.satisfies("+windfarm"):
            raise SkipTest("ERF +windfarm cannot run a case without a wind farm")
        erf_exec = self._erf_exe()
        with working_dir(self.test_suite.current_test_cache_dir, create=True):
            erf_exec(inputs, "max_step=2", "erf.plot_int_1=-1", "erf.check_int=-1")
