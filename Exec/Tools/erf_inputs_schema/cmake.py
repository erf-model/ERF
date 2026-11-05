"""Which ERF_USE_* macros a build defines: ERF's CMake option defaults, the build's own -D options over them, a set()
of an ERF option in a taken if() branch (the ERF_ENABLE_SHOC alias, ERF_ENABLE_EKAT), applied in file order with
CMakeLists.txt first, and the if() around each target_compile_definitions in ERF's CMake.

Mirrored from other testing infrastructure. The build options are an argument (a -D list or a CMakeCache.txt),
not a container recipe.
"""
import re

TRUE = {"ON", "TRUE", "YES", "Y", "1"}
FALSE = {"OFF", "FALSE", "NO", "N", "0", ""}


def as_bool(value):
    if isinstance(value, bool):
        return value
    up = str(value).upper()
    return True if up in TRUE else False if up in FALSE else value


def parse_d_options(text):
    """ERF_* options from a cmake command line or a CMakeCache.txt ('-DERF_X=ON' or 'ERF_X:BOOL=ON')."""
    out = {}
    for name, value in re.findall(r"-D(ERF_\w+)(?::\w+)?=(\S*)", text):
        out[name] = as_bool(value)
    for name, value in re.findall(r"^(ERF_\w+):\w+=(.*)$", text, re.M):
        out[name] = as_bool(value.strip())
    return dict(sorted(out.items()))


def _options(erf_root, overrides):
    """Every ERF option: a -D override, else its CMakeLists default, where a default of ${OTHER} follows OTHER (after
    OTHER's own override), as in option(ERF_ENABLE_HDF5 "..." ${ERF_ENABLE_NETCDF})."""
    opts = dict(re.findall(r'option\s*\(\s*(ERF_\w+)\s+"[^"]*"\s+([^)\s]+)\s*\)',
                           (erf_root / "CMakeLists.txt").read_text()))

    def value(name, seen=()):
        if isinstance(overrides.get(name), bool):
            return overrides[name]
        raw = opts.get(name, "OFF")
        ref = re.fullmatch(r"\$\{(\w+)\}", raw)
        if ref:
            return value(ref.group(1), seen + (name,)) if ref.group(1) not in seen else False
        return as_bool(raw) is True

    resolved = {name: value(name) for name in opts}
    resolved.update({k: v for k, v in overrides.items() if isinstance(v, bool)})
    return resolved


def _truth(expr, resolved):
    expr = expr.strip()
    if " OR " in expr:
        return any(_truth(e, resolved) for e in expr.split(" OR "))
    if " AND " in expr:
        return all(_truth(e, resolved) for e in expr.split(" AND "))
    if expr.startswith("NOT "):
        return not _truth(expr[4:], resolved)
    return resolved.get(expr, False)


def _scan_file(path, rel, resolved, features, sources):
    stack = []                                       # per open if(): [active, a branch of this chain was taken]
    for i, line in enumerate(path.read_text().split("\n"), 1):
        s = re.sub(r"\s+#.*$", "", line.strip())    # if(X) # note
        m = re.match(r"(if|elseif)\s*\((.*)\)\s*$", s)
        if m and m.group(1) == "if":
            on = _truth(m.group(2), resolved)
            stack.append([on, on])
        elif m and stack:
            on = not stack[-1][1] and _truth(m.group(2), resolved)
            stack[-1] = [on, stack[-1][1] or on]
        elif re.match(r"else\s*\(", s) and stack:
            stack[-1] = [not stack[-1][1], True]
        elif re.match(r"endif\s*\(", s) and stack:
            stack.pop()
        m_set = re.match(r"set\s*\(\s*(ERF_\w+)\s+(\w+)", s)   # set(ERF_X ON [CACHE BOOL "..." FORCE]) in a taken branch
        if m_set and all(active for active, _ in stack):
            resolved[m_set.group(1)] = as_bool(m_set.group(2)) is True
        for macro in re.findall(r"\b(ERF_USE_\w+)\b", s) if "compile_definitions" in s else []:
            on = all(active for active, _ in stack)
            if on or macro not in features:
                features[macro] = features.get(macro, False) or on
                sources[macro] = f"{rel}:{i}"


def cmake_features(erf_root, cmake_options):
    """(features, sources): ERF_USE_X -> on?, and the CMake line that defines it."""
    resolved = _options(erf_root, cmake_options)
    features, sources = {}, {}
    for path in sorted([erf_root / "CMakeLists.txt", *(erf_root / "CMake").glob("*.cmake")]):
        _scan_file(path, path.relative_to(erf_root).as_posix(), resolved, features, sources)
    return features, sources
