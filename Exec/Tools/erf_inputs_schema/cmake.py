"""Which ERF_USE_* macros a build defines: ERF's CMake option defaults, the build's own -D options over them, and the
if() around each target_compile_definitions in ERF's CMake.

Ported from the ASPIRE harness's gen_erf_schema.py. The build options are an argument (a -D list or a CMakeCache.txt),
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


def _option_defaults(erf_root):
    opts = dict(re.findall(r'option\s*\(\s*(ERF_\w+)\s+"[^"]*"\s+([^)\s]+)\s*\)',
                           (erf_root / "CMakeLists.txt").read_text()))
    resolved = {}
    for name, value in opts.items():
        ref = re.fullmatch(r"\$\{(\w+)\}", value)
        resolved[name] = (opts.get(ref.group(1), "OFF") if ref else value).upper() == "ON"
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
    stack = []
    for i, line in enumerate(path.read_text().split("\n"), 1):
        s = line.strip()
        m = re.match(r"(if|elseif)\s*\((.*)\)\s*$", s)
        if m and m.group(1) == "if":
            stack.append(_truth(m.group(2), resolved))
        elif m:
            stack[-1] = (not any(stack[-1:]) and _truth(m.group(2), resolved)) if stack else False
        elif re.match(r"else\s*\(", s) and stack:
            stack[-1] = not stack[-1]
        elif re.match(r"endif\s*\(", s) and stack:
            stack.pop()
        for macro in re.findall(r"\b(ERF_USE_\w+)\b", s) if "compile_definitions" in s else []:
            on = all(stack)
            if on or macro not in features:
                features[macro] = features.get(macro, False) or on
                sources[macro] = f"{rel}:{i}"


def cmake_features(erf_root, cmake_options):
    """(features, sources): ERF_USE_X -> on?, and the CMake line that defines it."""
    resolved = _option_defaults(erf_root)
    resolved.update({k: v for k, v in cmake_options.items() if isinstance(v, bool)})
    features, sources = {}, {}
    for path in sorted([erf_root / "CMakeLists.txt", *(erf_root / "CMake").glob("*.cmake")]):
        _scan_file(path, str(path.relative_to(erf_root)), resolved, features, sources)
    return features, sources
