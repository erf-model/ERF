"""Find every literal ParmParse query in a C++ tree, with the type of the variable it fills.

Mirrored from other testing infrastructure. The changes: conditions carry
the line of their #ifdef/#ifndef (the evidence a build requirement cites), and nothing is pinned to one ERF commit.
"""
import re

from .canonical import canonical_value

ERF_ALIASES = {"pp_prefix", "m_pp_prefix", "a_pp_prefix", "face_pp_prefix"}      # all "erf" (ERF.H: pp_prefix {"erf"})
# Prefixes built at run time, keyed by (file, the ParmParse constructor's argument): a name read as a family (the
# name is checked separately), or a member whose value is fixed (m_name). schema.generate stops when an entry whose
# file exists matches nothing, so an edit of the constructor argument is caught rather than dropping the keys.
FAMILIES = {("Source/Initialization/ERF_InitBCs.cpp", "pp_text"): (["<face>", "erf.<face>"], "face"),
            ("Source/Prob/ERF_CloudChamber.H", "face.c_str()"): (["<face>"], "face"),
            ("Source/IO/ERF_StationSampler.cpp", 'pp_prefix + "." + name'): (["erf.<station>"], "station"),
            ("Source/ObsNudging/ERF_ObsNudging.cpp", '"erf.obs_nudging." + name'):
                (["erf.obs_nudging.<station>"], "station"),
            ("Source/Refinement/ERF_Tagging.cpp", "ref_prefix"): (["erf.<indicator>"], "indicator"),
            ("Source/Refinement/ERF_RefineBox.cpp", "ref_prefix"): (["erf.<indicator>"], "indicator"),
            ("Source/Particles/ERFPCInitializations.cpp", "m_name"):             # ERFPC and SuperDropletPC (derived)
                (["tracer_particles", "super_droplets_moisture"], None),
            ("Source/Particles/ERF_SuperDropletPCInitializations.cpp", "m_name"): (["super_droplets_moisture"], None),
            ("Source/Microphysics/SuperDropletsMoist/ERF_SuperDropletsMoistInit.cpp", "m_name"):
                (["super_droplets_moisture"], None)}
# The read methods of AMReX_ParmParse.H: METHODS take (key, variable, ...); SKIP_ONE take an int before the variable
# (the occurrence k of the kth reads, or the value count of the array parser reads).
METHODS = ("query_enum_case_insensitive", "get_enum_case_insensitive", "query_enum_sloppy", "get_enum_sloppy",
           "query_enum", "get_enum", "queryAddWithParser", "queryWithParser", "getWithParser", "queryAsDouble",
           "getAsDouble", "queryarrAdd", "queryarr", "getarr", "queryline", "getline", "querytable", "gettable",
           "queryAdd", "query", "get")
SKIP_ONE = ("queryarrWithParser", "getarrWithParser", "queryarrAsDouble", "getarrAsDouble", "queryktharr",
            "getktharr", "querykth", "getkth")
CPP_TYPES = {"int": "int", "long": "int", "unsigned": "int", "Long": "int", "amrex::Long": "int", "bool": "bool",
             "Real": "real", "amrex::Real": "real", "double": "real", "float": "real", "std::string": "string",
             "string": "string", "IntVect": "int_list", "amrex::IntVect": "int_list", "RealVect": "real_list",
             "amrex::RealVect": "real_list"}
MODEL_BUILT = re.compile(r"SetModel\s*<|make_unique\s*<|make_shared\s*<|\bnew\s+\w")
LIST_METHODS = ("queryarr", "getarr", "queryarrAdd", "get_indexed", "query_indexed", "queryarrWithParser",
                "getarrWithParser", "queryarrAsDouble", "getarrAsDouble", "queryktharr", "getktharr")


def strip_comments(text):
    """C++ comments removed, strings kept, line numbers preserved."""
    out, i, n, state = [], 0, len(text), None
    while i < n:
        c = text[i]
        if state is None:
            if text.startswith("//", i):
                j = text.find("\n", i)
                i = n if j < 0 else j
                continue
            if text.startswith("/*", i):
                j = text.find("*/", i + 2)
                j = n if j < 0 else j + 2
                out.append("\n" * text.count("\n", i, j))
                i = j
                continue
            if c in "\"'":
                state = c
        elif c == "\\":
            out.append(text[i:i + 2])
            i += 2
            continue
        elif c == state:
            state = None
        out.append(c)
        i += 1
    return "".join(out)


def _directive(stack, kind, rest, line_no):
    if kind in ("ifdef", "ifndef"):
        stack.append((rest.split()[0] if rest else None, kind == "ifdef", line_no))
    elif kind == "if":
        simple = re.fullmatch(r"(!)?\s*(?:defined\s*\(?\s*)?(\w+)\s*\)?", rest)
        stack.append((simple.group(2), not simple.group(1), line_no) if simple else (None, True, line_no))
    elif kind == "elif" and stack:
        stack[-1] = (None, True, line_no)
    elif kind == "else" and stack:
        macro, positive, _ = stack[-1]
        stack[-1] = (macro, not positive, line_no)
    elif kind == "endif" and stack:
        stack.pop()


def conditions(lines):
    """For each line, the stack of (macro, positive, directive line) conditions on ERF_USE_* macros."""
    stack, out = [], []
    for line_no, line in enumerate(lines, 1):
        m = re.match(r"#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)", line.strip())
        if m:
            _directive(stack, m.group(1), m.group(2).strip(), line_no)
        out.append([c for c in stack if c[0] and c[0].startswith("ERF_USE_")])
    return out


def enums(files):
    found = {}
    for _, text in files:
        for name, body in re.findall(r"AMREX_ENUM\(\s*(\w+)\s*,([^)]*)\)", text):
            found.setdefault(name, [v.strip() for v in body.split(",") if v.strip()])
    return found


def declared(text, enum_names):
    """(name, schema type, 1-based line) for every declaration in this text, every declarator of a line included."""
    scalar = "|".join(sorted({re.escape(t) for t in CPP_TYPES} | {re.escape(e) for e in enum_names}, key=len,
                             reverse=True))
    rx = re.compile(r"(?<![\w:])(?:(?:amrex::)?(?:Vector|Array|GpuArray)\s*<\s*(" + scalar + r")\s*(?:,\s*\w+\s*)?>"
                    r"|std::vector\s*<\s*(" + scalar + r")\s*>|std::array\s*<\s*(" + scalar + r")\s*,\s*\w+\s*>"
                    r"|(" + scalar + r"))\s+&?\s*(\w+)\s*(?=[;=,\[{()])")
    for m in rx.finditer(text):
        element = m.group(1) or m.group(2) or m.group(3)
        base = element or m.group(4)
        kind = "enum:" + base if base in enum_names else CPP_TYPES[base]
        if element and not kind.endswith("_list"):
            kind += "_list"
        line = text.count("\n", 0, m.start()) + 1
        for name in [m.group(5)] + _more_declarators(text, m.end()):
            yield name, kind, line


def declarations(text, enum_names):
    """name -> set of schema types declared in this text."""
    out = {}
    for name, kind, _ in declared(text, enum_names):
        out.setdefault(name, set()).add(kind)
    return out


def pick_kind(candidates, method, per_level, line):
    """One name's kind from its (kind, declaration line) candidates, by the owner's rule of 10-30: keep the kinds the
    read call takes (a list for the array reads, a scalar otherwise, any for a per-level read), then the nearest
    declaration before the read; otherwise unknown."""
    kinds = {k for k, _ in candidates}
    if len(kinds) == 1:
        return next(iter(kinds))
    if not per_level:
        want_list = method in LIST_METHODS
        candidates = [(k, ln) for k, ln in candidates if k.endswith("_list") == want_list]
        if len({k for k, _ in candidates}) == 1:
            return candidates[0][0]
    before = [(ln, k) for k, ln in candidates if ln <= line]
    if before:
        nearest = max(ln for ln, _ in before)
        at = {k for ln, k in before if ln == nearest}
        if len(at) == 1:
            return next(iter(at))
    return "unknown"


def _more_declarators(text, pos):
    """The other names of `T a = x, b = y;`: top-level commas up to the statement's `;` (none for a call)."""
    depth, names = 0, []
    for i in range(pos, min(len(text), pos + 400)):
        c = text[i]
        if c in "([{<":
            depth += 1
        elif c in ")]}>":
            depth -= 1
            if depth < 0:
                return []                                # inside an argument list, not a declaration
        elif c == ";" and depth == 0:
            return names
        elif c == "," and depth == 0:
            m = re.match(r"\s*&?\s*(\w+)\s*(?=[=;,\[{(])", text[i + 1:])
            if not m:
                return names
            names.append(m.group(1))
    return []


def target_name(expr):
    expr = re.sub(r"\[[^\]]*\]", "", expr.strip().lstrip("&"))
    m = re.match(r"\w+", re.split(r"\.|->", expr)[-1].strip())
    return m.group(0) if m else None


def resolve_prefix(expr):
    expr = (expr or "").strip()
    if not expr:
        return ""
    out = []
    for p in (p.strip() for p in expr.split("+")):
        lit = re.fullmatch(r'"([^"]*)"', p)
        if lit:
            out.append(lit.group(1))
        elif p in ERF_ALIASES:
            out.append("erf")
        else:
            return None
    return "".join(out)


def literal_keys(expr):
    """Keys an expression can name: Concatenate("lit", n, w) or "lit" + std::to_string(n) is the family lit<n>;
    otherwise every string literal in it (e.g. cond ? "long" : "lon")."""
    m = re.search(r'Concatenate\s*\(\s*"([^"]+)"|"([^"]+)"\s*\+\s*std::to_string', expr)
    if m:
        return [(m.group(1) or m.group(2)) + "<n>"]
    return re.findall(r'"([^"]+)"', expr)


DECL = r"(?:std::string|char\s*\*|char\s+const\s*\*|auto)\s*&?\s*"


def identifier_keys(name, text, all_files):
    """Literal keys an identifier holds: its assignments in this file, or, for a function parameter, the literal
    (or literal-holding variable) passed at each call site of that function."""
    out = []
    for expr in re.findall(r"(?:const\s+)?" + DECL + name + r"\s*=(?!=)\s*([^;]+);", text):
        out += literal_keys(expr)
    if out:
        return sorted(set(out)) if len(set(out)) <= 4 else []
    sig = re.search(r"(\w+)\s*\(\s*const\s+std::string\s*&\s*" + name + r"\b", text)
    if not sig:
        return []
    for _, other in all_files:
        for arg in re.findall(r"\b" + sig.group(1) + r'\s*\(\s*("[^"]*"|\w+)\s*,', other):
            if arg.startswith('"'):
                out.append(arg.strip('"'))
                continue
            for expr in re.findall(r"(?:const\s+)?" + DECL + arg + r"\s*=(?!=)\s*([^;]+);", other):
                out += literal_keys(expr)
    return sorted(set(out)) if len(set(out)) <= 4 else []


class FileScan:
    """One source file's queries: (key, fact, conditions, method, target) hits and unresolved ParmParse reads."""

    def __init__(self, path, rel, enum_names, global_decls, erf_side, all_files):
        self.rel, self.enum_names, self.global_decls = rel, enum_names, global_decls
        self.erf_side, self.all_files = erf_side, all_files
        self.text = strip_comments(path.read_text(errors="replace"))
        self.lines = self.text.split("\n")
        self.conds = conditions(self.lines)
        self.local = {}                                    # name -> [(kind, declaration line)]
        for name, kind, line in declared(self.text, enum_names):
            self.local.setdefault(name, []).append((kind, line))
        header = path.with_suffix(".H")
        self.header_rel = self.rel[:-len(path.suffix)] + ".H"
        if path.suffix == ".cpp" and header.exists():          # a .cpp's own header is local scope too
            self.header_text = strip_comments(header.read_text(errors="replace"))
            for name, kind, _ in declared(self.header_text, enum_names):
                self.local.setdefault(name, []).append((kind, 0))     # declared before any line of the .cpp
        else:
            self.header_text = ""
        self.offsets = [0]
        for line in self.lines:
            self.offsets.append(self.offsets[-1] + len(line) + 1)
        self.decls = self._pp_decls()
        self.families_hit = set()                          # FAMILIES keys a declaration here matched

    def line_of(self, pos):
        return next(i for i in range(len(self.offsets) - 1) if self.offsets[i + 1] > pos)

    def _pp_decls(self):
        decls = []                                         # (line, var, prefix expression)
        for m in re.finditer(r"(?:amrex::)?ParmParse\s+(\w+)\s*(?:\(\s*(.*?)\s*\)|\{\s*(.*?)\s*\})?\s*;", self.text):
            decls.append((self.line_of(m.start()), m.group(1), m.group(2) if m.group(2) is not None else m.group(3)))
        for m in re.finditer(r"\b(\w+)\s*=\s*(?:amrex::)?ParmParse\s*\(\s*(.*?)\s*\)\s*;", self.text):
            decls.append((self.line_of(m.start()), m.group(1), m.group(2)))
        return sorted(decls)

    def prefixes(self, var, line):
        expr = next((e for ln, v, e in reversed(self.decls) if v == var and ln <= line), "__none__")
        if expr == "__none__":
            return self.parameter_prefixes(var)
        return self._prefixes_of(self.rel, expr)

    def _prefixes_of(self, rel, expr):
        key = (rel, (expr or "").strip())
        fam = FAMILIES.get(key)
        if fam:
            if rel == self.rel:
                self.families_hit.add(key)
            return fam
        p = resolve_prefix(expr)
        return ([p], None) if p is not None else (None, None)

    def parameter_prefixes(self, var):
        """A `ParmParse& var` parameter of a function here takes the prefixes of the ParmParse objects passed at the
        function's call sites in any file (ERF_WallScalarBC.H's parse_wall_face_scalars, called from ERF_InitBCs.cpp
        with the face object). None when var is not such a parameter or no call site resolves."""
        sig = re.search(r"(\w+)\s*\([^;{)]*?ParmParse\s*&\s*" + var + r"\b", self.text)
        if not sig:
            return None, None
        prefixes, family = [], None
        for rel, other in self.all_files:
            for call in re.finditer(r"\b" + sig.group(1) + r"\s*\(\s*(\w+)\s*[,)]", other):
                decl = None
                for m in re.finditer(r"(?:amrex::)?ParmParse\s+" + call.group(1) + r"\s*\(\s*(.*?)\s*\)\s*;|"
                                     r"\b" + call.group(1) + r"\s*=\s*(?:amrex::)?ParmParse\s*\(\s*(.*?)\s*\)\s*;",
                                     other):
                    if m.start() < call.start():
                        decl = m.group(1) if m.group(1) is not None else m.group(2)
                if decl is None:
                    continue
                pfx, fam = self._prefixes_of(rel, decl)
                if pfx:
                    prefixes += [p for p in pfx if p not in prefixes]
                    family = family or fam
        return (prefixes, family) if prefixes else (None, None)

    def kind_of(self, name, method="query", line=0, per_level=False):
        """The kind of the variable a read fills: its declarations in this file (pick_kind), else the one kind it has
        across ERF, else unknown."""
        if name in self.local:
            return pick_kind(self.local[name], method, per_level, line)
        kinds = self.global_decls.get(name, set())
        return next(iter(kinds)) if len(kinds) == 1 else "unknown"

    def names(self, arg):
        if arg.startswith('"'):
            return re.findall(r'"([^"]+)"', arg)
        return identifier_keys(arg.split(".")[0], self.text, self.all_files)

    def calls(self):
        """(position, ParmParse variable, method, key, target expression, per_level, extra) for every read."""
        key_rx = r'("[^"]+"(?:\s*,\s*"[^"]+")?|\w+(?:\.c_str\(\))?)'
        out = []
        for methods, skip in ((METHODS, ""), (SKIP_ONE, r"[^,()]+,\s*")):
            for m in re.finditer(r"\b(\w+)\s*\.\s*(" + "|".join(methods) + r")\s*\(\s*" + key_rx +
                                 r"\s*,\s*" + skip + r"([^,)]+)\s*([,)])", self.text):
                method = m.group(2) + ("_indexed" if m.group(5) == "," and m.group(2) in ("get", "query") else "")
                extra = {}
                if "sloppy" in method:                   # query_enum_sloppy(key, var, "_"): the ignored characters
                    ignores = re.match(r'\s*"([^"]*)"', self.text[m.end():])
                    extra["ignores"] = ignores.group(1) if ignores else None
                out += [(m.start(), m.group(1), method, key, m.group(4), False, extra) for key in self.names(m.group(3))]
        for m in re.finditer(r"\b(\w+)\s*\.\s*(contains|countval)\s*\(\s*" + key_rx + r"\s*\)", self.text):
            out += [(m.start(), m.group(1), m.group(2), key, None, False, {}) for key in self.names(m.group(3))]
        for m in re.finditer(r'\bquery_one_or_per_level(\w*)\s*\(\s*(\w+)\s*,\s*"([^"]+)"\s*,\s*([^,)]+)', self.text):
            method = "query_enum_case_insensitive" if "enum_case_insensitive" in m.group(1) else "query"
            out.append((m.start(), m.group(2), method, m.group(3), m.group(4), True, {}))
        if self.erf_side:
            for m in re.finditer(r'\bread_int_string\s*\(\s*\w+\s*,\s*"([^"]+)"', self.text):
                out.append((m.start(), "__erf__", "read_int_string", m.group(1), None, True, {}))
        return out

    def adds(self):
        """(key, raw literal, source) for each pp.add("k", literal) in ERF's add_par, which main.cpp passes to
        amrex::Initialize to run before the inputs are read: the inputs override it, so it is the default."""
        start = re.search(r"\bvoid\s+add_par\s*\(\s*\)\s*", self.text)
        if not self.erf_side or not start:
            return []
        body, at = _if_body(self.text, start.end())
        out = []
        for m in re.finditer(r'\b(\w+)\s*\.\s*add\s*\(\s*"([^"]+)"\s*,\s*(-?[\w.+-]+|"[^"]*")\s*\)\s*;', body):
            line = self.line_of(at + m.start())
            pfx, _ = self.prefixes(m.group(1), line)
            for p in pfx or []:
                out.append((f"{p}.{m.group(2)}" if p else m.group(2), m.group(3), f"{self.rel}:{line + 1}"))
        return out

    def kind(self, method, target, per_level, line=0):
        if method == "read_int_string":
            return "int"
        if method in ("contains", "countval", "querytable", "gettable"):
            return "unknown"           # a presence read (a typed read of the same key wins), or a {{...} {...}} table
        kind = self.kind_of(target_name(target), method, line, per_level)
        indexed = "[" in target                              # fixed_dt[0]: one element of a Vector
        if (per_level or indexed) and kind.endswith("_list") and method not in LIST_METHODS:
            kind = kind[:-5]
        if method in LIST_METHODS and not kind.endswith("_list"):
            kind += "_list"                                  # unknown_list: several words, each of unknown type
        return kind

    def default(self, target, kind, line, fact):
        """(value, origin, source): the variable's initializer as a canonical literal. Origins: initializer, none,
        expression (an initializer that is not a plain literal of the type)."""
        name = target_name(target) if target else None
        decl = max((ln for _, ln in self.local.get(name, []) if ln <= line + 1), default=None) if name else None
        found = _initializer(name, self.text, line) if name and decl != 0 else None
        if found is not None and decl is not None and found[1] != decl:
            found = None                             # that initializer belongs to an earlier, shadowed declaration
        rel = self.rel
        if found is None and name and self.header_text and decl in (None, 0):
            found, rel = _initializer(name, self.header_text, None), self.header_rel
        if found is None or kind.endswith("_list") or kind == "unknown" or not found[0]:
            return None, "none", None
        if kind == "unknown_list":
            return None, "none", None
        raw, at = found
        if kind == "string" and not re.fullmatch(r'"[^"]*"', raw):
            return None, "expression", f"{rel}:{at}"     # std::string p = prefix; is not a literal
        if kind.startswith("enum:"):                     # values attach later; schema._enum_facts canonicalizes
            return _literal(raw), "initializer", f"{rel}:{at}"
        try:
            value = canonical_value(kind, _literal(raw), fact)
        except ValueError:
            return None, "expression", f"{rel}:{at}"
        return value, "initializer", f"{rel}:{at}"

    def hits(self):
        hits, unresolved = [], []
        for pos, var, method, key, target, per_level, extra in self.calls():
            line = self.line_of(pos)
            pfx, family = (["erf"], None) if var == "__erf__" else self.prefixes(var, line)
            if pfx is None:
                is_pp = any(v == var for _, v, _ in self.decls) or re.search(r"ParmParse\s*&\s*" + var + r"\b",
                                                                              self.text)
                if is_pp:                                    # only ParmParse objects; not maps, NetCDF handles, etc.
                    unresolved.append(f"{self.rel}:{line + 1} {var}.{key}")
                continue
            fact = {"type": self.kind(method, target, per_level, line + 1), "source": f"{self.rel}:{line + 1}"}
            if per_level:
                fact["per_level"] = True
            if "enum" in method:
                fact["case_insensitive"] = "case_insensitive" in method or "sloppy" in method
            if "ignores" in extra:
                fact["ignores"] = extra["ignores"]
            if "WithParser" in method or "AsDouble" in method:
                fact["parser"] = True                        # an expression is meant; AsDouble also rounds to an int
                if "AsDouble" in method:
                    fact["as_double"] = True
            if method in ("queryline", "getline"):
                fact["line"] = True                          # the whole line, any number of words
            if method in ("querytable", "gettable"):
                fact["table"] = True
            if family:
                fact["family"] = family
            fact["required"] = method.startswith("get")
            value, origin, src = self.default(target, fact["type"], line, fact)
            fact.update(default=value, default_origin=origin)
            if src:
                fact["default_source"] = src
            retired = self.retired.get((var, key, line))
            if retired:
                fact = {"type": "removed", "source": fact["source"], "removed_shape": retired[0],
                        "removed_evidence": f"{self.rel}:{retired[1]}", "required": False, "default": None,
                        "default_origin": "none"}
            for p in pfx:
                hits.append((f"{p}.{key}" if p else key, fact, self.conds[line], method, target))
        return hits, unresolved

    @property
    def retired(self):
        """(var, key, 0-based line) -> (shape, abort line) for ERF's three ways of retiring a key: the read sits in an
        if() whose body calls Abort or Error.  contains(k) -> Abort; query(k, x) > 0 -> Error; countval(k) > 0 -> Error."""
        if not hasattr(self, "_retired"):
            self._retired = {}
            for shape, rx in RETIRED_SHAPES:
                for m in rx.finditer(self.text):
                    if _body_reads(self.text, m.end(), m.group(2)):
                        continue                         # read and mapped (AMReX's nosub): accepted, not removed
                    line = _abort_line(self.text, m.end())
                    if line:
                        self._retired[(m.group(1), m.group(2), self.line_of(m.start()))] = (shape, line)
            for m in RETIRED_VIA_VARIABLE.finditer(self.text):
                line = None if _body_reads(self.text, m.end(), m.group(4)) else _abort_line(self.text, m.end())
                if line:
                    shape = f"{m.group(3)}_var_gt0_error"
                    self._retired[(m.group(2), m.group(4), self.line_of(m.start()))] = (shape, line)
        return self._retired


RETIRED_SHAPES = (
    ("contains_abort", re.compile(r'\bif\s*\(\s*(\w+)\s*\.\s*contains\s*\(\s*"([^"]+)"\s*\)\s*\)')),
    ("query_gt0_error", re.compile(r'\bif\s*\(\s*(\w+)\s*\.\s*query\s*\(\s*"([^"]+)"\s*,[^()]*\)\s*>\s*0\s*\)')),
    ("countval_gt0_error", re.compile(r'\bif\s*\(\s*(\w+)\s*\.\s*countval\s*\(\s*"([^"]+)"\s*\)\s*>\s*0\s*\)')),
)
# int n = pp.countval(k); if (n > 0) Error(...)   (ERF_DataStruct.H's no_substepping)
RETIRED_VIA_VARIABLE = re.compile(r'\b(?:int|long|auto)\s+(\w+)\s*=\s*(\w+)\s*\.\s*(countval|query)\s*\(\s*"([^"]+)"'
                                  r'[^;]*;\s*if\s*\(\s*\1\s*>\s*0\s*\)')
ABORT = re.compile(r"\b(?:amrex::)?(?:Abort|Error)\s*\(")


def _body_reads(text, pos, key):
    """Does the if() body at pos read key with a typed query (so the key is mapped, not refused)?"""
    body, _ = _if_body(text, pos)
    return bool(re.search(r"\.\s*(?:" + "|".join(METHODS + SKIP_ONE) + r')\s*\(\s*"' + re.escape(key) + '"', body))


def _if_body(text, pos):
    """The statement an if() controls, starting at pos: a braced block, or one statement up to its semicolon."""
    i = pos + len(text[pos:]) - len(text[pos:].lstrip())
    if text[i:i + 1] != "{":
        end = text.find(";", i)
        return text[i:end if end >= 0 else len(text)], i
    depth = 0
    for j in range(i, len(text)):
        depth += {"{": 1, "}": -1}.get(text[j], 0)
        if depth == 0:
            return text[i:j + 1], i
    return text[i:], i


def _abort_line(text, pos):
    """The 1-based line of an Abort/Error the if() body always calls when entered, or None. An Abort nested in an
    inner if/else (a size or value check, as in ERF_InitZLevels.cpp) is not a removal."""
    body, start = _if_body(text, pos)
    top = 1 if body.startswith("{") else 0
    for m in ABORT.finditer(body):
        before = body[:m.start()]
        depth = before.count("{") - before.count("}")
        statement = re.split(r"[;{}]", before)[-1]
        if depth == top and not re.search(r"\b(?:if|else)\b", statement):
            return text.count("\n", 0, start + m.start()) + 1
    return None


def _in_string(text, pos):
    """Is pos inside a "..." literal on its line? (strip_comments keeps strings, so an Error("Set x = y instead.")
    text can look like an initializer of x.)"""
    line = text[text.rfind("\n", 0, pos) + 1:pos]
    return len(re.findall(r'(?<!\\)"', line)) % 2 == 1


def _initializer(name, text, before_line):
    """(raw initializer, 1-based line) of the last declaration of name before before_line (any line if None)."""
    rx = re.compile(r"(?:[\w>]\s+|,\s*)&?\s*" + re.escape(name) +
                    r"\s*(?:=(?!=)\s*([^;,{}]+?)|\{\s*([^{}]*?)\s*\}|\(\s*([^()]*?)\s*\))\s*[;,]")   # not ==
    best = None
    for m in rx.finditer(text):
        if _in_string(text, m.start() + 1):
            continue
        at = text.count("\n", 0, m.start()) + 1
        if before_line is None or at <= before_line + 1:
            best = (next(g for g in m.groups() if g is not None).strip(), at)
    return best


def _literal(raw):
    """A C++ literal as a ParmParse word: enum Name::Value -> Value, "text" -> text, 0.8_rt -> 0.8."""
    raw = raw.strip()
    m = re.fullmatch(r"(?:\w+::)*\w+::(\w+)", raw)
    if m and not re.fullmatch(r"[0-9.eE+-]+", raw):
        return m.group(1)
    m = re.fullmatch(r'"([^"]*)"', raw)
    if m:
        return m.group(1)
    return re.sub(r"(?<=[0-9.])(?:_rt|_prt|[fFlL])$", "", raw)


def posix_rel(path, root):
    """path relative to root with / separators on every platform (the schema's paths are byte-identical anywhere)."""
    return path.relative_to(root).as_posix()


def scan(root, rel_base, enum_names, global_decls, erf_side, all_files=()):
    """Every literal query in the tree: (key, fact, conditions, method, target), unresolved reads, the FAMILIES
    entries some declaration matched, and add_par's defaults (FileScan.adds)."""
    hits, unresolved, families_hit, adds = [], [], set(), []
    for path in sorted(p for p in root.rglob("*") if p.suffix in (".H", ".cpp")):
        rel = f"{rel_base}/{posix_rel(path, root)}"
        fs = FileScan(path, rel, enum_names, global_decls, erf_side, all_files)
        h, u = fs.hits()
        hits += h
        unresolved += u
        families_hit |= fs.families_hit
        adds += fs.adds()
    return hits, sorted(set(unresolved)), families_hit, adds


def value_requirements(files, enum_names):
    """Enum values that need a feature, with the line as evidence. Enforced: the value is refused at run time inside
    #ifndef X (RRTMGP), or its model is constructed only inside #ifdef X (Noah-MP: SetModel<NOAHMP>). Candidates, for
    review and not enforced: any other comparison inside #ifdef X."""
    out, candidates = {}, {}
    for rel, text in files:
        lines = text.split("\n")
        conds = conditions(lines)
        for i, line in enumerate(lines):
            if not conds[i] or not re.search(r"==|\bcase\b", line):
                continue
            macro, positive, _ = conds[i][-1]
            for enum, value in re.findall(r"\b(\w+)::(\w+)\b", line):
                if enum not in enum_names or value not in enum_names[enum]:
                    continue
                after = "\n".join(lines[i:i + 4])
                refused = not positive and re.search(r"\b(Error|Abort)\s*\(", after)
                built = positive and MODEL_BUILT.search(after)
                fact = {"feature": macro, "evidence": f"{rel}:{i + 1}"}
                if refused or built:
                    out.setdefault((enum, value), fact)
                elif positive:
                    candidates.setdefault((enum, value), fact)
    return out, {k: v for k, v in candidates.items() if k not in out}
