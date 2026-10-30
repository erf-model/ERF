"""Find every literal ParmParse query in a C++ tree, with the type of the variable it fills.

Ported from the ASPIRE harness's compute/tools/gen_erf_schema.py (Build B's generator). The changes: conditions carry
the line of their #ifdef/#ifndef (the evidence a build requirement cites), and nothing is pinned to one ERF commit.
"""
import re

from .canonical import canonical_value

ERF_ALIASES = {"pp_prefix", "m_pp_prefix", "a_pp_prefix", "face_pp_prefix"}      # all "erf" (ERF.H: pp_prefix {"erf"})
# Prefixes built at run time from a name, read as families (the name is checked separately).
FAMILIES = {("Source/Initialization/ERF_InitBCs.cpp", "pp_text"): (["<face>", "erf.<face>"], "face"),
            ("Source/IO/ERF_StationSampler.cpp", 'pp_prefix + "." + name'): (["erf.<station>"], "station")}
METHODS = ("query_enum_case_insensitive", "get_enum_case_insensitive", "query_enum", "get_enum", "queryarrAdd",
           "queryarr", "getarr", "queryAdd", "query", "get")
CPP_TYPES = {"int": "int", "long": "int", "unsigned": "int", "Long": "int", "amrex::Long": "int", "bool": "bool",
             "Real": "real", "amrex::Real": "real", "double": "real", "float": "real", "std::string": "string",
             "string": "string", "IntVect": "int_list", "amrex::IntVect": "int_list", "RealVect": "real_list",
             "amrex::RealVect": "real_list"}
MODEL_BUILT = re.compile(r"SetModel\s*<|make_unique\s*<|make_shared\s*<|\bnew\s+\w")
LIST_METHODS = ("queryarr", "getarr", "queryarrAdd", "get_indexed", "query_indexed")


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


def declarations(text, enum_names):
    """name -> set of schema types declared in this text."""
    scalar = "|".join(sorted({re.escape(t) for t in CPP_TYPES} | {re.escape(e) for e in enum_names}, key=len,
                             reverse=True))
    rx = re.compile(r"(?<![\w:])(?:(?:amrex::)?(?:Vector|Array|GpuArray)\s*<\s*(" + scalar + r")\s*(?:,\s*\w+\s*)?>"
                    r"|std::vector\s*<\s*(" + scalar + r")\s*>|std::array\s*<\s*(" + scalar + r")\s*,\s*\w+\s*>"
                    r"|(" + scalar + r"))\s+&?\s*(\w+)\s*(?=[;=,\[{()])")
    out = {}
    for m in rx.finditer(text):
        element = m.group(1) or m.group(2) or m.group(3)
        base = element or m.group(4)
        kind = "enum:" + base if base in enum_names else CPP_TYPES[base]
        if element and not kind.endswith("_list"):
            kind += "_list"
        out.setdefault(m.group(5), set()).add(kind)
    return out


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
        self.local = declarations(self.text, enum_names)
        header = path.with_suffix(".H")
        self.header_rel = self.rel[:-len(path.suffix)] + ".H"
        if path.suffix == ".cpp" and header.exists():          # a .cpp's own header is local scope too
            self.header_text = strip_comments(header.read_text(errors="replace"))
            for name, kinds in declarations(self.header_text, enum_names).items():
                self.local.setdefault(name, set()).update(kinds)
        else:
            self.header_text = ""
        self.offsets = [0]
        for line in self.lines:
            self.offsets.append(self.offsets[-1] + len(line) + 1)
        self.decls = self._pp_decls()

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
            return None, None
        fam = FAMILIES.get((self.rel, (expr or "").strip()))
        if fam:
            return fam
        p = resolve_prefix(expr)
        return ([p], None) if p is not None else (None, None)

    def kind_of(self, name):
        for table in (self.local, self.global_decls):
            kinds = table.get(name, set())
            if len(kinds) == 1:
                return next(iter(kinds))
        return "unknown"

    def names(self, arg):
        if arg.startswith('"'):
            return re.findall(r'"([^"]+)"', arg)
        return identifier_keys(arg.split(".")[0], self.text, self.all_files)

    def calls(self):
        key_rx = r'("[^"]+"(?:\s*,\s*"[^"]+")?|\w+(?:\.c_str\(\))?)'
        out = []
        for m in re.finditer(r"\b(\w+)\s*\.\s*(" + "|".join(METHODS) + r")\s*\(\s*" + key_rx +
                             r"\s*,\s*([^,)]+)\s*([,)])", self.text):
            method = m.group(2) + ("_indexed" if m.group(5) == "," and m.group(2) in ("get", "query") else "")
            out += [(m.start(), m.group(1), method, key, m.group(4), False) for key in self.names(m.group(3))]
        for m in re.finditer(r"\b(\w+)\s*\.\s*(contains|countval)\s*\(\s*" + key_rx + r"\s*\)", self.text):
            out += [(m.start(), m.group(1), m.group(2), key, None, False) for key in self.names(m.group(3))]
        for m in re.finditer(r'\bquery_one_or_per_level(\w*)\s*\(\s*(\w+)\s*,\s*"([^"]+)"\s*,\s*([^,)]+)', self.text):
            method = "query_enum_case_insensitive" if "enum_case_insensitive" in m.group(1) else "query"
            out.append((m.start(), m.group(2), method, m.group(3), m.group(4), True))
        if self.erf_side:
            for m in re.finditer(r'\bread_int_string\s*\(\s*\w+\s*,\s*"([^"]+)"', self.text):
                out.append((m.start(), "__erf__", "read_int_string", m.group(1), None, True))
        return out

    def kind(self, method, target, per_level):
        if method == "read_int_string":
            return "int"
        if method in ("contains", "countval"):
            return "unknown"                                 # a presence read; a typed read of the same key wins
        kind = self.kind_of(target_name(target))
        if per_level and kind.endswith("_list"):
            kind = kind[:-5]
        if method in LIST_METHODS and not kind.endswith("_list") and kind != "unknown":
            kind += "_list"
        return kind

    def default(self, target, kind, line, fact):
        """(value, origin, source): the variable's initializer as a canonical literal. Origins: initializer, none,
        expression (an initializer that is not a plain literal of the type)."""
        name = target_name(target) if target else None
        found = _initializer(name, self.text, line) if name else None
        rel = self.rel
        if found is None and name and self.header_text:
            found, rel = _initializer(name, self.header_text, None), self.header_rel
        if found is None or kind.endswith("_list") or kind == "unknown":
            return None, "none", None
        raw, at = found
        if kind.startswith("enum:"):                     # values attach later; schema._enum_facts canonicalizes
            return _literal(raw), "initializer", f"{rel}:{at}"
        try:
            value = canonical_value(kind, _literal(raw), fact)
        except ValueError:
            return None, "expression", f"{rel}:{at}"
        return value, "initializer", f"{rel}:{at}"

    def hits(self):
        hits, unresolved = [], []
        for pos, var, method, key, target, per_level in self.calls():
            line = self.line_of(pos)
            pfx, family = (["erf"], None) if var == "__erf__" else self.prefixes(var, line)
            if pfx is None:
                is_pp = any(v == var for _, v, _ in self.decls) or re.search(r"ParmParse\s*&\s*" + var + r"\b",
                                                                              self.text)
                if is_pp:                                    # only ParmParse objects; not maps, NetCDF handles, etc.
                    unresolved.append(f"{self.rel}:{line + 1} {var}.{key}")
                continue
            fact = {"type": self.kind(method, target, per_level), "source": f"{self.rel}:{line + 1}"}
            if per_level:
                fact["per_level"] = True
            if "enum" in method:
                fact["case_insensitive"] = "case_insensitive" in method
            if family:
                fact["family"] = family
            fact["required"] = method.startswith("get")
            value, origin, src = self.default(target, fact["type"], line, fact)
            fact.update(default=value, default_origin=origin)
            if src:
                fact["default_source"] = src
            for p in pfx:
                hits.append((f"{p}.{key}" if p else key, fact, self.conds[line], method, target))
        return hits, unresolved


def _initializer(name, text, before_line):
    """(raw initializer, 1-based line) of the last declaration of name before before_line (any line if None)."""
    rx = re.compile(r"[\w>]\s+&?\s*" + re.escape(name) + r"\s*(?:=\s*([^;,{}]+?)|\{\s*([^{}]*?)\s*\})\s*[;,]")
    best = None
    for m in rx.finditer(text):
        at = text.count("\n", 0, m.start()) + 1
        if before_line is None or at <= before_line + 1:
            best = ((m.group(1) if m.group(1) is not None else m.group(2)).strip(), at)
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


def scan(root, rel_base, enum_names, global_decls, erf_side, all_files=()):
    """Every literal query in the tree: (key, fact, conditions, method, target), and unresolved reads."""
    hits, unresolved = [], []
    for path in sorted(p for p in root.rglob("*") if p.suffix in (".H", ".cpp")):
        rel = f"{rel_base}/{path.relative_to(root)}"
        h, u = FileScan(path, rel, enum_names, global_decls, erf_side, all_files).hits()
        hits += h
        unresolved += u
    return hits, sorted(set(unresolved))


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
