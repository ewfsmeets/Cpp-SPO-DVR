"""The tool's structure: every module has a layer and imports only from lower layers (or its own package) and
the third-party packages its layer may use, and every `Shadows:` line of a module docstring names existing src/
files and functions."""

import ast
import re
import sys
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parent.parent
REPOSITORY = TOOL_DIR.parent.parent

# The layer of each package or module (README, "Architecture"): a module imports only from strictly lower
# layers or from its own package.
LAYERS = {
    "copies": 0,
    "readers": 1,
    "compute": 2,
    "style": 3,
    "explanations": 3,
    "panels": 4,
    "figures": 5,
    "verdict": 3,
    "report": 6,
    "viewer": 6,
    "spodvr_plot": 7,
    "view": 7,
    "tests": 8,
}
NUMPY = {"numpy"}
PLOTTING = {"numpy", "matplotlib"}
THIRD_PARTY = {
    0: NUMPY,
    1: NUMPY,
    2: NUMPY,
    3: PLOTTING,
    4: PLOTTING,
    5: PLOTTING,
    6: PLOTTING,
    7: PLOTTING,
    8: PLOTTING | {"pytest"},
}


def tool_modules() -> list[Path]:
    """Every Python file of the tool."""
    return sorted(p for p in TOOL_DIR.rglob("*.py") if "__pycache__" not in p.parts)


def layer_name(path: Path) -> str:
    """The package or top-level module a file belongs to."""
    return path.relative_to(TOOL_DIR).parts[0].removesuffix(".py")


def imported_names(path: Path) -> set[str]:
    """The top-level names a file imports."""
    names = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names |= {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0:
            names.add(node.module.split(".")[0])
    return names


def layer_breaches(path: Path, imports: set[str]) -> list[str]:
    """The imports of one module that break the layers: from its own or a higher layer (other than its
    package), or a package its layer may not use."""
    importer = layer_name(path)
    if importer not in LAYERS:
        return [f"{path.relative_to(TOOL_DIR)} has no layer in LAYERS"]
    rank = LAYERS[importer]
    wrong = []
    for name in sorted(imports):
        if name in LAYERS and name != importer and LAYERS[name] >= rank:
            wrong.append(f"{path.relative_to(TOOL_DIR)} (layer {rank}) imports {name} (layer {LAYERS[name]})")
        elif name not in LAYERS and name not in sys.stdlib_module_names and name not in THIRD_PARTY[rank]:
            wrong.append(f"{path.relative_to(TOOL_DIR)} imports {name}")
    return wrong


def test_layers() -> None:
    """No module imports from its own or a higher layer (other than its package) or an unlisted package."""
    wrong = [w for path in tool_modules() for w in layer_breaches(path, imported_names(path))]
    assert not wrong, "\n".join(wrong)


def test_layer_check_finds_a_breach() -> None:
    """Negative test: readers importing compute (a higher layer), and compute importing matplotlib."""
    assert layer_breaches(TOOL_DIR / "readers" / "tables.py", {"compute", "numpy"}) == [
        "readers/tables.py (layer 1) imports compute (layer 2)"
    ]
    assert layer_breaches(TOOL_DIR / "compute" / "results.py", {"matplotlib"}) == [
        "compute/results.py imports matplotlib"
    ]


def citation_breaches(module: str, docstring: str) -> tuple[int, list[str]]:
    """The functions a `Shadows:` line names (as `file (name, Class::name)`) that the file does not define or
    call by that name; and how many names were checked."""
    match = re.search(r"Shadows:(.*?)(?:\n\n|$)", docstring, re.DOTALL)
    if not match:
        return 0, []
    checked, missing = 0, []
    for file, names in re.findall(r"(src/[\w/.*]+)(?:\s*\(([^)]*)\))?", match.group(1)):
        files = sorted(REPOSITORY.glob(file))
        if not files:
            missing.append(f"{module}: {file} does not exist")
            continue
        text = "".join(f.read_text() for f in files)
        for name in filter(None, (n.strip() for n in names.split(","))):
            checked += 1
            if not re.search(rf"\b{re.escape(name)}\s*\(", text):
                missing.append(f"{module}: {name} not in {file}")
    return checked, missing


def test_cited_sources_exist() -> None:
    """Every file named in a `Shadows:` line exists and defines or calls every function named after it."""
    checked, missing = 0, []
    for path in tool_modules():
        count, wrong = citation_breaches(path.name, ast.get_docstring(ast.parse(path.read_text())) or "")
        checked, missing = checked + count, missing + wrong
    assert not missing, "\n".join(missing)
    assert checked >= 8, f"only {checked} citations found: the pattern no longer matches"


def test_citation_check_finds_a_renamed_function() -> None:
    """Negative test: a function the C++ file does not have (a rename) is reported."""
    docstring = "Shadows: src/analysis/s_matrix.cpp (SMatrixAccumulator::add, SMatrixAccumulator::accumulate)"
    assert citation_breaches("x.py", docstring) == (
        2,
        ["x.py: SMatrixAccumulator::accumulate not in src/analysis/s_matrix.cpp"],
    )


MAX_FUNCTION_LINES = 40  # code lines (docstring not counted); copies/ follow the C++ and are exempt
ALLOWED_DECORATORS = {"dataclass", "pytest.fixture", "pytest.mark.parametrize"}
ALLOWED_CLASSES = {
    "ReadError",
    "FrameViewer",
}  # besides @dataclass classes; the viewer's window (README, "Conventions")


def functions(tree: ast.AST) -> list[ast.FunctionDef]:
    """Every function definition of a module, nested ones included."""
    return [node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]


def code_lines(function: ast.FunctionDef) -> int:
    """The lines of a function from its `def` to its end, without its docstring."""
    lines = function.end_lineno - function.lineno + 1
    first = function.body[0]
    if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
        lines -= first.end_lineno - first.lineno + 1
    return lines


def decorator_name(decorator: ast.expr) -> str:
    """ "dataclass", "pytest.fixture", …: a decorator without its arguments."""
    target = decorator.func if isinstance(decorator, ast.Call) else decorator
    return ast.unparse(target)


def convention_breaches(source: str, path: Path) -> list[str]:
    """The conventions of one module (its source, at `path` in the tool) that ruff does not check: function
    length, docstrings of every function (nested and private ones too), decorators, classes, lambdas and
    partial in figures.py, matplotlib.use, Test* imports."""
    tree = ast.parse(source)
    name = path.relative_to(TOOL_DIR)
    breaches = []
    for function in functions(tree):
        if layer_name(path) != "copies" and code_lines(function) > MAX_FUNCTION_LINES:
            breaches.append(f"{name}: {function.name} has {code_lines(function)} lines (> {MAX_FUNCTION_LINES})")
        if ast.get_docstring(function) is None:
            breaches.append(f"{name}: {function.name} has no docstring")
    for node in ast.walk(tree):
        decorators = node.decorator_list if isinstance(node, (ast.FunctionDef, ast.ClassDef)) else []
        breaches += [
            f"{name}: decorator @{decorator_name(d)}" for d in decorators if decorator_name(d) not in ALLOWED_DECORATORS
        ]
        if isinstance(node, ast.ClassDef) and node.name not in ALLOWED_CLASSES and not decorators:
            breaches.append(f"{name}: class {node.name} is neither a @dataclass nor an allowed class")
        if isinstance(node, ast.Lambda) and path.name == "figures.py":
            breaches.append(f"{name}: a lambda in figures.py")
        if isinstance(node, ast.Attribute) and ast.unparse(node) == "matplotlib.use":
            breaches.append(f"{name}: matplotlib.use (the backend comes from the environment)")
        if isinstance(node, ast.ImportFrom) and any(a.name.startswith("Test") for a in node.names):
            breaches.append(f"{name}: imports a name starting with Test (pytest would collect it)")
        if isinstance(node, ast.Name) and node.id == "partial" and path.name == "figures.py":
            breaches.append(f"{name}: partial in figures.py")
    return breaches


def test_conventions() -> None:
    """Every module keeps the conventions of README "Conventions" that ruff does not check."""
    breaches = [b for path in tool_modules() for b in convention_breaches(path.read_text(), path)]
    assert not breaches, "\n".join(breaches)


def test_conventions_check_finds_a_breach() -> None:
    """Negative test: a long undocumented function, a @property and a lambda in figures.py are each found."""
    body = "\n".join(f"    x{i} = {i}" for i in range(MAX_FUNCTION_LINES + 1))
    source = f"from dataclasses import dataclass\n\n\ndef long():\n{body}\n\n\nclass Plain:\n    @property\n"
    source += "    def value(self):\n        '''Doc.'''\n        return 1\n\n\nf = lambda: 0\n"
    breaches = convention_breaches(source, TOOL_DIR / "figures.py")
    assert sorted(b.split(": ", 1)[1] for b in breaches) == sorted(
        [
            f"long has {MAX_FUNCTION_LINES + 2} lines (> {MAX_FUNCTION_LINES})",
            "long has no docstring",
            "decorator @property",
            "class Plain is neither a @dataclass nor an allowed class",
            "a lambda in figures.py",
        ]
    )
