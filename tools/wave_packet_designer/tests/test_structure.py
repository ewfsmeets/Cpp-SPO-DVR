"""Case structure: the `Shadows:` citations still exist, and the modules keep their layers (README,
"Architecture")."""

import ast
import re
import sys
from pathlib import Path

from tests.helpers import REPOSITORY, TOOL_DIR, check

# Layer of each module or package of the tool; a module may import only from lower layers, and from its own
# package (program imports program). Front ends (layer 4) never import each other; the tests sit on top.
LAYERS = {
    "program": 0,
    "design": 1,
    "settings": 1,
    "input_file": 1,
    "dz_gas_phase": 1,
    "evaluation": 2,
    "report": 3,
    "wave_packet_designer": 4,
    "tests": 5,
    "test_designer": 6,
}
# Third-party packages each layer may use (the standard library is always allowed).
THIRD_PARTY = {0: {"numpy"}, 1: {"numpy"}, 2: {"numpy"}, 3: {"numpy"}, 4: {"numpy"}, 5: {"numpy"}, 6: set()}


def tool_modules() -> list[Path]:
    return sorted(p for p in TOOL_DIR.rglob("*.py") if "__pycache__" not in p.parts)


def module_name(path: Path) -> str:
    """The top-level name of a tool file: its package (program, tests) or its file name."""
    return path.relative_to(TOOL_DIR).parts[0].removesuffix(".py")


def test_cited_sources_exist(work: Path) -> None:
    """Every file named in a `Shadows:` line exists, and every function named after it occurs in that file."""
    citations = 0
    for path in tool_modules():
        docstring = ast.get_docstring(ast.parse(path.read_text())) or ""
        match = re.search(r"Shadows:(.*?)(?:\n\n|$)", docstring, re.DOTALL)
        if not match:
            continue
        for file, names in re.findall(r"(src/[\w/.*]+)(?:\s*\(([^)]*)\))?", match.group(1)):
            files = sorted(REPOSITORY.glob(file))
            check(files != [], f"{path.name}: {file} does not exist")
            text = "".join(f.read_text() for f in files)
            for name in filter(None, (n.strip() for n in names.split(","))):
                check(name.split("::")[-1] in text, f"{path.name}: {name} not found in {file}")
                citations += 1
    check(citations > 20, f"only {citations} citations found: the pattern no longer matches")


def imported_names(path: Path) -> set[str]:
    """The top-level names a file imports."""
    names = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names |= {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0:
            names.add(node.module.split(".")[0])
    return names


def test_layers(work: Path) -> None:
    """No module imports from a higher layer, a front end, or a third-party package its layer does not allow."""
    for path in tool_modules():
        importer = module_name(path)
        check(importer in LAYERS, f"{path.relative_to(TOOL_DIR)} has no layer in LAYERS")
        rank = LAYERS[importer]
        for name in imported_names(path):
            if name in LAYERS:
                allowed = name == importer or LAYERS[name] < rank
                check(allowed, f"{path.relative_to(TOOL_DIR)} (layer {rank}) imports {name} (layer {LAYERS[name]})")
            elif name not in sys.stdlib_module_names:
                check(name in THIRD_PARTY[rank], f"{path.relative_to(TOOL_DIR)} imports {name}")


TESTS = [test_cited_sources_exist, test_layers]
