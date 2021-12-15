import tempfile
from pathlib import Path

import BioSimSpace as BSS


def test_read_write():
    methane = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/methane*"))
    methanol = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/methanol*"))
    water = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/water*"))
    merged_mol = BSS.Align.merge(methane[0], methanol[0]).toSystem()

    slow_system = merged_mol + water
    fast_system = BSS._SireWrappers._fast_system._fastSystemInit(merged_mol, water)

    with tempfile.TemporaryDirectory() as tempdir:
        base = Path(tempdir)
        BSS.IO.saveMolecules(f"{base}/slow", slow_system, "gro87,grotop")
        BSS.IO.saveMolecules(f"{base}/fast", fast_system, "gro87,grotop")

        gro_lines_slow = set(open(f"{base}/slow.gro").readlines())
        gro_lines_fast = set(open(f"{base}/fast.gro").readlines())
        gro_diff_slow = gro_lines_slow - gro_lines_fast
        gro_diff_fast = gro_lines_fast - gro_lines_slow

        top_lines_slow = set(open(f"{base}/slow.top").readlines())
        top_lines_fast = set(open(f"{base}/fast.top").readlines())
        top_diff_slow = top_lines_slow - top_lines_fast
        top_diff_fast = top_lines_fast - top_lines_slow

        assert gro_diff_fast == gro_diff_slow == set()
        assert top_diff_fast == top_diff_slow == set()
