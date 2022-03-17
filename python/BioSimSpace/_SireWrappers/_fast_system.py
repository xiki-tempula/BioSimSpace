# TODO: remove this file when BSS performance is sufficient for large systems
from collections import defaultdict as _defaultdict

from Sire import IO as _SireIO
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace._SireWrappers import Molecules as _Molecules
from BioSimSpace._SireWrappers import System as _System


def _fastSystemInit(left_molecules, right_molecules=None):
        """Add a molecule, or list of molecules to the system.

           Parameters
           ----------

           molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                       :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                       [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`], \
                       :class:`System <BioSimSpace._SireWrappers.System>`
              A Molecule, Molecules object, a list of Molecule objects, or a System containing molecules.
        """

        system = _SireSystem.System("BioSimSpace System")
        molgrp = _SireMol.MoleculeGroup("all")

        # Add the molecule groups in a quick way.
        for molecules in [left_molecules, right_molecules]:
            if molecules is None:
                continue
            # A BioSimSpace System object.
            if type(molecules) is _System:
                molecules = molecules.getMolecules()

            # Convert tuple to a list.
            if type(molecules) is tuple:
                molecules = list(molecules)

            # A Molecule object.
            if type(molecules) is _Molecule:
                molecules = [molecules]

            if type(molecules) is _SireMol.Molecule:
                molecules = [_Molecule(molecules)]

            if type(molecules) is not _Molecules:
                molecules = _Molecules(molecules)

            molgrp.add(molecules._sire_object)

        # Initialise the new system and renumber its molecules in a quick way.
        system.add(molgrp)
        _SireIO.renumberConstituents(system, system.nMolecules())
        system = _System(system)

        # Reset the name because it gets broken for some reason.
        system._sire_object.setName("BioSimSpace System")

        if type(left_molecules) is _System:
            # Copy all the properties from the old system.
            for prop in left_molecules._sire_object.propertyKeys():
                val = system._sire_object.property(prop)
                system._sire_object.setProperty(prop, val)

        return system


def _fastFixSomdPDB(filename):
    residues = _defaultdict(list)

    # Extract each molecule into a dictionary.
    with open(filename) as f:
        key = None
        for line in f.readlines():
            if line.strip() not in ["TER", "END"]:
                if key is None:
                    # The insertion code, followed by the molecule number.
                    key = line[26], int(line[22:26])
                residues[key].append(line)
            else:
                key = None

    # Write the sorted dictionary.
    sorted_keys = sorted(residues)
    with open(filename, "w") as f:
        for key in sorted_keys:
            residue = residues[key]
            for line in residue:
                f.write(line)
            f.write("TER\n")
        f.write("END\n")
