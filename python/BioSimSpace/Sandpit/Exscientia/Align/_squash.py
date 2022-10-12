import copy as _copy
import itertools as _it
import shutil as _shutil
import tempfile

import numpy as _np
import parmed as _pmd
from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from ._merge import _removeDummies
from ..IO import readMolecules as _readMolecules, saveMolecules as _saveMolecules
from .._SireWrappers import Molecule as _Molecule


def _squash(system):
    """Internal function which converts a merged BioSimSpace system into an AMBER-compatible format, where all perturbed
       molecules are represented sequentially, instead of in a mixed topology, like in GROMACS. In the current
       implementation, all perturbed molecules are moved at the end of the squashed system. For example, if we have an
       input system, containing regular molecules (M) and perturbed molecules (P):

       M0 - M1 - P0 - M2 - P1 - M3

       This function will return the following squashed system:

       M0 - M1 - M2 - M3 - P0_A - PO_B - P1_A - P1_B

       Where A and B denote the dummyless lambda=0 and lambda=1 states. In addition, we also
       return a mapping between the old unperturbed molecule indices and the new ones. This
       mapping can be used during coordinate update. Updating the coordinates of the perturbed
       molecules, however, has to be done manually through the Python layer.

       Parameters
       ----------

       system : BioSimSpace._SireWrappers.System
           The system.
    """
    # Create a copy of the original system.
    new_system = system.copy()

    # Get the perturbable molecules and their corresponding indices.
    pertmol_idxs = [i for i, molecule in enumerate(system.getMolecules()) if molecule.isPerturbable()]
    pert_mols = system.getPerturbableMolecules()

    # Remove the perturbable molecules from the system.
    new_system.removeMolecules(pert_mols)

    # Add them back at the end of the system. This is generally faster than keeping their order the same.
    new_indices = list(range(system.nMolecules()))
    for pertmol_idx, pert_mol in zip(pertmol_idxs, pert_mols):
        new_indices.remove(pertmol_idx)
        new_indices.append(pertmol_idx)
        new_system += _squash_molecule(pert_mol)

    # Create the old molecule index to new molecule index mapping.
    mapping = {_SireMol.MolIdx(idx): _SireMol.MolIdx(i) for i, idx in enumerate(new_indices)}

    return new_system, mapping


def _squash_molecule(molecule):
    if not molecule.isPerturbable():
        return molecule

    # Generate a "system" from the molecule at lambda = 0 and another copy at lambda = 1.
    mol0 = _squash_properties(_removeDummies(molecule, False), keep_lambda1=False)
    mol1 = _squash_properties(_removeDummies(molecule, True), keep_lambda0=False)
    system = (mol0 + mol1).toSystem()

    # We only need to call tiMerge for multi-residue molecules
    if molecule.nResidues() == 1:
        return system

    # Perform the multi-residue squashing with ParmEd as it is much easier and faster.
    with tempfile.TemporaryDirectory() as tempdir:
        # Load in ParmEd.
        _saveMolecules(f"{tempdir}/temp", mol0 + mol1, "prm7,rst7")
        _shutil.move(f"{tempdir}/temp.prm7", f"{tempdir}/temp.parm7")
        parm = _pmd.load_file(f"{tempdir}/temp.parm7", xyz=f"{tempdir}/temp.rst7")

        # Determine the molecule masks.
        mol_mask0 = f"@1-{mol0.nAtoms()}"
        mol_mask1 = f"@{mol0.nAtoms() + 1}-{system.nAtoms()}"

        # Determine the residue masks.
        atom0_offset, atom1_offset = 0, mol0.nAtoms()
        res_atoms0, res_atoms1 = [], []
        for res0, res1, res01 in zip(mol0.getResidues(), mol1.getResidues(), molecule.getResidues()):
            # We assume the residue is only perturbable if any of the matching element symbols differ
            atom_symbols0 = [atom._sire_object.property("element0").symbol() for atom in res01.getAtoms()]
            atom_symbols1 = [atom._sire_object.property("element1").symbol() for atom in res01.getAtoms()]
            if atom_symbols0 != atom_symbols1:
                res_atoms0 += list(range(atom0_offset, atom0_offset + res0.nAtoms()))
                res_atoms1 += list(range(atom1_offset, atom1_offset + res1.nAtoms()))
            atom0_offset += res0.nAtoms()
            atom1_offset += res1.nAtoms()
        res_mask0 = _amber_mask_from_indices(res_atoms0)
        res_mask1 = _amber_mask_from_indices(res_atoms1)

        # Merge the residues.
        action = _pmd.tools.tiMerge(parm, mol_mask0, mol_mask1, res_mask0, res_mask1)
        action.execute()

        # Reload into BioSimSpace.
        # TODO: prm7/rst7 doesn't work for some reason so we need to use gro/top
        parm.save(f"{tempdir}/squashed.gro", overwrite=True)
        parm.save(f"{tempdir}/squashed.top", overwrite=True)
        squashed_mol = _readMolecules([f"{tempdir}/squashed.gro", f"{tempdir}/squashed.top"])

    return squashed_mol


def _squash_properties(molecule, keep_lambda0=True, keep_lambda1=True):
    if not keep_lambda0 and not keep_lambda1:
        raise ValueError("Need to keep at least one of the endpoint lambda properties")

    sire_mol = molecule._sire_object.edit()
    types0 = sire_mol.property("ambertype0").toVector()
    dummies0 = [i for i, x in enumerate(types0) if "du" in x]

    # We start off by removing the "is_perturbable" attribute.
    if "is_perturbable" in sire_mol.propertyKeys():
        sire_mol = sire_mol.removeProperty("is_perturbable").molecule()

    alch_properties = {x[:-1] for x in sire_mol.propertyKeys() if x[-1] in ("0", "1")}
    for alch_property in alch_properties:
        prop0 = sire_mol.property(f"{alch_property}0")
        prop1 = sire_mol.property(f"{alch_property}1")

        # Delete the old properties.
        sire_mol = sire_mol.removeProperty(f"{alch_property}0").molecule()
        sire_mol = sire_mol.removeProperty(f"{alch_property}1").molecule()

        if prop0 == prop1 or (keep_lambda0 and not keep_lambda1):
            sire_mol.setProperty(alch_property, prop0)
        elif not keep_lambda0 and keep_lambda1:
            sire_mol.setProperty(alch_property, prop1)
        elif alch_property == "molecule":
            continue
        else:
            # We assume the bonded terms are equal (they should be) but we still do a basic sanity check.
            # TODO: add proper per-term checks?
            if hasattr(prop0, "nFunctions"):
                if prop0.nFunctions() != prop1.nFunctions():
                    raise ValueError("The number of functions between the two molecules doesn't match")
                sire_mol.setProperty(alch_property, prop0)
                continue

            prop0 = prop0.toVector()
            prop1 = prop1.toVector()
            prop = _copy.copy(prop0)
            for i in dummies0:
                prop[i] = prop1[i]

            # Set the properties (not yet possible to do it on a per-molecule basis, so we have to iterate).
            for i, atom_prop in enumerate(prop):
                idx = _SireMol.AtomIdx(i)
                sire_mol = sire_mol.atom(idx).setProperty(alch_property, atom_prop).molecule()

    return _Molecule(sire_mol.commit())


def _unsquash(system, squashed_system, mapping):
    """Internal function which converts an alchemical AMBER system where the perturbed molecules are
       defined sequentially and updates the coordinates and velocities of an input unsquashed system.
       Refer to the _squash() function documentation to see the structure of the squashed system
       relative to the unsquashed one.

       Parameters
       ----------

       system : BioSimSpace._SireWrappers.System
           The regular unsquashed system.

       squashed_system : BioSimSpace._SireWrappers.System
           The corresponding squashed system.

       mapping : dict(Sire.Mol.MolIdx, Sire.Mol.MolIdx)
           The molecule-molecule mapping generated by _squash().
    """
    # Create a copy of the original new_system.
    new_system = system.copy()

    # Update the unperturbed molecule coordinates in the original new_system using the mapping.
    if mapping:
        # Squash() puts all perturbable molecules at the end of the system so we prune them from the mapping now.
        pertmol_offset = len(new_system) - new_system.nPerturbableMolecules()
        nonpertmol_mapping = {k: v for k, v in mapping.items() if v.value() < pertmol_offset}
        if nonpertmol_mapping:
            new_system._sire_object, _ = _SireIO.updateCoordinatesAndVelocities(
                new_system._sire_object,
                squashed_system._sire_object,
                nonpertmol_mapping)

    # From now on we handle all perturbed molecules.
    pertmol_idxs = [i for i, molecule in enumerate(new_system.getMolecules()) if molecule.isPerturbable()]
    pertmols = new_system.getPerturbableMolecules()
    squashed_pertmols = squashed_system[-new_system.nPerturbableMolecules():]

    for i, pertmol, squashed_pertmol in zip(pertmol_idxs, pertmols, squashed_pertmols):
        new_pertmol = _unsquash_molecule(pertmol, squashed_pertmol)
        new_system.updateMolecule(i, new_pertmol)

    return new_system


def _unsquash_molecule(molecule, squashed_molecule):
    if molecule.nResidues() != squashed_molecule.nResidues():
        raise ValueError("The number of residues between the squashed and unsquashed molecules does not match")

    # Get the atom mapping and combine it with the lambda=0 molecule being prioritised
    atom_mapping0 = _squashed_atom_mapping(molecule, is_lambda1=False)
    atom_mapping1 = _squashed_atom_mapping(molecule, is_lambda1=True)
    atom_mapping = {**atom_mapping1, **atom_mapping0}
    atom_index_mask = [atom_mapping[i] for i in range(len(atom_mapping))]

    # Assign the coordinates based on the atom index mask
    siremol = molecule.copy()._sire_object.edit()
    squashed_coordinates = squashed_molecule._sire_object.property("coordinates").toVector()
    coordinates = _SireMol.AtomCoords([squashed_coordinates[i] for i in atom_index_mask])
    siremol = siremol.setProperty("coordinates0", coordinates).molecule()
    siremol = siremol.setProperty("coordinates1", coordinates).molecule()

    # Optionally update the velocities
    if squashed_molecule._sire_object.hasProperty("velocity"):
        squashed_velocities = squashed_molecule._sire_object.property("velocity").toVector()
        velocities = _SireMol.AtomVelocities([squashed_velocities[i] for i in atom_index_mask])
        siremol = siremol.setProperty("velocity0", velocities).molecule()
        siremol = siremol.setProperty("velocity1", velocities).molecule()

    return _Molecule(siremol.commit())


def _squashed_molecule_mapping(system):
    # TODO: incorporate into squash()
    # TODO: MolIdx vs int?
    # Get the perturbable molecules and their corresponding indices.
    pertmol_idxs = [i for i, molecule in enumerate(system.getMolecules()) if molecule.isPerturbable()]

    # Add them back at the end of the system. This is generally faster than keeping their order the same.
    new_indices = list(range(system.nMolecules()))
    for pertmol_idx in pertmol_idxs:
        new_indices.remove(pertmol_idx)
        new_indices.append(pertmol_idx)

    # Create the old molecule index to new molecule index mapping.
    mapping = {_SireMol.MolIdx(idx): _SireMol.MolIdx(i) for i, idx in enumerate(new_indices)}

    return mapping


def _squashed_atom_mapping(system, is_lambda1=False):
    if isinstance(system, _Molecule):
        return _squashed_atom_mapping(system.toSystem(), is_lambda1=is_lambda1)

    molecule_mapping = _squashed_molecule_mapping(system)
    molecule_mapping_rev = {v: k for k, v in molecule_mapping.items()}

    atom_mapping = {}
    atom_idx, squashed_atom_idx = 0, 0
    for i in range(len(system)):
        mol_idx = molecule_mapping_rev[_SireMol.MolIdx(i)].value()
        molecule = system[mol_idx]

        if not molecule.isPerturbable():
            atom_indices = _np.arange(atom_idx, atom_idx + molecule.nAtoms())
            squashed_atom_indices = _np.arange(squashed_atom_idx, squashed_atom_idx + molecule.nAtoms())
            atom_mapping.update(dict(zip(atom_indices, squashed_atom_indices)))
            atom_idx += molecule.nAtoms()
            squashed_atom_idx += molecule.nAtoms()
        else:
            for residue in molecule.getResidues():
                in_mol0 = ["du" not in atom._sire_object.property("ambertype0") for atom in residue.getAtoms()]
                in_mol1 = ["du" not in atom._sire_object.property("ambertype1") for atom in residue.getAtoms()]
                ndummy0 = residue.nAtoms() - sum(in_mol1)
                ndummy1 = residue.nAtoms() - sum(in_mol0)
                ncommon = residue.nAtoms() - ndummy0 - ndummy1
                natoms0 = ncommon + ndummy0
                natoms1 = ncommon + ndummy1

                if ndummy0 == ndummy1 == 0:
                    atom_indices = _np.arange(atom_idx, atom_idx + residue.nAtoms())
                    squashed_atom_indices = _np.arange(squashed_atom_idx, squashed_atom_idx + residue.nAtoms())
                    squashed_atom_idx += residue.nAtoms()
                else:
                    if not is_lambda1:
                        atom_indices = _np.arange(atom_idx, atom_idx + residue.nAtoms())[in_mol0]
                        squashed_atom_indices = _np.arange(squashed_atom_idx, squashed_atom_idx + natoms0)
                    else:
                        atom_indices = _np.arange(atom_idx, atom_idx + residue.nAtoms())[in_mol1]
                        offset = squashed_atom_idx + natoms0
                        squashed_atom_indices = _np.arange(offset, offset + natoms1)
                    squashed_atom_idx += natoms0 + natoms1
                atom_mapping.update(dict(zip(atom_indices, squashed_atom_indices)))
                atom_idx += residue.nAtoms()

    # Convert from NumPy integers to Python integers.
    return {int(k): int(v) for k, v in atom_mapping.items()}


def _amber_mask_from_indices(atom_idxs):
    """Internal helper function to create an AMBER mask from a list of atom indices.

       Parameters
       ----------

       atom_idxs : [int]
           A list of atom indices.

       Returns
       -------

       mask : str
           The AMBER mask.
    """
    # AMBER has a restriction on the number of characters in the restraint
    # mask (not documented) so we can't just use comma-separated atom
    # indices. Instead we loop through the indices and use hyphens to
    # separate contiguous blocks of indices, e.g. 1-23,34-47,...

    if atom_idxs:
        # AMBER masks are 1-indexed, while BioSimSpace indices are 0-indexed.
        atom_idxs = [x + 1 for x in sorted(list(set(atom_idxs)))]
        if not all(isinstance(x, int) for x in atom_idxs):
            raise TypeError("'atom_idxs' must be a list of 'int' types.")
        groups = []
        initial_idx = atom_idxs[0]
        for prev_idx, curr_idx in _it.zip_longest(atom_idxs, atom_idxs[1:]):
            if curr_idx != prev_idx + 1 or curr_idx is None:
                if initial_idx == prev_idx:
                    groups += [str(initial_idx)]
                else:
                    groups += [f"{initial_idx}-{prev_idx}"]
                initial_idx = curr_idx
        mask = "@" + ",".join(groups)
    else:
        mask = ""

    return mask