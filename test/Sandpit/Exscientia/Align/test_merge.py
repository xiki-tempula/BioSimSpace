import pytest

import BioSimSpace.Sandpit.Exscientia as BSS


@pytest.fixture
def perturbed_system():
    # N_atoms are: 12, 15, 18, 21, 24, 27 and 30.
    mol_smiles = ["c1ccccc1", "c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccccc1CCCC", "c1ccccc1CCCCC",
                  "c1ccccc1CCCCCC"]
    mols = [BSS.Parameters.gaff(smi).getMolecule() for smi in mol_smiles]
    pert_mols = [mols[0], BSS.Align.merge(mols[1], mols[2]), mols[3], mols[4], BSS.Align.merge(mols[5], mols[6])]
    system = BSS._SireWrappers.System(pert_mols)
    return system


def test_squash(perturbed_system):
    squashed_system, mapping = BSS.Align._merge._squash(perturbed_system)
    assert len(squashed_system) == 5
    n_atoms = [mol.nAtoms() for mol in squashed_system]
    assert squashed_system[-2].getResidues()[0].name() == "LIG"
    assert squashed_system[-1].getResidues()[0].name() == "LIG"
    # First we must have the unperturbed molecules, and then the perturbed ones.
    assert n_atoms == [12, 21, 24, 15 + 18, 27 + 30]
    python_mapping = {k.value(): v.value() for k, v in mapping.items()}
    assert python_mapping == {0: 0, 2: 1, 3: 2, 1: 3, 4: 4}


def test_unsquash(perturbed_system):
    squashed_system, mapping = BSS.Align._merge._squash(perturbed_system)
    new_perturbed_system = BSS.Align._merge._unsquash(perturbed_system, squashed_system, mapping)
    assert [mol0.nAtoms() == mol1.nAtoms() for mol0, mol1 in zip(perturbed_system, new_perturbed_system)]
    assert [mol0.isPerturbable() == mol1.isPerturbable() for mol0, mol1 in zip(perturbed_system, new_perturbed_system)]
