#!/usr/bin/env python

import numpy as np
import sys
from rdkit.Chem import MolToSmiles as mol2smi
from rdkit import RDLogger
from cell2mol.xyz2mol import xyz2mol
import ase
from ase.io.extxyz import read_xyz
from ase import neighborlist

RDLogger.DisableLog("rdApp.*")


def get_cutoffs(z, radii=ase.data.covalent_radii, mult=1):
    return [radii[zi] * mult for zi in z]


def check_symmetric(am, tol=1e-8):
    return np.linalg.norm(am - am.T, np.Inf) < tol


def check_connected(am, tol=1e-8):
    sums = am.sum(axis=1)
    lap = np.diag(sums) - am
    eigvals, eigvects = np.linalg.eig(lap)
    return len(np.where(abs(eigvals) < tol)[0]) < 2


def smilify(filename, charge=0):
    covalent_factors = [1.10, 1.15, 1.05, 1.20, 1.00]
    for covalent_factor in covalent_factors:
        assert filename[-4:] == ".xyz"
        name = filename.strip().split(".")[-2].split("/")[-1]
        mol = next(read_xyz(open(filename)))
        atoms = mol.get_chemical_symbols()
        z = [int(zi) for zi in mol.get_atomic_numbers()]
        coordinates = mol.get_positions()
        cutoff = get_cutoffs(z, radii=ase.data.covalent_radii, mult=covalent_factor)
        nl = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
        nl.update(mol)
        AC = nl.get_connectivity_matrix(sparse=False)
        try:
            assert check_connected(AC) and check_symmetric(AC)
            mol = xyz2mol(
                z,
                coordinates,
                AC,
                covalent_factor,
                charge=0,
                use_graph=True,
                allow_charged_fragments=True,
                embed_chiral=True,
                use_huckel=False,
            )
            if isinstance(mol, list):
                mol = mol[0]
            smiles = mol2smi(mol)
            if isinstance(smiles, list):
                smiles = smiles[0]
            print(smiles)
        except Exception as m:
            #print(m)
            continue


if __name__ == "__main__":
    smilify(sys.argv[1])
