"""

{{header}}

Meta
----
depends:
  - meta.pandas.pickl
  - trajs
  - top.pdb
"""
import mdtraj as md
from msmbuilder.preprocessing import RobustScaler
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic
from msmadapter.utils import get_sctrajs
from multiprocessing import Pool
import numpy as np


# Featurize logic
def feat(irow):
    i, row = irow
    print('Loading traj {}'.format(row['traj_fn']))
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    top = traj.topology

    ctni_atoms = []
    ctnt_atoms = []
    for res in [292, 329, 337, 383]:
        ca = top.select('resid {} and name CA'.format(res - 1))
        ctni_atoms.append(ca[0])
    for res in [162, 184, 185, 237]:
        ca = top.select('resid {} and name CA'.format(res - 1))
        ctnt_atoms.append(ca[0])

    atom_indices = np.array([
        ctni_atoms,
        ctnt_atoms
    ])

    diheds = md.compute_dihedrals(traj, atom_indices)

    return i, diheds


if __name__ == '__main__':
    meta = load_meta()
    tops = preload_tops(meta)

    with Pool() as pool:
        dtrajs = dict(
            pool.imap_unordered(
                feat,
                meta.iterrows()
            )
        )
    save_trajs(dtrajs, 'dtrajs', meta)
