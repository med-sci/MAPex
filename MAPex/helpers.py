from rdkit import Chem
from rdkit.Chem import AllChem, PyMol, SDWriter
from tqdm.notebook import tqdm
from pathlib import Path

colors = {
            "Donors": (0, 0.9, 0),  # Green
            "Acceptors": (0.9, 0, 0),  # Red
            "Hydrophobics": (1, 0.9, 0),  # Yellow
        }


def generate_conformations(smiles_dict, num_confs):
    mol_list = []
    for k, v in tqdm(smiles_dict.items()):
        m = Chem.MolFromSmiles(v)
        m.SetProp('_Name', k)
        m = AllChem.AddHs(m)
        c = AllChem.EmbedMultipleConfs(m, num_confs)
        assert len(c) == num_confs
        mol_list.append(AllChem.RemoveHs(m))
    return mol_list


def py_mol_viz(mols, chromosome, pharmacophore_coords=None, cols=None, render = True):
    if cols is None:
        cols = colors
    try:
        v = PyMol.MolViewer()
        for i, m in enumerate(mols):
            v.ShowMol(m, name=m.GetProp('_Name'), confId=chromosome[i], showOnly=False)
        if pharmacophore_coords:
            for k, w in pharmacophore_coords.items():
                v.AddPharmacophore(
                    w,
                    label=f"P-{k}",
                    colors=[cols[k] for _ in range(len(w))],
                    sphereRad=0.6,
                )
    except ConnectionRefusedError:
        print('Run your pyMol in server mode')
    if render:
        v.SaveFile(str(Path('examples', 'molecule_complex.png')))
