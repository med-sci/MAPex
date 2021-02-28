from sklearn import cluster
from rdkit.Chem import ChemicalFeatures
import math
from pathlib import Path
import numpy as np



class Pharmacophore:
    def __init__(self, mol, gene):
        self._fdef = str(Path('MAPex', 'BaseFeatures.fdef'))
        self._mol = mol
        self._gene = gene
        self._factory = ChemicalFeatures.BuildFeatureFactory(self._fdef)
        self.acceptors, self.donors, self.hydrophobics = self._features

    @property
    def _features(self):
        mol = self._mol
        acceptors = [
            list(f.GetPos(confId = self._gene))
            for f in self._factory.GetFeaturesForMol(mol, includeOnly="Acceptor", confId=self._gene)
        ]
        donors = [
            list(f.GetPos(confId = self._gene))
            for f in self._factory.GetFeaturesForMol(mol, includeOnly="Donor", confId=self._gene)
        ]
        hydrophobics = [
            list(f.GetPos(confId = self._gene))
            for f in self._factory.GetFeaturesForMol(mol, includeOnly="Hydrophobe", confId=self._gene)
        ]
        return acceptors, donors, hydrophobics


class PharmComplex:
    def __init__(self, molecules, chromosome):
        self.molecules = molecules
        self.chromosome = chromosome
        self._p_list = [Pharmacophore(m, g) for m, g in zip(self.molecules, self.chromosome)]
        self._feat_dict = self._features
        self._coords = {}

    def get_coords(self):
        return self._coords

    @property
    def _features(self):
        accs = []
        dons = []
        hyds = []
        for pharmacophore in self._p_list:
            accs.extend(pharmacophore.acceptors)
            dons.extend(pharmacophore.donors)
            hyds.extend(pharmacophore.hydrophobics)
        feat_dict = {"Acceptors": accs, "Donors": dons, "Hydrophobics": hyds}
        return feat_dict

    @property
    def _num_mols(self):
        return math.ceil(len(self._p_list) * 0.75)

    @staticmethod
    def _cluster_feats(feats, eps=1, min_samples=4):
        dbscan = cluster.DBSCAN(eps, min_samples=min_samples)
        dbscan.fit(feats)
        lbls = dbscan.labels_
        n_clus = len(np.unique(dbscan.labels_)) - 1
        clus_dict = {}
        for k in range(n_clus):
            clus_dict[k] = []

        for i in range(len(feats)):
            for j in range(n_clus):
                if lbls[i] == j:
                    clus_dict[j].append(feats[i])
        return clus_dict

    @staticmethod
    def _get_centers(clus_dict):
        centr_list = []
        for k in clus_dict.keys():
            x, y, z = (0, 0, 0)
            for i in clus_dict[k]:
                x = x + i[0]
                y = y + i[1]
                z = z + i[2]
            den = len(clus_dict[k])
            centr_list.append(
                [round((x / den), 4), round((y / den), 4), round((z / den), 4)]
            )
        return centr_list

    def create(self, distance=1, num_mols=None):
        if num_mols:
            num_mols = num_mols
        else:
            num_mols = self._num_mols
        feat_dict = self._feat_dict
        coords = {}
        for k in feat_dict.keys():
            clus_dict = self._cluster_feats(feat_dict[k], distance, num_mols)
            centr_list = self._get_centers(clus_dict)
            coords[k] = centr_list
        self._coords = coords

    @staticmethod
    def __add_spaces(length, base_space=6):
        diff = base_space - length
        total = base_space + diff
        return " " * total


    def write(self):
        coord_dict = self._coords
        l = 0
        for k in coord_dict.keys():
            l = l + len(coord_dict[k])
        f = open(str(Path('examples', 'pharmacophore.mol2')), "w")
        f.write("@<TRIPOS>MOLECULE \npharmacophore\n")
        f.write("    " + str(l) + "    0\nSMALL \nNO_CHARGES \n@<TRIPOS>ATOM")
        n = 1
        for k in coord_dict.keys():
            prop = None
            if k == "Acceptors":
                prop = "ACC"
            elif k == "Donors":
                prop = "DON"
            else:
                prop = "HYD"
            for v in coord_dict[k]:
                f.write(
                    f"\n      {n} {prop}   {self.__add_spaces(len(str(v[0])))}{v[0]}{self.__add_spaces(len(str(v[1])))}{v[1]}{self.__add_spaces(len(str(v[2])))}{v[2]} HB    {l - 1}    HB      0.0000"
                )
        f.write("\n@<TRIPOS>BOND")
        f.close()
