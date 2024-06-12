from pybmrb import ChemicalShiftStatistics



ATOM_DICT = {
    'A': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C2", "H2", "C8", "H8", "N6", "H61", "H62", "N1", "H1"],
    },

    'DA': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C2", "H2", "C8", "H8", "N6", "H61", "H62", "N1", "H1"],
    },

    'G': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["N1", "H1", "N2", "H21", "H22", "N3", "H3", "C8", "H8", "N7", "H7"],
    },

    'DG': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["N1", "H1", "N2", "H21", "H22", "N3", "H3", "C8", "H8", "N7", "H7"],
    },

    'C': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N4", "H41", "H42", "N3", "H3"],
    },

    'DC': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N4", "H41", "H42", "N3", "H3"],
    },

    'T': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N3", "H3"],
    },

    'DT': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N3", "H3","C7","H71","H72","H73"],
    },

    'U': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C7", "H71", "N3", "H3"],
    },

    'DU': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C7", "H71", "N3", "H3"],
    },

    'GLY': {
        'backbone': ['H', 'CA', 'HA2', 'HA3', 'N', 'C'],
    },

    'ALA': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB1', 'HB2', 'HB3'],
    },

    'VAL': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG1', 'CG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23'],
    },

    'LEU': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'CD2', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
    },

    'ILE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'CG2', 'HG21', 'HG22',
                      'HG23', ],
    },

    'MET': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CE', 'HE1', 'HE2', 'HE3'],
    },

    'PHE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ','CZ'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'TYR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2','CG','HH'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'TRP': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HD1', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1', 'NE1'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'CSE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'CYS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'PRO': {
        'backbone': ['HA', 'CA', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3'],
    },

    'SER': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'THR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23'],
    },

    'ASN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'ND2', 'HD21', 'HD22'],
    },

    'GLN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'NE2', 'HE21', 'HE22'],
    },

    'ASP': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG'],
    },

    'GLU': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD'],
    },

    'LYS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ'],
    },

    'HIS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD2', 'HD2', 'ND1', 'HD1', 'NE2', 'HE2', 'CE1', 'HE1'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'ARG': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'NH2',
                      'HH11', 'HH12', 'HH21', 'HH22'],
    },

}


def get_full_stat():
    atoms = ['C*','H*','N*','P*']
    #data = ChemicalShiftStatistics.get_data_from_bmrb('ALA', 'H*',)

    s={}
    for atm in atoms:
        c, data = ChemicalShiftStatistics.get_data('*', atm, filtered=False, standard_amino_acids=True)
        print (atm,len(data))
        for row in data:
            #print (row)
            if row[4] not in s:
                s[row[4]]={}
            if row[5] not in s[row[4]]:
                s[row[4]][row[5]]=0
            s[row[4]][row[5]]+=1
    f=open('stat_out.csv','w')
    for k in s:
        for k1 in s[k]:
            f.write(f'{k},{k1},{s[k][k1]}\n')
            #print (f'{k},{k1},{s[k][k1]}\n')
    f.close()
    print(data)

if __name__ == "__main__":
    get_full_stat()