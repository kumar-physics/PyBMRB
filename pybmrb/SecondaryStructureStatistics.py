#!/usr/bin/env python3
"""
Short description TBD
"""
import json
import logging
from urllib.request import urlopen, Request
import numpy
import gzip
from typing import Union, List, Optional
from mmcif.io.PdbxReader import PdbxReader
# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)

# _API_URL = "http://dev-api.bmrb.io/v2"
_API_URL = "http://api.bmrb.io/v2"
_PDB_BMRB_MAPPING = "/mappings/bmrb/pdb?format=json&match_type=exact"
_FTP_BMRB_PATH = "/projects/BMRB/public/ftp/pub/bmrb/entry_directories/bmr19109/bmr19109_3.str"
_FTP_PDB_PATH = "/projects/BMRB/public/ftp/pub/pdb/data/structures/divided/mmCIF"
three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


def _get_bmrb_pdb_mapping():
    url = Request(_API_URL+_PDB_BMRB_MAPPING)
    url.add_header('Application', 'PyBMRB')
    r = urlopen(url)
    dump = json.loads(r.read())
    for i in dump:
        print (i['bmrb_id'],i['pdb_ids'])
    return dump


def get_dssp_ss(cif_file):
    cif_data = []
    if cif_file.endswith(".gz"):
        ifh = gzip.open(cif_file,'rt')
    else:
        ifh = open(cif_file, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    struct_conf = c0.getObj('struct_conf')
    col_names = struct_conf.getAttributeList()
    print (col_names)
    conf_type_idx = col_names.index('conf_type_id')
    beg_comp_idx = col_names.index('beg_label_comp_id')
    beg_asym_idx = col_names.index('beg_label_asym_id')
    beg_seq_idx = col_names.index('beg_label_seq_id')
    beg_auth_comp_idx = col_names.index('beg_auth_comp_id')
    beg_auth_asym_idx = col_names.index('beg_auth_asym_id')
    beg_auth_seq_idx = col_names.index('beg_auth_seq_id')
    end_comp_idx = col_names.index('end_label_comp_id')
    end_asym_idx = col_names.index('end_label_asym_id')
    end_seq_idx = col_names.index('end_label_seq_id')
    end_auth_comp_idx = col_names.index('end_auth_comp_id')
    end_auth_asym_idx = col_names.index('end_auth_asym_id')
    end_auth_seq_idx = col_names.index('end_auth_seq_id')
    for dat in struct_conf.getRowList():
        conf_type = dat[conf_type_idx]
        beg_seq_id = dat[beg_auth_seq_idx]
        beg_asym_id = dat[beg_auth_asym_idx]
        beg_comp_id = dat[beg_auth_comp_idx]
        end_seq_id = dat[end_auth_seq_idx]
        end_asym_id = dat[end_auth_asym_idx]
        end_comp_id = dat[end_auth_comp_idx]
        print (conf_type,beg_asym_id,beg_comp_id,beg_seq_id,end_asym_id,end_comp_id,end_seq_id)
    # model_id = col_names.index('pdbx_PDB_model_num')
    # x_id = col_names.index('Cartn_x')
    # y_id = col_names.index('Cartn_y')
    # z_id = col_names.index('Cartn_z')
    # atom_id = col_names.index('label_atom_id')
    # comp_id = col_names.index('label_comp_id')
    # asym_id = col_names.index('label_asym_id')
    # entity_id = col_names.index('label_entity_id')
    # seq_id = col_names.index('label_seq_id')
    # icode_id = col_names.index('pdbx_PDB_ins_code')
    # alt_id = col_names.index('label_alt_id')
    # aut_seq_id = col_names.index('auth_seq_id')
    # aut_asym_id = col_names.index('auth_asym_id')
    # aut_atom_id = col_names.index('auth_atom_id')
    # aut_comp_id = col_names.index('auth_comp_id')
    # pdb_models = {}
    # atom_ids = {}
    # for model in range(1, max_models + 1):
    #     pdb = {}
    #     aid = {}
    #     for dat in atom_site.getRowList():
    #         if int(dat[model_id]) == model:
    #             if use_auth_tag:
    #                 aid[(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
    #                     (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
    #                      dat[alt_id], dat[icode_id], dat[aut_asym_id])
    #                 pdb[(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
    #                     numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
    #             else:
    #                 aid[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
    #                     (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
    #                      dat[alt_id], dat[icode_id], dat[aut_asym_id])
    #                 pdb[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
    #                     numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
    #     pdb_models[model] = pdb
    #     atom_ids[model] = aid
    # return pdb_models

if __name__ == "__main__":
    get_dssp_ss('/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/12gs.cif.gz')