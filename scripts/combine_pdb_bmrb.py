import json
import sys
import logging
from urllib.request import urlopen, Request
import gzip
from typing import Union, List, Optional
from mmcif.io.PdbxReader import PdbxReader
# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)
import pynmrstar
import os.path
import csv
#import plotly.express as px

# _API_URL = "http://dev-api.bmrb.io/v2"
_API_URL = "http://api.bmrb.io/v2"
_PDB_BMRB_MAPPING = "/mappings/bmrb/pdb?format=json&match_type=exact"
_FTP_BMRB_PATH = "/projects/BMRB/public/ftp/pub/bmrb/entry_directories"
_FTP_PDB_PATH = "/projects/BMRB/public/ftp/pub/pdb/data/structures/divided/mmCIF"
_FTP_VAL_PATH = "/projects/BMRB/public/ftp/pub/pdb/validation_reports"
_REBOXITORY_CIF = "/reboxitory/2024/04/PDB/data/structures/all/mmCIF"
_REBOXITORY_STR = "reboxitory/2024/04/BMRB/macromolecules"
three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


def _get_bmrb_pdb_mapping():
    url = Request(_API_URL+_PDB_BMRB_MAPPING)
    #url = "https://bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/adit_nmr_matched_pdb_bmrb_entry_ids.csv"
    url.add_header('Application', 'PyBMRB')
    r = urlopen(url)
    dump = json.loads(r.read())
    # for i in dump:
    #     print (i['bmrb_id'],i['pdb_ids'])
    return dump


def get_cs_data(str_file):
    try:
        ent = pynmrstar.Entry.from_file(str_file)
        cs_loop=ent.get_loops_by_category('Atom_chem_shift')
        cs_data={}
        for cs in cs_loop:
            col_names=cs.get_tag_names()
            auth_asym_idx = col_names.index('_Atom_chem_shift.Auth_asym_ID')
            auth_seq_idx = col_names.index('_Atom_chem_shift.Auth_seq_ID')
            seq_idx = col_names.index('_Atom_chem_shift.Comp_index_ID')
            auth_comp_idx = col_names.index('_Atom_chem_shift.Auth_comp_ID')
            comp_idx = col_names.index('_Atom_chem_shift.Comp_ID')
            atom_idx = col_names.index('_Atom_chem_shift.Atom_ID')
            cs_idx = col_names.index('_Atom_chem_shift.Val')
            list_idx = col_names.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')
            for row in cs.data:
                if row[list_idx] not in cs_data:
                    cs_data[row[list_idx]]={}
                if row[auth_seq_idx] == '.':
                    auth_seq_id = int(row[seq_idx])
                else:
                    try:
                        auth_seq_id = int(row[auth_seq_idx])
                    except ValueError:
                        auth_seq_id = int(row[seq_idx])
                if row[auth_asym_idx] == '.':
                    chain_id = 'A'
                else:
                    chain_id = row[auth_asym_idx]
                kk = (chain_id,int(row[seq_idx]),auth_seq_id,row[comp_idx],row[atom_idx])
                cs_data[row[list_idx]][kk]= float(row[cs_idx])
    except FileNotFoundError:
        cs_data={}
    return cs_data


def get_dssp_ss(cif_file):
    cif_data = []
    try:
        if cif_file.endswith(".gz"):
            ifh = gzip.open(cif_file,'rt')
        else:
            ifh = open(cif_file, 'r')
        pRd = PdbxReader(ifh)
        pRd.read(cif_data)
        ifh.close()
        c0 = cif_data[0]
        struct_conf = c0.getObj('struct_conf')
        entity_poly_seq = c0.getObj('entity_poly_seq')
        entity_poly = c0.getObj('entity_poly')
        try:
            col_names3 = entity_poly.getAttributeList()
            entity_id_idx = col_names3.index('entity_id',0)
            strand_id_idx = col_names3.index('pdbx_strand_id',0)
            entities = {}
            for dat in entity_poly.getRowList():
                entities[dat[entity_id_idx]]=dat[strand_id_idx].split(",")
            try:
                col_names = struct_conf.getAttributeList()
                col_names2 = entity_poly_seq.getAttributeList()
                poly_seq_num_idx = col_names2.index('num')
                poly_seq_mono_idx = col_names2.index('mon_id')
                poly_seq_entity_idx = col_names2.index('entity_id')
                sequence = {}
                for dat in entity_poly_seq.getRowList():
                    for c in entities[dat[poly_seq_entity_idx]]:
                        if c not in sequence:
                            sequence[c]={}
                    for c in entities[dat[poly_seq_entity_idx]]:
                        sequence[c][dat[poly_seq_num_idx]]=dat[poly_seq_mono_idx]
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
                ss_info={}
                ss_info['seq_id']={}
                ss_info['auth_seq_id']={}
                for dat in struct_conf.getRowList():
                    conf_type = dat[conf_type_idx]
                    beg_auth_seq_id = dat[beg_auth_seq_idx]
                    beg_auth_asym_id = dat[beg_auth_asym_idx]
                    beg_auth_comp_id = dat[beg_auth_comp_idx]
                    end_auth_seq_id = dat[end_auth_seq_idx]
                    end_auth_asym_id = dat[end_auth_asym_idx]
                    end_auth_comp_id = dat[end_auth_comp_idx]
                    beg_seq_id = dat[beg_seq_idx]
                    beg_asym_id = dat[beg_asym_idx]
                    beg_comp_id = dat[beg_comp_idx]
                    end_seq_id = dat[end_seq_idx]
                    end_asym_id = dat[end_asym_idx]
                    end_comp_id = dat[end_comp_idx]
                    ss_info['seq_id'][(beg_asym_id,beg_seq_id,beg_comp_id,end_asym_id,end_seq_id,end_asym_id)]=conf_type
                    ss_info['auth_seq_id'][(beg_auth_asym_id,beg_auth_seq_id,beg_auth_comp_id,end_auth_asym_id,end_auth_seq_id,end_auth_comp_id)]=conf_type
                ss={}
                for k1,k2 in zip(ss_info['seq_id'],ss_info['auth_seq_id']):
                    if int(k1[1]) != int(k2[1]):
                        offset = int(k1[1])-int(k2[1])
                    else:
                        offset = 0
                    for i in range(int(k1[1]),int(k1[4])+1):
                        try:
                            ss[(k1[0],i,i-offset,sequence[k1[0]][str(i)])] = ss_info['seq_id'][k1]
                        except KeyError:
                            print (k1)
                for k in sequence:
                    for i in sequence[k]:
                        kk = (k,int(i),int(i)-offset,sequence[k][i])
                        if kk not in ss:
                            ss[kk] = 'COIL'
            except AttributeError:
                ss={}
        except AttributeError:
            ss={}
    except FileNotFoundError:
        ss={}
    return ss
