#!/usr/bin/env python3
"""
Short description TBD
"""
import json
import logging
from urllib.request import urlopen, Request
import gzip
from typing import Union, List, Optional
from mmcif.io.PdbxReader import PdbxReader
# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)
import pynmrstar

# _API_URL = "http://dev-api.bmrb.io/v2"
_API_URL = "http://api.bmrb.io/v2"
_PDB_BMRB_MAPPING = "/mappings/bmrb/pdb?format=json&match_type=exact"
_FTP_BMRB_PATH = "/projects/BMRB/public/ftp/pub/bmrb/entry_directories"
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
    # for i in dump:
    #     print (i['bmrb_id'],i['pdb_ids'])
    return dump

def merge_cs_ss(pdb,bmrb):
    #pdb_file = _FTP_PDB_PATH+f'/{pdb[1]}{pdb[2]}/{pdb}.cif.gz'
    #bmrb_file = _FTP_BMRB_PATH+f'/bmr{bmrb}/bmr{bmrb}_3.str'
    pdb_file = f'/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/{pdb}.cif'
    bmrb_file = f'/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/bmr{bmrb}_3.str'
    ss_data = get_dssp_ss(pdb_file)
    err=''
    if len(ss_data)>0:
        cs_data = get_cs_data(bmrb_file)
        #fo=open(f'./cs_ss_out/{bmrb}_{pdb}.csv','w')
        fo = open(f'{bmrb}_{pdb}.csv', 'w')
        for cs_list in cs_data:
            for row in cs_data[cs_list]:
                #print (cs_data[cs_list][row],ss_data[(row[0],row[1],row[2],row[3])])
                try:
                    fo.write(f'{row[3]},{row[4]},{cs_data[cs_list][row]},{ss_data[(row[0],row[1],row[2],row[3])]},{row[1]},{row[2]},{row[0]},{pdb},{bmrb}\n')
                except KeyError:
                    err+=f'Missing atom {(row[0],row[1],row[2],row[3])} {pdb},{bmrb}\n'
        fo.close()
        msg = f"Success {pdb},{bmrb}"
    else:
        msg = f"No DSSP SS information found in CIF file {pdb},{bmrb}"
    return msg,err


def get_cs_data(str_file):
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
                auth_seq_id = int(row[auth_seq_idx])
            if row[auth_asym_idx] == '.':
                chain_id = 'A'
            else:
                chain_id = row[auth_asym_idx]
            kk = (chain_id,int(row[seq_idx]),auth_seq_id,row[comp_idx],row[atom_idx])
            cs_data[row[list_idx]][kk]= float(row[cs_idx])
    return cs_data



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
    entity_poly_seq = c0.getObj('entity_poly_seq')
    entity_poly = c0.getObj('entity_poly')
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
    return ss


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
    msg,err = merge_cs_ss('1nk2','4141')
    print (msg)
    print (err)
    # pair_list = _get_bmrb_pdb_mapping()
    # f=open('running_log.txt','w')
    # f1=open('running_err.txt','w')
    # for k in pair_list:
    #     bmrb = k['bmrb_id']
    #     for pdb in k['pdb_ids']:
    #         msg,err = merge_cs_ss(pdb.lower(),bmrb)
    #         f.write(f'{msg}\n')
    #         f1.write(f'{err}\n')
    # f.close()
    # f1.close()
    #####################
    #merge_cs_ss('1k8j','5716')
    #cs_data = get_cs_data('/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/bmr36445_3.str')
    #_get_bmrb_pdb_mapping()
    #ss_data=get_dssp_ss('/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/7VH9.cif')
    #4141 1nk2
