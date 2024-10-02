import json
import logging
import multiprocessing
from urllib.request import  urlopen,Request
import gzip
from mmcif.io.PdbxReader import PdbxReader
logging.getLogger().setLevel(logging.ERROR)
import pynmrstar
import sys

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
    url.add_header('Application', 'PyBMRB')
    r = urlopen(url)
    dump = json.loads(r.read())
    return dump

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
        exptl = c0.getObj('exptl')
        if 'NMR' not in exptl.getValue('method'):
            return {} #only NMR structures
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
                        sequence[c][int(dat[poly_seq_num_idx])]=dat[poly_seq_mono_idx]
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
                    ss_info['seq_id'][(beg_asym_id,beg_seq_id,beg_comp_id,end_asym_id,end_seq_id,end_comp_id)]=conf_type
                    ss_info['auth_seq_id'][(beg_auth_asym_id,beg_auth_seq_id,beg_auth_comp_id,end_auth_asym_id,end_auth_seq_id,end_auth_comp_id)]=conf_type
                try:
                    struct_sheet = c0.getObj('struct_sheet_range')
                    col_names = struct_sheet.getAttributeList()
                    #conf_type_idx = col_names.index('conf_type_id')
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
                    for dat in struct_sheet.getRowList():
                        conf_type = "SHEET"
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
                        ss_info['seq_id'][(beg_asym_id,beg_seq_id,beg_comp_id,end_asym_id,end_seq_id,end_comp_id)]=conf_type
                        ss_info['auth_seq_id'][(beg_auth_asym_id,beg_auth_seq_id,beg_auth_comp_id,end_auth_asym_id,end_auth_seq_id,end_auth_comp_id)]=conf_type
                except AttributeError:
                    pass
                ss={}
                for k1,k2 in zip(ss_info['seq_id'],ss_info['auth_seq_id']):
                    if int(k1[1]) != int(k2[1]):
                        offset = int(k1[1])-int(k2[1])
                    else:
                        offset = 0
                    for i in range(int(k1[1]),int(k1[4])+1):
                        try:
                            ss[(k2[0],i,i-offset,sequence[k2[0]][i])] = ss_info['seq_id'][k1]
                        except KeyError:
                            logging.warning(f'Key not found at instant a {k1}')
                for k in sequence:
                    for i in sequence[k]:
                        kk = (k,int(i),int(i)-offset,sequence[k][i])
                        if kk not in ss:
                            ss[kk] = 'COIL'
            except AttributeError:
                logging.info(f'No Struct_conf info found in {cif_file}; probably only beta sheets ')
                try:
                    struct_sheet = c0.getObj('struct_sheet_range')
                    entity_poly_seq = c0.getObj('entity_poly_seq')
                    entity_poly = c0.getObj('entity_poly')
                    try:
                        col_names3 = entity_poly.getAttributeList()
                        entity_id_idx = col_names3.index('entity_id', 0)
                        strand_id_idx = col_names3.index('pdbx_strand_id', 0)
                        entities = {}
                        for dat in entity_poly.getRowList():
                            entities[dat[entity_id_idx]] = dat[strand_id_idx].split(",")
                        try:
                            col_names = struct_sheet.getAttributeList()
                            col_names2 = entity_poly_seq.getAttributeList()
                            poly_seq_num_idx = col_names2.index('num')
                            poly_seq_mono_idx = col_names2.index('mon_id')
                            poly_seq_entity_idx = col_names2.index('entity_id')
                            sequence = {}
                            for dat in entity_poly_seq.getRowList():
                                for c in entities[dat[poly_seq_entity_idx]]:
                                    if c not in sequence:
                                        sequence[c] = {}
                                for c in entities[dat[poly_seq_entity_idx]]:
                                    sequence[c][int(dat[poly_seq_num_idx])] = dat[poly_seq_mono_idx]
                            col_names = struct_sheet.getAttributeList()
                            # conf_type_idx = col_names.index('conf_type_id')
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
                            ss_info = {}
                            ss_info['seq_id'] = {}
                            ss_info['auth_seq_id'] = {}
                            for dat in struct_sheet.getRowList():
                                conf_type = "SHEET"
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
                                ss_info['seq_id'][
                                    (beg_asym_id, beg_seq_id, beg_comp_id, end_asym_id, end_seq_id,
                                     end_comp_id)] = conf_type
                                ss_info['auth_seq_id'][(
                                    beg_auth_asym_id, beg_auth_seq_id, beg_auth_comp_id, end_auth_asym_id,
                                    end_auth_seq_id,
                                    end_auth_comp_id)] = conf_type
                            ss = {}
                            for k1, k2 in zip(ss_info['seq_id'], ss_info['auth_seq_id']):
                                if int(k1[1]) != int(k2[1]):
                                    offset = int(k1[1]) - int(k2[1])
                                else:
                                    offset = 0
                                for i in range(int(k1[1]), int(k1[4]) + 1):
                                    try:
                                        ss[(k2[0], i, i - offset, sequence[k2[0]][i])] = ss_info['seq_id'][k1]
                                    except KeyError:
                                        logging.warning(f'Key not found at instant b {k1}')
                            for k in sequence:
                                for i in sequence[k]:
                                    kk = (k, int(i), int(i) - offset, sequence[k][i])
                                    if kk not in ss:
                                        ss[kk] = 'COIL'
                        except AttributeError:
                            ss = {}
                    except AttributeError:
                        logging.warning(f'No DSSP information found {cif_file}')
                        ss = {}
                except AttributeError:
                    logging.warning(f'No DSSP information found {cif_file}')
                    ss = {}
        except AttributeError:
            logging.warning(f'Entity information missing in file {cif_file}')
            ss={}
    except FileNotFoundError:
        logging.warning(f'File not found {cif_file}')
        ss={}
    ss2={}
    for k in ss:
        if k[0] not in ss2:
            ss2[k[0]]={}
        ss2[k[0]][(k[1],k[3])]=ss[k]
    return ss2

def get_cs_data(str_file):
    try:
        ent = pynmrstar.Entry.from_file(str_file)
        cs_loop=ent.get_loops_by_category('Atom_chem_shift')
        cs_data={}
        for cs in cs_loop:
            col_names=cs.get_tag_names()
            seq_idx = col_names.index('_Atom_chem_shift.Comp_index_ID')
            entity_assembly_idx = col_names.index('_Atom_chem_shift.Entity_assembly_ID')
            comp_idx = col_names.index('_Atom_chem_shift.Comp_ID')
            atom_idx = col_names.index('_Atom_chem_shift.Atom_ID')
            cs_idx = col_names.index('_Atom_chem_shift.Val')
            list_idx = col_names.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')
            for row in cs.data:
                if row[list_idx] not in cs_data:
                    cs_data[row[list_idx]]={}
                atom_identifier = (row[entity_assembly_idx],int(row[seq_idx]),row[comp_idx],row[atom_idx])
                cs_data[row[list_idx]][atom_identifier]= float(row[cs_idx])
    except FileNotFoundError:
        logging.warning(f'FIle not found {str_file}')
        cs_data={}
    cs_data2={}
    for cs_list in cs_data:
        cs_data2[cs_list]={}
        for k in cs_data[cs_list]:
            if k[0] not in cs_data2[cs_list]:
                cs_data2[cs_list][k[0]]={}
            cs_data2[cs_list][k[0]][(k[1],k[2],k[3])]=cs_data[cs_list][k]
    return cs_data2

def match_chains(cs_data,ss_data):
    map={}
    for ss_chain in ss_data:
        ss_seq = list(ss_data[ss_chain].keys())
        for cs_list in cs_data:
            if ss_chain not in map: map[ss_chain]=[]
            for cs_chain in cs_data[cs_list]:
                cs_seq = list(set([(i[0],i[1]) for i in list(cs_data[cs_list][cs_chain].keys())]))
                match_value, offset = find_matching_and_offset(ss_seq,cs_seq)
                if match_value > 0.7:
                    map[ss_chain].append([cs_list,cs_chain,match_value,offset])
    return map

def find_matching_and_offset(seq1,seq2):
    offset = 0
    union = len(seq1)+len(seq2)
    n= [i[0] for i in seq1]
    for i in range(min(n)-len(n),max(n)):
        shifted_seq = [(j[0]+i,j[1]) for j in seq2]
        common_elements = list(set(shifted_seq+seq1))
        if len(common_elements)<= union:
            union = len(common_elements)
            offset = i
            matched_seq = shifted_seq
    match_count = 0
    for k in matched_seq:
        if k in seq1:
            match_count+=1
    match_value = float(match_count)/float(len(seq1))
    match_value2 = float(match_count)/float(len(seq2))
    m_v = max(match_value,match_value2)
    return m_v,offset

def merge_cs_ss(input_args):
    pdb = input_args[0]
    bmrb = input_args[1]
    out_dir = input_args[2]
    logging.info(f'Merging {pdb},{bmrb}')
    #pdb= pair.split("-")[0]
    #bmrb = f'bmr{pair.split("-")[1]}'
    #pdb_file = _REBOXITORY_CIF+f'/{pdb}.cif.gz'
    #bmrb_file = _REBOXITORY_STR+f'/{bmrb}/{bmrb}_3.str'
    pdb_file = _FTP_PDB_PATH+f'/{pdb[1]}{pdb[2]}/{pdb}.cif.gz'
    bmrb_file = _FTP_BMRB_PATH+f'/bmr{bmrb}/bmr{bmrb}_3.str'
    #pdb_file = f'/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/{pdb}.cif'
    #bmrb_file = f'/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/bmr{bmrb}_3.str'
    ss_data = get_dssp_ss(pdb_file)
    err=''
    msg=''
    if len(ss_data)>0:
        cs_data = get_cs_data(bmrb_file)
        if len(cs_data)>0:
            unique_cs_rows=[]
            fo = open(f'{out_dir}/{bmrb}_{pdb}.csv', 'w')
            map=match_chains(cs_data,ss_data)
            for k in map:
                if len(map[k])==0:
                    logging.warning(f'No matching chain found {pdb} {bmrb}')
                for cs_list in map[k]:
                    for atm in cs_data[cs_list[0]][cs_list[1]]:
                        try:
                            if f'{atm[0]},{atm[1]},{atm[2]},{cs_list[0]},{cs_list[1]},{cs_data[cs_list[0]][cs_list[1]][atm]},{ss_data[k][(atm[0]+cs_list[3],atm[1])]},{pdb},{bmrb}' not in unique_cs_rows:
                                fo.write(f'{atm[0]},{atm[1]},{atm[2]},{cs_list[0]},{cs_list[1]},{k},{cs_data[cs_list[0]][cs_list[1]][atm]},{ss_data[k][(atm[0]+cs_list[3],atm[1])]},{pdb},{bmrb}\n')
                                unique_cs_rows.append(f'{atm[0]},{atm[1]},{atm[2]},{cs_list[0]},{cs_list[1]},{cs_data[cs_list[0]][cs_list[1]][atm]},{ss_data[k][(atm[0]+cs_list[3],atm[1])]},{pdb},{bmrb}')
                        except KeyError:
                            logging.warning(f'No matching residue  found {bmrb},{pdb},{k},{cs_list},{atm}')

            fo.close()


if __name__ == "__main__":
    out_dir = sys.argv[1]
    pair_list = _get_bmrb_pdb_mapping()
    input_list = []
    for k in pair_list:
        bmrb = k['bmrb_id']
        for pdb in k['pdb_ids']:
            if pdb.lower() not in ['1dey','1ugt']:
                input_list.append((pdb.lower(),bmrb,out_dir))
    pool = multiprocessing.Pool()
    pool.map(merge_cs_ss,input_list)
