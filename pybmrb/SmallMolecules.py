#!/usr/bin/env python3
"""
This module extracts chemical shift information from BMRB small molecules entries, NMR-STAR files and PyNMRSTAR entry objects
"""
from __future__ import print_function
import logging
import os
import pynmrstar
from typing import Union, List, Optional

# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)


three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


def _from_pynmrstar_entry_object(entry_data: pynmrstar.Entry,
                                 data_set_id: str,
                                 auth_tag: Optional[bool] = False, ) -> dict:
    """
    Extracts chemical shift data as dictionary from PyNMRSTAR entry object

    :param entry_data: PyNMRSTAR entry object
    :type entry_data: pynmrstar.Entry
    :param data_set_id: Data set identifier
    :type data_set_id: str
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :return: Chemical shift dictionary {data_set_id:{chain_id:{seq_id:{atom_id:cs_value}},'seq_ids':[1,2,3,4..]}}
    :rtype: dict
    """
    cs_data = {}
    bond_info = {}
    # entry_data = pynmrstar.Entry.from_database(bmrb_id)
    cs_loops = entry_data.get_loops_by_category('Atom_chem_shift')
    logging.debug('Chemical shift loop read')
    cs_list_id = 1
    bond_loops = entry_data.get_loops_by_category('Chem_comp_bond')

    for bond_loop in bond_loops:
        bond_info[data_set_id]=[]
        col_names = bond_loop.get_tag_names()
        bond=col_names.index('_Chem_comp_bond.Value_order')
        atom1 = col_names.index('_Chem_comp_bond.Atom_ID_1')
        atom2 = col_names.index('_Chem_comp_bond.Atom_ID_2')
        chsqc_atoms=[]
        for row in bond_loop.data:
            if row[bond]=='SING' and (row[atom1][0] in ['C','H'] and row[atom2][0] in ['C','H']) and (row[atom1][0]!=row[atom2][0]):
                chsqc_atoms.append((row[atom1],row[atom2]))
        bond_info[data_set_id].append(chsqc_atoms)
    for cs_loop in cs_loops:
        cs_id = '{}-{}'.format(data_set_id, cs_list_id)
        cs_list_id += 1
        cs_data[cs_id] = {}
        column_name = cs_loop.get_tag_names()
        data = cs_loop.data
        if auth_tag:
            # chain = column_name.index('_Atom_chem_shift.Auth_asym_ID')
            seq_no = column_name.index('_Atom_chem_shift.Auth_seq_ID')
            # residue = column_name.index('_Atom_chem_shift.Auth_comp_ID')
            # atom = column_name.index('_Atom_chem_shift.Auth_atom_ID')
        else:
            # chain = column_name.index('_Atom_chem_shift.Entity_assembly_ID')
            seq_no = column_name.index('_Atom_chem_shift.Comp_index_ID')
            # residue = column_name.index('_Atom_chem_shift.Comp_ID')
            # atom = column_name.index('_Atom_chem_shift.Atom_ID')
        residue = column_name.index('_Atom_chem_shift.Comp_ID')
        atom = column_name.index('_Atom_chem_shift.Atom_ID')
        try:
            chain = column_name.index('_Atom_chem_shift.Entity_assembly_ID')
        except ValueError:
            chain = False
        cs_value = column_name.index('_Atom_chem_shift.Val')
        for cs_row in data:
            if chain == False:
                c=1
            else:
                c = cs_row[chain]
            s = int(cs_row[seq_no])
            r = cs_row[residue]
            a = cs_row[atom]
            v = float(cs_row[cs_value])
            if c not in cs_data[cs_id].keys():
                cs_data[cs_id][c] = {}
            if s not in cs_data[cs_id][c].keys():
                cs_data[cs_id][c][s] = {}
            if a not in cs_data[cs_id][c][s].keys():
                cs_data[cs_id][c][s][a] = []
            cs_data[cs_id][c][s][a].append(r)
            try:
                cs_data[cs_id][c][s][a].append(one_letter_code[r])
            except KeyError:
                cs_data[cs_id][c][s][a].append(r)
            # cs_data[cs_id][c][s][a].append(a)
            cs_data[cs_id][c][s][a].append(v)
    logging.debug('Getting sequence information from chemical shift loop')
    for cs_id in cs_data.keys():
        for chain in cs_data[cs_id].keys():
            seq_ids = sorted(set([i for i in cs_data[cs_id][chain].keys()]))
            cs_data[cs_id][chain]['seq_ids'] = seq_ids
    return cs_data,bond_info



def from_bmrb(bmrb_ids: Union[str, List[str], int, List[int]],
              auth_tag: Optional[bool] = False) -> dict:
    """
    Extracts chemical shift information directly from BMRB database for a given BMRB entry or list of entries

    :param bmrb_ids: List of BMRB entries ids or single BMRB ID
    :type bmrb_ids: str/int/list
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :return: Chemical shift dictionary {data_set_id:{chain_id:{seq_id:{atom_id:cs_value}},'seq_ids':[1,2,3,4..]}}
    :rtype: dict
    """

    if type(bmrb_ids) is list:
        all_cs_data = {}
        all_chsqc_atoms = {}
        for bmrb_id in bmrb_ids:
            try:
                logging.debug('Getting entry {} from BMRB'.format(bmrb_id))
                entry_data = pynmrstar.Entry.from_database(bmrb_id)
            except OSError:
                entry_data = None
            except KeyError:
                entry_data = None
            if entry_data is not None:
                cs_data,chsqc_atoms = _from_pynmrstar_entry_object(entry_data, bmrb_id, auth_tag)
            else:

                logging.error('Entry {} not found in public database'.format(bmrb_id))
                raise IOError('Entry not found in public database: {}'.format(bmrb_id))
            all_cs_data.update(cs_data)
            all_chsqc_atoms.update(chsqc_atoms)
    else:

        try:
            logging.debug('Getting entry {} from BMRB'.format(bmrb_ids))
            entry_data = pynmrstar.Entry.from_database(bmrb_ids)
        except OSError:
            entry_data = None
        except KeyError:
            entry_data = None
        if entry_data is not None:
            all_cs_data,chsqc_atoms = _from_pynmrstar_entry_object(entry_data, bmrb_ids, auth_tag)
        else:
            logging.error('Entry {} not found in public database'.format(bmrb_ids))
            raise IOError('Entry not found in public database: {}'.format(bmrb_ids))
    return all_cs_data,all_chsqc_atoms


def generate_spectra(bmrbids):
    cs_data,chsqc_atoms = from_bmrb(bmrbids)
    for dataset in cs_data:
        bid = dataset.split("-")[0]
        for atom in



if __name__ == "__main__":
    # x=from_bmrb(['bmse010240','bmse000462'])
    # print (x)
    generate_spectra(['bmse010240','bmse000462'])