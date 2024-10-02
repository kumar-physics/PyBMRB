import re

def parse_secondary_structure(cif_file):
    """
    Parse secondary structure information from a CIF file and assign it to each residue.

    Parameters:
    cif_file (str): Path to the CIF file.

    Returns:
    dict: A dictionary where each key is a tuple (chain_id, residue_number) and the value is the secondary structure type.
    """
    secondary_structure_residues = {}

    with open(cif_file, 'r') as file:
        lines = file.readlines()

    in_helix_section = False
    in_sheet_section = False

    for line in lines:
        # Identify helix information section
        if line.startswith('_struct_conf.'):
            in_helix_section = True
            in_sheet_section = False

        # Identify sheet information section
        if line.startswith('_struct_sheet_range.'):
            in_helix_section = False
            in_sheet_section = True

        # Process helix information
        if in_helix_section and not line.startswith('_struct_conf.'):
            parts = re.split(r'\s+', line.strip())
            if len(parts) >= 12:
                helix_info = {
                    'beg_label_asym_id': parts[3],
                    'beg_label_seq_id': int(parts[4]),
                    'end_label_asym_id': parts[7],
                    'end_label_seq_id': int(parts[8]),
                    'type': parts[2]  # helix type (alpha, 3-10, pi, etc.)
                }

                for res_num in range(helix_info['beg_label_seq_id'], helix_info['end_label_seq_id'] + 1):
                    key = (helix_info['beg_label_asym_id'], res_num)
                    secondary_structure_residues[key] = f"Helix ({helix_info['type']})"

        # Process sheet information
        if in_sheet_section and not line.startswith('_struct_sheet_range.'):
            parts = re.split(r'\s+', line.strip())
            if len(parts) >= 6:
                sheet_info = {
                    'beg_label_asym_id': parts[2],
                    'beg_label_seq_id': int(parts[3]),
                    'end_label_asym_id': parts[5],
                    'end_label_seq_id': int(parts[6]),
                }

                for res_num in range(sheet_info['beg_label_seq_id'], sheet_info['end_label_seq_id'] + 1):
                    key = (sheet_info['beg_label_asym_id'], res_num)
                    secondary_structure_residues[key] = "Sheet"

    return secondary_structure_residues

def write_secondary_structure_to_file(secondary_structure_residues, output_file):
    """
    Write secondary structure information for each residue to a file.

    Parameters:
    secondary_structure_residues (dict): A dictionary containing secondary structure information per residue.
    output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as file:
        file.write("Chain\tResidue\tSecondary Structure\n")
        for key, value in sorted(secondary_structure_residues.items()):
            chain_id, residue_number = key
            file.write(f"{chain_id}\t{residue_number}\t{value}\n")

# Example usage:
cif_file_path = '/Users/kumaranbaskaran/Projects/bmrb/PyBMRB/pybmrb/tests/test_data/1akk.cif'
output_file_path = 'secondary_structure_output.txt'

# Parse secondary structure
secondary_structure_info = parse_secondary_structure(cif_file_path)

# Write to file
write_secondary_structure_to_file(secondary_structure_info, output_file_path)

print(f"Secondary structure information written to {output_file_path}")
