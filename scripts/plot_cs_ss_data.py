import csv
import plotly.express as px
import os
from statistics import mean,stdev
import sys

standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                'ALA', 'MET', 'LEU', 'ARG', 'TYR',
               ]
def plot_csv(input_dir,out_dir):
    sd_filter=8
    full_cs={}
    for csv_file_name in os.listdir(input_dir):
        if csv_file_name.endswith(".csv"):
            csv_file=f'{input_dir}/{csv_file_name}'
            with open(csv_file, mode='r') as csvFile:
                cf = csv.reader(csv_file)
                for lines in cf:
                    if len(lines)>6:
                        if lines[1] in standard:
                            if lines[1] == "ALA" and lines[2] in ['HB1', 'HB2', 'HB3']:
                                if (lines[1], 'MB') not in full_cs:
                                    full_cs[(lines[1], 'MB')] = [[],[]]
                                full_cs[(lines[1], 'MB')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MB')][1].append(lines[7])
                            elif lines[1] == "ILE" and lines[2] in ['HD11', 'HD12', 'HD13']:
                                if (lines[1], 'MD') not in full_cs:
                                    full_cs[(lines[1], 'MD')] = [[],[]]
                                full_cs[(lines[1], 'MD')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MD')][1].append(lines[7])
                            elif lines[1] == "ILE" and lines[2] in ['HG21', 'HG22', 'HG23']:
                                if (lines[1], 'MG') not in full_cs:
                                    full_cs[(lines[1], 'MG')] = [[],[]]
                                full_cs[(lines[1], 'MG')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MG')][1].append(lines[7])
                            elif lines[1] == "LEU" and lines[2] in ['HD11', 'HD12', 'HD13']:
                                if (lines[1], 'MD1') not in full_cs:
                                    full_cs[(lines[1], 'MD1')] = [[],[]]
                                full_cs[(lines[1], 'MD1')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MD1')][1].append(lines[7])
                            elif lines[1] == "LEU" and lines[2] in ['HD21', 'HG22', 'HG23']:
                                if (lines[1], 'MD2') not in full_cs:
                                    full_cs[(lines[1], 'MD2')] = [[],[]]
                                full_cs[(lines[1], 'MD2')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MD2')][1].append(lines[7])
                            elif lines[1] == "MET" and lines[2] in ['HE1', 'HE2', 'HE3']:
                                if (lines[1], 'ME') not in full_cs:
                                    full_cs[(lines[1], 'ME')] = [[],[]]
                                full_cs[(lines[1], 'ME')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'ME')][1].append(lines[7])
                            elif lines[1] == "THR" and lines[2] in ['HG21', 'HG22', 'HG23']:
                                if (lines[1], 'MG') not in full_cs:
                                    full_cs[(lines[1], 'MG')] = [[],[]]
                                full_cs[(lines[1], 'MG')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MG')][1].append(lines[7])
                            elif lines[1] == "VAL" and lines[2] in ['HG11', 'HG12', 'HG13']:
                                if (lines[1], 'MG1') not in full_cs:
                                    full_cs[(lines[1], 'MG1')] = [[],[]]
                                full_cs[(lines[1], 'MG1')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MG1')][1].append(lines[7])
                            elif lines[1] == "VAL" and lines[2] in ['HG21', 'HG22', 'HG23']:
                                if (lines[1], 'MG2') not in full_cs:
                                    full_cs[(lines[1], 'MG2')] = [[],[]]
                                full_cs[(lines[1], 'MG2')][0].append(float(lines[6]))
                                full_cs[(lines[1], 'MG2')][1].append(lines[7])
                            else:
                                if (lines[1],lines[2]) not in full_cs:
                                    full_cs[(lines[1],lines[2])] = [[],[]]
                                else:
                                    full_cs[(lines[1],lines[2])][0].append(float(lines[6]))
                                    full_cs[(lines[1], lines[2])][1].append(lines[7])
    filtered_cs= {}
    for atom in full_cs:
        if len(full_cs[atom][0])>1:
            m = mean(full_cs[atom][0])
            sd = stdev(full_cs[atom][0])
            if atom not in filtered_cs:
                filtered_cs[atom] = [[],[]]
            for k in range(len(full_cs[atom][0])):
                if (m-(sd*sd_filter)) <= full_cs[atom][0][k] <= (m+(sd*sd_filter)):
                    filtered_cs[atom][0].append(full_cs[atom][0][k])
                    filtered_cs[atom][1].append(full_cs[atom][1][k])

    for k in full_cs:
        if len(full_cs[k][0]) != 0:
            foname = f'{out_dir}/{k[0]}-{k[1]}_ss_unfiltered.html'
            fig = px.histogram(full_cs[k][0],color=full_cs[k][1],barmode='overlay',opacity=0.7,labels={'value':'Chemical shift (ppm)'})
            fig.update_xaxes(autorange="reversed")
            fig.write_html(foname)
            if k in filtered_cs:
                foname2 = f'{out_dir}/{k[0]}-{k[1]}_ss.html'
                fig2 = px.histogram(filtered_cs[k][0], color=filtered_cs[k][1], barmode='overlay', opacity=0.7,
                                   labels={'value': 'Chemical shift (ppm)'})
                fig2.update_xaxes(autorange="reversed")
                fig2.write_html(foname2)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    plot_csv(input_dir,output_dir)
