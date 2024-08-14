import csv
import plotly.express as px

standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                'ALA', 'MET', 'LEU', 'ARG', 'TYR',
               ]
def read_csv(fname):
    with open('bmrb_stat.csv','r') as csfile:
        atm={}
        m={}
        sd={}
        cs = csv.reader(csfile)
        for l in cs:
            if l[0] in standard:
                atm[(l[0],l[1])]=[[],[]]
                m[(l[0],l[1])]=float(l[5])
                sd[(l[0], l[1])] = float(l[6])
    with open(fname, mode='r') as file:
        cs = []
        ss = []
        csvFile = csv.reader(file)
        plot_data = {}
        plot_data[('x','H')]= [[],[]]
        for lines in csvFile:
            if lines[0] in standard :
                if (lines[0],lines[1]) not in plot_data:
                    plot_data[(lines[0],lines[1])]=[[],[]]
                plot_data[(lines[0], lines[1])][0].append(float(lines[2]))
                plot_data[(lines[0], lines[1])][1].append(lines[3])
                if lines[1] == 'H':
                    plot_data[('x','H')][0].append(float(lines[2]))
                    plot_data[('x', 'H')][1].append(lines[3])
            try:
                atm[(lines[0],lines[1])][0].append(float(lines[2]))
                atm[(lines[0], lines[1])][1].append(lines[3])
            except KeyError:
                if lines[0] == "ALA" and lines[1] in ['HB1','HB2','HB3']:
                    atm[(lines[0], 'MB')][0].append(float(lines[2]))
                    atm[(lines[0], 'MB')][1].append(lines[3])
                if lines[0] == "ILE" and lines[1] in ['HD11','HD12','HD13']:
                    atm[(lines[0], 'MD')][0].append(float(lines[2]))
                    atm[(lines[0], 'MD')][1].append(lines[3])
                if lines[0] == "ILE" and lines[1] in ['HG21','HG22','HG23']:
                    atm[(lines[0], 'MG')][0].append(float(lines[2]))
                    atm[(lines[0], 'MG')][1].append(lines[3])
                if lines[0] == "LEU" and lines[1] in ['HD11','HD12','HD13']:
                    atm[(lines[0], 'MD1')][0].append(float(lines[2]))
                    atm[(lines[0], 'MD1')][1].append(lines[3])
                if lines[0] == "LEU" and lines[1] in ['HD21','HG22','HG23']:
                    atm[(lines[0], 'MD2')][0].append(float(lines[2]))
                    atm[(lines[0], 'MD2')][1].append(lines[3])
                if lines[0] == "MET" and lines[1] in ['HE1','HE2','HE3']:
                    atm[(lines[0], 'ME')][0].append(float(lines[2]))
                    atm[(lines[0], 'ME')][1].append(lines[3])
                if lines[0] == "THR" and lines[1] in ['HG21','HG22','HG23']:
                    atm[(lines[0], 'MG')][0].append(float(lines[2]))
                    atm[(lines[0], 'MG')][1].append(lines[3])
                if lines[0] == "VAL" and lines[1] in ['HG11','HG12','HG13']:
                    atm[(lines[0], 'MG1')][0].append(float(lines[2]))
                    atm[(lines[0], 'MG1')][1].append(lines[3])
                if lines[0] == "VAL" and lines[1] in ['HG21','HG22','HG23']:
                    atm[(lines[0], 'MG2')][0].append(float(lines[2]))
                    atm[(lines[0], 'MG2')][1].append(lines[3])


    for k in atm:
        if len(atm[k][0]) != 0:
            #print (k,len(atm[k][0]),m[k],sd[k])
            foname = f'./out/{k[0]}-{k[1]}_ss_unfiltered.html'
            fig = px.histogram(atm[k][0],color=atm[k][1],barmode='overlay',opacity=0.7,labels={'value':'Chemical shift (ppm)'})
            fig.update_xaxes(autorange="reversed")
            fig.write_html(foname)
            foname2 = f'./out/{k[0]}-{k[1]}_ss.html'
            filtered_data = cs_filter(atm[k],m[k],sd[k],8)
            fig2 = px.histogram(filtered_data[0], color=filtered_data[1], barmode='overlay', opacity=0.7,
                               labels={'value': 'Chemical shift (ppm)'})
            fig2.update_xaxes(autorange="reversed")
            fig2.write_html(foname2)

def cs_filter(d,m,sd,sd_limit=8):
    new_d=[[],[]]
    for i in range(len(d[0])):
        if m-(sd_limit*sd) <= d[0][i] <=  m+(sd_limit*sd):
            new_d[0].append(d[0][i])
            new_d[1].append(d[1][i])
    return new_d

if __name__ == "__main__":
    read_csv('./ss_cs.csv')
