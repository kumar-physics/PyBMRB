import csv
import plotly.express as px

standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                'ALA', 'MET', 'LEU', 'ARG', 'TYR',
                'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
def plot_cs_ss_data(fname,a):
    with open(fname, mode = 'r') as file:
        cs=[]
        res=[]
        atom=[]
        ss=[]
        csvFile = csv.reader(file)
        for lines in csvFile:
            if lines[1]==a and lines[0] in standard and 20.0 >= float(lines[2]) >0.0:
                res.append(lines[0])
                cs.append(float(lines[2]))
                ss.append(lines[3])
    print (len(cs))
    fig = px.histogram(x=cs,color=ss,barmode='overlay',opacity=0.3,facet_col=res,facet_col_wrap=5,labels={'x':'Chemical Shift [ppm]'},title=f'Chemical shift distribution of {a}')
    fig.update_yaxes(matches=None,)
    fig.update_xaxes( autorange='reversed')
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.show()
    fig2 = px.histogram(x=cs, color=ss, barmode='overlay', opacity=0.3, labels={'x':'Chemical Shift [ppm]'},title=f'Chemical shift distribution of {a}')
    fig2.update_xaxes(autorange='reversed')
    fig2.show()
    fig3 = px.histogram(x=cs, barmode='overlay', opacity=0.3, labels={'x': 'Chemical Shift [ppm]'},
                        title=f'Chemical shift distribution of {a}')
    fig3.update_xaxes(autorange='reversed')
    fig3.show()
    fig4 = px.histogram(x=cs, color=res,barmode='overlay', opacity=0.3, labels={'x': 'Chemical Shift [ppm]'},
                        title=f'Chemical shift distribution of {a}')
    fig4.update_xaxes(autorange='reversed')
    fig4.show()
    fig5 = px.histogram(x=res)
    fig5.show()

if __name__ == "__main__":
    plot_cs_ss_data('./ss_cs.csv','H')