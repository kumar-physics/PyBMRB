import plotly.express as px
import json
import logging
from urllib.request import urlopen, Request
import numpy
from typing import Union, List, Optional



# url = Request("http://api.bmrb.io/v2/search/chemical_shifts?atom_id=*&database=metabolomics")
# #url = Request("http://api.bmrb.io/v2/search/chemical_shifts?atom_id=*")#&database=macromolecule")
# url.add_header('Application', 'PyBMRB')
# r = urlopen(url)
# dump = json.loads(r.read())
# atom = dump['columns'].index('Atom_chem_shift.Atom_type')
# res = dump['columns'].index('Atom_chem_shift.Comp_ID')
# atom_type={}
#
# for i in dump['data']:
#     if i[atom] not in atom_type:
#         atom_type[i[atom]]=0
#     atom_type[i[atom]]+=1
# print (atom_type)
atom_type={}
atom_type['<sup>13</sup>C']=39622
atom_type['<sup>1</sup>H']=27285
atom_type['<sup>15</sup>N']=1
atom_type['<sup>31</sup>P']=0

proteins = {
    '<sup>13</sup>C' : 3801070,
    '<sup>1</sup>H' : 5836090,
    '<sup>15</sup>N' : 1230986,
    '<sup>31</sup>P' : 42,

}
na = {
    '<sup>13</sup>C':92273,
    '<sup>1</sup>H':240588,
    '<sup>15</sup>N':27606,
    '<sup>31</sup>P':3162,

}
fontsize=16
fontfamily="Arial"
bmrbdata={'Proteins':proteins,'Nucleic acids':na,'Metabolomics':atom_type}
print (bmrbdata)
bmrb=[]
mol_type=[]
atomtype=[]
value=[]
for k1 in bmrbdata:
    for k2 in bmrbdata[k1]:
        bmrb.append('BMRB')
        mol_type.append(k1)
        atomtype.append(k2)
        value.append(bmrbdata[k1][k2])
fig = px.sunburst(path=[bmrb,mol_type,atomtype],values=value,color=atomtype)
fig.update_traces(textinfo="label+value+percent parent", selector=dict(type='sunburst'))
#fig.update_traces(insidetextfont_size=20, selector=dict(type='sunburst'))
fig.update_traces(textfont_size=fontsize, textfont_family=fontfamily,selector=dict(type='sunburst'))
fig.write_image('All.pdf')
fig.write_html('all.html')


mol1=[]
atm1=[]
v1=[]
for k in proteins:
    mol1.append('Proteins')
    atm1.append(k)
    v1.append(proteins[k])
fig1 = px.sunburst(path=[mol1,atm1],values=v1,color=atm1)
fig1.update_traces(textinfo="label+value+percent parent", selector=dict(type='sunburst'))
#fig1.update_traces(insidetextfont_size=20, selector=dict(type='sunburst'))
fig1.update_traces(insidetextorientation='horizontal', selector=dict(type='sunburst'))
fig1.update_traces(textfont_size=fontsize, textfont_family=fontfamily,selector=dict(type='sunburst'))
fig1.write_image('proteins.pdf')
mol2=[]
atm2=[]
v2=[]
for k in na:
    mol2.append('Nucleic acids')
    atm2.append(k)
    v2.append(na[k])
fig2 = px.sunburst(path=[mol2,atm2],values=v2,color=atm2)
fig2.update_traces(textinfo="label+value+percent parent", selector=dict(type='sunburst'))
#fig2.update_traces(insidetextfont_size=20, selector=dict(type='sunburst'))
fig2.update_traces(textfont_size=fontsize, textfont_family=fontfamily,selector=dict(type='sunburst'))
fig2.update_traces(insidetextorientation='horizontal', selector=dict(type='sunburst'))
fig2.write_image('nucleicacids.pdf')
mol3=[]
atm3=[]
v3=[]
for k in atom_type:
    mol3.append('Metabolomics')
    atm3.append(k)
    v3.append(atom_type[k])
fig3 = px.sunburst(path=[mol3,atm3],values=v3,color=atm3)
fig3.update_traces(textinfo="label+value+percent parent", selector=dict(type='sunburst'))
fig3.update_traces(insidetextorientation='horizontal', selector=dict(type='sunburst'))
fig3.update_traces(rotation=160, selector=dict(type='sunburst'))
fig3.update_traces(textfont_size=fontsize, textfont_family=fontfamily,selector=dict(type='sunburst'))
fig3.write_image('metaboloites.pdf')

mol4=['BMRB','BMRB','BMRB']
atm4=['Proteins','Nucleic acids','Metabolomics']
v4=[10868188,363929,66908]
for k in atom_type:
    mol3.append('Metabolomics')
    atm3.append(k)
    v3.append(atom_type[k])
fig4 = px.sunburst(path=[mol4,atm4],values=v4,color=atm4)
fig4.update_traces(textinfo="label+value+percent parent", selector=dict(type='sunburst'))
fig4.update_traces(insidetextorientation='horizontal', selector=dict(type='sunburst'))
fig4.update_traces(rotation=160, selector=dict(type='sunburst'))
fig4.update_traces(textfont_size=fontsize, textfont_family=fontfamily,selector=dict(type='sunburst'))
fig4.write_image('bmrb1.pdf')