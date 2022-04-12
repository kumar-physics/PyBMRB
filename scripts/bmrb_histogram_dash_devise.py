import dash
from dash import html
from dash import dcc
import plotly.express as px
from dash.dependencies import Input, Output
import pandas as pd
import numpy
external_stylesheets = [
    'https://bmrb.io/stylesheets/bmrb.css',
]
app = dash.Dash(__name__,
                external_stylesheets=external_stylesheets)
#app = dash.Dash()
amino_acids_list = ['ALA','CYS','ASP','GLU','PHE',
                                              'GLY','HIS','ILE','LYS','LEU',
                                              'MET','ASN','PRO','GLN','ARG',
                                              'SER','THR','VAL','TRP','TYR']
methyl={
    'ALA':[['HB1','HB2','HB3']],
    'ILE':[['HG21','HG22','HG23'],['HD11','HD12','HD13']],
    'LEU':[['HD11','HD12','HD13'],['HD21','HD22','HD23']],
    'MET':[['HE1','HE2','HE3']],
    'THR':[['HG21','HG22','HG23']],
    'VAL':[['HG11','HG12','HG13'],['HG21','HG22','HG23']]
}
methyl_pseudo={
    'ALA':['MB'],
    'ILE':['MG2','MD1'],
    'LEU':['MD1','MD2'],
    'MET':['ME'],
    'THR':['MG2'],
    'VAL':['MG1','MG2']
}
def methyl_filter(df):
    for res in methyl:
        for i in range(len(methyl[res])):
            df.loc[(df['Atom_chem_shift.Comp_ID'] == res) & (df['Atom_chem_shift.Atom_ID'].isin(methyl[res][i])), 'Atom_chem_shift.Atom_ID'] = methyl_pseudo[res][i]
    df.drop_duplicates()
    return df

def filter2(df,c):
    df_filtered=None
    for aa in amino_acids_list:
        df2 = df[df['Atom_chem_shift.Comp_ID']==aa]
        atom_list = list(set(df2['Atom_chem_shift.Atom_ID']))
        for atom in atom_list:
            cs = df2[df2['Atom_chem_shift.Atom_ID'] == atom]
            m = numpy.mean(cs['Atom_chem_shift.Val'])
            sd = numpy.std(cs['Atom_chem_shift.Val'])
            if df_filtered is None:
                df_filtered = cs[(m-c*sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m+c*sd)]
            else:
                df_filtered = df_filtered.append(cs[(m-c*sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m+c*sd)])
    return df_filtered

api_url='http://api.bmrb.io/v2/search/chemical_shifts?'
for amino_acid in amino_acids_list:
    api_url = f'{api_url}comp_id={amino_acid}&'
api_url_h = f'{api_url}atom_type=H'
api_url_c = f'{api_url}atom_type=C'
api_url_n = f'{api_url}atom_type=N'
print ('Downloading Proton shifts')
print (api_url_h)
data_h = pd.read_json(api_url_h,orient='split')
data_h = methyl_filter(data_h)
print ('Downloading Carbon shifts')
print (api_url_c)
data_c = pd.read_json(api_url_c,orient='split')
print ('Downloading Nitrogen shifts')
print (api_url_n)
data_n = pd.read_json(api_url_n,orient='split')
print ('Filtering Proton shifts')
data_filtered_h = filter2(data_h,8)
print ('Filtering Carbon shifts')
data_filtered_c = filter2(data_c,8)
print ('Filtering Nitrogen shifts')
data_filtered_n = filter2(data_n,8)
data = data_h
data = data.append(data_c)
data = data.append(data_n)
data_filtered = data_filtered_h
data_filtered = data_filtered.append(data_filtered_c)
data_filtered = data_filtered.append(data_filtered_n)


app.layout = html.Div(id='parent', children=[
        html.H1(id='H1', children='BMRB Chemical Shift Statistics', style={'textAlign': 'center','marginTop': 40, 'marginBottom': 40}),
        html.Div( children=[
            html.Div( children=[
                html.H3('Amino acids'),
                dcc.Checklist(id='res',options= sorted(['ALA','CYS','ASP','GLU','PHE',
                                              'GLY','HIS','ILE','LYS','LEU',
                                              'MET','ASN','PRO','GLN','ARG',
                                              'SER','THR','VAL','TRP','TYR']),value=['ALA']),
                html.H3('Protons'),
                dcc.Checklist(id='atm_h',value=['H'])],style={'padding': 10, 'flex': 1}),#,options= sorted(alist_h),value=['H']),
            html.Div( children=[
                html.H3('Carbons'),
                dcc.Checklist(id='atm_c',value=['C']),
                html.H3('Nitrogens'),
                dcc.Checklist(id='atm_n',value=['N'])],style={'padding': 10, 'flex': 1}),
            html.Div( children=[
                html.H3('Plot type'),
                dcc.RadioItems(id='plot_type',options=[{'label':'Histogram','value':'hist'},
                                     {'label':'Violin','value':'violin'},
                                     {'label':'Box plot','value':'box'}],value='hist'),
                html.H3('Histogram norm'),
                dcc.RadioItems(id='hist_norm',options=[{'label':'Count','value':''},
                                     {'label':'Percent','value':'percent'},
                                     {'label':'Density','value':'probability'}],value='')],style={'padding': 10, 'flex': 1}),
            html.Div( children=[
                html.H3('Data filter'),
                dcc.RadioItems(id='filt',options=[{'label':'Full','value':'Full'},
                                     {'label':'Filtered','value':'Filtered'}],value='Filtered'),

                html.H3('Bin size'),
                dcc.Slider(0, 0.5, 0.05,value=0.1,id='bin_size')], style={'padding': 10, 'flex': 1}),
        ],style={'display': 'flex', 'flex-direction': 'row'}),
        dcc.Graph(id='hist')
    ])


@app.callback(Output('atm_h','options'),
              Input('res','value'))
def set_atomsh_for_amio_acids(aa):
    return sorted(list(set(data_filtered_h[data_filtered_h['Atom_chem_shift.Comp_ID'].isin(aa)]['Atom_chem_shift.Atom_ID'])))

@app.callback(Output('atm_c','options'),
              Input('res','value'))
def set_atomsc_for_amio_acids(aa):
    return sorted(list(set(data_filtered_c[data_filtered_c['Atom_chem_shift.Comp_ID'].isin(aa)]['Atom_chem_shift.Atom_ID'])))

@app.callback(Output('atm_n','options'),
              Input('res','value'))
def set_atomsn_for_amio_acids(aa):
    return sorted(list(set(data_filtered_n[data_filtered_n['Atom_chem_shift.Comp_ID'].isin(aa)]['Atom_chem_shift.Atom_ID'])))


@app.callback(Output(component_id='hist', component_property='figure'),
              [Input(component_id='res', component_property='value'),
                Input(component_id='atm_h', component_property='value'),
                Input(component_id='atm_c', component_property='value'),
               Input(component_id='atm_n', component_property='value'),
               Input(component_id='bin_size', component_property='value'),

               Input(component_id='filt', component_property='value'),

               Input(component_id='hist_norm', component_property='value'),
               Input(component_id='plot_type', component_property='value'),
               ])
def graph_update(res,
                 atm_h,atm_c,atm_n,
                 bin_size,
                 filt,
                 hist_norm,
                 plot_type):
    atm = atm_h+atm_c+atm_n
    if filt == 'Full':
        df = data[data['Atom_chem_shift.Comp_ID'].isin(res)]
        df = df[df['Atom_chem_shift.Atom_ID'].isin(atm)]
    else:
        df = data_filtered[data_filtered['Atom_chem_shift.Comp_ID'].isin(res)]
        df = df[df['Atom_chem_shift.Atom_ID'].isin(atm)]
    if len(res) == 1:
        if plot_type == 'hist':
            fig = px.histogram(df, x='Atom_chem_shift.Val', histnorm=hist_norm, color='Atom_chem_shift.Atom_ID',
                           barmode='overlay',facet_row='Atom_chem_shift.Atom_type',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type == 'violin':
            fig = px.violin(df, x='Atom_chem_shift.Val', color='Atom_chem_shift.Atom_ID',facet_row='Atom_chem_shift.Atom_type',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]'},
                                 )
        elif plot_type == 'box':
            fig = px.box(df, x='Atom_chem_shift.Val', color='Atom_chem_shift.Atom_ID',facet_row='Atom_chem_shift.Atom_type',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    elif len(atm) == 1:
        if plot_type == 'hist':
            fig = px.histogram(df, x='Atom_chem_shift.Val', histnorm=hist_norm, color='Atom_chem_shift.Comp_ID',
                            barmode='overlay',facet_row='Atom_chem_shift.Atom_type',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type == 'violin':
            fig = px.violin(df, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',facet_row='Atom_chem_shift.Atom_type',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type == 'box':
            fig = px.box(df, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',facet_row='Atom_chem_shift.Atom_type',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )

    else:
        if plot_type == 'hist':
            fig = px.histogram(df, x='Atom_chem_shift.Val', histnorm=hist_norm, color='Atom_chem_shift.Comp_ID',
                           pattern_shape='Atom_chem_shift.Atom_ID', barmode='overlay',facet_row='Atom_chem_shift.Atom_type',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                       )
        elif plot_type == 'violin':
            fig = px.violin(df, x='Atom_chem_shift.Val',  color='Atom_chem_shift.Comp_ID',facet_row='Atom_chem_shift.Atom_type',
                                 y='Atom_chem_shift.Atom_ID', violinmode='overlay',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type == 'box':
            fig = px.box(df, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                              y='Atom_chem_shift.Atom_ID', boxmode='overlay',facet_row='Atom_chem_shift.Atom_type',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.update_yaxes(matches=None)
    fig.update_xaxes(matches=None, showticklabels=True, autorange="reversed")
    fig.update_layout(legend_title_text='')
    #fig.update_xaxes(autorange="reversed")
    fig.update_layout(legend_title_text='')
    if plot_type == 'hist':
        fig.update_traces(xbins=dict(  size=bin_size))
    return fig

if __name__ == '__main__':
    app.run_server(port=8051)
