import dash
from dash import html
from dash import dcc
import plotly.express as px
from dash.dependencies import Input, Output
from multiprocessing import Pool
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
data=None
data_filtered=None
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


def filter(df,c):
    atom_list = list(set(df['Atom_chem_shift.Atom_ID']))
    s={}
    df2=None
    for atom in atom_list:
        cs = df[df['Atom_chem_shift.Atom_ID']==atom]
        m=numpy.mean(cs['Atom_chem_shift.Val'])
        sd=numpy.std(cs['Atom_chem_shift.Val'])
        s[atom]={'mean':m,'std':sd}
        if df2 is None:
            df2 = cs[(m-c*sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m+c*sd)]
        else:
            df2=df2.append(cs[(m-c*sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m+c*sd)])
    return df2

def atm_filter(p):
    df2= p[0]
    atom = p[1]
    c = p[2]
    cs = df2[df2['Atom_chem_shift.Atom_ID'] == atom]
    m = numpy.mean(cs['Atom_chem_shift.Val'])
    sd = numpy.std(cs['Atom_chem_shift.Val'])
    df_filtered = cs[(m - c * sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m + c * sd)]
    return df_filtered

def res_filter(p):
    df= p[0]
    res = p[1]
    c=p[2]
    df_filtered = None
    df2 = df[df['Atom_chem_shift.Comp_ID'] == res]
    atom_list = list(set(df2['Atom_chem_shift.Atom_ID']))
    # p = [(df2,i,c) for i in atom_list ]
    # with Pool() as pool:
    #     results = pool.map(atm_filter,p)
    # filtered = results[0]
    # for i in range(1, len(results)):
    #     filtered = filtered.append(results[i])
    # return filtered

    for atom in atom_list:
        cs = df2[df2['Atom_chem_shift.Atom_ID'] == atom]
        m = numpy.mean(cs['Atom_chem_shift.Val'])
        sd = numpy.std(cs['Atom_chem_shift.Val'])
        if df_filtered is None:
            df_filtered = cs[(m - c * sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m + c * sd)]
        else:
            df_filtered = df_filtered.append(
                cs[(m - c * sd < cs['Atom_chem_shift.Val']) & (cs['Atom_chem_shift.Val'] < m + c * sd)])
    return df_filtered


def filter_multiprocessing(df,c):
    p = [(df,i,c) for i in amino_acids_list]
    with Pool() as pool:
        results = pool.map(res_filter, p)
    filtered = results[0]
    for i in range(1,len(results)):
        filtered = filtered.append(results[i])
    return filtered

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


app.layout = html.Div(id='parent', children=[
    html.Div([
        html.H1(id='H1', children='Proton Chemical Shift Statistics', style={'textAlign': 'center', \
                                                                      'marginTop': 40, 'marginBottom': 40}),
        html.Div( children=[
            html.Div( children=[
        html.H3('Amino acids'),
        dcc.Checklist(id='res_h',options= sorted(['ALA','CYS','ASP','GLU','PHE',
                                              'GLY','HIS','ILE','LYS','LEU',
                                              'MET','ASN','PRO','GLN','ARG',
                                              'SER','THR','VAL','TRP','TYR']),value=['ALA']),
        html.H3('Atoms'),
        dcc.Checklist(id='atm_h',value=['H']),#,options= sorted(alist_h),value=['H']),
        ],style={'padding': 10, 'flex': 1}),
html.Div( children=[
        html.H3('Plot type'),
        dcc.RadioItems(id='plot_type_h',options=[{'label':'Histogram','value':'hist'},
                                     {'label':'Violin','value':'violin'},
                                     {'label':'Box plot','value':'box'}],value='hist'),

        html.H3('Histogram norm'),
        dcc.RadioItems(id='hist_norm_h',options=[{'label':'Count','value':''},
                                     {'label':'Percent','value':'percent'},
                                     {'label':'Density','value':'probability'}],value=''),
        ],style={'padding': 10, 'flex': 1}),
html.Div( children=[
        html.H3('Data filter'),
        dcc.RadioItems(id='filt_h',options=[{'label':'Full','value':'Full'},
                                     {'label':'Filtered','value':'Filtered'}],value='Filtered'),
        html.H3('Bin size'),
        dcc.Slider(0, 0.5, 0.05,
               value=0.05,
               id='bin_size_h'),
        ],style={'padding': 10, 'flex': 1}),],style={'display': 'flex', 'flex-direction': 'row'}),
    dcc.Graph(id='hist_h')
    ]),

html.Div([
        html.H1(id='H2', children='Carbon Chemical Shift Statistics', style={'textAlign': 'center', \
                                                                      'marginTop': 40, 'marginBottom': 40}),
        html.Div( children=[
            html.Div( children=[
        html.H3('Amino acids'),
        dcc.Checklist(id='res_c',options= sorted(['ALA','CYS','ASP','GLU','PHE',
                                              'GLY','HIS','ILE','LYS','LEU',
                                              'MET','ASN','PRO','GLN','ARG',
                                              'SER','THR','VAL','TRP','TYR']),value=['ALA']),
        html.H3('Atoms'),
        dcc.Checklist(id='atm_c',value=['C']),#,options= sorted(alist_c),value=['C']),
        ],style={'padding': 10, 'flex': 1}),
html.Div( children=[
        html.H3('Plot type'),
        dcc.RadioItems(id='plot_type_c',options=[{'label':'Histogram','value':'hist'},
                                     {'label':'Violin','value':'violin'},
                                     {'label':'Box plot','value':'box'}],value='hist'),
        html.H3('Histogram norm'),
        dcc.RadioItems(id='hist_norm_c',options=[{'label':'Count','value':''},
                                     {'label':'Percent','value':'percent'},
                                     {'label':'Density','value':'probability'}],value=''),
        ],style={'padding': 10, 'flex': 1}),

html.Div( children=[
html.H3('Data filter'),
        dcc.RadioItems(id='filt_c',options=[{'label':'Full','value':'Full'},
                                     {'label':'Filtered','value':'Filtered'}],value='Filtered'),
        html.H3('Bin size'),
        dcc.Slider(0, 2, 0.1,
               value=0.5,
               id='bin_size_c'),
        ],style={'padding': 10, 'flex': 1}),],style={'display': 'flex', 'flex-direction': 'row'}),
    dcc.Graph(id='hist_c')
    ]),

html.Div([
        html.H1(id='H3', children='Nitrogen Chemical Shift Statistics', style={'textAlign': 'center', \
                                                                      'marginTop': 40, 'marginBottom': 40}),
        html.Div( children=[
            html.Div( children=[
        html.H3('Amino acids'),
        dcc.Checklist(id='res_n',options= sorted(['ALA','CYS','ASP','GLU','PHE',
                                              'GLY','HIS','ILE','LYS','LEU',
                                              'MET','ASN','PRO','GLN','ARG',
                                              'SER','THR','VAL','TRP','TYR']),value=['ALA']),
        html.H3('Atoms'),
        dcc.Checklist(id='atm_n',value=['N']),#,options= sorted(alist_n),value=['N']),
            ],style={'padding': 10, 'flex': 1}),
html.Div( children=[
        html.H3('Plot type'),
        dcc.RadioItems(id='plot_type_n',options=[{'label':'Histogram','value':'hist'},
                                     {'label':'Violin','value':'violin'},
                                     {'label':'Box plot','value':'box'}],value='hist'),

        html.H3('Histogram norm'),
        dcc.RadioItems(id='hist_norm_n',options=[{'label':'Count','value':''},
                                     {'label':'Percent','value':'percent'},
                                     {'label':'Density','value':'probability'}],value=''),
        ],style={'padding': 10, 'flex': 1}),
html.Div( children=[
        html.H3('Data filter'),
        dcc.RadioItems(id='filt_n',options=[{'label':'Full','value':'Full'},
                                     {'label':'Filtered','value':'Filtered'}],value='Filtered'),
        html.H3('Bin size'),
        dcc.Slider(0, 2, 0.1,
               value=0.5,
               id='bin_size_n'),
        ],style={'padding': 10, 'flex': 1}),
],style={'display': 'flex', 'flex-direction': 'row'}),
    dcc.Graph(id='hist_n')
    ]),

])

@app.callback(Output('atm_h','options'),
              Input('res_h','value'))
def set_atomsh_for_amio_acids(aa):
    return sorted(list(set(data_filtered_h[data_filtered_h['Atom_chem_shift.Comp_ID'].isin(aa)]['Atom_chem_shift.Atom_ID'])))

@app.callback(Output('atm_c','options'),
              Input('res_c','value'))
def set_atomsc_for_amio_acids(aa):
    return sorted(list(set(data_filtered_c[data_filtered_c['Atom_chem_shift.Comp_ID'].isin(aa)]['Atom_chem_shift.Atom_ID'])))

@app.callback(Output('atm_n','options'),
              Input('res_n','value'))
def set_atomsn_for_amio_acids(aa):
    return sorted(list(set(data_filtered_n[data_filtered_n['Atom_chem_shift.Comp_ID'].isin(aa)]['Atom_chem_shift.Atom_ID'])))


@app.callback([Output(component_id='hist_h', component_property='figure'),
               Output(component_id='hist_c', component_property='figure'),
               Output(component_id='hist_n', component_property='figure')],
              [Input(component_id='res_h', component_property='value'),
               Input(component_id='res_c', component_property='value'),
               Input(component_id='res_n', component_property='value'),
                Input(component_id='atm_h', component_property='value'),
                Input(component_id='atm_c', component_property='value'),
               Input(component_id='atm_n', component_property='value'),
               Input(component_id='bin_size_h', component_property='value'),
               Input(component_id='bin_size_c', component_property='value'),
               Input(component_id='bin_size_n', component_property='value'),
               Input(component_id='filt_h', component_property='value'),
               Input(component_id='filt_c', component_property='value'),
               Input(component_id='filt_n', component_property='value'),
               Input(component_id='hist_norm_h', component_property='value'),
               Input(component_id='hist_norm_c', component_property='value'),
               Input(component_id='hist_norm_n', component_property='value'),
               Input(component_id='plot_type_h', component_property='value'),
               Input(component_id='plot_type_c', component_property='value'),
               Input(component_id='plot_type_n', component_property='value')])
def graph_update(res_h,res_c,res_n,
                 atm_h,atm_c,atm_n,
                 bin_size_h,bin_size_c,bin_size_n,
                 filt_h,filt_c,filt_n,
                 hist_norm_h,hist_norm_c,hist_norm_n,
                 plot_type_h,plot_type_c,plot_type_n):
    if filt_h == 'Full':
        df_h = data_h[data_h['Atom_chem_shift.Comp_ID'].isin(res_h)]
        df_h = df_h[df_h['Atom_chem_shift.Atom_ID'].isin(atm_h)]
    else:
        df_h = data_filtered_h[data_filtered_h['Atom_chem_shift.Comp_ID'].isin(res_h)]
        df_h = df_h[df_h['Atom_chem_shift.Atom_ID'].isin(atm_h)]
    if filt_c == 'Full':
        df_c = data_c[data_c['Atom_chem_shift.Comp_ID'].isin(res_c)]
        df_c = df_c[df_c['Atom_chem_shift.Atom_ID'].isin(atm_c)]
    else:
        df_c = data_filtered_c[data_filtered_c['Atom_chem_shift.Comp_ID'].isin(res_c)]
        df_c = df_c[df_c['Atom_chem_shift.Atom_ID'].isin(atm_c)]
    if filt_n == 'Full':
        df_n = data_n[data_n['Atom_chem_shift.Comp_ID'].isin(res_n)]
        df_n = df_n[df_n['Atom_chem_shift.Atom_ID'].isin(atm_n)]
    else:
        df_n = data_filtered_n[data_filtered_n['Atom_chem_shift.Comp_ID'].isin(res_n)]
        df_n = df_n[df_n['Atom_chem_shift.Atom_ID'].isin(atm_n)]


    if len(res_h) == 1:
        if plot_type_h == 'hist':
            fig_h = px.histogram(df_h, x='Atom_chem_shift.Val', histnorm=hist_norm_h, color='Atom_chem_shift.Atom_ID',
                           barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type_h == 'violin':
            fig_h = px.violin(df_h, x='Atom_chem_shift.Val', color='Atom_chem_shift.Atom_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_h == 'box':
            fig_h = px.box(df_h, x='Atom_chem_shift.Val', color='Atom_chem_shift.Atom_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    elif len(atm_h) == 1:
        if plot_type_h == 'hist':
            fig_h = px.histogram(df_h, x='Atom_chem_shift.Val', histnorm=hist_norm_h, color='Atom_chem_shift.Comp_ID',
                            barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type_h == 'violin':
            fig_h = px.violin(df_h, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_h == 'box':
            fig_h = px.box(df_h, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )

    else:
        if plot_type_h == 'hist':
            fig_h = px.histogram(df_h, x='Atom_chem_shift.Val', histnorm=hist_norm_h, color='Atom_chem_shift.Comp_ID',
                           pattern_shape='Atom_chem_shift.Atom_ID', barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                       )
        elif plot_type_h == 'violin':
            fig_h = px.violin(df_h, x='Atom_chem_shift.Val',  color='Atom_chem_shift.Comp_ID',
                                 y='Atom_chem_shift.Atom_ID', violinmode='overlay',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_h == 'box':
            fig_h = px.box(df_h, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                              y='Atom_chem_shift.Atom_ID', boxmode='overlay',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    fig_h.update_xaxes(autorange="reversed")
    fig_h.update_layout(legend_title_text='')
    if plot_type_h == 'hist':
        fig_h.update_traces(xbins=dict(  size=bin_size_h))

    if len(res_c) == 1:
        if plot_type_c == 'hist':
            fig_c = px.histogram(df_c, x='Atom_chem_shift.Val', histnorm=hist_norm_c, color='Atom_chem_shift.Atom_ID',
                           barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type_c == 'violin':
            fig_c = px.violin(df_c, x='Atom_chem_shift.Val',  color='Atom_chem_shift.Atom_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_c == 'box':
            fig_c = px.box(df_c, x='Atom_chem_shift.Val', color='Atom_chem_shift.Atom_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    elif len(atm_c) == 1:
        if plot_type_c == 'hist':
            fig_c = px.histogram(df_c, x='Atom_chem_shift.Val', histnorm=hist_norm_c, color='Atom_chem_shift.Comp_ID',
                            barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type_c == 'violin':
            fig_c = px.violin(df_c, x='Atom_chem_shift.Val',  color='Atom_chem_shift.Comp_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_c == 'box':
            fig_c = px.box(df_c, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    else:
        if plot_type_c == 'hist':
            fig_c = px.histogram(df_c, x='Atom_chem_shift.Val', histnorm=hist_norm_c, color='Atom_chem_shift.Comp_ID',
                           pattern_shape='Atom_chem_shift.Atom_ID', barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                       )
        elif plot_type_c == 'violin':
            fig_c = px.violin(df_c, x='Atom_chem_shift.Val',  color='Atom_chem_shift.Comp_ID',
                                 y='Atom_chem_shift.Atom_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_c == 'box':
            fig_c = px.box(df_c, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                              y='Atom_chem_shift.Atom_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
    fig_c.update_xaxes(autorange="reversed")
    fig_c.update_layout(legend_title_text='')
    if plot_type_c == 'hist':
        fig_c.update_traces(xbins=dict(size=bin_size_c))

    if len(res_n) == 1:
        if plot_type_n == 'hist':
            fig_n = px.histogram(df_n, x='Atom_chem_shift.Val', histnorm=hist_norm_n, color='Atom_chem_shift.Atom_ID',
                               barmode='overlay',
                               labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                       'Atom_chem_shift.Atom_ID': 'Atom'},
                               )
        elif plot_type_n == 'violin':
            fig_n = px.violin(df_n, x='Atom_chem_shift.Val' , color='Atom_chem_shift.Atom_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'})
        elif plot_type_n == 'box':
            fig_n = px.box(df_n, x='Atom_chem_shift.Val' , color='Atom_chem_shift.Atom_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'})
    elif len(atm_n) == 1:
        if plot_type_n == 'hist':
            fig_n = px.histogram(df_n, x='Atom_chem_shift.Val', histnorm=hist_norm_n, color='Atom_chem_shift.Comp_ID',
                            barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                           )
        elif plot_type_n == 'violin':
            fig_n = px.violin(df_n,x='Atom_chem_shift.Val',color='Atom_chem_shift.Comp_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )
        elif plot_type_n == 'box':
            fig_n = px.box(df_n, x='Atom_chem_shift.Val', color='Atom_chem_shift.Comp_ID',
                              labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                      'Atom_chem_shift.Atom_ID': 'Atom'},
                              )

    else:
        if plot_type_n == 'hist':
            fig_n = px.histogram(df_n, x='Atom_chem_shift.Val', histnorm=hist_norm_n, color='Atom_chem_shift.Comp_ID',
                           pattern_shape='Atom_chem_shift.Atom_ID', barmode='overlay',
                           labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                   'Atom_chem_shift.Atom_ID': 'Atom'},
                       )
        elif plot_type_n == 'violin':
            fig_n = px.violin(df_n, x='Atom_chem_shift.Val', y='Atom_chem_shift.Atom_ID', color='Atom_chem_shift.Comp_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )
        elif plot_type_n == 'box':
            fig_n = px.box(df_n, x='Atom_chem_shift.Val', y='Atom_chem_shift.Atom_ID', color = 'Atom_chem_shift.Comp_ID',
                                 labels={'Atom_chem_shift.Val': 'Chemical shift [ppm]',
                                         'Atom_chem_shift.Atom_ID': 'Atom'},
                                 )

    fig_n.update_xaxes(autorange="reversed")
    fig_n.update_layout(legend_title_text='')
    if plot_type_n == 'hist':
        fig_n.update_traces(xbins=dict(  size=bin_size_n))
    return [fig_h,fig_c,fig_n]


if __name__ == '__main__':
    app.run_server()