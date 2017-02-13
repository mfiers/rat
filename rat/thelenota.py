
import copy
from functools import lru_cache, partial
import io
import json
import logging
import pickle
import random

from IPython.display import display
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from path import Path
import requests
from scipy.stats import pearsonr

from bokeh.colors import RGB
import bokeh.charts as bokeh_chart
import bokeh.io  as bokeh_io
from bokeh.models import ColumnDataSource, CategoricalColorMapper, \
                         LinearColorMapper, ColorBar, LogTicker, \
                         BasicTicker, BasicTickFormatter, CustomJS
import bokeh.models as bmodels
import bokeh.palettes
from bokeh.plotting import figure as bokeh_figure

lg = logging.getLogger(__name__)

# CONSTANTS
CONFIG_FILE_NAME = 'config.yaml'
COUNT_TABLE_DIR = '55.normalized_counts'
COUNT_TABLE_EXTENSION = '.tsv'
CACHE_DIR = '99.thelenota_cache'
DIMRED_DIR = '60.dimred'
FIGSIZE=(600,500)
FIGDPI=100
FIGSIZE2 = FIGSIZE[0]/FIGDPI, FIGSIZE[1]/FIGDPI

#
# Helper functions
#


# ipywidget layout
def ilabel(lab, wid):
    if not isinstance(wid, (list, tuple)):
        wid = [wid]

    return widgets.HBox(
            [widgets.Label(lab, layout=widgets.Layout(width='120px',
                                                  max_width='120px'))] +
            wid )

# ipywidget layout
def ilabel2(lab, wid1, wid2):
    return widgets.HBox(
        [widgets.Label(lab, layout=widgets.Layout(width='120px',
                                                  max_width='120px')),
         wid1, wid2])

# save/load functions
def io_load_pkl(filename):
    with open(filename, 'rb') as F:
        return pickle.load(F)


def io_save_pkl(filename, blob):
    with open(filename, 'wb') as F:
        pickle.dump(blob, F)


def io_load_str(filename):
    with open(filename, 'r') as F:
        return pickle.load(F)


def io_save_str(filename, blob):
    with open(filename, 'w') as F:
        pickle.dump(blob, F)


def io_load_bin(filename):
    with open(filename, 'rb') as F:
        return F.read()


def io_save_bin(filename, blob):
    with open(filename, 'wb') as F:
        F.write(blob)


def io_load_tsv(filename):
    return pd.read_csv(filename, compression='gzip', sep="\t", index_col=0)


def io_save_tsv(filename, blob):
    blob.to_csv(filename, compression='gzip', sep="\t")


def io_load_meta_tsv(filename):
    metaname = filename.replace('.tsv.gz', '.meta.tsv')
    a = pd.read_csv(metaname, sep="\t", index_col=0)
    b = pd.read_csv(filename, compression='gzip', sep="\t", index_col=0)
    return a, b


def io_save_meta_tsv(filename, blob):
    a, b = blob
    metaname = filename.replace('.tsv.gz', '.meta.tsv')
    a.to_csv(metaname, sep="\t")
    b.to_csv(filename, compression='gzip', sep="\t")


IOFUNC = dict(
    tsv = (io_load_tsv, io_save_tsv),
    png = (io_load_bin, io_save_bin),
    str = (io_load_str, io_save_str),
    mtsv=(io_load_meta_tsv, io_save_meta_tsv))

#
# Dimred methods
#

def run_pca(mat):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=10)
    fit = pca.fit(mat.T)
    tra = fit.transform(mat.T)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    meta['variance_explained'] = fit.explained_variance_ratio_
    meta['method'] = 'pca'
    return meta, tra

def run_nmf(mat):
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=10)
    tmat = mat.T + abs(mat.min().min())
    fit = nmf.fit(tmat)
    tra = fit.transform(tmat)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    #meta['variance_explained'] = fit.explained_variance_ratio_
    meta['method'] = 'pca'
    return meta, tra

def run_tsne(mat, get_param):
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    param = get_param()
    tsne_param_names = '''early_exaggeration
                          perplexity
                          learning_rate'''.split()
    tsne_param = {a:b for a,b in param.items() if a in tsne_param_names}
    pca = PCA(n_components=param.get('pca_components', 20))
    ptra = pca.fit_transform(mat.T)
    tsne = TSNE(random_state=42, **tsne_param)
    ttra = pd.DataFrame(tsne.fit_transform(ptra), index=mat.columns)

    meta = pd.DataFrame(index=ttra.columns)
    meta['method'] = 'tsne'
    return meta, ttra


def setcache(thelenota, category, fmt, namer, value):
    cache_dir = thelenota.cachedir / category
    cache_dir.makedirs_p()
    objname = namer()
    extension = dict(mtsv='tsv.gz',
                     tsv='tsv.gz').get(fmt, fmt)
    cache_file = cache_dir / '{}.{}'.format(objname, extension)
    _, io_save = IOFUNC.get(fmt, (io_load_pkl, io_save_pkl))
    io_save(cache_file, value)


def dcache(thelenota, category, fmt, namer, force=lambda: False):
    """cache decorator"""
    def cfunc(func):
        def cfunc2(*args, **kwargs):
            cache_dir = thelenota.cachedir / category
            cache_dir.makedirs_p()

            objname = namer()
            extension = dict(mtsv='tsv.gz',
                             tsv='tsv.gz').get(fmt, fmt)

            cache_file = cache_dir / '{}.{}'.format(objname, extension)
            io_load, io_save = IOFUNC.get(fmt, (io_load_pkl, io_save_pkl))

            if cache_file.exists() and not force():
                # print('load from cache ({})'.format(force()))
                return io_load(cache_file)
            else:
                # print('redo')
                rv = func(*args, **kwargs)
                io_save(cache_file, rv)
                return rv
        return cfunc2
    return cfunc

#@lru_cache(32)
def get_project_metadata_info_2(metadata_dir):
    all_m = []
    for meta in metadata_dir.glob('*.meta.tsv'):
        #print(meta)
        m = pd.read_csv(meta, sep="\t", header=None, names=['column', 'datatype'],
                        index_col=0)
        m['filename'] = meta.basename()
        all_m.append(m)

    rv = pd.concat(all_m)
    rv.loc['thelenota'] = dict(datatype='categorical', filename='thelenota.tsv')
    return rv

@lru_cache(32)
def _tsv_loader(filename):
    rv = pd.read_csv(filename, index_col=0, sep="\t")
    return rv

def create_distribution_plot(data, title, ylog=False):
    """ Create a simple plot
    """
    fig = plt.figure(figsize=FIGSIZE2)
    plt.plot(range(len(data)), data)
    plt.title(title)
    if ylog:
        plt.yscale('log')
    img = io.BytesIO()
    plt.savefig(img, format='png', dpi=FIGDPI)
    plt.close()
    img.seek(0)
    return img.read()

@lru_cache(32)
def get_dr_meta(drdir, drset):
    drmeta = drdir / '{}.meta.tsv'.format(drset)
    return pd.read_csv(drmeta, sep="\t", index_col=0)

@lru_cache(32)
def get_dr_data(drdir, drset):
    drdata = drdir / '{}.tsv'.format(drset)
    rv = pd.read_csv(drdata, sep="\t", index_col=0)
    rv.columns = map(str, rv.columns)
    return rv


def create_scatter_plot(x, y, title):
    """ Create a simple plot
    """
    fig = plt.figure(figsize=FIGSIZE2)
    plt.scatter(x, y, s=5)
    plt.title(title)
    img = io.BytesIO()
    plt.savefig(img, format='png', dpi=FIGDPI)
    plt.close()
    img.seek(0)
    return img.read()


def create_count_table_stats(c):
    rv = {}
    rv['no_genes'] = c.shape[0]
    rv['no_cells'] = c.shape[1]
    rv['min_value'] = c.min().min()
    rv['max_value'] = c.max().max()
    return pd.Series(rv)

class Thelenota:

    def __init__(self,
                 basedir,
                 experiment_name=None,
                 count_table_dir=COUNT_TABLE_DIR,
                 ):

        #this will contain the indici of selected spots
        self._DR_INDICI_SEL_ = []
        self._DF_INDICI_SEL_ = []

        # WIDGET INSTANCES
        self.basedir = Path(basedir).expanduser().abspath()
        self.count_table_dir = count_table_dir
        if experiment_name is None:
            experiment_name = self.experiments[0]
        self.experiment_w = widgets.Dropdown(
            options=self.experiments, value=experiment_name)
        self.counttable_w = widgets.Dropdown(
            options=self.counttables, value=self.counttables[0])
        self.geneset_dir = None

    def __str__(self):
        return "<thelenota {}>".format(self.basedir)
    __repr__ = __str__


    #
    # Experiment
    #
    @property
    @lru_cache(12)
    def experiments(self):
        "Return a list of Single Cell experiments in the basedir"
        rv = []
        for d in self.basedir.dirs():
            if not (d / CONFIG_FILE_NAME).exists():
                continue
            if not (d / self.count_table_dir).exists():
                continue
            rv.append(str(d.basename()))
        return list(sorted(rv))

    @property
    def cachedir(self):
        """ Return a epxeriment specific cache directory
        """
        cachedir = self.basedir / self.experiment_name / CACHE_DIR
        cachedir.makedirs_p()
        return cachedir

    def get_experiment_name(self):
        return self.experiment_w.value
    def set_experiment_name(self, value):
        lg.warning("using set_experiment_name ({})".format(value))
        self.experiment_w.value = value
    experiment_name = property(get_experiment_name, set_experiment_name)


    def select_experiment(self):
        display(widgets.HBox([widgets.Label('Set Experiment: '),
                         self.experiment_w]))
    def show_experiment(self):
        expshow = widgets.Text(value=self.experiment_name, disabled=True)
        def _update_expshow(_):
            expshow.value = self.experiment_name
        self.experiment_w.observe(_update_expshow, 'value')
        display(widgets.HBox([widgets.Label('Experiment: '),
                              expshow]))

    @property
    def metadata_dir(self):
        return self.basedir / self.experiment_name \
                / '90.metadata'

    @property
    def metadata_info(self):
        return self.get_metadata_info()

    def get_metadata_info(self, mtype=None):
        rv = get_project_metadata_info_2(self.metadata_dir)

        if mtype is None:
            return rv
        else:
             return rv[rv['datatype'] == mtype]

    def get_metadata(self, name):
        minfo = self.metadata_info
        meta_filename = self.metadata_dir / minfo.loc[name]['filename']
        data_filename = Path(meta_filename.replace('.meta.tsv', '.tsv'))
        if not data_filename.exists() and name == 'thelenota':
            return pd.Series('na', index=self.counttable.columns)

        m = pd.read_csv(data_filename, sep="\t", index_col=0)
        m = m.loc[sorted(m.index)]
        return m[name]


    #
    # Count table
    #
    @property
    def countdir(self):
        return self.basedir / self.experiment_name \
            / self.count_table_dir

    @property
    def countfile(self):
        rv = self.basedir / self.experiment_name \
            / self.count_table_dir
        rv /= self.counttable_name + COUNT_TABLE_EXTENSION
        return rv

    @property
    def counttables(self):
        rv = []
        for f in self.countdir.glob('*' + COUNT_TABLE_EXTENSION):
            name = f.basename().replace(COUNT_TABLE_EXTENSION, '')
            rv.append(name)
        return list(sorted(rv))

    def select_counttable(self):
        def _change_counttable_w(_):
            self.counttable_w.options = self.counttables
            self.counttable_w.value = self.counttables[0]

        self.experiment_w.observe(_change_counttable_w, 'value')
        display(widgets.HBox([widgets.Label('Set Counttable: '),
                         self.counttable_w]))

    def show_counttable(self):
        expshow = widgets.Text(value=self.counttable_name, disabled=True)
        def _update_expshow(_):
            expshow.value = self.counttable_name
        self.counttable_w.observe(_update_expshow, 'value')
        display(widgets.HBox([widgets.Label('Counttable: '),
                              expshow]))

    def get_counttable_name(self):
        return self.counttable_w.value

    def set_counttable_name(self, value):
        self.counttable_w.value = value

    counttable_name = property(get_counttable_name, set_counttable_name)

    @property
    def counttable(self):
        rv = _tsv_loader(self.countfile)
        rv = rv[sorted(rv.columns)]
        rv = rv.loc[sorted(rv.index)]
        return rv

    #
    # basic count table stats
    #
    @property
    def counttable_stats(self):
        return tcache(
            self, 'counttable_stats', self.counttable_name, 'tsv',
            partial(create_count_table_stats, self.counttable))


    def show_counttable_stats(self):
        show_w = widgets.HTML()
        def _fill_stats(*args):
            stab = pd.DataFrame(self.counttable_stats).to_html()
            show_w.value = stab
        self.counttable_w.observe(_fill_stats, 'value')
        display(show_w)
        _fill_stats()

    #
    # mean gene expression plot
    #
    def mean_gene_expression(self):
        img_w = widgets.Image(format='png',
                width=FIGSIZE[0], height=FIGSIZE[1])
        log_w = widgets.ToggleButton(value=False, description='log y')

        def image_name():
            return '{}_{}'.format(self.counttable_name, log_w.value)

        @dcache(self, 'mean_gene_expression', 'png', image_name)
        def get_image():
            return create_distribution_plot(
                data = self.counttable.mean(1).sort_values(),
                title = 'mean gene expression {}'.format(self.counttable_name),
                ylog = log_w.value)

        def on_change(*_):
            img_w.value = get_image()

        log_w.observe(on_change, 'value')
        self.counttable_w.observe(on_change, 'value')
        display(log_w)
        display(img_w)
        on_change()

    #
    # mean cell expression plot
    #
    def mean_cell_expression(self):
        img_w = widgets.Image(format='png',
                width=FIGSIZE[0], height=FIGSIZE[1])
        log_w = widgets.ToggleButton(value=False, description='log y')

        def _makeplot(*args):
            _makeplot = partial(create_distribution_plot,
                data = self.counttable.mean().sort_values(),
                title = 'mean cell expression {}'.format(self.counttable_name),
                ylog=log_w.value)
            cat = 'mean_cell_expression' if not log_w.value else "mean_cell_expression_log"
            _makeplot_c = partial(tcache, self, cat, self.counttable_name,
                'png', _makeplot)
            img_w.value = _makeplot_c()

        log_w.observe(_makeplot, 'value')
        self.counttable_w.observe(_makeplot, 'value')
        display(log_w)
        display(img_w)
        _makeplot()

    def get_geneset_sets(self):
        if self.geneset_dir is None:
            return []

    @property
    def dimred_dir(self):
        drdir = self.basedir  / self.experiment_name / DIMRED_DIR
        drdir /= self.counttable_name
        return drdir

    @property
    def dimred_sets(self):
        rv = []
        drdir = self.dimred_dir

        for f in drdir.glob('*.tsv'):
            if f.endswith('.meta.tsv'):
                continue
            name = f.basename().replace('.tsv', '')
            rv.append(name)
        return list(sorted(rv))

    def DIFX(self):

        # define widgets
        cat_meta_grp = list(
            self.get_metadata_info('categorical').index)
        wstyle =  widgets.Layout(width='120', max_width='120px')
        sl_group_a = widgets.Dropdown(
            options=cat_meta_grp, layout=wstyle, value='thelenota')
        sl_set_a = widgets.Dropdown(options=[], layout=wstyle)
        sl_group_b = widgets.Dropdown(
            options=cat_meta_grp, layout=wstyle, value='thelenota')
        sl_set_b = widgets.Dropdown(options=[], layout=wstyle)
        sl_go = widgets.Button(description='Go')
        sl_enrichr_link = widgets.Button(description='Enrichr')
        sl_norm = widgets.Checkbox(value=False)
        html_w = widgets.HTML()
        html_link_w = widgets.HTML()

        nogenes = self.counttable.shape[0]
        colv = [1] * nogenes
        color_mapper = LinearColorMapper(palette="Viridis256", low=-0.3, high=2.5)
        pdata = ColumnDataSource(dict(
            x=[random.uniform(-10, 10) for x in range(nogenes)],
            y=[random.uniform(-10, 10) for x in range(nogenes)],
            size=[0.1] * nogenes,
            desc=list(self.counttable.index),
            score=colv))
        self._fdata = pdata

        select_callback = CustomJS(
            args = {'dsource': pdata},
            code = """
               var indici = dsource.selected['1d'].indices;
               console.log(indici);
               IPython.notebook.kernel.execute(
                    'T._DF_INDICI_SEL_ = ' + indici);
                """ )

        bokeh_tools = [
            bmodels.HoverTool(tooltips=[
                                    ("(mean,lfc)", "($x, $y)"),
                                    ("desc", "@desc"),
                                ]),
            bmodels.BoxSelectTool(callback = select_callback),
            bmodels.PanTool(),
            bmodels.WheelZoomTool(),
            bmodels.BoxZoomTool(),
            bmodels.LassoSelectTool(callback = select_callback),
            bmodels.SaveTool(),
            bmodels.ResetTool(),
            bmodels.HelpTool(),
        ]

        bfigure = bokeh_figure(plot_width=FIGSIZE[0],
                              plot_height=FIGSIZE[1],
                              tools = bokeh_tools,
                              toolbar_sticky = False,
                              toolbar_location='left',
                              title='diffexplot')
        bfigure.title.text_color = 'darkgrey'
        bfigure.title.text_font_style = 'normal'
        bfigure.title.text_font_size= "12px"
        self._ffigure = bfigure

        bplot = bfigure.circle(
                x='x', y='y', radius='size', source=pdata,
                legend='score',
                color=dict(field='score', transform=color_mapper))

        self._fplot = bplot
        bhandle= bokeh_io.show(bfigure, notebook_handle=True)

        def go_enrichr(*args):
            idx = self._DF_INDICI_SEL_
            if len(idx) == 0:
                return
            genes = list(self.counttable.index.to_series()\
                        .iloc[list(self._DR_INDICI_SEL_)])
            genes_str = "\n".join(genes)
            description = 'subset 1'
            payload = {
                'list': (None, genes_str),
                'description': (None, description)
            }

            ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
            response = requests.post(ENRICHR_URL, files=payload)

            if not response.ok:
                raise Exception('Error analyzing gene list')

            data = json.loads(response.text)
            shortid = data['shortId']
            newurl = 'http://amp.pharm.mssm.edu/Enrichr/enrich?dataset='
            newurl += shortid
            js = '<script>window.open("{}", "_blank")</script>'.format(newurl)
            html_link_w.value=js

        sl_enrichr_link.on_click(go_enrichr)
        # group defintion logic
        def update_groups(*args):
            #fill the menus with he correct values
            group_a = sl_group_a.value
            group_b = sl_group_b.value
            meta_a = self.get_metadata(group_a)
            meta_b = self.get_metadata(group_b) \
                if group_b != group_a else meta_a
            valc_a = meta_a.value_counts().sort_values(ascending=False)
            valc_b = meta_b.value_counts().sort_values(ascending=False)
            sl_set_a.options = \
                ['<all> -- {}'.format(len(meta_a))] + \
                ['{} -- {}'.format(a, b)
                 for (a,b) in valc_a.items()]
            sl_set_b.options = \
                ['<all> -- {}'.format(len(meta_a))] + \
                ['{} -- {}'.format(a, b)
                 for (a,b) in valc_b.items()]
            sl_set_a.value = sl_set_a.options[1]

        update_groups()
        sl_group_a.observe(update_groups, 'value')
        sl_group_b.observe(update_groups, 'value')

        def run(*args):
            all_samples_set = set(self.counttable.columns)
            logdata = pd.Series()
            normalize = sl_norm.value
            logdata['total cells'] = len(all_samples_set)
            group_a = sl_group_a.value
            group_b = sl_group_b.value

            set_a = sl_set_a.value.split('--')[0].strip()
            set_b = sl_set_b.value.split('--')[0].strip()

            meta_a = self.get_metadata(group_a)
            meta_b = self.get_metadata(group_b) \
                if group_b != group_a else meta_a
            logdata['group a'] = set_a
            logdata['group b'] = set_b

            sample_a = copy.copy(all_samples_set) \
                if set_a == '<all>' else set (meta_a[meta_a == set_a].index)
            sample_b = copy.copy(all_samples_set) \
                if set_b == '<all>' else set (meta_b[meta_b == set_b].index)

            sample_a &= all_samples_set #ensure overlap with this count table
            sample_b &= all_samples_set
            sample_b -= sample_a # so we don't have any duplicates

            logdata['cells in a'] = len(sample_a)
            logdata['cells in b'] = len(sample_b)

            cnts = self.counttable

            if normalize:
                cnts = 1e6 * cnts / cnts.sum()

            if cnts.min().min() < 0:
                #assume this is in log_space
                logdata['assuming log space'] = True
                cnts = 10 ** cnts

            cnts_a = cnts.loc[:,sample_a]
            cnts_b = cnts.loc[:,sample_b]
            if cnts_a.shape[1] < 1: return
            if cnts_b.shape[1] < 1: return

            html_w.value = pd.DataFrame(logdata).to_html(header=False)

            stats = pd.DataFrame(dict(
                mean_a = cnts_a.mean(1),
                mean_b = cnts_b.mean(1),
                mean_all = np.log10(cnts.mean(1)),
            ))
            stats['a/b'] = stats['mean_a'] / stats['mean_b']
            stats['lfc'] = np.log2(stats['a/b'])
            stats = stats.sort_values(by='lfc', ascending=False)
            bplot.data_source.data['x'] = stats['mean_all']
            bplot.data_source.data['y'] = stats['lfc']
            m = stats['mean_all'].max()
            bplot.data_source.data['size'] = [0.012 * m] * nogenes

            bokeh_io.push_notebook(handle=bhandle)

        sl_go.on_click(run)
        run()
        # display interface
        display(ilabel('Group A', [sl_group_a, sl_set_a]))
        display(ilabel('Group B (-A)', [sl_group_b, sl_set_b]))
        display(ilabel('TPM normalize', sl_norm))
        # tabs = widgets.Tab(children=tab_children)
        # tabs.set_title(0, 'Define sets')
        # display(tabs) sl_enrichr_link html_link_w

        display(widgets.HBox([sl_go, sl_enrichr_link]))
        display(html_w)
        display(html_link_w)




    def DRED(self):

        # DIMRED widgets

        def get_default_drmethod():
            return 'tsne'

        wstyle =  widgets.Layout(width='120', max_width='120px')
        drmethod_w = widgets.Dropdown(options='pca tsne nmf'.split(),
                                      value='tsne', layout=wstyle)
        drperplex_w = widgets.IntSlider(value=67, min=2, max=99, step=5)
        drlrate_w = widgets.IntSlider(value=920, min=20, max=10000, step=50)
        drearly_w = widgets.FloatSlider(value=3.5, min=1, max=20, step=0.5)
        drrun_w = widgets.Button(description='GO!')
        drforce_w = widgets.Checkbox(description='force',
                        layout = widgets.Layout(width='400px'))

        clmethod_w = widgets.Dropdown(layout=wstyle, options=[
            'gene', 'genes of interest', 'metadata', 'gene sigatures'])
        clgeneset_w = widgets.Dropdown(layout=wstyle, options=[
            'Top SD', 'Top expressing', 'Top DR Correlating'], disabled=True)
        clgenesetchoice_w = widgets.Dropdown(layout=wstyle, options=[], disabled=True)
        clmetadata_w = widgets.Dropdown(layout=wstyle, options=[], disabled=True)
        clgene_w = widgets.Text(layout=wstyle)

        sl_group_name_w = widgets.Text(layout=wstyle, value='thelenota')
        sl_group_set_w = widgets.Text(layout=wstyle)
        sl_group_set_go_w = widgets.Button(description='set', layout=wstyle)

        sl_groupextractname_w = widgets.Text(layout=wstyle)
        sl_group_extract_go_w = widgets.Button(description='extract')

        plotinfo = {}

        cl_w = {
            'gene': [clgene_w],
            'genes of interest': [clgeneset_w, clgenesetchoice_w],
            'metadata': [clmetadata_w],
            'gene sigatures': [],
            }

        html_w = widgets.HTML()
        dr_w = dict(method=drmethod_w,
                    perplexity=drperplex_w,
                    learning_rate=drlrate_w,
                    early_exaggeration=drearly_w)

        # data!

        samples = self.counttable.columns
        nosamples = len(samples)
        color_mapper = LinearColorMapper(palette="Viridis256", low=-0.3, high=2.5)
        topgene = self.counttable.std(1).sort_values().tail(1).index[0]
        colv = self.counttable.loc[topgene]
        clgene_w.value = topgene
        pdata = ColumnDataSource(dict(
            x=[random.uniform(-10, 10) for x in range(nosamples)],
            y=[random.uniform(-10, 10) for x in range(nosamples)],
            size=[0.3] * nosamples,
            score=colv))
        self._pdata = pdata

        # set up bokeh

        select_callback = CustomJS(
            args = {'dsource': pdata},
            code = """
               var indici = dsource.selected['1d'].indices;
               console.log(indici);
               IPython.notebook.kernel.execute(
                    'T._DR_INDICI_SEL_ = ' + indici);
                """ )

        bokeh_tools = [
            bmodels.HoverTool(), # bst,
            bmodels.BoxSelectTool(callback = select_callback),
            bmodels.PanTool(),
            bmodels.WheelZoomTool(),
            bmodels.BoxZoomTool(),
            bmodels.LassoSelectTool(callback = select_callback),
            bmodels.SaveTool(),
            bmodels.ResetTool(),
            bmodels.HelpTool(),
        ]

        bfigure = bokeh_figure(plot_width=FIGSIZE[0],
                              plot_height=FIGSIZE[1],
                              tools = bokeh_tools,
                              toolbar_sticky = False,
                              toolbar_location='left',
                              title='dimredplot')
        bfigure.title.text_color = 'darkgrey'
        bfigure.title.text_font_style = 'normal'
        bfigure.title.text_font_size= "12px"
        self._bfigure = bfigure

        bplot = bfigure.circle(
                x='x', y='y', radius='size', source=pdata,
                legend='score',
                color=dict(field='score', transform=color_mapper))

        self._bplot = bplot

        bcolorbar = ColorBar(
            color_mapper=color_mapper, ticker=BasicTicker(),
            formatter=BasicTickFormatter(precision=1), label_standoff=10,
            border_line_color=None, location=(0,0))

        bfigure.add_layout(bcolorbar, 'right')
        blegend = bfigure.legend[0].items[0]
        bhandle= bokeh_io.show(bfigure, notebook_handle=True)


        # GROUP SET
        def define_group(*args):
            groupset = sl_group_name_w.value
            groupname = sl_group_set_w.value
            self.selected = list(
                self.counttable.columns.to_series()\
                        .iloc[list(self._DR_INDICI_SEL_)])

            if groupset in self.metadata_info.index:
                tmeta = self.get_metadata(groupset)
            else:
                tmeta = pd.Series('na', index=self.counttable.columns)

            tmeta.loc[self.selected] = groupname

            outfile = self.metadata_dir / '{}.tsv'.format(groupset)
            moutfile = self.metadata_dir / '{}.meta.tsv'.format(groupset)

            tmeta = pd.DataFrame({groupset: tmeta})
            tmeta.to_csv(outfile, sep="\t")
            with open(moutfile, 'w') as F:
                F.write("{}\tcategorical\n".format(groupset))
            run_dr_color()

        sl_group_set_go_w.on_click(define_group)

        # GROUP EXTRACT
        #sl_groupextractname_w
        def extract_group(*args):
            groupname = sl_groupextractname_w.value
            self.selected = list(
                self.counttable.columns.to_series()\
                        .iloc[list(self._DR_INDICI_SEL_)])


            outfile =  Path(self.countdir) / '{}.tsv'.format(groupname)
            newcounts = self.counttable[self.selected]
            newcounts.to_csv(outfile, sep="\t")
            html_w.value = "save new count table to {}".format(outfile)

        sl_group_extract_go_w.on_click(extract_group)

        def get_dr_param():
            """ Get a dictionary with the relevant parameters
                given the method """
            method = dr_w['method'].value
            d = dict(method=method)

            if method == 'pca':
                return d

            for k, v in dr_w.items():
                d[k] = v.value
            return d

        def on_colormethod_change(*args):
            clmethod = clmethod_w.value
            for cl_w_name, cl_widgets in cl_w.items():
                for cl_widget in cl_widgets:
                    cl_widget.disabled = (cl_w_name != clmethod)

            if clmethod == 'metadata':
                clmetadata_w.options = list(self.metadata_info.index)
            elif clmethod == 'genes of interest':
                gset = clgeneset_w.value
                if gset == 'Top SD':
                    clgenesetchoice_w.options = list(
                        self.counttable.std(1).sort_values(ascending=False)\
                            .head(50).index)
                elif gset == 'Top expressing':
                    clgenesetchoice_w.options = list(
                        self.counttable.mean(1).sort_values(ascending=False)\
                            .head(50).index)
                elif gset == 'Top DR Correlating':
                    corr = get_top_correlating().max(1)
                    corr = corr.sort_values(ascending=False)[:50]
                    corr = ['{} -- {:.2g}'.format(a, b)
                            for (a, b) in zip(corr.index, corr)]
                    clgenesetchoice_w.options = corr

        clmethod_w.observe(on_colormethod_change, 'value')
        clgeneset_w.observe(on_colormethod_change, 'value')


        def dr_name():
            """ Make a unqiue name for this run - for caching purposes
            """
            rv = '{}_{}'.format(self.counttable_name,
                                      dr_w['method'].value)
            d = get_dr_param()
            for k, v in sorted(d.items()):
                if k == 'method': continue
                rv += '__{}_{}'.format(k,v)
            return rv

        def force_refresh():
            "Force to refresh calculations? Or can we use the cache??"
            return drforce_w.value


        @dcache(self, 'dimred_correlating_genes', 'tsv', dr_name, force_refresh)
        def get_top_correlating(*_):
            stats, trans = run_dr_2()
            c1 = np.abs(self.counttable.apply(
                lambda x: pearsonr(x, trans.iloc[:,0])[0], axis=1))
            c2 = np.abs(self.counttable.apply(
                lambda x: pearsonr(x, trans.iloc[:,1])[0], axis=1))
            return pd.DataFrame({0: c1, 1: c2})

        @dcache(self, 'dimred_table', 'mtsv', dr_name, force_refresh)
        def run_dr_2(*_):
            param = get_dr_param()
            method = param['method']

            del param['method']
            counts = self.counttable
            if method == 'pca':
                stats, trans = run_pca(counts)
            elif method == 'nmf':
                stats, trans = run_nmf(counts)
            elif method == 'tsne':
                stats, trans = run_tsne(counts, get_dr_param)
            else:
                raise NotImplemented()
            return stats, trans

        def _get_color_scale():
            clmethod = clmethod_w.value
            plotinfo['color_on'] = clmethod
            #html_w.value = clmethod
            ctype = 'continuous'
            if clmethod == 'gene':
                plotinfo['color_gene'] = '<invalid>'
                gene = clgene_w.value
                if "Invalid gene" in clgene_w.value:
                    return
                if not gene in self.counttable.index:
                    clgene_w.value = 'Invalid gene: {}'.format(gene)
                    return
                plotinfo['color_gene'] = gene
                colv = self.counttable.loc[gene]
            elif clmethod == 'genes of interest':
                gene = clgenesetchoice_w.value
                if ' -- ' in gene:
                    gene = gene.split('--')[0].strip()
                plotinfo['color_gene'] = gene
                plotinfo['color_set'] = clgeneset_w.value
                colv = self.counttable.loc[gene]

            elif clmethod == 'metadata':
                mkey = clmetadata_w.value
                plotinfo['metadata_key'] = mkey

                minfo = self.metadata_info.loc[mkey]['datatype']
                ctype = minfo

                colv = self.get_metadata(mkey)
                colv = colv[sorted(self.counttable.columns)]

            if ctype == 'categorical':
                color_mapper = CategoricalColorMapper(
                    palette=bokeh.palettes.Category10[10] ,
                    factors=list(colv.value_counts().index))
                colorbar_mapper = LinearColorMapper(palette='Viridis256',
                                                    low=0, high=0)
                bcolorbar.color_mapper = colorbar_mapper
                if not bfigure.legend[0].items:
                    bfigure.legend[0].items.append(blegend)

                #print(len(bfigure.legend[0].items))
#                print(bplot.legend)
            else:
                color_mapper = LinearColorMapper(palette="Viridis256",
                                                  low=colv.min(),
                                                  high=colv.max())

                bcolorbar.color_mapper = color_mapper
                if bfigure.legend[0].items:
                    bfigure.legend[0].items.pop()


            bplot.data_source.data['score'] = colv
            bplot.glyph.fill_color['transform'] = color_mapper
            bplot.glyph.line_width=0
            bplot.glyph.line_alpha=0

        def _create_title():
            cmethod = clmethod_w.value
            if cmethod == 'gene':
                cmethod += ':{}'.format(clgene_w.value)
            elif cmethod == 'metadata':
                cmethod += ':{}'.format(clmetadata_w.value)
            elif cmethod == 'genes of interest':
                cmethod += ':{} {}'.format(
                    clgeneset_w.value, clgenesetchoice_w.value)
            return '{}/{}/{}'.format(
                self.experiment_name, self.counttable_name,
                cmethod)

        def _create_scatter_plot(only_color=False):
            if not only_color:
                meta, trans = run_dr_2()
                col1, col2 = trans.columns[:2]
                d1 = trans[col1]; d2=trans[col2]
                title = dr_name().replace('_', ' ')

                m = max(d2.max() - d2.min(), d1.max()-d1.min())
                bplot.data_source.data['size'] = [0.008 * m] * nosamples
                bplot.data_source.data['x'] = d1
                bplot.data_source.data['y'] = d2

            bfigure.title.text = _create_title()
            _get_color_scale()
            #html_w.value = '<pre>{}</pre>'.format(
            #        str(pd.Series(plotinfo)).rsplit("\n", 1)[0])
            bokeh_io.push_notebook(handle=bhandle)

        def run_dr_color(*_):
            """ refresh dr scatter plot - color only """
            drrun_w.button_style = 'info'
            try:
                _create_scatter_plot(only_color=True)
            except:
                drrun_w.button_style = 'danger'
                raise
            drrun_w.button_style = ''

        def run_dr(*_):
            drrun_w.button_style = 'info'
            try:
                _create_scatter_plot()
            except:
                drrun_w.button_style = 'danger'
                raise
            drrun_w.button_style = ''


        drrun_w.on_click(run_dr)
        clgenesetchoice_w.observe(run_dr_color, 'value')
        clmetadata_w.observe(run_dr_color, 'value')

        tab_children = []
        tab_children.append(widgets.VBox([
            ilabel('method', drmethod_w),
            ilabel('perplexity', drperplex_w),
            ilabel('learning rate', drlrate_w),
            ilabel('early exagg.', drearly_w)]))

        tab_children.append(widgets.VBox([
            ilabel('method', clmethod_w),
            ilabel('geneset', [clgeneset_w, clgenesetchoice_w]),
            ilabel('metadata', clmetadata_w),
            ilabel('gene', clgene_w),
        ]))

        tab_children.append(widgets.VBox([
            ilabel('group define', [sl_group_name_w, sl_group_set_w,
                                    sl_group_set_go_w]),
            ilabel('count extract', [sl_groupextractname_w,
                                     sl_group_extract_go_w]),
        ]))

        tabs = widgets.Tab(children=tab_children)
        tabs.set_title(0, 'DimRed')
        tabs.set_title(1, 'Color')
        tabs.set_title(2, 'Selected')
        tabs.selected_index=1
        display(tabs)

        display(widgets.HBox([drrun_w, drforce_w]))
        display(html_w)

        #run a few on_change functions so that all is in sync
        on_colormethod_change()
        run_dr()
