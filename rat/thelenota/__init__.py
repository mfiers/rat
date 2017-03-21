
import collections
import copy
from functools import lru_cache, partial
import io
import itertools
import json
import logging
from multiprocessing.dummy import Pool as mpPool
import pickle
import random
import time

from IPython.display import display
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
    

from path import Path
import requests
from scipy.stats import pearsonr

from rat import celery_core, ctask
from bokeh.colors import RGB
import bokeh.charts as bokeh_chart
import bokeh.io  as bokeh_io
from bokeh.models import ColumnDataSource, CategoricalColorMapper, \
                         LinearColorMapper, ColorBar, LogTicker, \
                         BasicTicker, BasicTickFormatter, CustomJS
import bokeh.models as bmodels
import bokeh.palettes
from bokeh.plotting import figure as bokeh_figure

from rat.thelenota.utils import dcache, cache_widget_value, get_cachefile_name, setcache

USE_CELERY = False

lg = logging.getLogger(__name__)

# CONSTANTS
CONFIG_FILE_NAME = 'config.yaml'
COUNT_TABLE_DIR = '55.normalized_counts'
COUNT_TABLE_EXTENSION = '.tsv'
CACHE_DIR = '99.thelenota_cache'
DIMRED_DIR = '60.dimred'
FIGSIZE=(800,600)
FIGDPI=100
FIGSIZE2 = FIGSIZE[0]/FIGDPI, FIGSIZE[1]/FIGDPI

#
# Helper functions
#

#@lru_cache(32)
def get_project_metadata_info_2(metadata_dir):

    if metadata_dir is None:
        return  pd.DataFrame(columns = ['datatype', 'filename'])
    else:
        metadata_dir = Path(metadata_dir)
        
    all_m = []
    for meta in metadata_dir.glob('*.meta.tsv'):
        #print(meta)
        m = pd.read_csv(meta, sep="\t", header=None, names=['column', 'datatype'],
                        index_col=0)
        m['filename'] = meta.basename()
        all_m.append(m)

    if len(all_m) > 0:
        rv = pd.concat(all_m)
    else:
        rv = pd.DataFrame(columns = ['datatype', 'filename'])

    rv.loc['thelenota'] = dict(datatype='categorical', filename='thelenota.tsv')
    return rv


@lru_cache(32)
def _tsv_loader(filename):
    pklfile = Path(filename + '.pkl')
    if pklfile.exists():
        return pd.read_pickle(pklfile)
    
    rv = pd.read_csv(filename, index_col=0, sep="\t")
    rv = rv[sorted(rv.columns)]
    rv = rv.loc[sorted(rv.index)]

    #no NA's allowed - but to rpevern crashes, fill this
    if rv.shape != rv.dropna().shape:
        minv = rv.min().min()
        repv = 0
        if minv < 0:
            repv = 1.1 * minv
        lg.warning(
            "Counttable {} contains NANs -(NAN->{:.2f})!"\
                .format(filename, repv))
        rv = rv.fillna(minv)

    rv.to_pickle(pklfile)
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
                 counttable,
                 metadata_dir='.',
                 geneset_dir='.',
                 cachedir=None):
        
        self.counttable_file = Path(counttable)
        self.counttable_name = Path(counttable).basename().replace('.tsv','')
        self.metadata_dir = Path(metadata_dir)
        self.geneset_dir = Path(geneset_dir)
        if cachedir is None:
            cachedir = self.counttable_file.dirname() \
                       / '_cache_{}.d'.format(self.counttable_name)
        cachedir.makedirs_p()
        self.cachedir = cachedir
        self.counttable = _tsv_loader(self.counttable_file)
        

    def dr_name_simple(self):
        """Basename for caching"""
        return self.counttable_name


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
    


class ThelenotaOld:

    def __init__(self,
                 basedir = '.',
                 experiment_name=None,
                 geneset_dir=None,
                 count_table_dir=COUNT_TABLE_DIR ):

        #this will contain the indici of selected spots
        self._DR_INDICI_SEL_ = []
        self._DF_INDICI_SEL_ = []
        self._CLUSTER_LABELS = {}

        # WIDGET INSTANCES
        self.basedir = Path(basedir).expanduser().abspath()
        self.count_table_dir = count_table_dir
        if experiment_name is None:
            experiment_name = self.experiments[0]
        self.experiment_w = widgets.Dropdown(
            options=self.experiments, value=experiment_name)
        self.counttable_w = widgets.Dropdown(
            options=self.counttables, value=self.counttables[0])

        if geneset_dir is None:
            self.geneset_dir = self.basedir / 'signatures'
        else:
            self.geneset_dir = Path(geneset_dir)

        self._tsne_params = """angle early_exaggeration learning_rate
            pca_var_cutoff perplexity""".split()

    def __str__(self):
        return "<thelenota {}>".format(self.basedir)
    __repr__ = __str__

    def save_geneset(self, name, gset):
        outfile = self.geneset_dir / '{}.grp'.format(name)
        with open(outfile, 'w') as F:
            F.write("\n".join(sorted(list(gset))))

    @property
    def genesets(self):
        rv = list(self.geneset_dir.glob('*.grp'))
        rv = [x.replace('.grp', '') for x in rv]
        return rv

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
        lg.debug("using set_experiment_name ({})".format(value))
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
    # Precalc DimRed
    #
    def precalc(self, **param):

        assert set(param.keys()) == set(self._tsne_params)

        keys, vals = zip(*param.items())
        vals = [ x if isinstance(x, collections.Iterable) else [x]
                for x in vals]

        allmat = []

        celery_app = celery_core.get_celery_app()
        counts = self.counttable

        cache_dir = self.cachedir /  'dimred_table'
        cache_dir.makedirs_p()

        def run_one(v):

            tparam = dict(zip(keys, v))
            def dr_name():
                """ Make a unqiue name for this run -
                for caching purposes - including tsne parameters
                """
                rv = '{}_tsne'.format(self.counttable_name)
                for k, v in sorted(tparam.items()):
                    rv += '__{}_{}'.format(k,v)
                return rv

            cache_file = get_cachefile_name(self, 'dimred_table', 'mtsv', dr_name)
            if cache_file.exists() and cache_file.getsize() > 0:
                #already done
                return

            stats, trans = run_tsne(counts, lambda: tparam)

            setcache(self, 'dimred_table', 'mtsv', dr_name, (stats, trans))

        pool = mpPool(50)
        allmat = pool.map(run_one, itertools.product(*vals))

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
        wstyle =  widgets.Layout(width='200')

        sl_group_a = self.ccwidget(
            "diffx_group_a", "dropdown", lambda: 'difx', 'thelenota',
            options=cat_meta_grp)

        sl_set_a = self.ccwidget(
            "diffx_group_a_set", "dropdown", lambda: 'difx', None)

        sl_group_b = self.ccwidget(
            "diffx_group_b", "dropdown", lambda: 'difx', 'thelenota',
            options=cat_meta_grp)

        sl_set_b = self.ccwidget(
            "diffx_group_b_set", "dropdown", lambda: 'difx', None)

        self.DIFX_stats = None

        butlay = widgets.Layout(width="120px")
        sl_go = widgets.Button(description='Go', layout=butlay)
        sl_check = widgets.Button(description='Check', layout=butlay)
        sl_save_set = widgets.Button(description='Save', layout=butlay)
        sl_enrichr_link = widgets.Button(description='S&Enrichr', layout=butlay)
        sl_set_name = self.ccwidget('diffx_setname', 'text', lambda: 'difx',
                                    "set_name")

        sl_norm = widgets.Checkbox(value=False)
        html_w = widgets.HTML()
        html_link_w = widgets.HTML()

        nogenes = self.counttable.shape[0]
        colv = [1] * nogenes
        color_mapper = LinearColorMapper(palette="Inferno256", low=-0.3, high=2.5)
        pdata = ColumnDataSource(dict(
            x=[random.uniform(-10, 10) for x in range(nogenes)],
            y=[random.uniform(-10, 10) for x in range(nogenes)],
            mean_a = [3.1] * nogenes,
            mean_b = [-3.1] * nogenes,
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
                                    ('mean a', "@mean_a"),
                                    ('mean b', "@mean_b"),
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
        bhandle = bokeh_io.show(bfigure, notebook_handle=True)

        #def save_geneset
        def go_enrichr(*args):
            idx = list(self._DF_INDICI_SEL_)
            if len(idx) == 0:
                return

            genes = list(self.counttable.index.to_series()\
                        .iloc[idx])

            if len(genes) == 0:
                return

            setname = sl_set_name.value
            self.save_geneset(setname, genes)

            genes_str = "\n".join(genes)
            description = setname
            if not description:
                description = 'a set'
            payload = {
                'list': (None, genes_str),
                'description': (None, description)
            }

            ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
            response = requests.post(ENRICHR_URL, files=payload)
            if not response.ok:
                #print(response)
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
            title = 'QDDE'
            all_samples_set = set(self.counttable.columns)
            logdata = pd.Series()
            normalize = sl_norm.value
            logdata['total cells'] = len(all_samples_set)
            group_a = sl_group_a.value
            group_b = sl_group_b.value

            set_a = sl_set_a.value.split('--')[0].strip()
            set_b = sl_set_b.value.split('--')[0].strip()

            title += ' A:{}/{} B:{}/{}'.format(group_a, set_a, group_b, set_b)

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
            # lfc_l = stats['lfc']

            stats['no_cells_in_a'] = len(sample_a)
            stats['no_cells_in_b'] = len(sample_b)

            stats['name_a'] =  '{}/{}'.format(group_a, set_a)
            stats['name_b'] =  '{}/{}'.format(group_b, set_b)

            #print(stats.head())
            #stats = stats.sort_values(by='lfc', ascending=False)
            #bplot.data_source.data['x'] = stats['mean_a']
            #bplot.data_source.data['y'] = stats['mean_b']
            bplot.data_source.data['x'] = stats['mean_all']
            bplot.data_source.data['y'] = stats['lfc']
            bplot.data_source.data['mean_a'] = stats['mean_a']
            bplot.data_source.data['mean_b'] = stats['mean_b']
            m = stats['mean_all'].max()
            bfigure.title.text = title
            self.DIFX_stats = stats
            bplot.data_source.data['size'] = [0.01 * m] * nogenes
            bokeh_io.push_notebook(handle=bhandle)

        sl_go.on_click(run)
        run()
        # display interface
        display(ilabel('Group A', sl_group_a, sl_set_a))
        display(ilabel('Group B (-A)', sl_group_b, sl_set_b))
        display(ilabel('TPM normalize', sl_norm))
        display(ilabel('Set Name', sl_set_name))
        # tabs = widgets.Tab(children=tab_children)
        # tabs.set_title(0, 'Define sets')
        # display(tabs) sl_enrichr_link html_link_w

        display(widgets.HBox([sl_go, sl_check, sl_save_set, sl_enrichr_link]))
        display(html_w)
        display(html_link_w)


    def tsne_parameter_space(self, **kwargs):
        dcache = self.cachedir / 'dimred_table'
        def fn2dat(f):
            f = f.basename()[:-7]
            d = f.split('__')
            args = dict([x.rsplit('_', 1) for x in d[1:]])
            args['method'] = d[0].split('_')[-1]
            return args
        
        rv = pd.DataFrame(list(map(fn2dat, dcache.glob('*.tsv.gz'))))
        rv = rv[rv['method'] == 'tsne']
        del rv['method']
        rv['perplexity'] = rv['perplexity'].astype(int)
        rv['learning_rate'] = rv['learning_rate'].astype(int)
        rv['angle'] = rv['angle'].astype(float)
        rv['early_exaggeration'] = rv['early_exaggeration'].astype(float)
        rv['pca_var_cutoff'] = rv['pca_var_cutoff'].astype(float)
        for k,v in kwargs.items():
            rv = rv[np.isclose(rv[k], v)]
        return rv
    
    def PEXP(self, **kwargs):
        params = 'angle early_exaggeration perplexity learning_rate pca_var_cutoff'.split()
        for p in params:
            print(p, kwargs[p])

            

    def GENE(self):

        def gcw(*args, **kwargs):
            return self.ccwidget('geneview_' + args[0], *args[1:],
                                 **kwargs, setnamer='gene')
        #widgets
        w = dict(
            gene = gcw('gene_name', 'text'),
            warn = widgets.HTML(),
        )

        # data
        samples = self.counttable.columns
        nosamples = len(samples)

        #gene mapper
        gmapper = CMAPPER(self, name='gene')
        gmapper.method = 'gene'
        defgene = self.counttable.std(1).sort_values().tail(1).index[0]
        gmapper.value = defgene
        w['gene'].value = defgene

        #color mapper
        cmapper = CMAPPER(self, name='color')
        cmapper.method = 'metadata'
        cmapper.value = 'plate'

        #sort order
        smapper = CMAPPER(self, name='sort')
        smapper2 = CMAPPER(self, name='sort 2')

        def get_sort_order():
            sdf = pd.DataFrame(
                {1: smapper.score.sort_values(),
                 2: smapper2.score.sort_values()})
            sdf = sdf.sort_values(by=[1,2])
            return sdf.index

        sort_order = get_sort_order()

        pdata = ColumnDataSource(dict(
            x=pd.Series(range(nosamples), index=sort_order),
            y=gmapper.score.loc[sort_order],
            score=cmapper.score.loc[sort_order],
            ))

        pdata2 = pd.DataFrame(dict(
            x=pd.Series(range(nosamples), index=sort_order),
            y=gmapper.score.loc[sort_order],
            score=cmapper.score.loc[sort_order],
        ))

        bfigure = bokeh_figure(
            plot_width=FIGSIZE[0],
            plot_height=int(FIGSIZE[1]*0.8),
            # tools = bokeh_tools,
            y_range=bmodels.Range1d(gmapper.min(), gmapper.max()),
            toolbar_sticky = False,
            toolbar_location='left',
            title='geneplot')

        #bbar = bokeh_chart.Bar(pdata2, 'x', values='y', group='plate')
        bfigure.title.text_color = 'darkgrey'
        bfigure.title.text_font_style = 'normal'
        bfigure.title.text_font_size= "12px"
        bplot = bfigure.vbar(
                x='x', width=0.5, bottom=0, top='y', source=pdata,
                legend='score',
                color=dict(field='score',
                           transform=cmapper.colormapper))

        blegend = bfigure.legend[0].items[0]
        bcolorbar = ColorBar(
            color_mapper=gmapper.colormapper, ticker=BasicTicker(),
            formatter=BasicTickFormatter(precision=1), label_standoff=10,
            border_line_color=None, location=(0,0))

        null_colorbar_mapper = LinearColorMapper(palette='Inferno256',
                                                 low=0, high=0)

        if cmapper.discrete:
            #remove ColorBar
            bcolorbar.color_mapper = null_colorbar_mapper
        else:
            #remove legend
            bfigure.legend[0].items.pop() #remove legend - we can add this later again

        bfigure.add_layout(bcolorbar, 'right')

        # # display widgets
        display(ilabel('gene', w['gene']))
        cmapper.display()
        smapper.display()
        smapper2.display()
        display(w['warn'])
        #for k, v in bplot.data_source.data.items():
        #        print(k, v.shape, v.dropna().shape)
        bhandle= bokeh_io.show(bfigure, notebook_handle=True)

        #bhandle = bokeh_io.show(bbar, notebook_handle=True)

        def warn(message):
            w['warn'].value = '<b>{}</b>'.format(message)


        def on_gene_change(*args):
            gene = w['gene'].value
            if not gene in self.counttable.index:
                warn("gene {} is not in current counttable".format(gene))
                return

            sortorder = get_sort_order()
            gmapper.value = gene

            yval = gmapper.score.loc[sortorder]
            bplot.data_source.data['y'] = yval
            bokeh_io.push_notebook(handle=bhandle)

        def on_sort_change(*args):
            order = get_sort_order()
            d = bplot.data_source.data
            d['x'].index = order
            d['y'] = d['y'].loc[order]
            d['score'] = d['score'].loc[order]
            bokeh_io.push_notebook(handle=bhandle)

        def on_color_change(*args):
            order = get_sort_order()
            score = cmapper.score
            score = score.loc[order]
            bplot.data_source.data['score'] = score
            bplot.glyph.fill_color['transform'] = cmapper.colormapper
            cm = cmapper.colormapper
            self._cm =cm

            if cmapper.discrete:
                warn('discrete')
                bcolorbar.color_mapper = null_colorbar_mapper
                if not bfigure.legend[0].items:
                    bfigure.legend[0].items.append(blegend)
            else:
                warn('cont')
                bcolorbar.color_mapper = cmapper.colormapper
                if bfigure.legend[0].items:
                    bfigure.legend[0].items.pop()

            bokeh_io.push_notebook(handle=bhandle)

        smapper.on_change = on_sort_change
        smapper2.on_change = on_sort_change
        cmapper.on_change = on_color_change
        w['gene'].on_submit(on_gene_change)
        on_gene_change
        on_color_change
        on_sort_change



