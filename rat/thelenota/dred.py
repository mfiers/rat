
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
import uuid

from IPython.display import display
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from path import Path
import requests
from scipy.stats import pearsonr

from rat import celery_core, ctask
import rat.thelenota.utils
from bokeh.colors import RGB
import bokeh.charts as bokeh_chart
import bokeh.io  as bokeh_io
from bokeh.models import ColumnDataSource, CategoricalColorMapper, \
                         LinearColorMapper, ColorBar, LogTicker, \
                         BasicTicker, BasicTickFormatter, CustomJS
import bokeh.models as bmodels
import bokeh.palettes
from bokeh.plotting import figure as bokeh_figure

from rat.thelenota.utils import dcache, cache_widget_value, ccwidget, ilabel
from rat.thelenota import cmapper

from rat import rqtask

USE_CELERY = False

lg = logging.getLogger(__name__)
#
# Dimred methods
#

def run_dbscan(mat, eps):
    from sklearn.cluster import DBSCAN
    db = DBSCAN(eps=eps).fit(mat)
    rv = pd.Series(db.labels_, index=mat.index)
    return rv

def run_dr(mat, name, drfunc, **kwargs):
    dr = drfunc(**kwargs)
    fit = dr.fit(mat.T)
    tra = dr.transform(mat.T)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    meta['method'] = name
    return meta, tra

def celery_wait(job):
    import time
    while not job.ready():
        time.sleep(0.5)
    return job.result

def run_scorpius_dimred(mat):
    import rat.ctask
    
    if USE_CELERY:
        job = rat.ctask.run_scorpius_dimred.delay(mat)
        tra = celery_wait(job)
        if isinstance(tra, Exception):
            raise tra
    else:
        tra = rat.ctask.run_scorpius_dimred(mat)
        
    meta = pd.DataFrame(index=tra.columns)
    meta['method'] = 'scorpius'
    return meta, tra
    
#def run_ica(mat):
#    from sklearn.decomposition import FastICA
#    return run_dr(mat, 'ica', FastICA, n_components=10)


# def run_KernelPCA(mat):
#     from sklearn.decomposition import KernelPCA
#     dr = KernelPCA()
#     fit = dr.fit(mat.T)
#     tra = dr.transform(mat.T)
#     tra = pd.DataFrame(tra, index=mat.columns)
#     meta = pd.DataFrame(index=tra.columns)
#     meta['method'] = name
#     return meta, tra
#     return run_dr(mat, 'kernel_pca', KernelPCA, n_components=10)
# #    from rat.rq_tasks import kernelPCA
# #    from rat.rq_base import get_queue, async
# #    q1 = get_queue('t1')
# #    async(q1, kernelPCA, mat)
    
#     from sklearn.decomposition import KernelPCA
#     return run_dr(mat, 'kernel_pca', KernelPCA, n_components=10)

                       
def run_pca(mat):
    from sklearn.decomposition import PCA
    return run_dr(mat, 'pca', PCA, n_components=20)

def run_nmf(mat):
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=10, solver='cd')
    tmat = mat.T + abs(mat.min().min())
    fit = nmf.fit(tmat)
    tra = fit.transform(tmat)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    #meta['variance_explained'] = fit.explained_variance_ratio_
    meta['method'] = 'pca'
    return meta, tra




class DUMMY:
    pass


class DRED:


    def warn(self, message):
        """
        Show a warning message
        """
        self.w.html.value = '<b>{}</b>'.format(message)
        


    def get_dr_param(self):
        """Get a dictionary with the relevant dimred parameters
        """
        
        method = self.w.drmethod.value
        d = dict(method=method)

        if method in 'pca ica scorpius KernelPCA nmf'.split():
            return d

        d.update(dict(method=self.w.drmethod.value,
                      perplexity=self.w.drperplex.value,
                      angle=self.w.drangle.value,
                      learning_rate=self.w.drlrate.value,
                      early_exaggeration=self.w.drearly.value,
                      pca_var_cutoff=self.w.dr_tsne_pcavar_cutoff.value))
        return d
    
    def dr_name_simple(self):
        """Basename for caching"""
        return '{}'.format(self.thelenota.counttable_name)

                               
    def dr_name_filter(self):
        """ Make a unqiue name for this run - 
        for caching purposes - including filter parameters
        """
        rv = self.dr_name_simple()
        if self.w.drfilter_lower_perc_switch.value:
            rv += '_fltperc_{}'.format(self.w.drfilter_lower_percentage.value)
        if self.w.drfilter_log_switch.value:
            rv += '_log'
        #lg.warning('dr name filter: {}'.format(rv))
        return rv

                           
    def dr_name(self):
        """ Make a unqiue name for this run -
        for caching purposes - including tsne parameters
        """
        rv = self.dr_name_filter()

        method = self.w.drmethod.value
            
        rv += '_{}'.format(method)
        if method == 'tsne':
            d = self.get_dr_param()
            for k, v in sorted(d.items()):
                if k == 'method': continue
                rv += '__{}_{}'.format(k,v)
        #lg.warning('dr name: {}'.format(rv))
        return rv
    

    def __init__(self, thelenota, figsize=(800,600)):

        # DIMRED widgets

        self.thelenota = thelenota
        self.figsize = figsize
        self.counts = None
        
        #widgets
        self.w = DUMMY()

        current_cluster_labels = None

        wstyle =  widgets.Layout(width='200px') #, max_width='120px',
                                 # min_width='120px')

        self.w.drfilter_lower_perc_switch = ccwidget(self.thelenota, 'drfilter_lower_perc_switch',
                                                     'bool', self.dr_name_simple, False, lwidth='40px')
        self.w.drfilter_log_switch = ccwidget(self.thelenota, 'drfilter_log_switch',
                                              'bool', self.dr_name_simple, False, lwidth='40px')
        self.w.drfilter_lower_percentage = ccwidget(self.thelenota, 'drfilter_lower_percentage', 'float',
                                          self.dr_name_simple, 0.2, min=0.01, max=1, value=0.2)

        self.w.drmethod =  ccwidget(self.thelenota,
            'dimred_method', 'dropdown', self.dr_name_simple, 'tsne',
            options='pca ica scorpius KernelPCA tsne nmf'.split())

        self.w.drperplex = cache_widget_value(
            widgets.IntSlider(value=67, min=2, max=99, step=5),
            67, self.thelenota, 'dimred_perplexity', self.dr_name_simple, fmt='int')

        self.w.drangle = cache_widget_value(
            widgets.FloatSlider(value=0.5, min=0.05, max=0.95, step=0.05),
            0.5, self.thelenota, 'dimred_angle', self.dr_name_simple, fmt='float')

        self.w.drlrate = cache_widget_value(
            widgets.IntSlider(value=920, min=20, max=10000, step=50),
            920, self.thelenota, 'dimred_learning_rate', self.dr_name_simple, fmt='int')

        self.w.drearly = cache_widget_value(
            widgets.FloatSlider(value=3.5, min=1, max=20, step=0.5),
            3.5, self.thelenota, 'dimred_early_exag', self.dr_name_simple, fmt='float')

        self.w.dr_tsne_pcavar_cutoff = ccwidget(self.thelenota,
            'tsne_pca_var_cutoff', 'float', self.dr_name_simple, 0.05,
            min=0.01, max=0.2, step=0.01)


        self.w.drrun = widgets.Button(description='GO!')

        self.w.drforce = widgets.Checkbox(description='force',
                        layout = widgets.Layout(width='300px'))


        @dcache(self.thelenota, 'dimred_correlating_genes', 'tsv', self.dr_name)
        def get_dr_correlating(*_):
            stats, trans = run_dr_2()
            c1 = self.thelenota.counttable.apply(
                lambda x: pearsonr(x, trans.iloc[:,0])[0], axis=1)
            c2 = self.thelenota.counttable.apply(
                lambda x: pearsonr(x, trans.iloc[:,1])[0], axis=1)
            return pd.DataFrame({0: c1, 1: c2})


        def get_top_correlating(*_):
            d = get_dr_correlating().abs().sum(1).sort_values(ascending=False)
            d = d.head(40)
            d = d.reset_index()
            d = d.apply(lambda x: '{}  ({:.3g})'.format(x.iloc[0], x.iloc[1]), axis=1)
            return list(d)

        tcolormap = cmapper.CMAPPER(
            self.thelenota, extra_intrinsic_methods={
                "top DR correlate": (get_top_correlating, cmapper.intrinsic_gene_score)
                })

        self.w.sl_group_name = cache_widget_value(
            widgets.Text(layout=wstyle), 'self.thelenota', self.thelenota,
            'group_define_name', self.dr_name_simple)

        self.w.sl_group_set = widgets.Text(layout=wstyle)
        self.w.sl_group_set_go = widgets.Button(description='set', layout=wstyle)

        self.w.sl_groupextractname = widgets.Text(layout=wstyle)
        self.w.sl_group_extract_go = widgets.Button(description='extract')

        self.w.clu_method = ccwidget(self.thelenota,
            'cluster_method', 'dropdown', self.dr_name_simple, 'dbscan',
            options=['dbscan'])
        self.w.clu_dbscan_eps = ccwidget(self.thelenota,
            'clu_dbscan_eps_w', 'float', self.dr_name_simple, 2.5,
            min=0.1, max=10.0, step=0.1)
        self.w.clu_go = widgets.Button(description='Cluster!')
        self.w.clu_name = ccwidget(self.thelenota,
            "cluster_name", "text", self.dr_name_simple, 'cluster')
        self.w.clu_store_go = widgets.Button(description='save')

        plotinfo = {}

        self.w.html = widgets.HTML()

        # data!

        samples = self.thelenota.counttable.columns
        nosamples = len(samples)
        color_mapper = LinearColorMapper(palette="Inferno256", low=-0.3, high=2.5)
        topgene = self.thelenota.counttable.std(1).sort_values().tail(1).index[0]
        colv = list(self.thelenota.counttable.loc[topgene])
        pdata = ColumnDataSource(dict(
            x=[random.uniform(-10, 10) for x in range(nosamples)],
            y=[random.uniform(-10, 10) for x in range(nosamples)],
            desc=list(samples),
            size=[0.3] * nosamples,
            score=colv))
        self.thelenota._pdata = pdata

        # Clustering

        def run_clustering(*args):
            method = self.w.clu_method.value
            stats, trans = run_dr_2()

            if method == 'dbscan':
                labels = run_dbscan(trans, self.w.clu_dbscan_eps.value)
            else:
                lg.warning('Not implemented cluster method: {}'.format(method))
                return

            self.thelenota._CLUSTER_LABELS[self.dr_name()] = labels
            colv = labels
            color_mapper = CategoricalColorMapper(
                palette=bokeh.palettes.Category20[20] ,
                factors=list(colv.value_counts().index))
            colorbar_mapper = LinearColorMapper(palette='Inferno256',
                                                low=0, high=0)
            bcolorbar.color_mapper = colorbar_mapper
            if not bfigure.legend[0].items:
                bfigure.legend[0].items.append(blegend)
            bplot.data_source.data['score'] = colv
            bplot.glyph.fill_color['transform'] = color_mapper
            bplot.glyph.line_width=0
            bplot.glyph.line_alpha=0
            bokeh_io.push_notebook(handle=bhandle)

        def store_current_cluster(*args):
            labels = self.thelenota._CLUSTER_LABELS.get(self.dr_name())

            if labels is None:
                self.w.html.value = "no cluster defined"
                return

            cname = self.w.clu_name.value
            labels = labels.apply(lambda x: '{}_{}'.format(cname, x))
            outfile = self.thelenota.metadata_dir / '{}.tsv'.format(cname)
            moutfile = self.thelenota.metadata_dir / '{}.meta.tsv'.format(cname)

            tmeta = pd.DataFrame({cname: labels})
            tmeta.to_csv(outfile, sep="\t")
            with open(moutfile, 'w') as F:
                F.write("{}\tcategorical\n".format(cname))
            self.w.html.value = 'saved cluster to {}'.format(outfile)

        self.w.clu_store_go.on_click(store_current_cluster)
        self.w.clu_go.on_click(run_clustering)

        # GROUP SET
        def define_group(*args):
            groupset = self.w.sl_group_name.value
            groupname = self.w.sl_group_set.value
            self.thelenota.selected = list(
                self.thelenota.counttable.columns.to_series()\
                        .iloc[list(self.thelenota._DR_INDICI_SEL_)])

            if groupset in self.thelenota.metadata_info.index:
                tmeta = self.thelenota.get_metadata(groupset)
            else:
                tmeta = pd.Series('na', index=self.thelenota.counttable.columns)

            tmeta.loc[self.thelenota.selected] = groupname

            self.thelenota.metadata_dir.makedirs_p()
            outfile = self.thelenota.metadata_dir / '{}.tsv'.format(groupset)
            moutfile = self.thelenota.metadata_dir / '{}.meta.tsv'.format(groupset)

            tmeta = pd.DataFrame({groupset: tmeta})
            tmeta.to_csv(outfile, sep="\t")
            with open(moutfile, 'w') as F:
                F.write("{}\tcategorical\n".format(groupset))
            run_dr_color()

        self.w.sl_group_set_go.on_click(define_group)

        # GROUP EXTRACT
        #sl_groupextractname_w
        def extract_group(*args):
            groupname = self.w.sl_groupextractname.value
            self.thelenota.selected = list(
                self.thelenota.counttable.columns.to_series()\
                        .iloc[list(self.thelenota._DR_INDICI_SEL_)])


            outfile =  self.thelenota.counttable_file.dirname() \
                       / '{}__{}.tsv'.format(self.thelenota.counttable_name, groupname)

            newcounts = self.thelenota.counttable[self.thelenota.selected]
            newcounts.to_csv(outfile, sep="\t")
            self.w.html.value = "save new count table to {}".format(outfile)

        self.w.sl_group_extract_go.on_click(extract_group)


        def force_refresh():
            "Force to refresh calculations? Or can we use the cache??"
            return self.w.drforce.value
            
        @dcache(self.thelenota, 'filtered', 'pickle', self.dr_name_filter, force_refresh)
        def filter_counts(*_):
            counts = self.thelenota.counttable

            if self.w.drfilter_log_switch.value:
                minval = 0.5 * counts[counts>0].min().min()
                lg.warning('log tranform log10(x+{:.4g})'.format(minval))
                counts = np.log10(counts + minval)
                
            
            if self.w.drfilter_lower_perc_switch.value:
                perc = self.w.drfilter_lower_percentage.value
                csum = counts.sum(1)
                counts = counts[csum >= csum.quantile(perc)]
                
            return counts
        
        
        @dcache(self.thelenota, 'dimred_table', 'mtsv', self.dr_name, force_refresh)
        def run_dr_2(*_):
            param = self.get_dr_param()
            method = param['method']

            del param['method']

            if method == 'pca':
                stats, trans = qrtask.pca.sync(self.counts)
            elif method == 'ica':
                stats, trans = rqtask.ica.sync(self.counts)
            elif method == 'KernelPCA':                
                stats, trans = rqtask.kernelPCA.sync(self.counts)
            elif method == 'scorpius':
                stats, trans = rqtask.scorpius_dr.sync(self.counts)
            elif method == 'nmf':
                stats, trans = rqtask.nmf.sync(self.counts)
            elif method == 'tsne':
                stats, trans = rqtask.tsne.sync(self.counts, self.get_dr_param())
            else:
                raise NotImplemented
            return stats, trans

        def set_color_scale():
            score = tcolormap.score
            method = tcolormap.method

            self.warn('set color scale {}/{}'.format(method, tcolormap.value))

            color_mapper = tcolormap.colormapper
            bplot.glyph.fill_color['transform'] = color_mapper
            bplot.data_source.data['score'] = score

            if not tcolormap.discrete:
                bcolorbar.color_mapper = color_mapper

            if tcolormap.discrete:
                if not bfigure.legend[0].items:
                    bfigure.legend[0].items.append(blegend)
            else:
                if bfigure.legend[0].items:
                    bfigure.legend[0].items.pop()

            bplot.glyph.line_width=0
            bplot.glyph.line_alpha=0

        def _create_title():
            cmethod = '{}/{}'.format(tcolormap.method,
                                      tcolormap.value)
            
            if self.w.drfilter_lower_perc_switch.value:
                cmethod += '/low%{}'.format(self.w.drfilter_lower_percentage.value)
            if self.w.drfilter_log_switch.value:
                cmethod += '/log'
            cmethod += '/{}/{}'.format(*self.counts.shape)
            
            return '{}/{}'.format(
                self.thelenota.counttable_name,
                cmethod)

        def _create_scatter_plot(only_color=False):
            if (self.w.drfilter_lower_perc_switch.value) or (self.w.drfilter_log_switch.value):
                self.counts = filter_counts()
            else:
                self.counts = self.thelenota.counttable

            self.warn("count table: {}".format(str(self.counts.shape)))
            
            if not only_color:
                meta, trans = run_dr_2()
                self.thelenota.dr = trans

                self.thelenota._dr_meta = meta
                self.thelenota._dr_trans = trans
                col1, col2 = trans.columns[:2]
                d1 = trans[col1]; d2=trans[col2]
                title = self.dr_name().replace('_', ' ')

                m = max(d2.max() - d2.min(), d1.max()-d1.min())
                bplot.data_source.data['size'] = [0.008 * m] * nosamples
                bplot.data_source.data['x'] = d1
                bplot.data_source.data['y'] = d2

            bfigure.title.text = _create_title()
            set_color_scale()
            bokeh_io.push_notebook(handle=bhandle)


        def run_dr_color(*_):
            """ refresh dr scatter plot - color only """
            self.w.drrun.button_style = 'info'
            try:
                _create_scatter_plot(only_color=True)
            except:
                self.w.drrun.button_style = 'danger'
                raise
            self.w.drrun.button_style = ''

        def run_dr(*_):
            self.w.drrun.button_style = 'info'
            try:
                _create_scatter_plot()
            except:
                self.w.drrun.button_style = 'danger'
                raise
            self.w.drrun.button_style = ''

        tcolormap.on_change = run_dr_color
        self.w.drrun.on_click(run_dr)
        # clgenesetchoice_w.observe(run_dr_color, 'value')
        # clmetadata_w.observe(run_dr_color, 'value')


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
            bmodels.HoverTool(tooltips=[
                        ("(x,y)", "($x, $y)"),
                        ("score", "$score"),
                        ("desc", "@desc"),
                    ]),            bmodels.BoxSelectTool(callback = select_callback),
            bmodels.PanTool(),
            bmodels.WheelZoomTool(),
            bmodels.BoxZoomTool(),
            bmodels.LassoSelectTool(callback = select_callback),
            bmodels.SaveTool(),
            bmodels.ResetTool(),
            bmodels.HelpTool(),
        ]

        bfigure = bokeh_figure(plot_width=figsize[0],
                               plot_height=figsize[1],
                               tools = bokeh_tools,
                               toolbar_sticky = False,
                               toolbar_location='left',
                               title='dimredplot')
        bfigure.title.text_color = 'darkgrey'
        bfigure.title.text_font_style = 'normal'
        bfigure.title.text_font_size= "12px"
        self.thelenota._bfigure = bfigure

        bplot = bfigure.circle(
                x='x', y='y', radius='size', source=pdata,
                legend='score',
                color=dict(field='score', transform=color_mapper))

        self.thelenota._bplot = bplot

        bcolorbar = ColorBar(
            color_mapper=tcolormap.colormapper, ticker=BasicTicker(),
            formatter=BasicTickFormatter(precision=1), label_standoff=10,
            border_line_color=None, location=(0,0))



        bfigure.add_layout(bcolorbar, 'right')
        blegend = bfigure.legend[0].items[0]

        bhandle= bokeh_io.show(bfigure, notebook_handle=True)

        tab_children = []
        tab_children.append(widgets.VBox([
            ilabel('filter sum%',
                   self.w.drfilter_lower_perc_switch,
                   self.w.drfilter_lower_percentage),
            ilabel('log tranform', self.w.drfilter_log_switch),
            ilabel('method', self.w.drmethod),
            ilabel('perplexity', self.w.drperplex),
            ilabel('angle', self.w.drangle),
            ilabel('learning rate', self.w.drlrate),
            ilabel('early exagg.', self.w.drearly),
            ilabel('PCA var. cutoff', self.w.dr_tsne_pcavar_cutoff),
            widgets.HBox([self.w.drrun, self.w.drforce]),
            ]))

        tab_children.append(widgets.VBox([
            tcolormap.prepare_display()
        ]))

        tab_children.append(widgets.VBox([
            ilabel('method', self.w.clu_method),
            ilabel('dbscan:eps', self.w.clu_dbscan_eps),
            ilabel('store', self.w.clu_name, self.w.clu_store_go),
            self.w.clu_go
        ]))

        tab_children.append(widgets.VBox([
            ilabel('group define', self.w.sl_group_name, self.w.sl_group_set,
                                    self.w.sl_group_set_go),
            ilabel('count extract', self.w.sl_groupextractname,
                                     self.w.sl_group_extract_go),
        ]))

        tabs = widgets.Tab(children=tab_children)
        tabs.set_title(0, 'DimRed')
        tabs.set_title(1, 'Color')
        tabs.set_title(2, 'Cluster')
        tabs.set_title(3, 'Select')
        tabs.selected_index=0
        display(tabs)

        display(self.w.html)

        #run a few on_change functions so that all is in sync
        #on_colormethod_change()
        run_dr()

