
import numpy as np
from path import Path
import ipywidgets as widgets
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

from functools import partial

def find_signature_sets(thelenota):
    if not thelenota.geneset_dir is None:
        return [x.basename()
                for x in Path(thelenota.geneset_dir).dirs()]
    else:
        return []


def find_signatures(setname, thelenota):
    gsets = [x.basename().replace('.grp', '')
             for x
             in (Path(thelenota.geneset_dir) / setname).glob('*.grp')]
    return gsets


def signature_score(signature_set, thelenota, select):
    gset =  Path(thelenota.geneset_dir) \
                / signature_set \
                / '{}.grp'.format(select)
    gset = gset.open().read().split()

    # ## use AUCell
    # def namer():
    #     return '{}__{}'.format(thelenota.counttable_name,
    #                         select)

    ct = thelenota.counttable # shortcut
    rnks = ct.rank(method='min', ascending=True,
                   na_option='bottom')
    #normalize to the number of genes
    rnks /= ct.shape[0]
    notgset = set(ct.index) - set(gset)
    q = rnks.loc[gset].dropna().mean()
    r = rnks.loc[notgset].dropna().mean()
#    print( (q-r).min(), (q-r).max())
    score = rnks.loc[gset].dropna().mean() - rnks.loc[notgset].mean()
    return score

    # @dcache(thelenota, 'aucell_score', 'tsv', namer)
    # def run_aucell():
    #     import warnings
    #     import rpy2.robjects as robjects
    #     from rpy2.robjects.packages import importr
    #     from rpy2.robjects import pandas2ri
    #     with warnings.catch_warnings():
    #         warnings.simplefilter("ignore")
    #         pandas2ri.activate()
    #         aucell = importr('AUCell')
    #         genesets = robjects.r.list(gset=gset)
    #         aurank = aucell.AUCell_buildRankings(thelenota.counttable, plotStats=False)
    #         try:
    #             cauc = aucell.AUCell_calcAUC(genesets, aurank)
    #             cauc = pd.DataFrame({select:np.array(cauc)[:,0]},
    #                                 index=thelenota.counttable.columns)
    #         except RRuntimeError:
    #             lg.warning("AUCell fail!")
    #             cauc = pd.DataFrame({select:[0]*len(thelenota.counttable.columns)},
    #                                 index=thelenota.counttable.columns)
    #     cauc = pd.DataFrame({select:np.array(cauc)[:,0]},
    #                         index=thelenota.counttable.columns)
    #     return cauc
    #
    # return run_aucell()[select]


# generic functions retrieving intrinsic counttable & dimred metrics
def intrinsic_generic_options(func, thelenota):
    #appy `func` across the counttable & return the highest scoring values
    rv = thelenota.counttable.apply(func, axis=1).sort_values(ascending=False).head(40)
    rv = rv.reset_index()
    rv = rv.apply(lambda x: '{}  ({:.3g})'.format(x.iloc[0], x.iloc[1]), axis=1)
    return list(rv)


def intrinsic_gene_score(thelenota, select):
    if '(' in select:
        select = select.split('(')[0].strip()
    return thelenota.counttable.loc[select]


class CMAPPER:

    def __init__(self,
                 thelenota,
                 name='map',
                 extra_intrinsic_methods = None):

        #some default settings
        self.cpalette = 'Viridis256'
        self._default_method = 'gene'
        self.thelenota = thelenota
        self.name = name
        self._methods = ['gene', 'metadata']
        self.on_change = lambda *args: None

#            "signatures": (find_signatures, signature_score),

        self.intrinsic_methods = {
            "high exp": (partial(intrinsic_generic_options, np.sum),
                         intrinsic_gene_score),
            "high stdev": (partial(intrinsic_generic_options, np.std),
                           intrinsic_gene_score),
            "high variance": (partial(intrinsic_generic_options, np.var),
                           intrinsic_gene_score),
        }

        for sigset in find_signature_sets(self.thelenota):
            self.intrinsic_methods['S:{}'.format(sigset)] = (
                partial(find_signatures, sigset),
                partial(signature_score, sigset)
            )

        if not extra_intrinsic_methods is None:
            self.intrinsic_methods.update(extra_intrinsic_methods)

        self._methods.extend(self.intrinsic_methods.keys())

        # widgets
        self.w_method = widgets.Dropdown(
            value=self._default_method, options=self._methods,
            layout=widgets.Layout(width='140px'))
        self.w_select_txt = widgets.Text(
            layout=widgets.Layout(width='140px'))
        self.w_select_dd = widgets.Dropdown(
            layout=widgets.Layout(width='140px'))
        self.w_binnify = widgets.Checkbox(
            value=False, layout=widgets.Layout(width='30px'))
        self.w_bins = widgets.IntSlider(
            min=2, max=10, value=4,
            layout=widgets.Layout(display='none', width='150px'))

        self.method = self._default_method

        #select gene with highest sd as default
        tgene = self.thelenota.counttable.std(1).sort_values().tail(1).index[0]
        self.value = tgene


    @property
    def options(self):
        return self._get_options(self.method)

    def _get_options(self, method):
        if method == 'gene':
            return []
        elif method == 'metadata':
            return list(self.thelenota.metadata_info.index)
        elif method in self.intrinsic_methods.keys():
            option_func = self.intrinsic_methods[method][0]
            return option_func(self.thelenota)

        raise Exception("Invalid method", method)

    def prepare_display(self):
        def on_method_change(*args):
            method = self.w_method.value
            # only 'gene' has a free entry text field, make this visible
            if method == 'gene':
                self.w_select_txt.layout.display = 'flex'
                self.w_select_dd.layout.display = 'none'
            else:
                self.w_select_txt.layout.display = 'none'
                self.w_select_dd.layout.display = 'flex'
            self.w_select_dd.options = self._get_options(method)

        def on_bin_change(*args):
            self.w_bins.layout.display = 'flex' \
                if self.w_binnify.value else 'none'

        def on_select_change(*args):
            #I hope late binding
            return self.on_change()

        self.w_binnify.observe(on_bin_change, 'value') # show binslider
        self.w_method.observe(on_method_change, 'value') # show method vars
        self.w_binnify.observe(on_select_change, 'value')
        self.w_bins.observe(on_select_change, 'value')
        self.w_select_txt.on_submit(on_select_change)
        self.w_select_dd.observe(on_select_change, 'value')

        on_method_change()

        return widgets.HBox([
            widgets.Label(self.name,
                          layout=widgets.Layout(width='100px', padding="5px 0px 0px 0px")),
            self.w_method, self.w_select_txt, self.w_select_dd,
             widgets.Label('bin', layout=widgets.Layout(padding="5px 0px 0px 0px")),
             self.w_binnify, self.w_bins])


    def display(self):
        display(self.prepare_display())


    def get_method(self):
        return self.w_method.value

    def set_method(self, method):
        self.w_method.value = method
        self.w_select_dd.options = self._get_options(method)


    method = property(get_method, set_method)

    def set_value(self, value):
        self._value = value
        if self.method == 'gene':
            self.w_select_txt.value = value
        else:
            self.w_select_dd.value = value

    def get_value(self):
        if  self.method == 'gene':
            return self.w_select_txt.value
        return self.w_select_dd.value

    value = property(get_value, set_value)

    @property
    def discrete_intrinsic(self):
        """ Is the underlying data discrete - independent of
            forced binning """
        if self.method == 'gene':
            return False
        elif self.method == 'metadata':
            mtype = self.value
            minfo = self.thelenota.metadata_info.loc[mtype]
            return minfo['datatype'].lower().startswith('cat')
        else:
            return False

    @property
    def discrete(self):
        if self.w_binnify.value:
            return True
        else:
            return self.discrete_intrinsic

    def max(self):
        return self.score.max()

    def min(self):
        return self.score.min()

    @property
    def colormapper(self):
        s = self.score
        if self.discrete:
            return self.categorical_colormapper
        else:
            return self.linear_colormapper

    @property
    def linear_colormapper(self):
        return LinearColorMapper(
            palette=self.cpalette, low=self.min(), high=self.max())

    @property
    def categorical_colormapper(self):
        if self.discrete_intrinsic:
            return CategoricalColorMapper(
                palette=bokeh.palettes.Category20[20] ,
                factors=list(self.score.value_counts().index))
        else:
            return self.linear_colormapper

    @property
    def score(self):

        if self.method in self.intrinsic_methods.keys():
            value_func = self.intrinsic_methods[self.method][1]
            rv = value_func(self.thelenota, self.value)
        else:
            rv = dict(
                gene= self._score_gene,
                metadata= self._score_metadata,
            )[self.method]()

        #ensure the returning score has the same index as the current counttable
        if self.discrete:
            rv = rv.fillna('<na>')
        else:
            rv = rv.fillna(0)

        rv = rv.loc[self.thelenota.counttable.columns]

        if self.w_binnify.value:
            nobins = self.w_bins.value
            nobins = 2 if nobins < 2 else nobins
            nobins = 10 if nobins > 10 else nobins
            return pd.cut(rv, nobins, labels=range(1, nobins+1)).astype('object')
        else:
            return rv


    def _score_metadata(self):
        md = self.thelenota.get_metadata(self.value)
        return md

    def _score_gene(self):
        return self.thelenota.counttable.loc[self.value]
