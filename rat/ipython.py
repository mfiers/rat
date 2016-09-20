# -*- coding: utf-8 -*-
"""Ipython helper tools

a collection of utililties to use from an ipython notebook

"""

import glob
import os
import pickle
import hashlib
import zipfile
import sys
from jinja2 import Template
from IPython.core.display import HTML
import pandas as pd
import seaborn as sns
from path import Path
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
# from sklearn.cluster import MeanShift, estimate_bandwidth

from ipywidgets import interact, IntSlider, FloatSlider  # NOQA
from functools import partial

import rat.celltypes


def serr(t):
    sys.stderr.write(str(t) + "\n")


class MemoizeMutableDisk:

    def __init__(self, fn):
        self.fn = fn
        self.memo = {}
        self.cachedir = Path('~/.cache/rat/memoize').expanduser()
        self.cachedir.makedirs_p()

    def __call__(self, *args, **kwds):
        strick = pickle.dumps(args, 1) + pickle.dumps(kwds, 1)
        sid = hashlib.sha1()
        sid.update(strick)
        sid = sid.hexdigest()
        picklefile = self.cachedir / ('memoize_%s.pkl' % sid)

        # serr('len strick %d' % len(strick))
        # serr('unique id %s' % sid)

        if sid in self.memo:
            # serr('from memo %s' % sid)
            return self.memo[sid]

        if picklefile.exists():
            # serr('load from pickle: %s' % picklefile)
            try:
                with open(picklefile, 'rb') as F:
                    self.memo[sid] = pickle.load(F)
                # serr('from pickle %s' % sid)
                return self.memo[sid]
            except:
                serr("error loading pickled results")

        # serr('rerun')
        self.memo[sid] = self.fn(*args, **kwds)
        with open(picklefile, 'wb') as F:
            pickle.dump(self.memo[sid], F)
        return self.memo[sid]


@MemoizeMutableDisk
def pca(m, **kwargs):
    """ Run a PCA """
    _pca = PCA(**kwargs).fit(m).transform(m)
    return pd.DataFrame(_pca, index=m.index)


@MemoizeMutableDisk
def tsne(m, **kwargs):
    t = TSNE(**kwargs).fit_transform(m)
    return pd.DataFrame(t, index=m.index)


@MemoizeMutableDisk
def signature(sgn, mat):
    return rat.celltypes.signature(sgn, mat, plot=False)['signature_score']


COLORS = ["#3366cc", "#dc3912", "#ff9900", "#109618", "#990099",
          "#0099c6", "#dd4477", "#66aa00", "#b82e2e", "#316395",
          "#994499", "#22aa99", "#aaaa11", "#6633cc", "#e67300",
          "#8b0707", "#651067", "#329262", "#5574a6", "#3b3eac"]

PALETTE = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
           "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
           "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
           "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]


def tsneer(c, markers, ax=None, palette=PALETTE,
           pcafunc=pca, tsnefunc=tsne,
           cmap=plt.cm.inferno, cax=None,
           colors=None):

    if ax is None:
        if colors is None:
            plt.figure(figsize=(14, 10))
            ax = plt.gca()
        else:
            fig, (ax, cax) = plt.subplots(
                1, 2, sharex=False, sharey=False,
                figsize=(14, 10),
                gridspec_kw=dict(width_ratios=(10, 1)))

    cpca = pcafunc(c.T)
    tsne = tsnefunc(cpca)
    plotted = set()
    scatargs = dict(edgecolor='none')
    norm = plt.matplotlib.colors.Normalize(vmin=colors.quantile(0.1),
                                           vmax=colors.quantile(0.9))

    for mlabel, marker, msize, mcolor, mselect in markers[:-1]:
        t2 = tsne.loc[mselect]
        plotted |= set(t2.index)
        if colors is not None:
            scatargs['c'] = cmap(norm(colors.loc[mselect]))
        else:
            scatargs['c'] = palette[mcolor],
        ax.scatter(t2[0], t2[1], marker=marker,
                   s=msize, label=mlabel,
                   **scatargs)
    mlabel, marker, msize, mcolor = markers[-1]
    not_plotted = list(set(tsne.index) - plotted)
    if colors is not None:
        scatargs['c'] = cmap(norm(colors.loc[mselect]))
    else:
        scatargs['c'] = palette[mcolor],
    ax.scatter(tsne.loc[not_plotted][0],
               tsne.loc[not_plotted][1], label=mlabel,
               marker=marker, s=msize, **scatargs)
    ax.legend(loc='best')

    if colors is not None:
        mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)


def _find_palette(start, rot, gamma, dark, light):
    chp = partial(sns.cubehelix_palette, start=start, rot=rot,
                  gamma=gamma, dark=dark, light=light)
    print(("cmap = sns.cubehelix_palette(start=%(start)s, "
           "rot=%(rot)s, gamma=%(gamma)s, dark=%(dark)s, "
           "light=%(light)s)") % locals())
    sns.palplot(chp())
    return chp(as_cmap=True)


def choose_cmap():
    fs = partial(FloatSlider, min=0.0, step=0.05)
    return interact(
        _find_palette,
        start=fs(max=4.0, value=0.2),
        rot=fs(max=4., value=1.),
        gamma=fs(max=3, value=0.8),
        dark=fs(max=1., value=0.3),
        light=fs(max=1., value=0.6))


fqc_out = """
<table style="border:0px;"><tr style="border:0px; margin: 0px;">
{% for n in names %}
<td  style="border:0px; padding:0px; margin: 0px; padding-left:15px;">
<span style="font-family: monospace;">
<a title="{{n}}" href="{{n}}">FQ
{%- for fc in fqcols %}{% set result = summ[n][fc] -%}
{%- if result == "PASS" -%}
<span title="{{fc}}" style="background-color: palegreen;">&nbsp;</span>
{%- elif result == "WARN" -%}
<span title="{{fc}}" style="background-color: gold;">&nbsp;</span>
{%- elif result == "FAIL" -%}
<span title="{{fc}}" style="background-color: tomato;">&nbsp;</span>
{%- endif -%}
{%- endfor -%}</a></span></td>
{% if loop.index % 6 == 0 %}</tr><tr style="border:0px;">{% endif %}
{% endfor %}
</tr>
</table>

"""


def zip_parse(html_file):
    """Helper function to parse a fastqc zip output file

    :param html_file: the html file for which the associated zip
      file needs to be parsed
    :returns: two dictionaries: summary of stats, results of fqc tests
    """
    zpf = html_file.replace('.html', '.zip')
    summhead = []
    datahead = []
    rv_summary = []
    rv_data = []
    with zipfile.ZipFile(zpf, 'r') as Z:
        summfile = [x for x in Z.namelist()
                    if x.endswith("summary.txt")][0]
        datafile = summfile.replace('summary.txt', 'fastqc_data.txt')

        summ = Z.read(summfile).decode('ASCII')
        data = Z.read(datafile).decode('ASCII')

        for line in summ.split("\n"):
            ls = line.strip().split("\t")
            if len(ls) < 3:
                continue
            summhead.append(ls[1])
            rv_summary.append(ls[0])

        for line in data.split("\n"):
            if ">>END_MODULE" in line:
                break
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                continue
            ls = line.strip().split("\t")
            if ls[0] == 'Filename':
                continue
            datahead.append(ls[0])
            rv_data.append(ls[1])

    return dict(zip(summhead, rv_summary)), \
        dict(zip(datahead, rv_data))

FQCOLUMNS = """Per base sequence quality
Per tile sequence quality
Per sequence quality scores
Per base sequence content
Per sequence GC content
Per base N content
Sequence Length Distribution
Sequence Duplication Levels
Overrepresented sequences
Adapter Content
Kmer Content""".split("\n")


def fastqc_display_dir(path, ignore=[]):
    """Returns iPython-HTML summarizing a folder of fastqc outputs

    :param path: Path containing a number of fastqc output html/zips
    :returns: ipython.core.display.HTML object
    """

    html_files = set(glob.glob(os.path.join(path, '*.html')))
    html_files -= set(ignore)
    html_files = list(sorted(html_files))

    zsumm, zdata = zip(*[zip_parse(x) for x in html_files])
    zsumm = {a: b for (a, b) in zip(html_files, zsumm)}
    zdata = {a: b for (a, b) in zip(html_files, zdata)}

    rv = pd.DataFrame(zdata).T
    rvs = pd.DataFrame(zsumm).T

    rv = pd.concat([rv, rvs], axis=1)

    rv.columns = "fastqc_" + rv.columns.to_series().replace(' ', '_')

    return rv, HTML(Template(fqc_out).render(
        dict(names=html_files, data=zdata, summ=zsumm, fqcols=FQCOLUMNS)))
