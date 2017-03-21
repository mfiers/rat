
import pickle
import logging
import pandas as pd

from functools import lru_cache, partial
import ipywidgets as widgets
from traitlets import TraitError

lg = logging.getLogger(__name__)

# ipywidget layout
def ilabel(lab, *wid):
    vbw = '{}px'.format(150 * (1+len(wid)))
    return widgets.HBox(
            [widgets.Label(lab, layout=widgets.Layout(width='100px'))] +
            list(wid))

# save/load functions
def io_load_pkl(filename):
    with open(filename, 'rb') as F:
        return pickle.load(F)


def io_save_pkl(filename, blob):
    with open(filename, 'wb') as F:
        pickle.dump(blob, F)


def io_load_int(filename):
    with open(filename, 'r') as F:
        return int(F.read())


def io_save_int(filename, blob):
    with open(filename, 'w') as F:
        F.write(str(blob))

def io_load_float(filename):
    with open(filename, 'r') as F:
        return float(F.read())


def io_save_float(filename, blob):
    with open(filename, 'w') as F:
        F.write(str(blob))


def io_load_str(filename):
    with open(filename, 'r') as F:
        return F.read()


def io_save_str(filename, blob):
    with open(filename, 'w') as F:
        F.write(blob)

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
    int = (io_load_int, io_save_int),
    pickle = (io_load_pkl, io_save_pkl),
    bool = (io_load_pkl, io_save_pkl),
    float = (io_load_float, io_save_float),
    str = (io_load_str, io_save_str),
    mtsv=(io_load_meta_tsv, io_save_meta_tsv))



def get_cachefile_name(thelenota, category, fmt, namer):
    cache_dir = thelenota.cachedir / category
    cache_dir.makedirs_p()
    if isinstance(namer, str):
        objname = namer
    else:
        objname = namer()
    extension = dict(mtsv='tsv.gz',
                     tsv='tsv.gz').get(fmt, fmt)
    return cache_dir / '{}.{}'.format(objname, extension)

def setcache(thelenota, category, fmt, namer, value):
    cache_dir = thelenota.cachedir / category
    cache_dir.makedirs_p()
    if isinstance(namer, str):
        objname = namer
    else:
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

            if isinstance(namer, str):
                objname = namer
            else:
                objname = namer()

            extension = dict(mtsv='tsv.gz',
                             tsv='tsv.gz').get(fmt, fmt)

            cache_file = cache_dir / '{}.{}'.format(objname, extension)
            io_load, io_save = IOFUNC.get(fmt, (io_load_pkl, io_save_pkl))

            if cache_file.exists() and not force():
                # print('load from cache ({})'.format(force()))
                try:
                    rv = io_load(cache_file)
                    return rv
                except:
                    lg.warning("error loading from cache")
                    lg.warning(str(cache_file))
                    raise
            # print('redo')
            rv = func(*args, **kwargs)
            io_save(cache_file, rv)
            return rv
        return cfunc2
    return cfunc


def cache_widget_value(widget, default, thelenota, name, namer,
                       field_name='value', fmt='str'):
    """
    ensure the widget caches 'value' on disk & recreation
    """

    @dcache(thelenota, 'widget/{}_{}'.format(field_name, name), fmt, namer)
    def get_default_value():
        return default

    def set_default_value(*args):
        setcache(thelenota,  'widget/{}_{}'.format(field_name, name),
                 fmt, namer, getattr(widget, field_name))

    defval = get_default_value()
    # lg.warning('set widget {} to {}'.format(field_name, defval))
    setattr(widget, field_name, defval)
    # widget.value = defval
    widget.observe(set_default_value, field_name)
    return widget


def ccwidget(thelenota, name, wtype, setnamer, default=None, lwidth='200px', **kwargs):
    """ Create a widget with a disk-cached value for persistence
        inbetween instantiating the widget.
    """

    field_name = 'value'
    wmaker = {
        'dropdown' : partial(widgets.Dropdown, options=[]),
        'int' :  widgets.IntSlider,
        'bool' : widgets.Checkbox,
        'float' : widgets.FloatSlider}.get(wtype, widgets.Text)

    wstyle =  widgets.Layout(width=lwidth)
    widget = wmaker(layout=wstyle, **kwargs)

    fmt = dict(int='int', float='float', bool='bool').get(wtype, 'txt')

    @dcache(thelenota, 'widget/{}_{}'.format(field_name, name), fmt, setnamer)
    def get_default_value():
        return default

    def set_default_value(*args):
        setcache(thelenota,  'widget/{}_{}'.format(field_name, name),
                 fmt, setnamer, getattr(widget, field_name))

    try:
        defval = get_default_value()
        if not defval is None:
            setattr(widget, field_name, get_default_value())
    except TraitError:
        pass
        #setattr(widget, field_name, default)

    widget.observe(set_default_value, field_name)
    return widget


    return cache_widget_value(widget, default, thelenota, name)
    cache_widget_value(
        widgets.Text(layout=wstyle), 'thelenota', thelenota,
        'group_define_name', dr_name_simple)

