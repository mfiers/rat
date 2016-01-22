
from hashlib import sha1
import glob
import os
import pandas as pd
import numpy as np 
import sys
import math
    
from sh import java, shasum
import sh
from path import Path

GSEA_DB = dict(
    goall = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v5.0.symbols.gmt',
)

def _check_gsea(cachedir):
    gd = cachedir.glob('my_analysis.GseaPreranked.*')

    if len(gd) == 0:
        return False

    if not len(gd) == 1:
        print('error dir: %s' % gd)
        return False
    
    gd = gd[0]
    for rf in gd.glob('gsea_report_for_*.xls'):
        d = pd.read_csv(rf, sep="\t")
        if len(d) == 0:
            continue
        d.columns = '_ _ _ size es nes pval qval fwer rank_at_max _ _'.split()
        del d['_']
        d['cachedir'] = cachedir
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                d['slp'] = -np.sign(d['nes']) * np.log10(d['pval'])
            except RuntimeWarning:
                if d['pval'][0] == 0:
                    d['slp'] = d['nes'] * 100
            
        d.loc[d['slp'] < -10, 'slp'] = -10
        d.loc[d['slp'] > 10, 'slp'] = 10
        return d
    
    return False

def ipy_set2rank(*args, **kwargs):
    d = set2rank(*args, **kwargs)
    cd = Path(list(d['cachedir'])[0])
    ccd = cd.glob('my_analysis.*')[0]
    img = ccd.glob('enplot_*.png')[0]
    return img, d



def set2rank(
        rnk, gset,
        cachedir="~/data/rat/gsea_output",
        gseajar="~/bin/gsea2-2.2.1.jar",
        force=False):

    outpath = Path(outpath).expanduser()
    gseajar = Path(gseajar).expanduser()
    
    gset = frozenset(gset)

    rnk = rnk.sort_values()
    txtrnk = "\n".join(['%s\t%s' %(a,b) for a, b in rnk.iteritems()])

    uid = str(abs(hash(gset))) + "_" + str(abs(hash(txtrnk)))

    cachedir = outpath / uid
    
    if force:
        cachedir.rmtree()
        
    cachedir.makedirs_p()
    
    rv = _check_gsea(cachedir)
    if isinstance(rv, pd.DataFrame):
        return rv

    if cachedir.exists():
        cachedir.rmtree()
    cachedir.makedirs_p()
    
    rnkfile = cachedir / 'rank.rnk'
    if not rnkfile.exists():
        with open(rnkfile, 'w') as F:
            F.write(txtrnk)

    gsetfile = cachedir / 'gset.gmx'
    if not gsetfile.exists():
        with open(gsetfile, 'w') as F:
            F.write("gset\nna\n")
            F.write("\n".join(gset))
            F.write("\n")

    cl = ("""-cp %s 
            -Xmx2048m xtools.gsea.GseaPreranked 
            -gmx %s -collapse false 
            -mode Max_probe -norm meandiv 
            -nperm 1000 -rnk %s
            -scoring_scheme weighted -rpt_label my_analysis
            -include_only_symbols true
            -make_sets true -plot_top_x 1
            -rnd_seed timestamp -set_max 9999
            -set_min 5 -zip_report false
            -out %s -gui false """ % (
                gseajar, gsetfile, rnkfile, cachedir)).split()
    java(*cl)
    rv = _check_gsea(cachedir)
    if isinstance(rv, pd.DataFrame):
        return rv
    
    return set2rank(rnk, gset, outpath, gseajar, force=True)


def set2rank2(args):

    rnkfile = Path(args['rnkfile']).expanduser()
    gsetfile = Path(args['gsetfile']).expanduser()
    cachedir = Path(args['cachedir']).expanduser()
    gseajar = Path(args['gseajar']).expanduser()

    def fix_rv(rv):
        rv['gset'] = str(gsetfile.basename()).replace(".grp", '')
        rv['gset_type'] = str(gsetfile.dirname().basename())
        rv['rank'] = str(rnkfile.basename()).replace(".rnk", '')
        rv['rank_type'] = str(rnkfile.dirname().basename())
        del rv['cachedir']
        return rv
        
    if os.path.exists(cachedir):
        rv = _check_gsea(cachedir)
        if rv is False:
            cachedir.rmtree()
        else:
            assert isinstance(rv, pd.DataFrame)
            return fix_rv(rv)

    cachedir.makedirs_p()
        
    cl = ("""-cp %s 
            -Xmx2048m xtools.gsea.GseaPreranked 
            -gmx %s -collapse false 
            -mode Max_probe -norm meandiv 
            -nperm 1000 -rnk %s
            -scoring_scheme weighted -rpt_label my_analysis
            -include_only_symbols true
            -make_sets true -plot_top_x 1
            -rnd_seed timestamp -set_max 9999
            -set_min 5 -zip_report false
            -out %s -gui false """ % (
                gseajar, gsetfile, rnkfile, cachedir)).split()

    #print('run', " ".join(cl))
    java(*cl)#, _out = str(cachedir / 'gsea.out'),
         #_err = str(cachedir / 'gsea.err'))
    
    return fix_rv(_check_gsea(cachedir))



def run(rnk, database,
        outpath='~/data/rat/gsea_output',
        gseajar="~/Desktop/gsea2-2.2.0.jar"):
        
    ctx = dict(
        rnkfile = Path(rnkfile).expanduser().abspath(),        
        gseajar = Path(gseajar).expanduser().abspath(),
        outpath = outpath,    
        setdb=GSEA_DB[database] )

    uid = sha1()

    uid.update(str(hash(frozenset(rnk.head().apply(lambda x: str(x))))))

    rnkshasum = shasum(ctx['rnkfile'])
    uid.update(rnkshasum.encode('UTF-8'))
    uid.update(ctx['setdb'].encode('UTF-8'))
    uid = uid.hexdigest()[:9]

    outpath = Path(outpath).expanduser().abspath() / uid

    if not outpath.exists():
        os.makedirs(outpath)

    ctx['outpath'] = outpath
    if len(glob.glob(outpath / '*/*.xls')) == 0:

        print("# start gsea run (uid=%s)" % uid)
        sys.stdout.flush()

        cl = '''

            -cp {gseajar}           -Xmx2048M  
            xtools.gsea.GseaPreranked 
            -gmx {setdb}            -collapse false 
            -mode Max_probe         -norm meandiv 
            -nperm 1000             -scoring_scheme weighted 
            -rpt_label my_analysis  -include_only_symbols true
            -make_sets true         -plot_top_x 20 
            -rnd_seed timestamp     -set_max 1000 
            -set_min 10             -zip_report false 
            -out {outpath}          -gui false
            -rnk {rnkfile}   '''.strip().format(**ctx).split()
        java(*cl)


    #load results
    posfile = glob.glob(outpath / '*/gsea_report_for_na_pos_*.xls')[0]
    posdata = pd.read_csv(posfile, sep="\t")
    print(posdata.head())
    print(posfile)
    
