
from hashlib import sha1
import glob
import os
import pandas as pd
import sys
    
from sh import java, shasum
import sh
from path import Path

GSEA_DB = dict(
    goall = 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v5.0.symbols.gmt',
)

        
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
    
