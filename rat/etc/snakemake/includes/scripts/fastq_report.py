import pandas as pd
import zipfile

print("script")
DATA = {}
cols = []

for zpf in snakemake.input:
    zname = zpf.replace('_fastqc.zip', '')
    summrow = {}

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
            v, k = ls[:2]
            k = k.replace(' ', '_').lower()
            if not k in cols: cols.append(k)
            summrow[k] = v

        def modparser(data):
            mods = data.split(">>END_MODULE")
            for mod in mods:
                lines = mod.strip().split("\n")
                title = [l for l in lines if l.startswith('>>')]
                rmod =  [l for l in lines
                         if (not l.startswith('>>')) and (not l.startswith('#'))]
                if not title: return
                title = title[0]\
                  .lower()\
                  .replace('fail', '')\
                  .replace('warn', '')\
                  .replace('pass', '')\
                  .strip().strip('>')
                yield title, rmod

        def ignore(mod):
            return

        def basic_stats(mod):
            for l in mod:
                if 'Filename' in l: continue
                k, v = l.split('\t')
                k = k.lower().strip().replace(' ', '_')
                k = k.replace('%', 'percent_')
                
                if not k in cols: cols.append(k)
                summrow[k] = v.strip().strip()

        def per_base_seq_qual(mod):
            for l in mod:
                ls = l.split()
                pos = int(ls[0].split('-')[0])
                if pos in [1, 5] or pos % 10 == 0:
                    k = 'mean_qual_at_%d' % pos
                    
                    if not k in cols: cols.append(k)
                    summrow[k] = ls[1].strip()

        for title, module in modparser(data):
            {'basic statistics': basic_stats,
             'per base sequence quality': per_base_seq_qual,
            }.get(title, ignore)(module)

    DATA[zname] = summrow
            
DATA = pd.DataFrame(DATA).T
DATA = DATA[cols]

DATA.to_csv(snakemake.output[0], sep="\t")
