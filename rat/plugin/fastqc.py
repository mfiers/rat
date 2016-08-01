
import zipfile
import leip

@leip.arg('zipfile')
@leip.command
def fastqc_zip(app, args):
    """
    Scan a fastqc output zip file
    """
    zpf = args.zipfile
    with zipfile.ZipFile(zpf, 'r') as Z:
        summfile = [x for x in Z.namelist()
                    if x.endswith("summary.txt")][0]
        datafile = summfile.replace('summary.txt', 'fastqc_data.txt')
        summ = Z.read(summfile).decode('ASCII')
        data = Z.read(datafile).decode('ASCII')
        datahead = []
        for line in summ.split("\n"):
            ls = line.strip().split("\t")
            if len(ls) < 3:
                continue
            v, k = ls[:2]
            k = k.replace(' ', '_').lower()
            print("%s\t%s" % (k, v))
            if v in ['FAIL', 'WARN', 'PASS']:
                print("%s_score\t%s" %
                      (k, dict(FAIL=0, WARN=1, PASS=2)[v]))

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
                print('%s\t%s' % (k, v.strip()))

        def per_base_seq_qual(mod):
            for l in mod:
                ls = l.split()
                pos = int(ls[0].split('-')[0])
                if pos in [1, 5] or pos % 10 == 0:
                    key = 'mean_qual_at_%d' % pos
                    print('%s\t%s' % (key, ls[1]))

        def overrep_seq(mod):
            for i, l in enumerate(mod):
                ls = l.split()
                key = 'overrepresented_seq_%02d' % (i+1)
                print('%s\t%s\n%s_count\t%s\n%s_percent\t%s' % (
                    key, ls[0], key, ls[1], key, ls[1]))
                if i > 5:
                    break
        for title, module in modparser(data):
            {'basic statistics': basic_stats,
             'per base sequence quality': per_base_seq_qual,
             'overrepresented sequences': overrep_seq,
            }.get(title, ignore)(module)
              
