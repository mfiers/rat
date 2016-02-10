
import zipfile
import leip

@leip.arg('bamstat')
@leip.command
def flagstat(app, args):
    """
    parse a samtools bamstat file
    """
    bsf = args.bamstat
    with open(bsf) as F:
        for line in F:
            line = line.strip()
            if not line:
                continue

            
            ls = line.split(None, 3)
            key = ls[3].split('(')[0].strip().replace(' ', '_').lower()
            val = ls[0]
            print('%s\t%s' % (key, val))


@leip.arg('idxfile')
@leip.command
def samtools_idxstats(app, args):
    ixf = args.idxfile
    with open(ixf) as F:
        for line in F:
            line = line.strip()
            if not line: continue

            chrom, length, mapped, unmapped = line.split("\t")
            length, mapped = int(length), int(mapped)

            if chrom == '*':
                chrom = 'star'
            
            print("mapped_%s\t%d" % (chrom, mapped))
            if not chrom == 'star':
                print("mapped_per_mb_%s\t%.4g" % (chrom, (1e6 * (mapped / length))))


            
@leip.arg('statfile')
@leip.command
def samtools_stats(app, args):
    """
    parse a samtools bamstat file
    """
    bsf = args.statfile
    with open(bsf) as F:
        for line in F:
            line = line.strip()
            
            if not line:
                continue

            if line[0] == '#':
                continue

            lid, key, val = line.split('\t', 2)
            key = key.replace(':', '').lower().replace(' ', '_')
            key = key.replace('1st', 'first')
            val = val.split('#')[0].strip()
            if line.startswith('CHK'):
                continue

            if line.startswith('SN\t'):
                print('%s\t%s' % (key, val))
                
            if line.startswith('FFQ'):
                break

