
import logging
import re
import zipfile
import leip

lg = logging.getLogger(__name__)

@leip.arg('statfile')
@leip.command
def cutadapt_stats(app, args):
    """
    parse cutadapt stats
    """
    sf = args.statfile

    def blockparser(handle):
        block = []
        title = 'head'
        
        for line in handle:
            line = line.strip()
            if not line: continue
            if line.startswith('===') and line.endswith('==='):
                yield title, block
                block = []
                title = line.strip('=').strip()
            else:
                block.append(line.strip())
                
    with open(sf) as F:
        for title, block in blockparser(F):
            if title == 'head':
                sec = " ".join(block).split("Finished in")[1].strip().split()[0]
                print("runtime_sec\t" + sec)
            elif title == 'Summary':
                for line in block:
                    k, v = line.split(':', 1)
                    k = k.strip().lower().replace(' ', '_')
                    k = k.replace('(', '').replace(')', '')
                    v = v.strip().split('(')[0].replace(',', '')
                    k = k.replace('_reads', '_reads_count')
                    k = k.replace('_bp', '_bp_count')
                    k = k.replace('_basepairs', '_bp_count')
                    if 'bp' in v:
                        k += '_bp_count'
                        v = v.replace('bp', '').strip()
                    print("%s\t%s" % (k, v))
            else:
                exit()
                print()
                print()
                print(title)
                print("\n".join(block))
                exit()

@leip.arg('statfile')
@leip.command
def fastq_mcf_stats(app, args):
    """
    parse the fastq-mcf stats output
    """
    sf = args.statfile
    with open(sf) as F:
        for line in F:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('Command Line'): continue
            if line.startswith('No adapters found'):
                print("no_adapters_found\ttrue")
                continue
            
            if line.startswith('Scale used'): continue
            if line.startswith('Files'): continue

            if line.startswith('Adapter'):
                a, b = line.split(':', 1)
                a = a.replace('Adapter', '')
                adapter_name, adapter_seq = a.strip().rsplit(' ', 1)
                adapter_name = adapter_name.lower().replace(' ', '_')

                count = re.sub(r'\s*counted (\d+).*$', r'\1', b)
                print('adapter_%s\t%s' % (adapter_name, count))
            elif line.startswith('Clipped'):
                k, r = line.split(':', 1)
                k = k.strip().replace("'", '').replace(' ', '_').lower()
                count, mean, sd = map(lambda x: x.split(':')[1].strip(), r.split(','))
                print('%s_count\t%s' % (k, count))
                print('%s_mean\t%s' % (k, mean))
                print('%s_sd\t%s' % (k, sd))
            elif line.startswith('Trimmed'):
                ls = line.split()
                print('trimmed\t%s' % ls[1])
                print('average_trimmed\t%s' % ls[7])

                
            else:
                #lg.warning('-' * 80)
                #lg.warning(line)
                k, v = map(lambda x: x.strip(), line.split(':', 1))
                k = k.lower().replace(' ', '_')
                print('%s\t%s' % (k, v))


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

