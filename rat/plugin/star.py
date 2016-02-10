
import zipfile
import leip

@leip.arg('starlog')
@leip.command
def star_final_log(app, args):
    """
    parse a samtools bamstat file
    """
    slf = args.starlog
    with open(slf) as F:
        for line in F:
            line = line.strip()
            if not line:
                continue

            ls = line.split('|', 1)
            if len(ls) != 2:
                continue
            
            key = ls[0].strip().replace(' ', '_')\
              .replace('%', 'perc').lower()\
              .replace(',', '').replace('/', '_').replace(':', '')
            val = ls[1].strip()
            print('%s\t%s' % (key, val))



