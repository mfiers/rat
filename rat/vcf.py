
from functools import partial

import pandas as pd
import multiprocessing.dummy as mp

# simpel vcf parser


def read(filename):
    ifields = {}
    with open(filename) as F:
        for line in F:
            line = line.strip()
            if not line:
                continue
            if not line.startswith('##'):                
                break
            if line.startswith('##INFO=<ID=') or \
              line.startswith('##FORMAT=<ID='):
                line = line.split('=', 2)[-1]
                iname, inum, itype = line.split(',')[:3]
                inum = inum.replace('Number=', '')
                try:
                    inum = int(inum)
                except ValueError:
                    pass
                itype = itype.replace('Type=', '')
                ifields[iname] = (inum, itype)
        header = line
        assert header.startswith('#')
        header = header.lstrip('#').split("\t")

    def _fix_ifield(iv, inumber, itype):
        if isinstance(inumber, int) or inumber in 'ARG.'.split():
            if itype == 'Integer':
                iv = list(map(int, iv.split(',')))
            elif itype == 'Float':
                iv = list(map(float, iv.split(',')))
            if inumber == 1:
                iv = iv[0]
        return iv

    def _process_info(row, ifields):
        rv = {}
        for ifield in row['INFO'].split(';'):
            if not '=' in ifield:
                ik = ifield
                inum, itype = ifields[ik]
                assert inum == 0 and itype == 'Flag'
                rv['i_' + ik] = True
                continue
            ik, iv = ifield.split('=', 1)
            iv = _fix_ifield(iv, *ifields[ik])
            rv['i_'+ik] = iv
                
        return pd.Series(rv)

    d = pd.read_csv(filename, sep="\t", comment='#', names=header)
#    s = d.iloc[:,9:]
#    d = d.iloc[:,:9]
    
    def _process_sample(nr, ifields):
        n, r = nr
        frmt = r['FORMAT'].split(':')
        srv = []
        for s in r[9:]:
            rv = {}
            for ik, iv in zip(frmt, s.split(':')):
                rv[ik] = _fix_ifield(iv, *ifields[ik])
            srv.append(pd.Series(rv))
        srv = pd.DataFrame(srv, index=r[9:].index).T
        return srv

    _psp = partial(_process_sample, ifields=ifields)
    with mp.Pool() as P:
       samples = P.map(_psp, d.iterrows())
       
    samples = pd.Panel(dict(zip(d.index, samples)))

    info = d.apply(_process_info, ifields=ifields, axis=1)
    d =  pd.concat([d, info], axis=1)
    for ifield, (inum, itype) in ifields.items():
        if inum == 0 and itype == 'Flag' and (('i_' + ifield) in d):
            d['i_' + ifield] = d['i_' + ifield].fillna(False)
            
    return d, samples




