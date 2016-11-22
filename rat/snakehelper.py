from path import Path


def pathfinder(glob):
    rv = Path('.').glob(glob)
    rv = [x.abspath() for x in rv]
    return rv

def pathmapper(input, newdir, extfrom, extto):
    newdir = Path(newdir)
    def fix_ext(a, fr, to):
        assert a.endswith(fr)
        return a[:-len(fr)] + to    
    return [newdir / fix_ext(x.basename(), extfrom, extto) for x in input]


