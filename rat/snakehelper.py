from path import Path
import logging

lg = logging.getLogger(__name__)

def pathfinder(glob):
    rv = Path('.').glob(glob)
    rv = [x.abspath() for x in rv]
    return rv

def pathmapper(input, newdir, extfrom, extto):
    newdir = Path(newdir)
    def fix_ext(a, fr, to):
        if not a.endswith(fr):
            lg.critical("{} does not have extensions {} ({})"\
                .format(a, fr, extfrom))
            assert False
        return a[:-len(fr)] + to

    if newdir:
        return [newdir / fix_ext(x.basename(), extfrom, extto) for x in input]
    else:
        return [fix_ext(x.basename(), extfrom, extto) for x in input]
