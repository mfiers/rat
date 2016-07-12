
import hashlib
import json
import os

from path import Path
import requests


TMPDIR = Path('~/data/rat/gorilla/').expanduser()
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)


def run(rank):
    rank = list(rank)
    sha1 = hashlib.sha1()
    for r in rank:
        sha1.update(r.encode())
    uid = sha1.hexdigest()
    rundir = TMPDIR / uid
    if not rundir.exists():
        os.makedirs(rundir)
        
    fwd = rundir / 'forward.rnk'
    rev = rundir / 'reverse.rnk'
    with open(fwd, 'w') as F:
        F.write("\n".join(rank))
    with open(rev, 'w') as F:
        F.write("\n".join(reversed(rank)))

    fwd_files = {'target_file_name': open(fwd, 'rb') }
    rev_files = {'target_file_name': open(rev, 'rb') }
    
    payload = {
        'species': 'HOMO_SAPIENS',
        'run_mode': 'mhg',
        'db': 'all'}

    url = 'http://cbl-gorilla.cs.technion.ac.il/servlet/GOrilla'

    fwd_html_response = rundir / 'forward.response.html'
    rev_html_response = rundir / 'revese.response.html'
    
    r_fwd = requests.post(url, data=payload,
                          files=fwd_files)
    
    r_fwd_raw = r_fwd.text
    
    with open(fwd_html_response, 'w') as F:
        F.write(r_fwd_raw)

    print(uid)
    return r
    
#    print(len(rank))request
