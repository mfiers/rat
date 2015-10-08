"""
Ipython helper tools
"""

import glob
import os
import zipfile

from jinja2 import Template
from IPython.core.display import HTML
import pandas as pd

fqc_out = """
<table style="border:0px;"><tr style="border:0px; margin: 0px;">
{% for n in names %}
<td  style="border:0px; padding:0px; margin: 0px; padding-left:15px;"><span style="font-family: monospace;">
<a title="{{n}}" href="{{n}}">FQ
{%- for fc in fqcols %}{% set result = summ[n][fc] -%}
{%- if result == "PASS" -%}
<span title="{{fc}}" style="background-color: palegreen;">&nbsp;</span>
{%- elif result == "WARN" -%}
<span title="{{fc}}" style="background-color: gold;">&nbsp;</span>
{%- elif result == "FAIL" -%}
<span title="{{fc}}" style="background-color: tomato;">&nbsp;</span>
{%- endif -%}
{%- endfor -%}</a></span></td>
{% if loop.index % 6 == 0 %}</tr><tr style="border:0px;">{% endif %}
{% endfor %}
</tr>
</table>

"""

def zip_parse(html_file):
    """Helper function to parse a fastqc zip output file

    :param html_file: the html file for which the associated zip
      file needs to be parsed
    :returns: two dictionaries: summary of stats, results of fqc tests
    """
    zpf = html_file.replace('.html', '.zip')
    summhead = []
    datahead = []
    rv_summary = []
    rv_data = []
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
            summhead.append(ls[1])
            rv_summary.append(ls[0])

        for line in data.split("\n"):
            if ">>END_MODULE" in line:
                break
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                continue
            ls = line.strip().split("\t")
            if ls[0] == 'Filename':
                continue
            datahead.append(ls[0])
            rv_data.append(ls[1])

    return dict(zip(summhead, rv_summary)), \
           dict(zip(datahead, rv_data))

FQCOLUMNS = """Per base sequence quality
Per tile sequence quality
Per sequence quality scores
Per base sequence content
Per sequence GC content
Per base N content
Sequence Length Distribution
Sequence Duplication Levels
Overrepresented sequences
Adapter Content
Kmer Content""".split("\n")

def fastqc_display_dir(path, ignore = []):
    """Returns iPython-HTML summarizing a folder of fastqc outputs

    :param path: Path containing a number of fastqc output html/zips
    :returns: ipython.core.display.HTML object
    """
    
    html_files = set(glob.glob(os.path.join(path, '*.html')))
    html_files -= set(ignore)
    html_files = list(sorted(html_files))

    zsumm, zdata = zip(*[zip_parse(x) for x in html_files])
    zsumm = {a: b for (a, b) in zip(html_files, zsumm)}
    zdata = {a: b for (a, b) in zip(html_files, zdata)}
    return pd.DataFrame(zdata).T, HTML(Template(fqc_out).render(
        dict(names=html_files, data=zdata, summ=zsumm, fqcols = FQCOLUMNS)))
