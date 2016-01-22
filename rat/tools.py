
import hashlib
import os
import matplotlib.pyplot as plt

import xml.etree.cElementTree as ET

import requests
import pandas as pd

CACHEDIR = os.path.expanduser('~/.cache/rat/biomart_query')


def threewayplot(matrix):
    """
          A (1, 1)



     B (0,0)    C (0, 2)

    

    """
    color = matrix['color']
    datacols = list(sorted(list(set(matrix.columns) - set(['color']))))
    data = matrix[datacols]
    a, b, c = datacols
    print(data.head())
    data = data.subtract(data.min(0), axis=1)
    data = data.divide(data.max(0), axis=1)
    print(data.head())


def _get_cachedir():
    if not os.path.exists(CACHEDIR):
        os.makedirs(CACHEDIR)
    return CACHEDIR


class BiomartQuery():

    def __init__(self, dataset):
        self.root = ET.Element(
            "Query",
            virtualSchemaName="default",
            formatter="TSV",
            header="1",
            uniqueRows="1",
            count="",
            datasetConfigVersion="0.6")
        self.dataset = ET.SubElement(self.root,
                                     "Dataset",
                                     name=dataset,
                                     interface="default")
        self.attribs = []

    def add_filter(self, name, value):
        ET.SubElement(self.dataset, "Filter", name=name, value=value)

    def add_attribute(self, name):
        if name in self.attribs:
            return
        self.attribs.append(name)
        ET.SubElement(self.dataset, "Attribute", name=name)

    def __str__(self):
        return ET.tostring(self.root).decode('UTF-8')


def add_biomart_data(qset, df, column,
                     index_column='external_gene_name',
                     aggfunc=lambda x: ';'.join(map(str, set(x)))):

    if qset.lower() == 'mouse':
        qset = 'mmusculus_gene_ensembl'
    elif qset.lower() == 'human':
        qset = 'hsapiens_gene_ensembl'

    gdq = get_default_biomartquery(qset, attributes=[column])
    gd = biomart_get_query(gdq)
    prep_rv_raw = gd[[index_column, column]].drop_duplicates()
    prep_rv = prep_rv_raw.set_index(index_column).copy()

    if not prep_rv.index.is_unique:
        def _join_fields(r):
            return aggfunc(r)

        prep_rv_g = prep_rv_raw.groupby(index_column)
        prep_rv = prep_rv_g.agg(_join_fields)

    rv = prep_rv.loc[list(df.index)]
    df[column] = rv

    
def get_default_biomartquery(qset, attributes=[]):
    bq = BiomartQuery(qset)
    bq.add_attribute("ensembl_gene_id")
    bq.add_attribute("external_gene_name")
    bq.add_attribute("ensembl_transcript_id")
    for a in attributes:
        bq.add_attribute(a)
    return bq


def get_mouse_biomartquery(attributes=[]):
    bq = BiomartQuery("mmusculus_gene_ensembl")
    bq.add_attribute("ensembl_gene_id")
    bq.add_attribute("external_gene_name")
    bq.add_attribute("ensembl_transcript_id")
    for a in attributes:
        bq.add_attribute(a)
    return bq


def biomart_get_query(query):

    orig_query = query
    colnames = None
    if isinstance(query, str):
        query = " ".join(query.strip().split())
    elif isinstance(query, BiomartQuery):
        colnames = orig_query.attribs
        query = str(query)

    sha1 = hashlib.sha1()
    sha1.update(query.encode('UTF-8'))
    digest = sha1.hexdigest()
    cachefile = os.path.join(_get_cachedir(),
                             '{}.pickle'.format(digest))

    if os.path.exists(cachefile):
        return pd.read_pickle(cachefile)

    url = 'http://ensembl.org/biomart/martservice?query=' + query
    r = requests.get(url, stream=True)
    r.raw.decode_content = True

    
    rv = pd.read_csv(r.raw, sep="\t")
    if not colnames is None:
        #print(rv.columns, colnames)
        rv.columns = colnames
    rv.to_pickle(cachefile)
    return rv
