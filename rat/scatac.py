"""
Module for working with scATAC-seq data.

Includes functions for making and counting aggregates...
"""
import os

from celery import Celery
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth

import rat.celery_core
app = rat.celery_core.get_celery_app()

import os
import pandas as pd
import pysam
from scipy import stats



@app.task
def makeSimulatedCell(numReads, bulkFile, name, directory, debug=False):
    """
    Create a simulated cell.

    Make a simulated cell, with a (roughly) equal number of reads as another
    cell. Sort the resulting bam file. RNG seed is a concatination of the
    number of reads and the name of the cell, therefore using the same
    parameters should return the same result.

    Parameters
    ----------
    numReads : int or list of int
        Number of reads to generate for the simulated cell.
    bulkFile : string
        String of path to BAM/SAM file of bulk sample to generate simulated
        cell from.
    name : string or list of strings
        Name to give simulated cell. Automatically recieves ".sorted.bam"
        suffix and will generate index.
    directory : string
        Directory to put simulated cell BAM file in.
    """
    import random

    bulk = pysam.AlignmentFile(bulkFile, "rb")
    bulkReads = [read for read in bulk]

    if type(numReads) == int:
        numReads = [numReads]

    if type(name) == str:
        name = [name]

    for num, sample in enumerate(name):
        random.seed(str(numReads[num]) + sample)
        with pysam.AlignmentFile(os.path.join(directory,
                                 sample + ".bam"),
                                 "wb", header=bulk.header) as simCell:
            for read in random.sample(bulkReads, numReads[num]):
                simCell.write(read)
        pysam.sort(os.path.join(directory, sample + ".bam"),
                   '-T tmp',
                   '-o',
                   os.path.join(directory, sample + ".sorted.bam"),
                   catch_stdout=False)
        os.remove(os.path.join(directory, sample + ".bam"))
        pysam.index(os.path.join(directory, sample + ".sorted.bam"),
                    catch_stdout=False)

    # from numpy import random
    # bulk = pysam.AlignmentFile(bulkFile, "rb")
    # if type(numReads) == int:
    #     numReads = [numReads]
    # if type(name) == str:
    #     name = [name]
    # for num, cell in enumerate(name):
    #     random.seed(abs(hash(str(numReads[num]) + cell) % (10**9)))
    #     readsToGet.append(set(random.randint(1, bulk.mapped+1, numReads[num])))
    # for num, sample in enumerate(name):
    #     bulk = pysam.AlignmentFile(bulkFile, "rb")
    #     with pysam.AlignmentFile(os.path.join(directory,
    #                              sample + ".sorted.bam"),
    #                              "wb", header=bulk.header) as simCell:
    #         count = 0
    #         for number, read in enumerate(bulk):
    #             count += 1
    #             if debug:
    #                 if count % 1000000 == 0:
    #                     print("Finished %i reads" % (count), flush=True)
    #             if number in readsToGet[num]:
    #                 simCell.write(read)
    #             # if random.randint(0, bulk.mapped) <= numReads:
    #             #     simCell.write(read)
    #     pysam.index(os.path.join(directory, sample + ".sorted.bam"),
    #                 catch_stdout=False)

def makeAggregate(cells, directory, suffix, output):
    """
    Create aggregate sample.

    Make an aggregate bam file from a list of cells, sorts and indexes
    the file for easy use in IGV. Suffix is required to prevent non
    0-padded numbers matching the wrong files. Return final file name.

    Parameters
    ----------
    cells : list
        List of cell names to create aggregate from.
    directory : string
        Directory path with the bam files from each cell.
    suffix : string
        String to match the end of the bam file, use to add file extension
        and to anchor the extension after file numbers - this will prevent
        cell_4 matching cell_4*.
    output : string
        String containing output file location.
    """
    from glob import glob
    cells = set(cells)
    fileList = []
    for cell in cells:
        fileList.append(glob(os.path.join(directory, "*" + cell + suffix))[0])
    pysam.cat("-o", output + ".bam", *fileList, catch_stdout=False)
    pysam.sort(output + ".bam", output + ".sorted", catch_stdout=False)
    pysam.index(output + ".sorted.bam", catch_stdout=False)

    return output + ".sorted.bam"


def makeCounts(aggregates, output, bed):
    """
    Create counts matrix.

    Create a count matrix where each column is a sample and each row is a
    region from a bed file. Unstranded. Return countsTable and write it to
    disk.

    Parameters
    ----------
    aggregates : list
        List of file names of aggregates to count.
    output : string
        String containing output file location.
    bed : string
        String contining location of bed file with regions to count within.
    """
    regions = []
    data = []
    header = []
    with open(bed) as f:
        for line in f:
            line = line.split('\t')
            if int(line[1]) < 1:
                line[1] = "1"
            regions.append(line[0] + ":" + line[1] + "-" + line[2])
    for file in aggregates:
        alignments = pysam.AlignmentFile(file)
        header.append(os.path.basename(file).rstrip(".sorted.bam"))
        counts = []
        for region in regions:
            count = 0
            try:
                for read in alignments.fetch(region=region):
                    count += 1
            except ValueError:
                print(region)
                count = 0
            counts.append(count)
        data.append(counts)
    countsTable = pd.DataFrame(columns=regions, data=data, index=header).T
    countsTable.to_csv(output, header=True, index=True)

    return countsTable


def ptm(template, match, mask=[]):
    """Perform a PTM.

    Perform a PTM, checking that the vector sizes match and return the
    pValue and r^2 value.

    Parameters
    __________
    template : list
        List of values that make up the template.
    match : list
        List of values to compare to the template.
    mask : list of bools : Optional
        List of bools to decide whether to exclude any values.
    """
    if len(template) != len(match):
        raise IndexError("Lengths differ")
    if len(mask) > 0:
        newTemplate = []
        newMatch = []
        for n, i in enumerate(mask):
            if i:
                newTemplate.append(template[n])
                newMatch.append(match[n])
        template = newTemplate
        match = newMatch
    return stats.pearsonr(template, match)

@app.task
def anyTask(fun, *args):
    import dill as pickle
    f = pickle.loads(fun)
    return f(*args)
