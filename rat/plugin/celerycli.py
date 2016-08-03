
import glob
import re

from path import Path
import leip

@leip.subparser
def c(app, conf):
    """
    Celery scripts :)
    """
    pass


@leip.subcommand(c, "shutdown")
def c_shutdown(app, args):
    import rat.ctask
    rat.ctask.app.control.broadcast('shutdown')


@leip.subcommand(c, "ping")
def c_ping(app, args):
    import rat.ctask
    print(rat.ctask.app.control.ping())


STAR_SE_SCRIPT = """
module load STAR/2.4.2a-foss-2014a
module load SAMtools/1.2-foss-2014a
cd {outputdir}
echo $PWD

STAR \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 10 \
  --genomeDir {stardb} \
  --outFileNamePrefix {sample} \
  --chimSegmentMin 18 \
  --readFilesCommand zcat \
  --readFilesIn {inputfile}

"""

@leip.arg('--sample_re', default=r'(.*)\.fastq\.gz')
@leip.arg('outputdir')
@leip.arg('stardb')
@leip.arg('inputglob')
@leip.subcommand(c, "star_se")
def c_star_se(app, args):
    import rat.ctask
    sample_re = re.compile(args.sample_re)
    scripts = []
    for infile in Path('.').glob(args.inputglob):
        rematch = sample_re.search(infile.basename())
        if not rematch:
            print('warning: no sample found for %s' % infile)
            exit()
        sample = rematch.groups()[0]
        script = STAR_SE_SCRIPT.format(
            inputfile = Path(infile).expand().abspath(),
            stardb = Path(args.stardb).expand().abspath(),
            outputdir = Path(args.outputdir).expand().abspath(),
            sample = sample )
        scripts.append(script)

    for s in scripts:
        rat.ctask.runscript.delay(s)


FASTQCSCRIPT = """
module load FastQC/0.11.4
cd {outputdir}
lnfile="{sample}.fastq.gz"
[[ -e $lnfile ]] || ln -s {inputfile} $lnfile
echo "Processing {sample}"
fastqc -o . $lnfile
"""

@leip.arg('sample_re', default=r'(.*)\.fastq\.gz', nargs='?')
@leip.arg('outputdir')
@leip.arg('inputglob')
@leip.subcommand(c, "fastqc")
def c_fastqc(app, args):
    import rat.ctask
    sample_re = re.compile(args.sample_re)
    scripts = []
    for infile in Path('.').glob(args.inputglob):
        rematch = sample_re.search(infile.basename())
        if not rematch:
            print('warning: no sample found for %s' % infile)
            exit()

        script = FASTQCSCRIPT.format(
            inputfile = Path(infile).expand().abspath(),
            outputdir=Path(args.outputdir).expand().abspath(),
            sample=rematch.groups()[0])
        scripts.append(script)

    for s in scripts:
        rat.ctask.runscript.delay(s)
