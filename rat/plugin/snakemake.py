import logging

from path import Path
import leip

lg = logging.getLogger(__name__)

@leip.subparser
def sm(app, args):
    pass


@leip.flag('-F', '--full')
@leip.subcommand(sm, 'list')
def smlist(app, args):
    """ listaa known makefiles """
    mkpath = Path(__file__).dirname().dirname() / 'etc' / 'snakemake'
    mkfiles = mkpath.glob('*.sm')
    for f in mkfiles:
        if args.full:
            print(f)
        else:
            print(f.basename().replace('.sm', ''))

@leip.command
def smhelp(app, args):
    print('hi')


@leip.arg('-j', '--threads', type=int, default=8)
@leip.arg('snakemake_template')
@leip.subcommand(sm, 'prep')
def sminstall(app, args):
    
    """ list known makefiles """
    smpath = Path(__file__).dirname().dirname() / 'etc' / 'snakemake'
    smfile = smpath / ('%s.sm' % args.snakemake_template)
    smconf = smpath / ('%s.yaml' % args.snakemake_template)
    
    if not smfile.exists():
        print('cannot find snakemake file')
        exit(-1)

    with open('snake.sh', 'w') as F:
        F.write("#/bin/bash\n")
        F.write('snakemake -j %d' % (args.threads))
        F.write(' -s %s' % (smfile))
        #if args.target:
        #    F.write(' %s' % args.target)
        F.write(" $@\n")
    Path('snake.sh').chmod('ug+x')
    
    if not smconf.exists():
        return

    if Path('config.yaml').exists():
        lg.warning("config file exists, creating backup in:")
        lg.warning("  config.%s.yaml" % args.snakemake_template)
        smconf.copy('config.%s.yaml' % args.snakemake_template)
    else:
        smconf.copy('config.yaml')
    return

        
    
    
    local_config_file = Path('./config.yaml')
    
    if local_config_file.exists():
        lg.warning("local config exists - not overwriting")
        fix_config = False
        with open(config_file) as F:
            for line in F:
                if line.startswith('RAT_TEMPLATE'):
                    ct = line.split('=')[1].strip()
                    if not ct == args.makefile_name:
                        lg.warning("Switching template?")
                        lg.warning("Backing up config file to config.mk.backup")
                        config_file.move("config.mk.backup")
                    else:
                        break
            else:
                fix_config = True
        if fix_config:
            lg.warning("Config exists, but not template is not defined - fixing")
            config_file.move("config.mk.backup")
            with open("config.mk.backup") as F, open(config_file, 'w') as G:
                G.write("RAT_TEMPLATE := %s\n" % args.makefile_name)
                for line in F:
                    G.write(line)

    if not config_file.exists():
        with open(config_file, 'w') as F:
            F.write('RAT_TEMPLATE := %s\n' % args.makefile_name)
            with open(mkfile) as G:
                for line in G:
                    line = line.strip()
                    if not line: continue
                    if "end parameters" in line: break
                    if line.startswith('#?'):
                        F.write("\n# %s\n" % line[2:].strip())
                        F.write("# " + G.readline())

        
    if Path('Makefile').exists() and not args.force:
        print("Makefile exists - maybe try -f?")
        exit(-1)
    else:
        if Path('Makefile').exists():
            Path('Makefile').remove()
        if args.link:
            mkfile.symlink('Makefile')
        else:
            mkfile.copy('Makefile')

    if Path("rat_core.mk").exists() and not args.force:
        # ignore
        pass
    else:
        if Path("rat_core.mk").exists():
            Path("rat_core.mk").remove()
        if args.link:
            mkcore.symlink("rat_core.mk")
        else:
            mkcore.copy("rat_core.mk")
        
    for fn in smpath.glob('%s__*' % args.makefile_name):
        newfn = Path(fn.rsplit(args.makefile_name + '__')[-1])
        if newfn.exists():
            if args.force:
                newfn.remove()
            else:
                print("%s exists - maybe try -f?")
                exit(-1)
        if args.link:
            fn.symlink(newfn)
        else:
            fn.copy(newfn)
