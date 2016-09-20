import logging

from path import Path
import leip

lg = logging.getLogger(__name__)

@leip.subparser
def mf(app, args):
    pass


@leip.subcommand(mf, 'list')
def mf_list(app, args):
    """ list known makefiles """
    mkpath = Path(__file__).dirname().dirname() / 'etc' / 'makefiles'
    mkfiles = mkpath.glob('*.mk')
    mkfiles = [x.basename().replace('.mk', '') for x in mkfiles]
    for f in mkfiles:
        print(f)


@leip.command
def mf_help(app, args):
    print('hi')


@leip.flag('-f', '--force')
@leip.flag('-l', '--link')
@leip.arg('makefile_name')
@leip.subcommand(mf, 'install')
def mf_install(app, args):
    
    """ list known makefiles """
    mkpath = Path(__file__).dirname().dirname() / 'etc' / 'makefiles'
    mkfile = mkpath / ('%s.mk' % args.makefile_name)
    mkcore = mkpath / ('_rat_core.mk')
    config_file = Path('./config.mk')
    
    if config_file.exists():
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

    if not mkfile.exists():
        print('cannot find makefile')
        
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
        
    for fn in mkpath.glob('%s__*' % args.makefile_name):
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
