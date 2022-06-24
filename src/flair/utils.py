#! /usr/bin/env python3
import sys
import traceback
import os
from shutil import which

def handle_prog_errors(ex, debug):
    """Prints error messages without call stack and exit. For expected exceptions """
    print("Error: " + str(ex), file=sys.stderr)
    if debug:
        traceback.print_tb(ex.__traceback__, file=sys.stderr)
    exc = ex.__cause__
    while exc is not None:
        print("caused by: " + str(exc), file=sys.stderr)
        if debug:
            traceback.print_tb(exc.__traceback__, file=sys.stderr)
        exc = exc.__cause__
    exit(1)

def is_in_path(inf):
    """Returns True if file is in current PATH"""
    if which(inf):
        return True
    return False

def file_exists(inf):
    """Returns True if file exists"""
    if os.path.exists(inf):
        return True
    return False

def minimap_ok(name):
    """Make sure minimap2 is usable"""
    if name == 'minimap2':
        if not is_in_path(name):
            return False
        return name
    elif name[-8:] != 'minimap2':
        if name[-1] == '/':
                name += 'minimap2'
        else:
                name += '/minimap2'
    if not os.path.exists(name):
        return False
    return name


# call this as
#    try:
#        if args.subcommand == "check":
#            check_subcommand(args)
#        elif args.subcommand == "design":
#            design_subcommand(args)
#    except Exception as ex:
#        handle_prog_errors(ex, args.debug)
