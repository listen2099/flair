#! /usr/bin/env python3
import sys
import traceback

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

# call this as
#    try:
#        if args.subcommand == "check":
#            check_subcommand(args)
#        elif args.subcommand == "design":
#            design_subcommand(args)
#    except Exception as ex:
#        handle_prog_errors(ex, args.debug)
