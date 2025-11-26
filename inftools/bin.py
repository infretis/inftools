"""The functions to be used to run infretis via the terminal."""
import pathlib
import sys

import typer

from inftools.misc.binhelper import get_mapper
from inftools.analysis.wham import wham

# define constants
MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
FOLDERS = ["exercises", "tistools", "xyz", "report"]
MAPPER = get_mapper(FOLDERS, MOD_PATH)

# add individual functions
MAPPER["wham"] = wham

app = typer.Typer(
    no_args_is_help=True,
    help="inftools CLI",
    context_settings={"help_option_names": ["-h", "--help"]},
)

# decorating imported mapper functions
for func in MAPPER.values():
    app.command()(func)

# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `inft` call from the command line
