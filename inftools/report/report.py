from typing import Annotated
import typer

"""https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package"""
from importlib import resources as impresources
# from . import report
from . import __package__ as report_pkg
JUP_REPORT = "infretis_report.ipynb"
# INP_FILE = impresources.files(report) / JUP_REPORT
INP_FILE = impresources.files(report_pkg) / JUP_REPORT

def report(
    folder: Annotated[str, typer.Option("-folder", help="Output folder")] = "report",
    ):
    """Copies a jupyter notebook report template to the designated folder.

    Run this command in the root infretis folder containing the data and log files."""
    import shutil
    import os

    if not os.path.isdir(folder):
        os.mkdir(folder)

    dest = os.path.join(folder, JUP_REPORT)
    if os.path.exists(dest):
        print(f"The file {dest} already exists!!")
        print("Exit without copying file.")
        exit()
    else:
        shutil.copy(INP_FILE, dest)
        print(f"File copied to {dest}.")

    print(f"run: jupyter notebook {folder}/infretis_report.ipynb")
