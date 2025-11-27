from typing import Annotated
import typer


def report(
    folder: Annotated[str, typer.Option("-folder", help="Output folder")] = "report",
    ):
    """Copies a jupyter notebook report template to the designated folder.

    Run this command in the root infretis folder containing the data and log files."""
    import shutil
    import os

    """https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package"""
    from importlib import resources as impresources
    # from . import report
    from . import __package__ as report_pkg
    jup_report = "infretis_report.ipynb"
    inp_file = impresources.files(report_pkg) / jup_report

    if not os.path.isdir(folder):
        os.mkdir(folder)

    dest = os.path.join(folder, jup_report)
    if os.path.exists(dest):
        print(f"The file {dest} already exists!!")
        print("Exit without copying file.")
        exit()
    else:
        shutil.copy(inp_file, dest)
        print(f"File copied to {dest}.")

    print(f"run: jupyter notebook {folder}/infretis_report.ipynb")
