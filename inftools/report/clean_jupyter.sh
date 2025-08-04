# This is a way to clean the jupyter metadata, st. it does not have to be added to git.

jupyter nbconvert infretis_report.ipynb --to notebook --ClearMetadataPreprocessor.enabled=True --ClearOutputPreprocessor.enabled=True
