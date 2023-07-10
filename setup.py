#!/usr/bin/env python
from setuptools import setup
setup(name="rsv-typer",
      version="0.1",
      description="Genotyping RSV samples from nanopore sequencing data",
      url="https://github.com/DiltheyLab/RSVTyper",
      packages = find_packages(),
      #packages=["rsv_typer"],
      #package_dir={"rsv-typer": "rsv_typer"},
      #package_data={"rsv-typer": ["references/*", "primer_scheme/*"]},
      entry_points={
          "console_scripts": [
              "rsv-typer=rsv_typer.process_RSV_sample:main"
              ]
          }
      )


