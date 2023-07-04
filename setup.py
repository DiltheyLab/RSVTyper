#!/usr/bin/env python
from setuptools import setup
setup(name="rsv-typer",
      version="0.1",
      description="Genotyping RSV samples from nanopore sequencing data",
      url="https://github.com/DiltheyLab/RSVTyper",
      packages=["rsv-typer"],
      package_dir={"rsv-typer": "rsv-typer"},
      package_data={"rsv-typer": ["references/*", "primer_scheme/*"],
      entry_points={
          "console_scripts": [
              "rsv-typer=process_RSV_sample:main"
              ]
          }
      )


