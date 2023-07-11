#!/usr/bin/env python
import setuptools
from setuptools import setup
setup(name="rsv-typer",
      version="0.1",
      description="Genotyping RSV samples from nanopore sequencing data",
      url="https://github.com/DiltheyLab/RSVTyper",
      packages = setuptools.find_packages(),
      #packages=["rsv_typer"],
      #package_dir={"rsv-typer": "rsv_typer"},
      #package_data={"rsv_typer": ["references/*", "primer_scheme/*"]},
      include_package_data=True,
      entry_points={
          "console_scripts": [
              "rsv-typer=rsv_typer.process_RSV_sample:main"
              ]
          }
      )

