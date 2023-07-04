#!/usr/bin/env python
from setuptools import setup
setup(name="rsv-typer",
      version="0.1",
      description="Genotyping RSV samples from nanopore sequencing data",
      url="https://github.com/DiltheyLab/RSVTyper",
      entry_points={
          "console_scripts": [
              "rsvtyper=process_RSV_sample:main"
              ]
          }
      )

