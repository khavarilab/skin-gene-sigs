from glob import glob
from setuptools import setup
from setuptools import find_packages

if __name__== "__main__":
    setup(name="skin-gene-sigs",
          version="0.1.0",
          description="get gene signatures",
          author="Daniel Kim",
          author_email='danielskim@stanford.edu',
          url="https://github.com/khavarilab/skin-gene-sigs",
          license="MIT",
          install_requires=["numpy", "pandas", "scipy"],
          packages=find_packages(),
          package_data={"src":["data/*.json"]},
          scripts=["bin/runall"] + glob("R/*.R")
    )
