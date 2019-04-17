# skin-gene-sigs
Gene sets from skin datasets

To use, `git clone` this repo, go into repo, and install:

```
python setup.py develop
```

to run all code, just run

```
runall -o $OUT_DIR
```

all outputs will be in the **OUT_DIR**.

Note that the mapping/quant of the RNA is not in the runall script, the scripts used to align/quant the RNA-seq files is in the src dir.

Note that there are R requirements, check R scripts and install as necessary.