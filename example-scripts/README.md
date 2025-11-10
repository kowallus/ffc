# FFC usage examples

The data used in the example is a Listeria monocytogenes bacterial genome in FASTA format:
* GCA_000585755.1_Lm1823_genomic.fna

The genome was downloaded from [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/pathogens).

The folder contains:
* [ex1.sh](ex1.sh) -  a script performing an FFC (de)compression of a FASTA file
* [ffc_decompressor.py](ffc_decompressor.py) - an LLM-generated decompressor prompted by the FFC File Format Specification v1.0 and some examples, followed with minor corrections.

### Executing the examples on Linux

For the examples to work, the `ffc` binary needs to be located in a programs' folder 
included in the system path, e.g. `/usr/local/bin`.
```
cp ffc /usr/local/bin 
```

Then you can start the test by running:
```
bash ex1.sh
```
