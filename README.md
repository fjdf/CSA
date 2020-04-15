# Cyclic DNA Sequence Aligner

### Description

This tool finds the optimal rotation for a set of *circular DNA sequences* that are going to be aligned. It is very well suited to apply to _mitochondrial genome alignments_, and also for small plasmids, chloroplasts, circular viruses and bacterial chromosomes.

The best rotation is calculated based on the longest chain of non-repeated blocks that belongs to all the sequences simultaneously. These maximum common blocks are obtained with the help of a _generalized cyclic suffix tree_ data structure.

As of yet, the current version of this tool does not perform the multiple sequence alignment itself. But other existing multiple sequence alignment tools can be used to perform this final task.

### Downloads

[Source code](https://github.com/fjdf/CSA/blob/master/website/SourceCode.zip?raw=true)

[User manual](https://github.com/fjdf/CSA/blob/master/website/UserManual.pdf?raw=true)

[Example sequences](https://github.com/fjdf/CSA/blob/master/website/Examples.zip?raw=true)

### Manual

Use `unzip` to unpack and `make` to compile.

Run with:
```bash
./CSA R <multi-fasta-file>
```

### Reference

If you use *CSA*, please cite:

[Fernandes, F., Pereira, L. & Freitas, A.T. **CSA: An efficient algorithm to improve circular DNA multiple alignment**. BMC Bioinformatics 10, 230 (2009).](https://doi.org/10.1186/1471-2105-10-230)

### Contact

Please send your comments, suggestions, bug reports or questions to:

csatool@kdbio.inesc-id.pt
