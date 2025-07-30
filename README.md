# dna_processing

Holds two functions:

`encode_fasta` encodes specified chromosomes from a FASTA file into a quickly accessible format at the specified directory.

`mutate` mutates a reference chromosome from encoded FASTA data (needs to be specified) and writes it to a txt file wraped at 60 characters that you specify.

## Install:
`pip install .`

## Use:
```Python
import dna_processing

dna_processing.encode_fasta(set(["1", "X"]), "Homo_sapiens.GRCh38.dna.primary_assembly.fa", "test_out")
dna_processing.mutate("1", [(0, "SNV", "A"), (5, "INS", 'G')], 'test_out', 'out.txt')
```
Function signatures can be found in `dna_processing/dna_processing.pyi`