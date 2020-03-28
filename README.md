# trinotateExtractor

### Easy extraction of transcripts or peptide sequences based on annotated keywords from Trinotate annotation.

#### This Python script is meant to allow you to extract a subset of nucleotide transcripts or peptide sequences generated by Transdecoder based on the Blastx annotation used by Trinotate.  You can feel free to modify this script to extract transscripts based on other Trinotate columns.  

#### I use this script if I want to pull out all of sequences annotated with a specific keyword in the Blastx column of Trinotate.  Based on the keyword, the Trinity isoform names are compiled in a list which is then used to extract the sequences from the Trinity fasta or Transdecoder pep file.  The sequence headers along with sequences are written to a new file.

#### Basic usage is as follows and is written in Python 2. 

```bash
python trinotateExtractor.py -k [keyword] -t [trinotate file] -f [Trinity Fasta or Transdecoder PEP file] -o [output file to write to]
```
#### Flags

	-k  --keyword 	A keyword that will be searched in the Trinotate file.
	-t  --trinotate The Trinotate.xls file generated by Trinotate
	-f  --file      Trinity.fasta or Trinity/Transdecoder.pep file. These files contain the desired sequences.
	-o  --output    The name of the file to write the sequences.  Will retain fasta format.

Here is an example:

```bash
python trinotateExtractor.py -k "Auxin-responsive protein IAA25" -t Trinotate.xls -f Paspalum_urvillei_GENES.fasta -o test.fasta
```
**** Links to [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki), and TransDecoder (https://github.com/TransDecoder/TransDecoder/wiki)
