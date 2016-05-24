Example of usage findltr and ltrfamilies tools on part of chr6 of human assembly:

Firstly you should run `findltr` program with appropriate parameters to find putative LTR retrotransposons:

`findltr chr6_example.fasta -o chr6_example.gff3`

As output you will get GFF3 file with annotated retrotransposons. For further analysis will be used LTRdigest program from [genome tools](https://github.com/genometools/genometools). Processing of output of `findltr` in LTRdigest required few additional steps:

`gt suffixerator -db chr6_example.fasta`

`sed 's/\.\t.\t./.\t?\t./' chr6_example.gff3 > chr6_example_converted.gff3`  

Running of LTRdigest:

`gt ltrdigest -hmms ./GyDB_collection/profiles/* -trnas ./hg38-tRNAs.fa  chr6_example_converted.gff3 chr6_example.fasta > chr6_example_ltrdigest.gff3`

File chr6_example_ltrdigest.gff3 contains the same retrotransposons with mapped proteins(protein_match), PBS(primer_binding_site) and RR_tracts.
To include LTR-families information into GFF File used `ltrfamilies` script:

`ltrfamilies chr6_example_ltrdigest.gff ../family_annotation.yaml -o chr6_example_ltrfamilies.gff3`

As output this script will create chr6_example_ltrfamilies.gff3 with family info in Retrotransposon features.
