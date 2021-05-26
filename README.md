# coronavirus_recombination
Scripts for "Modular evolution of coronavirus genomes" (in review)


## parser_gb.py

```
parser_gb.py [-h] -input INPUT_FILE -min MIN_LENGTH -max MAX_LENGTH -f
                    FEATURES

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in genbank format
  -min MIN_LENGTH, --min_length MIN_LENGTH
                        Minimal length of sequence. Sequences shorter than min
                        length will not be included in the final dataset
  -max MAX_LENGTH, --max_length MAX_LENGTH
                        Maximal length of sequence. Sequences longer than max
                        length will not be included in the final dataset
  -f FEATURES, --features FEATURES
                        string with qualifiers to retrieve from GenBank
                        annotation, e.g. 'country,host,collection_date'
```

## get_orfs.py

```
get_orfs.py [-h] -input INPUT_FILE -orf_map ORF_MAP_FILE [-r]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in genbank format
  -orf_map ORF_MAP_FILE, --orf_map_file ORF_MAP_FILE
                        Csv-file with short codes for ORFs
  -r, --remove_exceptions
                        Remove exceptions
```


## split_genome.py

```
split_genome.py [-h] -input INPUT_FILE -coord COORD_FILE

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in fasta-format
  -coord COORD_FILE, --coord_file COORD_FILE
                        Csv-file with coordinates of ORFs
```

## trans_alignment.py

```
trans_alignment.py [-h] -input INPUT_FILE -path_mafft PATH_TO_MAFFT

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in fasta-format
  -path_mafft PATH_TO_MAFFT, --path_to_mafft PATH_TO_MAFFT
                        Path to mafft
```

## gradient_color.py

```
gradient_color.py [-h] -taxa TAXA -t1 TREE1 -t2 TREE2

optional arguments:
  -h, --help            show this help message and exit
  -taxa TAXA, --taxa TAXA
                        Text file with taxa names separated by
  -t1 TREE1, --tree1 TREE1
                        Tree file in nexus format
  -t2 TREE2, --tree2 TREE2
                        Tree file in nexus format
```