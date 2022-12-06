# AAQuant
Annotation-Agnostic RNA-seq Quantification

`AAQuant` performs annotation-agnostic quantification of RNA-seq reads for a cohort of individuals. The method is described [here](https://www.biorxiv.org/content/10.1101/2022.12.02.518787v1). Each individual is identified with a seven letter alphanumeric string. `AAQuant` parses sequences from a FASTA-file, containing one read per entry, where the first seven characters of the header for each read consist of the identifier of the individual the read belongs to (followed by a space): 

```
>XXXXXXX [auxiliary information]
RNA-SEQ READ 1...
>YYYYYYY [auxiliary information]
RNA-SEQ READ 2...
```

## Installation

`AAQuant` requires [Bifrost](https://github.com/pmelsted/bifrost.git) to be installed at a location where your compiler can find it. After that, it can be compiled using the supplied `Makefile`.


## Running
`AAQuant` can be run using the command

```
./AAQuant -f path/to/fasta.fa -o path/to/output/directory
```

It will yield a number of files:

```
graph.gfa
kmer2unitig
pn2count
pn2count.normalized
```

Of these, `pn2count.normalized` is of most note, as it contains the normalized abundances of each vertex in the graph for each individual in the cohort. Zero abundances are not reported. Each line contains the identifier for the indeividual, the first `k` basepairs in the vertex, and the normalized abundance. The file `kmer2unitig` maps the first `k` basepairs in a vertex to the full sequence of the vertex. `graph.gfa` contains the _de Bruijn_ graph itself, which can be visualized using [Bandage](https://rrwick.github.io/Bandage/).

## Other options
```
  -h,--help                   Print a help message and exit
  -v,--verbose                Print information about the pruning process to stdout.
  -d,--dry_run                Constructs and prunes a DBG, and counts abundances without saving to disk.
  -t,--threads INT            Maximum number of threads to be used.
  -g,--graph TEXT             An existing de Bruijn graph to be pruned.
  -f,--fasta TEXT             A single FASTX file containing RNA reads for association.
  -s,--subtract TEXT ...      DBGs to be subtracted from the target DBG.
  -r,--remove-kmers TEXT ...  Newline separated kmers whose encapsulating unitigs are to be removed from the target DBG.
  -l,--transcriptomic-filter TEXT
                              Assume unitigs that have coverage < max(1, 0.5% of median individual transcriptomic coverage) for an individual are sequencing errors and prune them. arg is a newline separated list of the kmers in all annotated transcripts of the gene.
  -c,--low-coverage-filter    Filter out unitigs with per-individual coverage low w.r.t. the per-individual coverage of the H1 neighborhood around the unitig.
  --ref-fasta-filter TEXT     Fasta file containing one sequence per entry, for example all known cDNA sequences for the organism. We find all kmers within H1 of the constituent kmers and prune the low coverage ones.
  -o,--output TEXT            Directory to which the output files are to be saved.
```
