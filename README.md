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

Of these, `pn2count.normalized` is of most note, as it contains the normalized abundances of each vertex in the graph for each individual in the cohort. Zero abundances are not reported. Each line contains the identifier for the indeividual, the first `k` basepairs in the vertex, and the normalized abundance. The file `kmer2unitig` maps the first `k` basepairs in a vertex to the full sequence of the vertex`. `graph.gfa` contains the _de Bruijn_ graph itself, which can be visualized using [Bandage](https://rrwick.github.io/Bandage/).
