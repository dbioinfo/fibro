This repository contains an RNA-Seq analysis pipeline designed for comparative genomics across multiple species. It handles alignment, ID translation between species, transcript quantification, and differential expression analysis.

```
rnaseq-pipeline
├── main workflow
│   ├── alignment
│   │   └── hisat.sh                      # Alignment with HISAT2
│   ├── id-translation
│   │   ├── translate_ids.R               # Main ID translation script
│   │   ├── translate_ids.helper.sh       # Helper script for ID translation
│   │   └── translate_squirrel_ids.sh     # Specific for squirrel IDs
│   └── quantification
│       ├── filter_gtf.sh                 # GTF ortholog filtering
│       ├── stringtie.sh                  # Transcript assembly/quantification
│       ├── gather_counts.R               # Collect count data
│       └── gather_human_transcripts.R    # Collection of human transcripts
├── analysis
│   └── DESeq_allspecs.R                  # Differential expression with DESeq2
├── docs
│   └── README.md                         # Project documentation
└── .gitignore                            # Ignore large data files, etc.
```
