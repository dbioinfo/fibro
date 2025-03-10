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

### Shiny app installation

1. Install Docker using snap:
   ```bash
   snap install docker
   ```

2. Pull the Docker image:
   ```bash
   docker pull dbioinfo/shinyfibro
   ```

3. Unzip the provided data:
   ```bash
   unzip shipping_container.zip -d tmp/
   ```
   This will extract the contents into a directory called `tmp/`

4. Run the Docker container:
   ```bash
   docker run -p 3838:3838 -v ./tmp:/data dbioinfo/shinyfibro
   ```
   This command:
   - Maps port 3838 from the container to your local machine
   - Mounts the local `tmp/` directory to `/data` in the container

5. Access the application:
   Open your web browser and navigate to:
   ```
   http://localhost:3838/shinyfibro/
   ```

### Troubleshooting

- If you encounter permission issues with Docker, you may need to run the commands with `su
- Ensure ports are not already in use by another application
 
