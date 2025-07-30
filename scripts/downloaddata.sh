#!/bin/bash

set -e  # Exit in case of error

# Colors for better visualization
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "${GREEN}>> Downloading files from repository...${NC}"
mkdir -p ../data/hg38_datasets
wget --recursive \
     --no-parent \
     --no-directories \
     --accept "BSS*.bed.gz" \
     --cut-dirs=5 \
     --directory-prefix=../data/hg38_datasets \
     https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/
echo -e "${GREEN}>> Data downloaded successfully!${NC}"