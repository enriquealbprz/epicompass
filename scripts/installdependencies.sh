#!/bin/bash

set -e  # Exit in case of error

# Colors for better visualization
GREEN='\033[0;32m'
NC='\033[0m'

sudo -v
( while true; do sudo -n true; sleep 60; done ) 2>/dev/null &
KEEP_ALIVE_PID=$!
trap 'kill $KEEP_ALIVE_PID' EXIT

echo -e "${GREEN}>> Updating system...${NC}"
sudo apt update -y && sudo apt upgrade -y
sudo apt install -y software-properties-common curl

echo -e "${GREEN}>> Adding deadsnakes repository...${NC}"
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt update -y

echo -e "${GREEN}>> Installing Python 3.11...${NC}"
sudo apt install -y python3.11 python3.11-venv python3.11-distutils
curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11

echo -e "${GREEN}>> Creating virtual environment...${NC}"
python3.11 -m venv epicompassvenv

echo -e "${GREEN}>> Activating virtual environment...${NC}"
source epiepicompassvenv/bin/activate

echo -e "${GREEN}>> Installing dependencies in venv...${NC}"
pip install --upgrade pip
pip install tqdm dask pandas pyarrow scipy seaborn matplotlib statsmodels

echo -e "${GREEN}>> Library check...${NC}"
python -c "import tqdm, dask, pandas, pyarrow; print('Installation process finished successfully.')"

echo -e "${GREEN}>> To activate the virtual environment run: source epicompassvenv/bin/activate${NC}"