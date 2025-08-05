# README
## EpiCompass 1.0

_EpiMap repository analysis pipeline_

### Description
[EpiMap](https://compbio.mit.edu/epimap/) (Epigenome Integration across Multiple Annotation Projects) is a public repository that provides genome-wide epigenetic state annotations across 859 biosamples, using an 18-state chromatin model.

This pipeline, **EpiCompass**, was developed to download, process, and interpret EpiMap data through classification, matrix generation, and plotting of chromatin state distributions.

### Requirements
- Linux
- Python 3.11 and virtualenv
- Python libraries: `dask`, `matplotlib`, `pandas`, `pyarrow`, `scipy`, `seaborn`, `statsmodels`, `tqdm`
- Git
- Command-line tools: `wget`

### Installation
1. Clone this repository:
```bash
git clone https://github.com/enriquealbprz/epicompass.git
```
2. Install dependencies by running the corresponding script:
```bash
chmod +x installdependencies.sh
./installdependencies.sh
```
> This script may prompt for your password to install system packages.

Alternatively, install everything manually:
<details>
<summary>Click here to expand</summary>

```bash
sudo apt update -y && sudo apt upgrade -y
sudo apt install -y software-properties-common curl
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt update -y
sudo apt install -y python3.11 python3.11-venv python3.11-distutils
curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11
python3.11 -m venv epicompassvenv
source epicompassvenv/bin/activate
pip install --upgrade pip
pip install tqdm dask pandas pyarrow scipy seaborn matplotlib statsmodels
```
</details>

3. Download the dataset from the EpiMap repository by running the corresponding script:
```bash
./datadownload.sh
```
or else, download it manually from this link: <https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/>
4. (OPTIONAL) You can download the metadata table in .xlsx format from this link:
<https://personal.broadinstitute.org/cboix/epimap/metadata/Imputation_Metadata.xlsx>

### Repository structure
```
epicompass/
├── data/
│   ├── hg38_datasets/
│   │   ├── BSS00001_18_CALLS_segments.bed.gz
│   │   ├── BSS00002_18_CALLS_segments.bed.gz
│   │   └── [...] (+800 files)
│   ├── main_metadata_table.tsv
│   └── statemap.tsv
├── docs/
│   └── EpiCompass Documentation.pdf
├── scripts/
│   ├── classmap.py
│   ├── epimatrix.py
│   ├── hypertest.py
│   ├── installdependencies.sh
│   └── downloaddata.sh
└── README.md
```

### Execution
1. Run `classmap.py` to create a sample classification map (saved as a '.tsv' file).
2. (OPTIONAL) You may modify the `statemap.tsv` file located in the `/data` directory if you wish to group the 18 chromatin states into a reduced set of custom categories. The first column contains the names of the 18 original states, and the second column, the custom group names. By default, this file defines a classification into 7 broader states.
3. Run `epimatrix.py` to generate a chromatin state matrix (saved as a '.tsv' file) and optionally, create a state plot (saved as a '.png' file).
4. Run `hypertest.py` to perform a hypergeometric test on the previously generated chromatin state matrix, with the option to create a heatmap and save it as a '.png' file.
> **NOTE:** It is recommended to store all data files (such as `classmap.tsv`, `hg38_datasets`, etc.) inside the `/data` directory at the project `root`.

For more details, refer to the [EpiCompass Documentation](docs/EpiCompass%20Documentation.pdf).

### Example
```bash
# Create a sample classification map, for example, sorting by tissue:
python3 classmap.py tissue --output sort_by_tissue.tsv

# Execute epimatrix.py to generate the chromatin state
# matrix and the state plot:
python3 epimatrix.py hg38_datasets chr7:140000-150000,chr10:100000-150000
5000 --classmap sort_by_tissue.tsv --statemap statemap.tsv --output
cancer_matrix.tsv --plot cancer_plot.png

# Execute hypertest.py to perform a hypergeometric test on the previously
# generated matrix and save test as a TSV file and heatmap as a PNG file:
python3 hypertest.py cancer_matrix.tsv --output htest.tsv --plot
--plot-output cancer_heatmap.png --rangecap 10 --fdr
```

### License
This project is licensed under the MIT License.

### Author
**Enrique Albarrán Pérez**

### Contact
Email address: [enriqueap@uma.es](mailto:enriqueap@uma.es)