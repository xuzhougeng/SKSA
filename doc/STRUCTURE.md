# CHIL Analyzer - Project Documentation

## Project Overview

CHIL Analyzer is a Flask-based web application for protein sequence classification. It uses BLASTP to compare user-submitted protein sequences against a reference database containing CHI (Chitinase) and CHIL (Chitinase-Like) sequences, automatically classifying input sequences as CHI, CHIL, or Unknown.

### Key Features

- **Direct Sequence Input**: Supports FASTA format or plain sequence text pasting
- **BLAST Alignment**: Uses BLASTP to compare against CHI/CHIL reference sequences
- **Intelligent Classification**: Automatically determines CHI/CHIL/Unknown based on bit score and similarity
- **Key Position Analysis**: Visualizes differences between user sequence and reference sequences at critical positions
- **Results Display**: Real-time display of classification results with different colors for different types
- **Detailed Reports**: Provides downloadable detailed alignment results in JSON format
- **Automated Testing**: Built-in test suite to verify classification accuracy

### Technology Stack

- **Backend**: Python 3.x, Flask
- **Sequence Alignment**: NCBI BLAST+ (BLASTP)
- **Frontend**: HTML5, CSS3, JavaScript
- **Web Server**: Nginx (reverse proxy)

### Online Demo

**Live Application**: [http://wanglab.sippe.ac.cn/chil_analyzer/](http://wanglab.sippe.ac.cn/chil_analyzer/)

Try the CHIL Analyzer online without installation! The live demo provides full functionality for protein sequence classification and analysis.

## Documentation Structure

This project documentation is split into multiple files for clarity:

- **STRUCTURE.md** (this file): Project overview, directory structure, installation, and quick start
- **[ALGORITHM.md](ALGORITHM.md)**: Classification algorithm, BLAST parameters, and position analysis logic
- **[WEB_INTERFACE.md](WEB_INTERFACE.md)**: Web interface usage, API documentation, and deployment configuration

## Directory Structure

```
chil_analyzer/
├── app.py                 # Core Flask application and backend logic
├── test_analyzer.py       # Automated testing script
├── nginx.conf.example     # Nginx reverse proxy configuration example
│
├── doc/                   # Documentation
│   ├── STRUCTURE.md       # This file - project overview
│   ├── ALGORITHM.md       # Algorithm and classification logic
│   └── WEB_INTERFACE.md   # Web interface and deployment
│
├── templates/             # Flask templates
│   └── index.html         # Main web interface
│
├── static/                # Static assets
│   └── styles.css         # CSS file (deprecated, now embedded in HTML)
│
├── data/                  # Reference data and BLAST database
│   ├── sequence.fasta     # CHI/CHIL reference sequences
│   └── sequence_db.*      # BLAST database files (generated)
│
├── tests/                 # Test sequences
│   ├── CHI_pep.fasta      # CHI test sequences (should classify as CHI)
│   ├── CHIL_pep.fasta     # CHIL test sequences (should classify as CHIL)
│   └── CHS_pep.fasta      # CHS test sequences (should classify as Unknown)
│
├── uploads/               # Temporary file storage (auto-generated)
└── processed/             # Analysis results storage (auto-generated)
```

## Installation

### Prerequisites

- Python 3.6 or higher
- NCBI BLAST+ toolkit
- pip (Python package installer)

### Step 1: Install Python Dependencies

```bash
pip install flask werkzeug requests
```

### Step 2: Install BLAST+ Tools

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install ncbi-blast+
```

#### macOS
```bash
brew install blast
```

#### Other Systems
Download and install from [NCBI BLAST+ downloads](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Step 3: Initialize BLAST Database

Before first run, create the BLAST database from reference sequences:

```bash
makeblastdb -in data/sequence.fasta -dbtype prot -out data/sequence_db
```

This command will generate several database files:
- `sequence_db.phr`
- `sequence_db.pin`
- `sequence_db.psq`

### Step 4: Create Required Directories

The application will automatically create these directories if they don't exist:

```bash
mkdir -p uploads processed
```

## Quick Start

### Running Locally

1. Navigate to the project directory:
   ```bash
   cd /path/to/chil_analyzer
   ```

2. Run the Flask application:
   ```bash
   python app.py
   ```

3. Open your browser and visit:
   ```
   http://127.0.0.1:50000/chil_analyzer/
   ```

### Basic Usage Example

1. Paste a protein sequence in FASTA format or plain text:
   ```
   >test_sequence
   MSSSNACASPSPFPAVTKLHVDSVTFVPSVKSPASSNPLFLGGAGVRGLDIQGKFV
   ```

2. Click "Analyze Sequence"

3. View the classification result and key position analysis

4. Download the detailed report (JSON + ZIP)

## Testing

### Running Automated Tests

The project includes automated tests to verify classification accuracy:

```bash
python test_analyzer.py
```

### Test Coverage

The test suite includes:

1. **CHI Classification Test**
   - Input: `tests/CHI_pep.fasta`
   - Expected: Classification as CHI

2. **CHIL Classification Test**
   - Input: `tests/CHIL_pep.fasta`
   - Expected: Classification as CHIL

3. **Unknown Classification Test**
   - Input: `tests/CHS_pep.fasta` (Chalcone Synthase)
   - Expected: Classification as Unknown

### Test Output

Successful test output should show:
```
Testing CHI sequence... ✓ Passed
Testing CHIL sequence... ✓ Passed
Testing Unknown sequence (CHS)... ✓ Passed
All tests passed!
```

## Configuration

### Application Configuration

Edit `app.py` to modify application settings:

```python
# Flask configuration
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['PROCESSED_FOLDER'] = 'processed'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB max file size

# BLAST configuration
BLAST_DB = 'data/sequence_db'
BLAST_EVALUE = 10
```

### Classification Thresholds

Modify classification thresholds in `app.py`:

```python
MIN_BITSCORE = 20      # Minimum bit score for CHI/CHIL classification
MIN_SIMILARITY = 25    # Minimum similarity percentage
```

## Troubleshooting

### Common Issues

1. **BLAST command not found**
   - Ensure BLAST+ is installed and in your PATH
   - Test with: `blastp -version`

2. **Database files not found**
   - Run `makeblastdb` command to create database files
   - Verify `data/sequence_db.*` files exist

3. **Permission denied on uploads/ or processed/**
   - Ensure directories exist and have write permissions
   - Run: `chmod 755 uploads processed`

4. **Port 50000 already in use**
   - Change the port in `app.py`: `app.run(port=<new_port>)`
   - Or kill the process using port 50000

## Further Documentation

For detailed information on specific topics, see:

- **Classification Algorithm**: [ALGORITHM.md](ALGORITHM.md)
  - BLAST parameters and workflow
  - Classification logic and thresholds
  - Key position analysis methodology
  - Future algorithm improvements

- **Web Interface & Deployment**: [WEB_INTERFACE.md](WEB_INTERFACE.md)
  - API endpoint documentation
  - Nginx deployment configuration
  - Frontend interface details
  - Production deployment guide

## Contributing

When contributing to this project:

1. Follow the existing code style
2. Add tests for new features
3. Update documentation accordingly
4. Test both CHI and CHIL classification accuracy

## License

[Specify your license here]

## Contact

[Specify contact information or repository URL here]