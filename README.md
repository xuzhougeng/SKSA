# CHIL Analyzer

A protein sequence classification tool that uses BLASTP to distinguish CHI (Chitinase) and CHIL (Chitinase-Like) proteins from other sequences.

## Quick Start

### Online Demo

Try the live web application without installation:

**🌐 [http://wanglab.sippe.ac.cn/chil_analyzer/](http://wanglab.sippe.ac.cn/chil_analyzer/)**

Simply paste your protein sequence and get instant classification results with key position analysis.

### Web Interface (Local)
For running your own instance with interactive analysis and visualization:
```bash
python app.py
# Visit http://127.0.0.1:50000/chil_analyzer/
```

### Command Line Usage (No Web Interface Required)

You can use CHIL Analyzer directly from the command line without running the web server:

#### Single Sequence Analysis

```bash
# Create a test sequence file
cat > my_sequence.fasta << 'EOF'
>test_seq
MSSSNACASPSPFPAVTKLHVDSVTFVPSVKSPASSNPLFLGGAGVRGLDIQGKFV
EOF

# Run BLASTP directly
blastp -query my_sequence.fasta -db data/sequence_db -outfmt 6 -out results.txt

# View results
cat results.txt
```

#### Batch Processing Multiple Sequences

```bash
# Process all FASTA files in a directory
for file in *.fasta; do
    echo "Processing $file..."
    blastp -query "$file" -db data/sequence_db -outfmt 6 -out "${file%.fasta}_results.txt"
done
```

#### Using Python Script for Batch Analysis

Create a Python script `batch_analyze.py`:

```python
#!/usr/bin/env python3
import subprocess
import os
import sys

def analyze_sequence(fasta_file, output_dir='batch_results'):
    """Analyze a single FASTA file using BLASTP"""
    os.makedirs(output_dir, exist_ok=True)

    basename = os.path.basename(fasta_file).replace('.fasta', '')
    output_file = os.path.join(output_dir, f"{basename}_blast.txt")

    # Run BLASTP
    cmd = [
        'blastp',
        '-query', fasta_file,
        '-db', 'data/sequence_db',
        '-outfmt', '6',
        '-out', output_file
    ]

    subprocess.run(cmd, check=True)
    print(f"Results saved to: {output_file}")

    # Parse and classify
    classify_results(output_file)

def classify_results(blast_output):
    """Parse BLAST results and classify sequences"""
    with open(blast_output, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                query_id = fields[0]
                subject_id = fields[1]
                identity = float(fields[2])
                bitscore = float(fields[11])

                # Classification logic
                if bitscore >= 20 and identity >= 25:
                    if 'CHI' in subject_id and 'CHIL' not in subject_id:
                        classification = 'CHI'
                    elif 'CHIL' in subject_id:
                        classification = 'CHIL'
                    else:
                        classification = 'Unknown'
                else:
                    classification = 'Unknown'

                print(f"{query_id}: {classification} (bitscore={bitscore:.1f}, identity={identity:.1f}%)")
                break  # Use best match only

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python batch_analyze.py <fasta_file> [output_dir]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'batch_results'

    analyze_sequence(fasta_file, output_dir)
```

Usage:
```bash
# Analyze single file
python batch_analyze.py my_sequence.fasta

# Analyze with custom output directory
python batch_analyze.py my_sequence.fasta custom_output/

# Analyze multiple files
for file in tests/*.fasta; do
    python batch_analyze.py "$file"
done
```

#### Advanced: Custom BLAST Parameters

```bash
# More stringent alignment
blastp -query my_sequence.fasta \
       -db data/sequence_db \
       -outfmt 6 \
       -evalue 1e-5 \
       -qcov_hsp_perc 50

# Different output format (human-readable)
blastp -query my_sequence.fasta \
       -db data/sequence_db \
       -outfmt 7

# Output with alignment details
blastp -query my_sequence.fasta \
       -db data/sequence_db \
       -outfmt "6 qseqid sseqid pident length bitscore evalue qseq sseq"
```

## Installation

### Prerequisites
- Python 3.6+
- NCBI BLAST+ toolkit

### Install Dependencies
```bash
# Python packages
pip install flask werkzeug requests

# BLAST+ (Ubuntu/Debian)
sudo apt-get install ncbi-blast+

# BLAST+ (macOS)
brew install blast
```

### Initialize Database
```bash
makeblastdb -in data/sequence.fasta -dbtype prot -out data/sequence_db
```

## Classification Logic

- **CHI/CHIL**: Bit score ≥ 20 AND Similarity ≥ 25%
- **Unknown**: Bit score < 20 OR Similarity < 25%

The algorithm prioritizes bit score as it's more reliable than raw similarity percentage.

## Documentation

Detailed documentation is available in the `doc/` directory:

- **[doc/STRUCTURE.md](doc/STRUCTURE.md)** - Project overview, installation, and quick start guide
- **[doc/ALGORITHM.md](doc/ALGORITHM.md)** - Classification algorithm, BLAST parameters, and position analysis
- **[doc/WEB_INTERFACE.md](doc/WEB_INTERFACE.md)** - Web interface usage, API documentation, and deployment

## Testing

Run automated tests:
```bash
python test_analyzer.py
```

Test individual sequences:
```bash
blastp -query tests/CHI_pep.fasta -db data/sequence_db -outfmt 6
blastp -query tests/CHIL_pep.fasta -db data/sequence_db -outfmt 6
blastp -query tests/CHS_pep.fasta -db data/sequence_db -outfmt 6
```

## Directory Structure

```
chil_analyzer/
├── README.md              # This file
├── app.py                 # Flask web application
├── test_analyzer.py       # Automated tests
│
├── doc/                   # Documentation
│   ├── STRUCTURE.md       # Project overview
│   ├── ALGORITHM.md       # Algorithm details
│   └── WEB_INTERFACE.md   # Web interface guide
│
├── data/                  # BLAST database
│   ├── sequence.fasta     # Reference sequences
│   └── sequence_db.*      # BLAST database files
│
└── tests/                 # Test sequences
    ├── CHI_pep.fasta
    ├── CHIL_pep.fasta
    └── CHS_pep.fasta
```

## Examples

### Example 1: Quick Classification
```bash
# Create test sequence
echo ">test" > test.fasta
echo "MSSSNACASPSPFPAVTKLHVDSVTFVPSVKSPASSNPLFLGGAGVRGLDIQGKFV" >> test.fasta

# Run BLAST
blastp -query test.fasta -db data/sequence_db -outfmt 6

# Output will show: query_id, subject_id, %identity, length, ..., bitscore
# Classify based on bitscore ≥ 20 and identity ≥ 25
```

### Example 2: Batch Processing with Detailed Output
```bash
#!/bin/bash
# batch_process.sh

mkdir -p batch_results

for fasta in input_sequences/*.fasta; do
    name=$(basename "$fasta" .fasta)

    # Run BLAST with detailed output
    blastp -query "$fasta" \
           -db data/sequence_db \
           -outfmt "7 qseqid sseqid pident bitscore" \
           -out "batch_results/${name}_results.txt"

    # Parse and classify
    awk '/^[^#]/ {
        if ($4 >= 20 && $3 >= 25) {
            if ($2 ~ /CHIL/) print $1 ": CHIL"
            else if ($2 ~ /CHI/) print $1 ": CHI"
            else print $1 ": Unknown"
        } else {
            print $1 ": Unknown"
        }
        exit
    }' "batch_results/${name}_results.txt"
done
```

## Troubleshooting

**BLAST not found:**
```bash
# Check BLAST installation
blastp -version

# Add to PATH if needed
export PATH=/path/to/blast/bin:$PATH
```

**Database not found:**
```bash
# Verify database files exist
ls data/sequence_db.*

# Recreate if needed
makeblastdb -in data/sequence.fasta -dbtype prot -out data/sequence_db
```

## License

[Specify your license]

## Contact

[Specify contact information]
