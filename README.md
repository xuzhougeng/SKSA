# SKSA pipeline

SKSA, Structure-based Key Amino acid residues clustering plus Sequence Alignment

## Quick Start

### Online Demo: 

Try the live web application [http://wanglab.sippe.ac.cn/chil_analyzer/](http://wanglab.sippe.ac.cn/chil_analyzer/) without installation.

Simply paste your protein sequence and get instant classification results with key position analysis.

### Command Line Usage:

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

## Installation

Prerequisites:

- Python 3.6+
- NCBI BLAST+ toolkit


Initialize Database:

```bash
makeblastdb -in data/sequence.fasta -dbtype prot -out data/sequence_db
```

## Documentation

Detailed documentation is available in the `doc/` directory:

- **[doc/STRUCTURE.md](doc/STRUCTURE.md)** - Project overview, installation, and quick start guide
- **[doc/ALGORITHM.md](doc/ALGORITHM.md)** - Classification algorithm, BLAST parameters, and position analysis
- **[doc/WEB_INTERFACE.md](doc/WEB_INTERFACE.md)** - Web interface usage, API documentation, and deployment
