# CHIL Analyzer - Classification Algorithm & Analysis

## Overview

The CHIL Analyzer uses BLASTP (Basic Local Alignment Search Tool for Proteins) to compare user-submitted protein sequences against a reference database containing CHI and CHIL sequences. The classification is based on alignment quality metrics, primarily bit score and sequence similarity.

## Key Workflow

### Step-by-Step Process

1. **User Input Processing**
   - User inputs protein sequence on the frontend (FASTA format or plain text)
   - Backend validates and processes the sequence

2. **Temporary FASTA File Creation**
   - Backend receives sequence and creates a temporary FASTA file
   - File is saved to the `uploads/` directory with a unique identifier

3. **BLASTP Alignment** (`run_blastp()`)
   - Calls BLASTP to align user sequence against reference database
   - Command: `blastp -query <user_sequence> -db data/sequence_db -outfmt 6`
   - Output format 6 provides tabular results with key alignment metrics

4. **Results Analysis** (`analyze_blast_results()`)
   - Parses BLASTP output
   - Analyzes alignment results based on bit score and similarity thresholds
   - Determines best match and classification

5. **Position Analysis** (`analyze_user_sequence_positions()`)
   - For CHI/CHIL classifications, analyzes amino acids at key positions
   - Maps user sequence positions to reference sequence positions using alignment data
   - Identifies matches and mismatches at critical positions

6. **Results Return**
   - Returns classification result (CHI/CHIL/Unknown)
   - Provides detailed position analysis
   - Includes complete BLAST alignment information

7. **Cleanup**
   - Automatically removes temporary files
   - Archives results in `processed/` directory

## Classification Logic

### Decision Criteria

The classification algorithm uses a hierarchical decision tree based on alignment quality metrics:

#### Primary Classification Rules

1. **High-Quality Match** (CHI or CHIL)
   - **Bit score ≥ 20** AND **Similarity ≥ 25%**
   - Classified as the type (CHI or CHIL) with the highest bit score
   - Bit score is prioritized over raw similarity percentage

2. **Low-Quality Match** (Unknown)
   - **Bit score < 20**
   - Classified as "Unknown" regardless of similarity
   - Indicates insufficient homology to confidently assign a type

### Why Bit Score Over Similarity?

- **Bit score** is a normalized metric that accounts for:
  - Alignment length
  - Database size
  - Scoring matrix parameters
  - Statistical significance

- **Raw similarity percentage** can be misleading:
  - Short, high-identity alignments may not be biologically meaningful
  - Does not account for gap penalties or alignment quality

- Bit score provides a more reliable measure of true homology

### Threshold Selection

The current thresholds (bit score ≥ 20, similarity ≥ 25%) are designed to:
- Minimize false positives (incorrectly assigning CHI/CHIL to non-homologous sequences)
- Maintain sensitivity for divergent CHI/CHIL family members
- Balance specificity and recall

**Note**: These thresholds may be adjusted based on validation testing and user requirements.

## Key Position Analysis

### Purpose

Key positions represent functionally or structurally important residues that differentiate CHI from CHIL proteins. These positions may be involved in:
- Substrate binding
- Catalytic activity
- Protein stability
- Functional divergence between CHI and CHIL

### Reference Key Positions

#### CHI Key Positions
- **Positions**: 47, 59, 124, 201, 202
- **Total**: 5 positions
- These positions are characteristic of CHI protein family members

#### CHIL Key Positions
- **Positions**: 33, 35, 36, 37, 38, 40, 42, 94, 97, 100, 104, 110, 138, 184
- **Total**: 14 positions
- These positions are characteristic of CHIL protein family members

### Position Mapping Algorithm

The position analysis uses **BLAST alignment-based mapping** to establish correspondence between user sequence and reference sequence positions:

1. **Parse BLAST Alignment**
   - Extract query (user) sequence alignment
   - Extract subject (reference) sequence alignment
   - Identify alignment start positions for both sequences

2. **Build Position Map**
   - Create mapping from reference positions to user positions
   - Account for gaps in both sequences
   - Handle insertions and deletions

3. **Extract Key Position Residues**
   - For each reference key position:
     - Find corresponding user position (if aligned)
     - Extract amino acid at that position
     - Compare with reference amino acid

4. **Categorize Matches**
   - **Match**: User amino acid = Reference amino acid (green)
   - **Mismatch**: User amino acid ≠ Reference amino acid (red)
   - **Unaligned**: Position not covered by alignment (gray)

### Position Analysis Output

For each key position, the analysis provides:

```json
{
  "position": 47,
  "reference_aa": "R",
  "user_aa": "R",
  "user_position": 45,
  "match": true,
  "aligned": true
}
```

- **position**: Reference sequence position number
- **reference_aa**: Amino acid in reference sequence
- **user_aa**: Amino acid in user sequence (or "-" if gap)
- **user_position**: Corresponding position in user sequence
- **match**: Boolean indicating if amino acids match
- **aligned**: Boolean indicating if position is covered by alignment

## Visualization

### Position Display

The web interface visualizes key positions with color coding:

| Color | Meaning | Description |
|-------|---------|-------------|
| Green | Match | User amino acid matches reference |
| Red | Mismatch | User amino acid differs from reference |
| Gray | Unaligned | Position not covered by BLAST alignment |

### Interactive Features

- **Hover over position number**: View detailed position information
  - Reference sequence position
  - User sequence position
  - Alignment coverage

## BLAST Parameters

### Current Configuration

```bash
blastp \
  -query <user_sequence.fasta> \
  -db data/sequence_db \
  -outfmt 6 \
  -evalue 10 \
  -max_target_seqs 10
```

### Parameter Explanations

- **-query**: User sequence file
- **-db**: Reference database path
- **-outfmt 6**: Tabular output format
  - Columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
- **-evalue**: Expectation value threshold (default: 10)
  - Controls statistical significance
  - Lower values = more stringent
- **-max_target_seqs**: Maximum number of aligned sequences to return

### Customization Options

For specific use cases, parameters can be adjusted:

- **Stricter alignment**: Lower e-value (e.g., 1e-5)
- **More permissive**: Higher e-value (e.g., 100)
- **Longer alignment coverage**: Add `-qcov_hsp_perc` parameter
- **Different scoring matrix**: Add `-matrix` parameter (e.g., BLOSUM45, BLOSUM80)

## Reference Database

### Database Structure

```
data/
├── sequence.fasta     # FASTA file with CHI and CHIL reference sequences
└── sequence_db.*      # BLAST database files (generated by makeblastdb)
```

### Database Creation

```bash
makeblastdb -in data/sequence.fasta -dbtype prot -out data/sequence_db
```

### Database Updates

To add new reference sequences:

1. Edit `data/sequence.fasta` to add new sequences
2. Regenerate BLAST database:
   ```bash
   makeblastdb -in data/sequence.fasta -dbtype prot -out data/sequence_db
   ```
3. Restart the application

## Performance Considerations

### Computational Complexity

- **BLAST alignment**: O(mn) where m = query length, n = database size
- **Position analysis**: O(k) where k = number of key positions
- **Total time**: Typically < 1 second for single sequence

### Optimization Strategies

- Pre-built BLAST database for fast lookups
- Efficient parsing of BLAST output
- Caching of reference sequence data
- Automatic cleanup of temporary files

## Future Improvements

### Algorithm Enhancements

1. **Machine Learning Classification**
   - Train ML model on CHI/CHIL features
   - Use alignment features, key position conservation, sequence composition
   - Potentially improve classification accuracy

2. **Multiple Sequence Alignment**
   - Align user sequence with multiple CHI and CHIL references
   - Build consensus classification from multiple alignments
   - More robust to individual sequence variations

3. **Phylogenetic Analysis**
   - Construct phylogenetic tree including user sequence
   - Classify based on tree topology
   - Provide evolutionary context

4. **Position Weight Matrices**
   - Build PWMs for CHI and CHIL key positions
   - Score user sequence against both PWMs
   - Quantitative conservation scores

5. **Structural Prediction**
   - Predict 3D structure of user sequence
   - Compare with CHI/CHIL structural templates
   - Identify functional residues in structural context

6. **Configurable Thresholds**
   - Allow users to adjust bit score and similarity thresholds
   - Provide sensitivity/specificity trade-off controls
   - Custom classification rules for specific applications

### Data Enhancements

1. **Expanded Reference Database**
   - Add more CHI and CHIL sequences from diverse organisms
   - Include experimentally validated sequences
   - Curate sequence annotations

2. **Key Position Validation**
   - Experimental validation of key positions
   - Structural and functional analysis
   - Literature-based curation

3. **Cross-validation**
   - Test classification accuracy on known sequences
   - Optimize thresholds based on validation data
   - Report confidence scores
