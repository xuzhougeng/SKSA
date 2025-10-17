# CHIL Analyzer - Web Interface & Deployment

## Online Demo

**Live Application**: [http://wanglab.sippe.ac.cn/chil_analyzer/](http://wanglab.sippe.ac.cn/chil_analyzer/)

Try the CHIL Analyzer without any installation! The online demo provides full functionality including:
- Instant protein sequence classification
- Interactive key position analysis visualization
- Downloadable detailed reports in JSON format
- Support for both FASTA format and plain text sequences

## Usage

### How to Use
1. Paste protein sequence in the text box (supports FASTA format or plain sequence)
2. Click "Analyze Sequence" to analyze
3. View classification results and key position analysis
4. Download detailed report

### Example Sequence Formats
```
>my_sequence
MSSSNACASPSPFPAVTKLHVDSVTFVPSVKSPASSNPLFLGGAGVRGLDIQGKFV

Or direct sequence input:
MSSSNACASPSPFPAVTKLHVDSVTFVPSVKSPASSNPLFLGGAGVRGLDIQGKFV
```

## Frontend Interface

### Results Display
- Real-time display of classification results
- Different colors for different classification types:
  - CHI: Specific color coding
  - CHIL: Specific color coding
  - Unknown: Specific color coding

### Position Analysis Visualization
- When a sequence is classified as CHI or CHIL, automatically displays key position analysis
- Visual amino acid comparison between user and reference sequences:
  - **Green background**: Match
  - **Red background**: Mismatch
  - **Gray background**: Unaligned
- Displays reference position number, user amino acid (top) and reference amino acid (bottom)
- Hover over position number to view detailed position information (reference position and corresponding user sequence position)

## API Interface

### POST /chil_analyzer/upload

**Endpoint**: `/chil_analyzer/upload`

**Method**: POST

**Request Format**:
- `sequence`: Directly submitted sequence text (required)

**Response Format**:
```json
{
  "classification": "CHI|CHIL|Unknown",
  "download_url": "/download/filename.zip",
  "position_analysis": {
    "type": "CHI|CHIL",
    "positions": [47, 59, 124, 201, 202],
    "comparison": {
      "47": {
        "user": "R",
        "reference": "R",
        "match": true
      },
      "59": {
        "user": "A",
        "reference": "T",
        "match": false
      }
    }
  }
}
```

### GET /chil_analyzer/download/<filename>

**Endpoint**: `/chil_analyzer/download/<filename>`

**Method**: GET

**Description**: Download analysis results as a compressed package

**Returns**: ZIP file containing:
- `result_<filename>.json`: Detailed classification results, BLAST output, and position analysis
- Temporarily generated FASTA file

## Deployment Configuration

### Local Development

```bash
python app.py
# Visit http://127.0.0.1:50000/chil_analyzer/
```

### Subpath Access

The application is configured to run under `/chil_analyzer` subpath for convenient nginx reverse proxy:

- **Application homepage**: `http://your-domain.com/chil_analyzer/`
- **API endpoint**: `http://your-domain.com/chil_analyzer/upload`
- **Download link**: `http://your-domain.com/chil_analyzer/download/<filename>`

### Nginx Reverse Proxy Setup

#### Step 1: Copy Configuration File
```bash
sudo cp nginx.conf.example /etc/nginx/sites-available/chil_analyzer
```

#### Step 2: Edit Configuration
```bash
sudo nano /etc/nginx/sites-available/chil_analyzer
```

Replace the domain name and adjust paths as needed.

#### Step 3: Enable Site
```bash
sudo ln -s /etc/nginx/sites-available/chil_analyzer /etc/nginx/sites-enabled/
```

#### Step 4: Test Configuration
```bash
sudo nginx -t
```

#### Step 5: Reload Nginx
```bash
sudo systemctl reload nginx
```

### Production Considerations

- Ensure BLAST+ is installed on the production server
- Configure proper file permissions for `uploads/` and `processed/` directories
- Set up automatic cleanup for temporary files
- Consider using a process manager like systemd or supervisor for the Flask application
- Configure proper logging for monitoring and debugging

## File Structure

### Frontend Files
```
templates/
└── index.html      # Frontend page (with embedded styles and classification results display)

static/
└── styles.css      # Style file (deprecated, CSS now embedded in HTML)
```

### Upload and Processing Directories
```
uploads/            # Temporary file storage directory
processed/          # Analysis results storage directory
```

## Output Files

After analysis completion, the system provides a compressed package containing:

- **`result_<filename>.json`**: Detailed classification results, BLAST output, and position analysis
- **Temporarily generated FASTA file**: User sequence in FASTA format

### Example Output JSON Structure

```json
{
  "classification": "CHI",
  "best_match": {
    "subject_id": "CHI_reference_001",
    "bitscore": 245.5,
    "identity_percent": 85.3,
    "alignment_length": 210
  },
  "blast_output": "...",
  "position_analysis": {
    "type": "CHI",
    "positions": [47, 59, 124, 201, 202],
    "comparison": {
      "47": {"user": "R", "reference": "R", "match": true},
      "59": {"user": "T", "reference": "T", "match": true},
      "124": {"user": "A", "reference": "S", "match": false},
      "201": {"user": "L", "reference": "L", "match": true},
      "202": {"user": "K", "reference": "K", "match": true}
    }
  }
}
```

## Security Considerations

- Input validation for sequence format
- File size limits for uploads
- Sanitization of user-provided sequence data
- Automatic cleanup of temporary files
- Rate limiting for API endpoints (recommended for production)
