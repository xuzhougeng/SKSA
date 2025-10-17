from flask import Flask, render_template, request, send_file, url_for, jsonify, Blueprint
import os
import zipfile
from werkzeug.utils import secure_filename
import subprocess
import tempfile
import json

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['PROCESSED_FOLDER'] = 'processed'
app.config['BLAST_DB'] = 'data/sequence_db'

# Create Blueprint for the analyzer
analyzer_bp = Blueprint('analyzer', __name__, url_prefix='/chil_analyzer')

# Key positions for visualization
CHI_POSITIONS = [47, 59,117, 124, 201, 202]
CHIL_POSITIONS = [33, 35, 36, 37, 38, 40, 42, 94, 97, 100, 104, 110, 138, 184]

# Reference sequences
REFERENCE_SEQUENCES = {
    'CHI': 'MSSSNACASPSPFPAVTKLHVDSVTFVPSVKSPASSNPLFLGGAGVRGLDIQGKFVIFTVIGVYLEGNAVPSLSVKWKGKTTEELTESIPFFREIVTGAFEKFIKVTMKLPLTGQQYSEKVTENCVAIWKQLGLYTDCEAKAVEKFLEIFKEETFPPGSSILFALSPTGSLTVAFSKDDSIPETGIAVIENKLLAEAVLESIIGKNGVSPGTRLSVAERLSQLMMKNKDEKEVSDHSVEEKLAKEN',
    'CHIL': 'MGTEMVMVHEVPFPPQIITSKPLSLLGQGITDIEIHFLQVKFTAIGVYLDPSDVKTHLNWKGKTGKELAGDDDFFDALASAEMEKVIRVVVIKEIKGAQYGVQLENTVRDRLAEEDKYEEEEETELEKVVGFFQSKYFKANSVITYHFSAKDGICEIGFETEGKEEEKLKVENANVVGMMQRWYLSGSRGVSPSTIVSIADSISAVLT'
}

def extract_key_positions(sequence, positions):
    """Extract amino acids at specified positions (1-indexed)"""
    key_residues = {}
    for pos in positions:
        if pos <= len(sequence):
            key_residues[pos] = sequence[pos - 1]  # Convert to 0-indexed
        else:
            key_residues[pos] = '-'  # Position beyond sequence length
    return key_residues

def get_reference_key_positions():
    """Get key positions for both CHI and CHIL reference sequences"""
    return {
        'CHI': {
            'positions': CHI_POSITIONS,
            'residues': extract_key_positions(REFERENCE_SEQUENCES['CHI'], CHI_POSITIONS)
        },
        'CHIL': {
            'positions': CHIL_POSITIONS,
            'residues': extract_key_positions(REFERENCE_SEQUENCES['CHIL'], CHIL_POSITIONS)
        }
    }

def run_blastp(query_file):
    """Run BLASTP against CHI/CHIL database"""
    try:
        # Run BLASTP with detailed alignment output
        result = subprocess.run([
            'blastp',
            '-query', query_file,
            '-db', app.config['BLAST_DB'],
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-max_target_seqs', '2'
        ], capture_output=True, text=True, check=True)
        
        # Also get alignment details for position mapping
        alignment_result = subprocess.run([
            'blastp',
            '-query', query_file,
            '-db', app.config['BLAST_DB'],
            '-outfmt', '0',  # pairwise alignment format
            '-max_target_seqs', '1'
        ], capture_output=True, text=True, check=True)
        
        return result.stdout.strip(), alignment_result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"BLAST error: {e}")
        return None, None

def parse_blast_alignment(alignment_text, classification):
    """Parse BLAST alignment to map query positions to reference positions"""
    if not alignment_text or classification not in ['CHI', 'CHIL']:
        return None
    
    lines = alignment_text.split('\n')
    position_map = {}
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Look for Query line
        if line.startswith('Query'):
            query_parts = line.split()
            if len(query_parts) >= 4:
                query_start = int(query_parts[1])
                query_seq = query_parts[2]
                query_end = int(query_parts[3])
                
                # Look for corresponding Sbjct line (should be 2 lines down)
                sbjct_line = None
                for j in range(i+1, min(i+4, len(lines))):
                    if lines[j].startswith('Sbjct'):
                        sbjct_line = lines[j]
                        break
                
                if sbjct_line:
                    sbjct_parts = sbjct_line.split()
                    if len(sbjct_parts) >= 4:
                        sbjct_start = int(sbjct_parts[1])
                        sbjct_seq = sbjct_parts[2]
                        sbjct_end = int(sbjct_parts[3])
                        
                        # Create position mapping for this alignment block
                        query_pos = query_start - 1  # Convert to 0-indexed
                        sbjct_pos = sbjct_start - 1  # Convert to 0-indexed
                        
                        for k in range(len(query_seq)):
                            if k < len(sbjct_seq):  # Ensure we don't go beyond sbjct_seq
                                query_char = query_seq[k]
                                sbjct_char = sbjct_seq[k]
                                
                                if query_char != '-':
                                    query_pos += 1
                                if sbjct_char != '-':
                                    sbjct_pos += 1
                                
                                # Map query position to subject position
                                if query_char != '-' and sbjct_char != '-':
                                    position_map[query_pos - 1] = sbjct_pos - 1  # Store as 0-indexed
        i += 1
    
    return position_map

def analyze_user_sequence_positions_aligned(user_sequence, classification, blast_alignment):
    """Analyze user sequence positions based on BLAST alignment"""
    if classification not in ['CHI', 'CHIL']:
        return None
    
    # Clean user sequence
    if '>' in user_sequence:
        lines = user_sequence.strip().split('\n')
        sequence_lines = [line for line in lines if not line.startswith('>')]
        clean_sequence = ''.join(sequence_lines)
    else:
        clean_sequence = user_sequence.replace('\n', '').replace(' ', '')
    
    # Get position mapping from BLAST alignment
    position_map = parse_blast_alignment(blast_alignment, classification)
    if not position_map:
        return None
    
    # Get key positions for the classification type
    key_positions = CHI_POSITIONS if classification == 'CHI' else CHIL_POSITIONS
    reference_sequence = REFERENCE_SEQUENCES[classification]
    
    comparison = {}
    
    for ref_pos in key_positions:
        ref_pos_0 = ref_pos - 1  # Convert to 0-indexed
        
        # Find corresponding query position
        query_pos_0 = None
        for q_pos, r_pos in position_map.items():
            if r_pos == ref_pos_0:
                query_pos_0 = q_pos
                break
        
        if query_pos_0 is not None and query_pos_0 < len(clean_sequence):
            user_aa = clean_sequence[query_pos_0]
            ref_aa = reference_sequence[ref_pos_0] if ref_pos_0 < len(reference_sequence) else '-'
            
            comparison[ref_pos] = {
                'user': user_aa,
                'reference': ref_aa,
                'match': user_aa == ref_aa,
                'aligned': True,
                'query_pos': query_pos_0 + 1  # Store 1-indexed for display
            }
        else:
            # Position not aligned or beyond sequence
            ref_aa = reference_sequence[ref_pos_0] if ref_pos_0 < len(reference_sequence) else '-'
            comparison[ref_pos] = {
                'user': '-',
                'reference': ref_aa,
                'match': False,
                'aligned': False,
                'query_pos': None
            }
    
    return {
        'type': classification,
        'positions': key_positions,
        'comparison': comparison,
        'alignment_info': {
            'total_mapped_positions': len(position_map),
            'sequence_length': len(clean_sequence)
        }
    }

def analyze_blast_results(blast_output):
    """Analyze BLAST results to determine CHI/CHIL/Unknown"""
    if not blast_output:
        return "Unknown", 0
    
    lines = blast_output.strip().split('\n')
    if not lines or lines == ['']:
        return "Unknown", 0
    
    best_hit = None
    best_bitscore = 0
    best_identity = 0
    
    for line in lines:
        fields = line.split('\t')
        if len(fields) >= 12:
            subject_id = fields[1]
            identity = float(fields[2])
            bitscore = float(fields[11])
            
            # Use bitscore as primary criterion, with minimum identity requirement
            if bitscore > best_bitscore and identity >= 25:  # Basic identity threshold
                best_bitscore = bitscore
                best_hit = subject_id
                best_identity = identity
    
    # Minimum bitscore threshold for classification
    if best_bitscore < 20:  # Adjust based on test results
        return "Unknown", best_identity
    
    classification = best_hit if best_hit in ['CHI', 'CHIL'] else "Unknown"
    return classification, best_identity

def process_fasta(file_path):
    """Process uploaded FASTA file with BLAST analysis"""
    try:
        # Read user sequence for position analysis
        with open(file_path, 'r') as f:
            user_sequence = f.read()
        
        # Run BLAST analysis
        blast_output, blast_alignment = run_blastp(file_path)
        classification, similarity_percentage = analyze_blast_results(blast_output)
        
        # Analyze key positions using alignment
        position_analysis = analyze_user_sequence_positions_aligned(user_sequence, classification, blast_alignment)
        
        # Additional evaluation criterion for CHIL: if more than 1/3 positions are missing, classify as Unknown
        if classification == 'CHIL' and position_analysis and position_analysis.get('comparison'):
            missing_count = 0
            total_positions = len(CHIL_POSITIONS)
            
            for pos in CHIL_POSITIONS:
                if pos in position_analysis['comparison']:
                    if position_analysis['comparison'][pos]['user'] == '-':
                        missing_count += 1
            
            # If more than 1/3 of positions are missing, reclassify as Unknown
            if missing_count > total_positions / 3:
                classification = "Unknown"
                # Update position analysis to reflect the new classification
                position_analysis = None
        
        # Get reference key positions
        reference_positions = get_reference_key_positions()
        
        # Create result file
        result_data = {
            'classification': classification,
            'similarity_percentage': similarity_percentage,
            'blast_output': blast_output,
            'blast_alignment': blast_alignment,
            'position_analysis': position_analysis,
            'reference_positions': reference_positions
        }
        
        result_filename = f"result_{os.path.basename(file_path)}.json"
        result_path = os.path.join(app.config['PROCESSED_FOLDER'], result_filename)
        
        with open(result_path, 'w') as f:
            json.dump(result_data, f, indent=2)
        
        # Create zip file with result
        zip_filename = f"analysis_{os.path.basename(file_path)}.zip"
        zip_path = os.path.join(app.config['PROCESSED_FOLDER'], zip_filename)
        
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            zipf.write(result_path, os.path.basename(result_path))
            zipf.write(file_path, os.path.basename(file_path))
        
        return zip_path, classification, position_analysis, similarity_percentage
        
    except Exception as e:
        print(f"Processing error: {e}")
        return None, "Unknown", None, 0

# Routes using Blueprint
@analyzer_bp.route('/')
def index():
    return render_template('index.html')

@analyzer_bp.route('/upload', methods=['POST'])
def upload_file():
    try:
        # Check if sequence is provided
        if 'sequence' not in request.form or not request.form['sequence'].strip():
            return jsonify({'error': 'No sequence provided'}), 400
            
        sequence = request.form['sequence'].strip()
        
        # Create temporary file for sequence
        import uuid
        temp_filename = f"temp_{uuid.uuid4().hex[:8]}.fasta"
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], temp_filename)
        
        # Format sequence if not in FASTA format
        if not sequence.startswith('>'):
            sequence = f">user_sequence\n{sequence}"
        
        with open(file_path, 'w') as f:
            f.write(sequence)
        
        # Process the sequence
        zip_path, classification, position_analysis, similarity_percentage = process_fasta(file_path)
        
        # Clean up temporary file
        if os.path.exists(file_path):
            os.remove(file_path)
        
        if zip_path:
            download_url = url_for('analyzer.download_file', filename=os.path.basename(zip_path))
            return jsonify({
                'download_url': download_url,
                'classification': classification,
                'position_analysis': position_analysis,
                'similarity_percentage': similarity_percentage
            })
        else:
            return jsonify({'error': 'Processing failed'}), 500
            
    except Exception as e:
        print(f"Processing error: {e}")
        return jsonify({'error': 'Processing failed'}), 500

@analyzer_bp.route('/download/<filename>')
def download_file(filename):
    return send_file(os.path.join(app.config['PROCESSED_FOLDER'], filename), as_attachment=True)

# Optional: Redirect root to the analyzer
@app.route('/')
def root():
    return f'<h1>CHIL Analyzer</h1><p>Access the application at <a href="/chil_analyzer/">/chil_analyzer/</a></p>'

# Register the blueprint AFTER all routes are defined
app.register_blueprint(analyzer_bp)

if __name__ == '__main__':
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    os.makedirs(app.config['PROCESSED_FOLDER'], exist_ok=True)
    app.run(debug=True, port=50000)
