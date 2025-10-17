#!/usr/bin/env python3
"""
Test script for CHIL Analyzer
Tests three protein sequences:
- CHI_pep.fasta (expected: CHI)
- CHIL_pep.fasta (expected: CHIL) 
- CHS_pep.fasta (expected: Unknown)
"""

import requests
import os
import sys

def read_first_sequence(filepath):
    """Read the first sequence from a FASTA file"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    sequence = ""
    header_found = False
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header_found:  # Start of second sequence, stop here
                break
            header_found = True
            sequence += line + '\n'
        elif header_found and line:
            sequence += line + '\n'
    
    return sequence.strip()

def test_sequence(sequence, expected_result, test_name):
    """Test a sequence against the analyzer"""
    url = "http://127.0.0.1:50000/chil_analyzer/upload"
    
    try:
        response = requests.post(url, data={'sequence': sequence})
        
        if response.status_code == 200:
            result = response.json()
            classification = result.get('classification', 'Error')
            
            status = "✓ PASS" if classification == expected_result else "✗ FAIL"
            print(f"{test_name:15} | Expected: {expected_result:7} | Got: {classification:7} | {status}")
            
            return classification == expected_result
        else:
            print(f"{test_name:15} | ERROR: HTTP {response.status_code}")
            return False
            
    except requests.exceptions.ConnectionError:
        print(f"{test_name:15} | ERROR: Cannot connect to server")
        return False
    except Exception as e:
        print(f"{test_name:15} | ERROR: {str(e)}")
        return False

def main():
    print("CHIL Analyzer Test Suite")
    print("=" * 50)
    
    # Check if server is running
    try:
        response = requests.get("http://127.0.0.1:50000/chil_analyzer/")
        if response.status_code != 200:
            print("ERROR: Server not responding correctly")
            sys.exit(1)
    except requests.exceptions.ConnectionError:
        print("ERROR: Server not running. Please start with 'python app.py'")
        sys.exit(1)
    
    test_files = [
        ('tests/CHI_pep.fasta', 'CHI', 'CHI Test'),
        ('tests/CHIL_pep.fasta', 'CHIL', 'CHIL Test'),
        ('tests/CHS_pep.fasta', 'Unknown', 'CHS Test')
    ]
    
    results = []
    
    print(f"{'Test Name':15} | {'Expected':9} | {'Result':9} | Status")
    print("-" * 50)
    
    for filepath, expected, name in test_files:
        if not os.path.exists(filepath):
            print(f"{name:15} | ERROR: File {filepath} not found")
            results.append(False)
            continue
            
        sequence = read_first_sequence(filepath)
        if not sequence:
            print(f"{name:15} | ERROR: No sequence found in file")
            results.append(False)
            continue
            
        result = test_sequence(sequence, expected, name)
        results.append(result)
    
    print("-" * 50)
    passed = sum(results)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests PASSED!")
        sys.exit(0)
    else:
        print("❌ Some tests FAILED!")
        sys.exit(1)

if __name__ == "__main__":
    main() 