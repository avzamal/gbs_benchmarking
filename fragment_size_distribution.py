#!/usr/bin/env python
'''
Calculate fragment size distribution from a SAM/BAM file and output a table with two columns:
fragment size and count.

Usage: 
    fragment_size_distribution.py <SAM file>
    samtools view <BAM file> | fragment_size_distribution.py -

Output:
    The script outputs a tab-delimited table to standard output with two columns:
    1. fragment size
    2. count (number of read pairs with this fragment size)
'''

import sys
import argparse
import re
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Calculate fragment size distribution from a SAM/BAM file.')
    parser.add_argument('samfile', nargs='?', type=argparse.FileType('r'), 
                        default=sys.stdin, help='Input SAM file (use - for stdin)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file (default: stdout)')
    parser.add_argument('-m', '--min-size', type=int, default=0,
                        help='Minimum fragment size to include (default: 0)')
    parser.add_argument('-M', '--max-size', type=int, default=float('inf'),
                        help='Maximum fragment size to include (default: no limit)')
    parser.add_argument('-s', '--summary', action='store_true',
                        help='Also output summary statistics (mean, median, mode, std)')
    args = parser.parse_args()

    # Store fragment sizes
    fragment_sizes = defaultdict(int)
    total_fragments = 0
    
    # Process the SAM file
    line_count = 0
    for line in args.samfile:
        line_count += 1
        if line_count % 1000000 == 0:
            print(f"Processed {line_count // 1000000}M lines...", file=sys.stderr)
            
        # Skip header lines
        if line.startswith('@'):
            continue
            
        fields = line.strip().split('\t')
        
        # Check if we have enough fields and if the read is properly paired
        if len(fields) < 11:
            continue
            
        flag = int(fields[1])
        
        # Check if read is paired and is the first in pair (to avoid counting twice)
        if not (flag & 1) or not (flag & 64):  # not paired or not first in pair
            continue
            
        # Check if read is properly paired and mapped to same chromosome
        if not (flag & 2) or fields[6] != '=':
            continue
            
        # Get the fragment size (absolute value of TLEN field)
        fragment_size = abs(int(fields[8]))
        
        # Skip if fragment size is outside the specified range
        if fragment_size < args.min_size or fragment_size > args.max_size:
            continue
            
        # Increment count for this fragment size
        fragment_sizes[fragment_size] += 1
        total_fragments += 1
    
    # Output the fragment size distribution
    print("Fragment_Size\tCount", file=args.output)
    for size in sorted(fragment_sizes.keys()):
        print(f"{size}\t{fragment_sizes[size]}", file=args.output)
    
    # Calculate and output summary statistics if requested
    if args.summary and total_fragments > 0:
        # Calculate mean
        mean = sum(size * count for size, count in fragment_sizes.items()) / total_fragments
        
        # Calculate mode (most common fragment size)
        mode = max(fragment_sizes.items(), key=lambda x: x[1])[0]
        
        # Calculate median
        sizes_sorted = []
        for size, count in fragment_sizes.items():
            sizes_sorted.extend([size] * count)
        sizes_sorted.sort()
        if total_fragments % 2 == 0:
            median = (sizes_sorted[total_fragments//2-1] + sizes_sorted[total_fragments//2]) / 2
        else:
            median = sizes_sorted[total_fragments//2]
        
        # Calculate standard deviation
        variance = sum(((size - mean) ** 2) * count for size, count in fragment_sizes.items()) / total_fragments
        std_dev = variance ** 0.5
        
        print("\nSummary Statistics:", file=args.output)
        print(f"Total fragments: {total_fragments}", file=args.output)
        print(f"Mean fragment size: {mean:.2f}", file=args.output)
        print(f"Median fragment size: {median:.2f}", file=args.output)
        print(f"Mode fragment size: {mode}", file=args.output)
        print(f"Standard deviation: {std_dev:.2f}", file=args.output)

if __name__ == "__main__":
    main()
