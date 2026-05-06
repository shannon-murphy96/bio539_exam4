#!/usr/bin/env python3

# kmer_analyzer.py
# Analyzes k-mer frequencies and next characters in the DNA sequence 
# Usage: python kmer_analyzer.py [sequence file] [k] [output file]
# Example: python kmer_analyzer.py example_sequences.txt 2 results.txt
# Output: a text file with each kmer, its total count, and next character frequencies

import sys
 
#Starting by establishing kmer size and vaidating for only nucleotides and the correct "k" length 
def validate_sequence(sequence, k):
    # nucleotide sequence must be longer than "k"" so there is at least one kmer with a next character. Otherwise it rejects it
    if len(sequence) <= k: #for example, "2" would mean 2 nucleotides
        return False
    # check every character is a valid nucletoide
    valid_bases = {'A', 'C', 'G', 'T'} #prevents any misstypes outsdie of the nucleotides I want, original script only accounted for numbers
    for base in sequence:
        if base not in valid_bases:
            return False #Stops there
    return True #keeps on moving along and will analyze the kmer
 

#Count the kmers into kmer_data dictionary. Fixing bug where there was an initial double count if we did not start at 0
def update_kmer_count(kmer_data, kmer, next_char):
    # if we haven't seen this kmer yet, add it with count 0
    if kmer not in kmer_data: 
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
 
    # if we have seeen this kmer before then count it
    kmer_data[kmer]['count'] += 1
 
    # follow up nucleotide, add it to the dictionary if we havent seen it or add 1 if we have, aka tracking how often we see this follow up
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1
 
    return kmer_data
 
 
#comb through sequence checking each kmer and tallying in "update_kmer_count" 
def count_kmers_with_context(sequence, k):
    kmer_data = {}
 
    # slide through sequence - stop k positions before the end so next_char always exists
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]       # get the kmer at position i in the length specified by k at the beginningg
        next_char = sequence[i+k]    # get the nucleotide after the kmer
 
        kmer_data = update_kmer_count(kmer_data, kmer, next_char) #put into the count dictionary
 
    return kmer_data
 
 
#takes finished kmer_data with counts and writes it to a file 
def write_results_to_file(kmer_data, output_filename):
    with open(output_filename, 'w') as f: #open writing file
        for kmer in sorted(kmer_data.keys()):  # sort kmers alphabetically
            total_count = kmer_data[kmer]['count'] # pull out total counts for each kmer and the follow up nucleotide from the dictionary
            next_chars = kmer_data[kmer]['next_chars']
 
            
            next_char_str = " ".join( #formatting my reported value of the nucleotide and how often
                f"{char}:{freq}"
                for char, freq in sorted(next_chars.items())
            )
 
            f.write(f"{kmer} {total_count} {next_char_str}\n") #prints the values in a line that shows the kmer, the follow up nucleotide, and the frequency it comes up
 
 
#Tying all of this together 
#The fix in here removes the "for dequence in f" portion of the loop and rewrites it so it does not overwrite and resent the last saved sequence
def main(): #go back to specified python prompt from the terminal command for the file, kmer length, and output file
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
 
    print(f"Reading sequences from {sequence_file}...")
 
    kmer_data = {}  # accumulate counts across all sequences
 
    with open(sequence_file, 'r') as f: #reading one file at a time 
        for line in f:
            sequence = line.strip()  # remove  charcater at end
 
            if not sequence:  # skips blank lines
                continue
 
            if not validate_sequence(sequence, k): #skipping invalid sequences and providing a warning message
                print(f"  Warning: Skipping invalid sequence: {sequence!r}")
                continue
 
            # count kmers in the sequence
            sequence_kmer_data = count_kmers_with_context(sequence, k) 
 
            # merge into running totals
            for kmer, data in sequence_kmer_data.items():
                if kmer not in kmer_data:
                    kmer_data[kmer] = {'count': 0, 'next_chars': {}}
                kmer_data[kmer]['count'] += data['count']
                for char, freq in data['next_chars'].items():
                    if char not in kmer_data[kmer]['next_chars']:
                        kmer_data[kmer]['next_chars'][char] = 0
                    kmer_data[kmer]['next_chars'][char] += freq
 
    #write up the final totals and put them together 
    write_results_to_file(kmer_data, output_file)
    print(f"Results written to {output_file}")
 
if __name__ == '__main__':
    main()

