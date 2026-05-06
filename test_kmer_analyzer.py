#Importing test tool and also print out my functions this way I can test to make sure they all work with the user running my code
import pytest
from kmer_analyzer import validate_sequence, update_kmer_count, count_kmers_with_context, write_results_to_file


#Making sure my validate_sequence command actually detects a normal valid sequence
def test_validate_sequence_basic():
    result = validate_sequence('ATGC', 2)
    assert result == True


#Makes sure the script returns a FALSE when running into a kmer with nothing after it
def test_validate_sequence_too_short():
    result = validate_sequence('AT', 2)
    assert result == False


#Makes sure we are only counting valid letters
def test_validate_sequence_invalid_base():
    result = validate_sequence('ATGX', 2)
    assert result == False


#makes sure we are only counting letters not numbers
def test_validate_sequence_digit():
    result = validate_sequence('ATG1', 2)
    assert result == False


#Only sees upper case no lower case
def test_validate_sequence_lowercase():
    result = validate_sequence('atgc', 2)
    assert result == False


#Only looking at non-empty strings
def test_validate_sequence_empty():
    result = validate_sequence('', 2)
    assert result == False


#Testing on ATG to make sure it is counting exactly 1 not 2 like before
def test_update_kmer_count_new_kmer():
    result = update_kmer_count({}, 'AT', 'G')
    obs = result['AT']['count']
    exp = 1
    assert obs == exp


#Making sure the follow up of G in ATG is being counted and reported properly
def test_update_kmer_count_next_char():
    result = update_kmer_count({}, 'AT', 'G')
    obs = result['AT']['next_chars']['G']
    exp = 1
    assert obs == exp


#Seeing the same kmer twice should give count of 2 - confirms the double count bug is fixed
def test_update_kmer_count_increments():
    data = update_kmer_count({}, 'AT', 'G')
    data = update_kmer_count(data, 'AT', 'G')
    obs = data['AT']['count']
    exp = 2
    assert obs == exp


#Two different follow up nucleotides after the same kmer should be tracked separately
def test_update_kmer_count_different_next_chars():
    data = update_kmer_count({}, 'TG', 'C')
    data = update_kmer_count(data, 'TG', 'A')
    assert data['TG']['next_chars']['C'] == 1
    assert data['TG']['next_chars']['A'] == 1


#Basic check that the sliding window finds a kmer correctly
def test_count_kmers_basic():
    result = count_kmers_with_context('ATGT', 2)
    obs = result['AT']['count']
    exp = 1
    assert obs == exp


#AT appears twice in ATATA so count should be 2
def test_count_kmers_repeated():
    result = count_kmers_with_context('ATATA', 2)
    obs = result['AT']['count']
    exp = 2
    assert obs == exp


#In ATGT, AT is followed by G - making sure next char is recorded correctly
def test_count_kmers_next_char_correct():
    result = count_kmers_with_context('ATGT', 2)
    obs = result['AT']['next_chars']['G']
    exp = 1
    assert obs == exp


#Making sure a kmer at the end with no follow up is not recorded
def test_count_kmers_last_kmer_excluded():
    result = count_kmers_with_context('ATGC', 3)
    assert 'TGC' not in result


#Making sure I'm not counting empty spaces
def test_count_kmers_empty_sequence():
    result = count_kmers_with_context('', 2)
    assert result == {}


#Making sure the output line looks exactly right - kmer, total count, then next char frequencies
def test_write_results_format(tmp_path):
    kmer_data = {'AT': {'count': 2, 'next_chars': {'G': 2}}}
    output_file = tmp_path / "out.txt"
    write_results_to_file(kmer_data, str(output_file))
    lines = output_file.read_text().strip().split('\n')
    assert lines[0] == 'AT 2 G:2'


#Kmers should come out in alphabetical order in the output file
def test_write_results_sorted(tmp_path):
    kmer_data = {
        'TG': {'count': 1, 'next_chars': {'A': 1}},
        'AT': {'count': 1, 'next_chars': {'G': 1}},
    }
    output_file = tmp_path / "out.txt"
    write_results_to_file(kmer_data, str(output_file))
    lines = output_file.read_text().strip().split('\n')
    assert lines[0].startswith('AT')
    assert lines[1].startswith('TG')


#Making sure bug was fixed where counts are translating to the output file
def test_write_results_includes_total_count(tmp_path):
    kmer_data = {'AT': {'count': 3, 'next_chars': {'G': 3}}}
    output_file = tmp_path / "out.txt"
    write_results_to_file(kmer_data, str(output_file))
    content = output_file.read_text()
    assert '3' in content
