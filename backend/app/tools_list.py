tools_list = [
    {
        'popularity': 0,
        'URL': "approximate_pattern_occurrences",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Approximate Pattern Occurrences",
        'textDescription': "Find all occurrences of a pattern with some mismatches",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds all approximate occurrences of a pattern in a text. This is defined by all patterns with d or less mismatches to the input pattern.',
        'inputFormat': 'distance (d) - the maximum number of mismatches for a string to be considered an approximate occurrence of the input pattern\npattern - the pattern being searched for\ntext - the full text that the pattern is being searched for in',
        'outputFormat': '[index, occurrence]\nindex - the location in text where the approximate pattern occurs\noccurrence - the pattern at that index',
        'inputExample': '1\nAGGA\nAGGATCGAGATTACGAATGCTAGATATCG',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "frequent_words_mismatch_reverse_complement",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Find Frequent Kmers",
        'textDescription': "Find the most frequent kmers with mismatches and reverse complement",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds the most frequent kmers in a text. The count for a kmer is incremented for all kmers within d mismatches from it or its reverse complement.',
        'inputFormat': 'distance (d) - the maximum number of mismatches for two patterns to be counted for the same kmer\nk - the length of the frequent kmers you are searching for\ntext - the text in which you are searching for frequent kmers',
        'outputFormat': 'pattern - the most frequent pattern with mismatches and reverse complement (can be multiple)',
        'inputExample': '1\n5\nAGGATCGAGATTACGATTGCTAGATATCG',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "reverse_complement",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Reverse Complement",
        'textDescription': "Generate the reverse complement of a string",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool generates the reverse complement of a dna string. This is created by taking the complement of each nucleotide, then reversing the entire string.',
        'inputFormat': 'string - the dna string',
        'outputFormat': 'reverse complement - the reverse complement of the dna string',
        'inputExample': 'AGGATCGAGATTACGAATGCTAGATATCG\n',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "find_clumps",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Find Clumps",
        'textDescription': "Find clumps of kmers in a text",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds kmers that occur very frequently in a small section of a text, creating a clump.',
        'inputFormat': 'k - the length of the kmers\nwindow size - the length of the area where kmers could clump\nminimum occurrences - the minimum number of occurrences of a kmer in a small section for it to be considered a clump\ntext - the text in which you are searching for clumps',
        'outputFormat': '(index, kmer)\nindex - the index in text of the start of the window where there is a clump\nkmer - the kmer that is clumped in a window starting at index',
        'inputExample': '3\n20\n4\nAGGATCGAGATTACGAATGGATCGTATAGCGTAGATAGCGTATAAGAGATCGCGTCTTAGAGATTAAGCTAGGCTACTAGATATCG',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "find_skew",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Find Skew",
        'textDescription': "Find the GC skew of a genome",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds the GC skew of a genome. At each point in the genome it adds 1 for a G and subtracts 1 for a C. This can be used to determine the ori region of some bacterial genomes. Many bacterial genomes as they have a uneven distribution of cytosine and guanine that makes their skews reach a minimum at the point where genome transcription begins.',
        'inputFormat': 'genome - the genome you are calculating the skew for',
        'outputFormat': 'skew matrix - a list of skews at every point in the genome',
        'inputExample': 'AGGATCGAGATTACGAATGGATCGTATAGCGTAGATAGCGTATAAGAGATCGCGTCTTAGAGATTAAGCTAGGCTACTAGATATCG',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "hamming_distance",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Hamming Distance",
        'textDescription': "Find the hamming distance between two strings",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds the hamming distance between two strings. This is the number of mismatched characters between them. The strings must be the same length.',
        'inputFormat': 'string 1 - the first string\nstring 2 - the second string',
        'outputFormat': 'hamming distance - the hamming distance between the strings',
        'inputExample': 'CTGATCGAGATGACGAATGATAGATATCG\nAGGATCGAGATTACGAATGCTAGATATCG\n',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "neighbors",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Find Neighbors",
        'textDescription': "Find the neighbors of a pattern",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds the neighbors of a pattern. This includes all possible kmers that are a distance d from the input pattern.',
        'inputFormat': 'distance (d) - the distance between the pattern and its neighbors\npattern - the pattern you want to find the neighbors of',
        'outputFormat': 'neighbors list - a list of all the neighbors to a pattern',
        'inputExample': '2\nAGGA',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': "find_potential_dnaa_boxes",
        'imagePath': "/images/Frequent Pattern Finding Logo.png",
        'toolTitle': "Find DnaA Box",
        'textDescription': "Find possible DnaA boxes for bacterial genomes",
        'dateAdded': "2025-05-20T03:12:00",
        'category': "Frequent Pattern Finding",
        'longDescription': 'This tool finds a list of possible DnaA boxes for bacterial genomes. A DnaA box is a sequence in bacterial genomes where the DnaA protein binds, initiating genome replication. This tool works by finding the region containing the minimum skew of the genome, then finding a frequent pattern with mismatches and reverse complement in that region.',
        'inputFormat': 'k - the length of the DnaA box\ndistance (d) - the maximum number of mismatches for a pattern to be counted as an occurrence of the DnaA box\nwindow size - the size of the region around the minimum skew where you want to search for a DnaA box\nwindow location - the location of the window around the minimum skew index, the only options are "before", "after", and "centered"\ngenome - the genome in which you want to find a DnaA box',
        'outputFormat': 'DnaA boxes - a list of potential DnaA boxes (may be multiple)',
        'inputExample': '3\n1\n30\nafter\nATATAGTAGATAAGGTATCGATTTAGCGATCGATTCGGCTATCGTTCGCCTCTCATATCTTCTATCTTAGCATGCTATGCTACTGCTAGGATCGTACGTCAGGCGACGCGAGCGGACGGCGTATCGACGTACTAGCGCTAGCTACTG',
        'listArgumentIndex': -1
    },













    {
        'popularity': 0,
        'URL': 'count_matrix',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Count Matrix',
        'textDescription': 'Generate a count matrix for a set of motifs',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool generates a count matrix for a set of motifs. The matrix counts the number of nucleotides from all the motifs at each position.',
        'inputFormat': 'motifs - a list of motifs, each on a new line. Motifs must be the same length',
        'outputFormat': 'count matrix - The rows count nucleotides in this order, from top to bottom: A,C,G,T. The columns represent each position of the motif.',
        'inputExample': 'CTTGATAGG\nCATCATAAG\nCTTGTTAGG\nCCTGATAGT\nCTTAATTGG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'profile_no_psuedocounts',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Profile Matrix',
        'textDescription': 'Generate the profile matrix for a set of motifs',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool generates the profile matrix for a set of motifs. This is similar to the count matrix, but instead of individual nucleotide counts, the matrix stores the fraction of each nucleotide in all the motifs at a specific position.',
        'inputFormat': 'motifs - a list of motifs, each on a new line. Motifs must be the same length',
        'outputFormat': 'profile matrix - The rows count nucleotides in this order, from top to bottom: A,C,G,T. The columns represent each position of the motif.',
        'inputExample': 'CTTGATAGG\nCATCATAAG\nCTTGTTAGG\nCCTGATAGT\nCTTAATTGG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'consensus',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Consensus Motif',
        'textDescription': 'Generate a consensus motif for a list of motifs',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool generates a consensus motif for a list of motifs. This is the string made up of the most common nucleotide at each position.',
        'inputFormat': 'motifs - a list of motifs, each on a new line. Motifs must be the same length',
        'outputFormat': 'consensus string - the consensus string',
        'inputExample': 'CTTGATAGG\nCATCATAAG\nCTTGTTAGG\nCCTGATAGT\nCTTAATTGG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'score',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Motif Score',
        'textDescription': 'Calculate the score for a list of motifs',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool calculates the score for a list of motifs. This is the number of mismatches between all motifs and the consensus motif. A lower score means the motifs are more similar to each other.',
        'inputFormat': 'motifs - a list of motifs, each on a new line. Motifs must be the same length',
        'outputFormat': 'score - the score of the motifs',
        'inputExample': 'CTTGATAGG\nCATCATAAG\nCTTGTTAGG\nCCTGATAGT\nCTTAATTGG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'median_string_all_kmers',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Brute Force Motif Search',
        'textDescription': 'Search all possible kmers for the best motifs',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool searches multiple dna sequences for the set of motifs with the lowest score. The tool searches all possible kmers for the consensus motif.',
        'inputFormat': 'k - the length of the motifs\ndna - a list of dna sequences that you are searching for motifs in',
        'outputFormat': 'best motifs - a set of motifs, one from each dna sequence',
        'inputExample': '4\nGATAGGTCGCGTATAGATTCGATACGCTAG\nGATCGATATGCGGCTATAGCTAGCTAGCGA\nTGAGTCGGCTAGTAGCGCTAGTATCTCCGA\nTAGGCTGAGATCTCTGGATTAGGCAGATCG',
        'listArgumentIndex': 1
    },
    {
        'popularity': 0,
        'URL': 'greedy_motif_search',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Greedy Motif Search',
        'textDescription': 'Find the best motifs using a greedy algorithm',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool uses a greedy algorithm to find the best motifs in a set of dna sequences. This is less accurate than the brute force motif search, but much faster to compute.',
        'inputFormat': 'k - the length of the motifs\ndna - a list of dna sequences that you are searching for motifs in',
        'outputFormat': 'best motifs - a set of motifs, one from each dna sequence',
        'inputExample': '4\nGATAGGTCGCGTATAGATTCGATACGCTAG\nGATCGATATGCGGCTATAGCTAGCTAGCGA\nTGAGTCGGCTAGTAGCGCTAGTATCTCCGA\nTAGGCTGAGATCTCTGGATTAGGCAGATCG',
        'listArgumentIndex': 1
    },
    {
        'popularity': 0,
        'URL': 'iterate_randomized_motif_search',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Randomized Motif Search',
        'textDescription': 'Find the best motifs using a randomized algorithm',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool finds motifs using a randomized algorithm, using the profile matrix of one iteration to determine the set of motifs for the next iteration. This is faster than the greedy algorithm, but less accurate and reliable, as it is randomized.',
        'inputFormat': 'k - the length of the motifs\niterations - the number of iterations for the algorithm to do. More iterations means a slower calculation but more accurate results.\ndna - a list of dna sequences that you are searching for motifs in',
        'outputFormat': 'best motifs - a set of motifs, one from each dna sequence',
        'inputExample': '4\n10\nGATAGGTCGCGTATAGATTCGATACGCTAG\nGATCGATATGCGGCTATAGCTAGCTAGCGA\nTGAGTCGGCTAGTAGCGCTAGTATCTCCGA\nTAGGCTGAGATCTCTGGATTAGGCAGATCG',
        'listArgumentIndex': 2
    },
    {
        'popularity': 0,
        'URL': 'iterate_gibbs_sampler',
        'imagePath': "/images/Motif Finding Logo.png",
        'toolTitle': 'Gibbs Sampler Motif Search',
        'textDescription': 'Find the best motifs using a gibbs sampler',
        'dateAdded': "2025-05-21T09:50:00",
        'category': "Motif Finding",
        'longDescription': 'This tool uses the gibbs sampling technique to find motifs in a set of dna sequences. This is more accurate than the randomized motif search, because it randomly selects motifs one at a time rather than all at once.',
        'inputFormat': 'k - the length of the motifs\nrepeats - the number of repeats the gibbs sampler does for one interation\niterations - the number of iterations of the gibbs sampler. Each iteration randomly selects new starting motifs.\ndna - a list of dna sequences that you are searching for motifs in',
        'outputFormat': 'best motifs - a set of motifs, one from each dna sequence',
        'inputExample': '4\n10\n20\nGATAGGTCGCGTATAGATTCGATACGCTAG\nGATCGATATGCGGCTATAGCTAGCTAGCGA\nTGAGTCGGCTAGTAGCGCTAGTATCTCCGA\nTAGGCTGAGATCTCTGGATTAGGCAGATCG',
        'listArgumentIndex': 3
    },


















    {
        'popularity': 0,
        'URL': 'kmer_composition',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'Kmer Composition',
        'textDescription': 'Generate the kmer composition of a read',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool generates all substrings of length k in a read.',
        'inputFormat': 'k - the length of the kmers\ntext - the text you want to get the kmer composition of',
        'outputFormat': 'kmers - the list of kmers',
        'inputExample': '3\nAGGATATAGGATCGA',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'path_to_genome',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'Kmer Path to Genome',
        'textDescription': 'Generate a genome from a path of kmers',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool generates a genome given the kmers in a pathway through a De Bruijn graph.',
        'inputFormat': 'kmers - the list of kmers',
        'outputFormat': 'genome - the completed genome',
        'inputExample': 'AGG\nGGA\nGAT\nATA\nTAT\nATA\nTAG\nAGG\nGGA\nGAT\nATC\nTCG\nCGA',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'read_pairs_path_to_genome',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'Read Pair Path to Genome',
        'textDescription': 'Generate a genome from a path of read pairs',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool generates a genome given the read pairs in a pathway through a paired De Bruijn graph.',
        'inputFormat': 'd - the distance in between the reads in the pair\nread pairs - the list of read pairs, with the read pair separated by a | symbol.',
        'outputFormat': 'genome - the completed genome',
        'inputExample': '3\nAGG|TAG\nGGA|AGG\nGAT|GGA\nATA|GAT\nTAT|ATC\nATA|TCG\nTAG|CGA\nAGG|GAG\nGGA|AGA\nGAT|GAT',
        'listArgumentIndex': 1
    },
    {
        'popularity': 0,
        'URL': 'debruijn_graph_from_kmers',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'De Bruijn Graph From Kmers',
        'textDescription': 'Create a De Bruijn graph from a list of kmers',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool generates an adjacency list representing the De Bruijn graph for a list of kmers.',
        'inputFormat': 'kmers - the list of kmers',
        'outputFormat': 'graph - the adjacency list',
        'inputExample': 'CGA\nAGG\nGAT\nATA\nTAG\nATA\nATC\nGGA\nGGA\nGAT\nTAT\nTCG\nAGG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'paired_debruijn_graph',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'Paired De Bruijn Graph',
        'textDescription': 'Generate a paired de bruijn graph from read pairs',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool generates an adjacency list representing the paired De Bruijn graph for a list of read pairs.',
        'inputFormat': 'read pairs - the list of read pairs',
        'outputFormat': 'graph - the adjacency list',
        'inputExample': 'ATA|GAT\nTAT|ATC\nGAT|GGA\nAGG|TAG\nGGA|AGG\nATA|TCG\nGGA|AGA\nGAT|GAT\nTAG|CGA\nAGG|GAG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'debruijn_assembler',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'Single Read Genome Assembler',
        'textDescription': 'Assemble a genome based on single reads',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool uses a de bruijn graph to assemble a genome from reads. This is a simplified assembler, and assumes total coverage with no errors.',
        'inputFormat': 'reads - the reads used to assemble the genome',
        'outputFormat': 'genome - the completed genome',
        'inputExample': 'CGA\nAGG\nGAT\nATA\nTAG\nATA\nATC\nGGA\nGGA\nGAT\nTAT\nTCG\nAGG',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'debruijn_assembler_read_pairs',
        'imagePath': "/images/Genome Assembly Logo.png",
        'toolTitle': 'Read Pair Genome Assembler',
        'textDescription': 'Assemble a genome based on read pairs',
        'dateAdded': "2025-05-23T20:50:00",
        'category': "Genome Assembly",
        'longDescription': 'This tool uses a de bruijn graph to assemble a genome from read pairs. This is a simplified assembler, and assumes total coverage with no errors.',
        'inputFormat': 'd - the distance in between the reads in the pair\nread pairs - the read pairs used to assemble the genome',
        'outputFormat': 'genome - the completed genome',
        'inputExample': '3\nATA|GAT\nTAT|ATC\nGAT|GGA\nAGG|TAG\nGGA|AGG\nATA|TCG\nGGA|AGA\nGAT|GAT\nTAG|CGA\nAGG|GAG',
        'listArgumentIndex': 1
    },

















    {
        'popularity': 0,
        'URL': 'transcribe',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Transcribe DNA',
        'textDescription': 'Transcribe a DNA sequence into RNA',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tool replaces all the T\'s with U\'s in the DNA sequence, transcribing it into RNA.',
        'inputFormat': 'DNA - the DNA sequence',
        'outputFormat': 'RNA -  the transcribed RNA sequence',
        'inputExample': 'ACCGATAGCTATAGCGCT',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'translate',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Translate RNA',
        'textDescription': 'Translate an RNA sequence into an amino acid sequence',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tool generates an amino acid sequence from the given RNA sequence, stopping either when the string ends or when it reaches a stop codon.',
        'inputFormat': 'RNA - the RNA sequence',
        'outputFormat': 'Protein - the translated amino acid sequence, in single letter amino acid code',
        'inputExample': 'ACCGAUAGCUAUAGCGCU',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'amino_acid_masses_to_letters',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Amino Acid Mass to Letter',
        'textDescription': 'Convert a list of amino acid masses to letters',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tool converts a list of amino acid masses to the letter codes corresponding to each amino acid.',
        'inputFormat': 'masses - a list of amino acid masses',
        'outputFormat': 'peptide - the corresponding letters based on the masses. For a mass that does not correspond to a main amino acid, a ? is shown.',
        'inputExample': '71 101 115 87 23 163 87',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'linear_and_cyclo_peptide_spectrum',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Peptide Spectrum',
        'textDescription': 'Generate the theoretical linear and cyclical spectra of a peptide',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tool generates two theoretical mass spectra of a peptide, one assuming the peptide is linear, and one assuming it is cyclical. This is a simplified mass spectra, so instead of a mass/charge ratio, this simply calculates the mass of each subpeptide.',
        'inputFormat': 'peptide - the peptide that the spectra are being generated for, in single letter amino acid code format.',
        'outputFormat': 'linear spectrum - the spectrum if the peptide is linear\ncyclical spectrum - the spectrum if the peptide is cyclical',
        'inputExample': 'TDSYSA',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'brute_force_cyclopeptide_sequencer',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Brute Force Cyclopeptide Sequencer',
        'textDescription': 'Sequence a spectrum by checking all possible cyclopeptides',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tool is a very inefficient algorithm that checks all possible cyclopeptides for a match with the given spectrum.',
        'inputFormat': 'spectrum - the spectrum you want to sequence. This is a list of masses separated by whitespace.',
        'outputFormat': 'peptide - the peptide that matches the spectrum',
        'inputExample': '0 87 101 115 188 202 216 303',
        'listArgumentIndex': 0
    },
    {
        'popularity': 0,
        'URL': 'leaderboard_cyclopeptide_sequencer',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Leaderboard Cyclopeptide Sequencer',
        'textDescription': 'Sequence a cyclopeptide with a leaderboard branch and bound algorithm',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tools sequences cyclopeptides with a branch and bound algorithm. The algorithm adds one amino acid at a time, generating all possible new peptides. Then, it scores each one and takes the top n scoring peptides for the next iteration. It outputs the highest scoring peptide from all iterations.',
        'inputFormat': 'n - the number of peptides in the leaderboard. A higher n means more accurate results, but a slower runtime.\nspectrum - the spectrum you want to sequence',
        'outputFormat': 'cyclopeptide - the highest scoring cyclopeptide from all iterations. Each number is the mass of an amino acid in the cyclopeptide, and it is represented this way because some amino acids have the same mass.',
        'inputExample': '100\n0 71 87 87 101 115 158 163 172 202 216 250 250 259 287 303 321 337 365 374 374 408 422 452 461 466 509 523 537 537 553 624',
        'listArgumentIndex': 1
    },
    {
        'popularity': 0,
        'URL': 'convolution_cyclopeptide_sequencing',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': 'Convolution Cyclopeptide Sequencer',
        'textDescription': 'Sequence a peptide with non-proteinogenic amino acids',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': 'This tool uses a spectral convolution to first determine what the amino acid composition of a peptide is before sequencing it. It is also useful for peptides with non-proteinogenic amino acids, because the algorithm treats any mass between 57 and 200 as a possible amino acid.',
        'inputFormat': 'm - the m most frequently occurring amino acids in the convolution will be considered as a possible composition for the peptide. All other amino acids will not be considered.\nn - the number of peptides in the leaderboard. A higher n means more accurate results, but a slower runtime.\nspectrum - the spectrum you want to sequence',
        'outputFormat': 'cyclopeptide - the highest scoring cyclopeptide from all iterations. Each number is the mass of an amino acid in the cyclopeptide, and it is represented this way because some amino acids have the same mass.',
        'inputExample': '6\n100\n0 71 87 87 101 115 158 163 172 202 216 250 250 259 287 303 321 337 365 374 374 408 422 452 461 466 509 523 537 537 553 624',
        'listArgumentIndex': 2
    },















    {
        'popularity': 0,
        'URL': 'scored_pairwise_global_alignment',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Global Alignment',
        'textDescription': 'Align two nucleotide or protein sequences',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise alignment using a dynamic programming algorithm. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\nindel score: the score for insertions or deletions. Usually a negative number.\nsequences: two sequences, each on their own line, to align.',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\nalignment: the aligned sequences',
        'inputExample': 'DNA_MATRIX\n-2\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nCGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGAC',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_local_alignment',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Local Alignment',
        'textDescription': 'Align two nucleotide or protein sequences locally',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise local alignment using a dynamic programming algorithm. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\nindel score: the score for insertions or deletions. Usually a negative number.\nsequences: two sequences, each on their own line, to align locally.',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\naligned sequence 1\nstarting index 1: the index of sequence 1 that the alignment starts at\naligned sequence 2\nstarting index 2: the index of sequence 2 that the alignment starts at',
        'inputExample': 'DNA_MATRIX\n-2\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nCGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGAC',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_fitting_alignment',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Fitting Alignment',
        'textDescription': 'Align a shorter sequence to a section of a longer sequence',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise fitting alignment using a dynamic programming algorithm. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\nindel score: the score for insertions or deletions. Usually a negative number.\ncontainer sequence: the longer sequence that contains the fitting sequence\nfitting sequence: the shorter sequence that will be aligned with a substring of the container sequence',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\naligned container sequence\ncontainer sequence index: the index of the container sequence that the alignment starts at\naligned fitting sequence',
        'inputExample': 'DNA_MATRIX\n-2\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nACCTCGGATC',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_overlap_alignment',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Overlap Alignment',
        'textDescription': 'Align the overlap of two nucleotide or protein sequences',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise overlap alignment using a dynamic programming algorithm. Overlap alignment means that instead of the entire string being aligned, the ending of the first sequence is aligned with the beginning of the second, creating an overlap. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\nindel score: the score for insertions or deletions. Usually a negative number.\nsequences: two sequences, each on their own line, to align.',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\naligned sequence 1: this is an alignment of the suffix of sequence 1. It starts from somewhere in the middle of the sequence and goes to the end\nstarting index 1: the index of sequence 1 that the alignment starts at\naligned sequence 2: this is an alignment of the prefix of sequence 2. It starts from the beginning of the sequence and goes to somewhere in the middle',
        'inputExample': 'DNA_MATRIX\n-2\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nGATCTAGCTAGCTAGCGACCGCCCTACCTTTCGTGCCTATCAAA',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_global_alignment_affine_gap_penalty',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Global Alignment With Affine Gap',
        'textDescription': 'Align two nucleotide or protein sequences with an affine gap penalty',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise alignment using a dynamic programming algorithm, and uses affine gap penalties to reduce the number of gaps. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\ninitial gap cost: the cost for starting a gap. Generally this is a negative number to avoid creating gaps in the alignment.\nadditional gap cost: this is the added cost for every space added to a gap. This is usually a higher value than the initial cost, to incentivise longer gaps over shorter ones.\nsequences: two sequences, each on their own line, to align.',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\nalignment: the aligned sequences',
        'inputExample': 'DNA_MATRIX\n-3\n-1\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nCGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGAC',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_local_alignment_affine_gap_penalty',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Local Alignment With Affine Gap',
        'textDescription': 'Align two nucleotide or protein sequences locally with an affine gap penalty',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise local alignment using a dynamic programming algorithm, and uses affine gap penalties to reduce the number of gaps. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\ninitial gap cost: the cost for starting a gap. Generally this is a negative number to avoid creating gaps in the alignment.\nadditional gap cost: this is the added cost for every space added to a gap. This is usually a higher value than the initial cost, to incentivise longer gaps over shorter ones.\nsequences: two sequences, each on their own line, to align locally.',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\naligned sequence 1\nstarting index 1: the index of sequence 1 that the alignment starts at\naligned sequence 2\nstarting index 2: the index of sequence 2 that the alignment starts at',
        'inputExample': 'DNA_MATRIX\n-3\n-1\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nCGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGAC',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_fitting_alignment_affine_gap_penalty',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Fitting Alignment With Affine Gap',
        'textDescription': 'Align a shorter sequence to a section of a longer sequence with an affine gap penalty',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise fitting alignment using a dynamic programming algorithm, and uses affine gap penalties to reduce the number of gaps. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\ninitial gap cost: the cost for starting a gap. Generally this is a negative number to avoid creating gaps in the alignment.\nadditional gap cost: this is the added cost for every space added to a gap. This is usually a higher value than the initial cost, to incentivise longer gaps over shorter ones.\ncontainer sequence: the longer sequence that contains the fitting sequence\nfitting sequence: the shorter sequence that will be aligned with a substring of the container sequence',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\naligned container sequence\ncontainer sequence index: the index of the container sequence that the alignment starts at\naligned fitting sequence',
        'inputExample': 'DNA_MATRIX\n-3\n-1\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nACCTCGGATC',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'scored_pairwise_overlap_alignment_affine_gap_penalty',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Pairwise Overlap Alignment With Affine Gap',
        'textDescription': 'Align the overlap of two nucleotide or protein sequences with an affine gap penalty',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs pairwise overlap alignment using a dynamic programming algorithm, and uses affine gap penalties to reduce the number of gaps. Overlap alignment means that instead of the entire string being aligned, the ending of the first sequence is aligned with the beginning of the second, creating an overlap. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\ninitial gap cost: the cost for starting a gap. Generally this is a negative number to avoid creating gaps in the alignment.\nadditional gap cost: this is the added cost for every space added to a gap. This is usually a higher value than the initial cost, to incentivise longer gaps over shorter ones.\nsequences: two sequences, each on their own line, to align.',
        'outputFormat': 'score: the score of the alignment. A higher score means a better alignment.\naligned sequence 1: this is an alignment of the suffix of sequence 1. It starts from somewhere in the middle of the sequence and goes to the end\nstarting index 1: the index of sequence 1 that the alignment starts at\naligned sequence 2: this is an alignment of the prefix of sequence 2. It starts from the beginning of the sequence and goes to somewhere in the middle',
        'inputExample': 'DNA_MATRIX\n-3\n-1\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nGATCTAGCTAGCTAGCGACCGCCCTACCTTTCGTGCCTATCAAA',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'greedy_multiple_alignment',
        'imagePath': "/images/Sequence Alignment Icon.png",
        'toolTitle': 'Greedy Multiple Alignment',
        'textDescription': 'Align more than two nucleotide or protein sequences',
        'dateAdded': "2025-07-16T07:52:00",
        'category': "Sequence Alignment",
        'longDescription': 'This tool performs multiple global alignment using a greedy algorithm. It repeatedly finds the most similar sequence to the multiple alignment, then performs pairwise alignment with the consensus string and adds the aligned sequence to the multiple alignment. It can be used for nucleotide or protein sequences, depending on the scoring matrix chosen.',
        'inputFormat': 'scoring matrix: the scoring matrix used for the alignment. The options for matrices are BLOSUM62 and PAM250 for protein sequences, and DNA_MATRIX for dna sequences. DNA_MATRIX is a simple scoring matrix with +1 for matches and -1 for mismatches.\nindel score: the score for insertions or deletions. Usually a negative number.\nsequences: many sequences, each on their own line, to align.',
        'outputFormat': 'alignment: the aligned sequences',
        'inputExample': 'DNA_MATRIX\n-2\nCAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTAT\nCGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGAC\nGCCCTCCTAGTCGTGCCCTTCAAATGCTCTGCTAGCGTAGGAC\nGCGCCTTCTTTGTGCCTCATCAACGATCTAGCTAAACTACGAC',
        'listArgumentIndex': 2
    },
]