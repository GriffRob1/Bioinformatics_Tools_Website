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
        'toolTitle': '',
        'textDescription': '',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': '',
        'inputFormat': '',
        'outputFormat': '',
        'inputExample': '',
        'listArgumentIndex': -1
    },
    {
        'popularity': 0,
        'URL': 'convolution_cyclopeptide_sequencing',
        'imagePath': "/images/Cyclopeptide Sequencing Icon.png",
        'toolTitle': '',
        'textDescription': '',
        'dateAdded': "2025-06-04T09:38:00",
        'category': "Cyclopeptide Sequencing",
        'longDescription': '',
        'inputFormat': '',
        'outputFormat': '',
        'inputExample': '',
        'listArgumentIndex': -1
    },
]