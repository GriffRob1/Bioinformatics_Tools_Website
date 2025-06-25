import copy
import itertools
import math

codon_table = {
    'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAU': 'N',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACU': 'T',
    'AGA': 'R',
    'AGC': 'S',
    'AGG': 'R',
    'AGU': 'S',
    'AUA': 'I',
    'AUC': 'I',
    'AUG': 'M',
    'AUU': 'I',
    'CAA': 'Q',
    'CAC': 'H',
    'CAG': 'Q',
    'CAU': 'H',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCU': 'P',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGU': 'R',
    'CUA': 'L',
    'CUC': 'L',
    'CUG': 'L',
    'CUU': 'L',
    'GAA': 'E',
    'GAC': 'D',
    'GAG': 'E',
    'GAU': 'D',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCU': 'A',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGU': 'G',
    'GUA': 'V',
    'GUC': 'V',
    'GUG': 'V',
    'GUU': 'V',
    'UAA': '*',
    'UAC': 'Y',
    'UAG': '*',
    'UAU': 'Y',
    'UCA': 'S',
    'UCC': 'S',
    'UCG': 'S',
    'UCU': 'S',
    'UGA': '*',
    'UGC': 'C',
    'UGG': 'W',
    'UGU': 'C',
    'UUA': 'L',
    'UUC': 'F',
    'UUG': 'L',
    'UUU': 'F'
}

amino_acid_mass = {
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'I': 113,
    'L': 113,
    'N': 114,
    'D': 115,
    'K': 128,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186
}




def transcribe(dna):
    return dna.replace('T','U')



def translate(rna):
    amino_acid_string = ''
    for i in range(0, len(rna), 3):
        amino_acid_string += codon_table[rna[i:i+3]]
        if amino_acid_string[-1] == '*':
            return amino_acid_string[0:-1]
    return amino_acid_string



def amino_acid_masses_to_letters(masses):
    masses = [int(item) for item in masses]
    peptide = ''
    mass_to_letter = {v: k for k, v in amino_acid_mass.items()}
    for mass in masses:
        if mass in mass_to_letter.keys():
            peptide += mass_to_letter[mass]
        else:
            peptide += '?'
    return peptide



def cyclopeptide_to_spectrum(cyclopeptide):
    masses = [0, peptide_mass(cyclopeptide)] # initialize with full cyclopeptide as one of its subpeptides
    length = len(cyclopeptide)
    cyclopeptide = cyclopeptide + cyclopeptide # makes it easier to get substrings that wrap around
    for i in range(length):
        for j in range(1, length):
            subpeptide = cyclopeptide[i:i+j]
            masses.append(peptide_mass(subpeptide))
    masses.sort()
    return masses



def linear_and_cyclo_peptide_spectrum_all_amino_acids(peptide):
    linear_masses = [0, peptide_mass_all_amino_acids(peptide)]
    cyclo_masses = [0, peptide_mass_all_amino_acids(peptide)]
    length = len(peptide)
    peptide = peptide + peptide  # makes it easier to get subpeptides that wrap around
    for i in range(length):
        for j in range(1, length):
            subpeptide = peptide[i:i + j]
            if i + j <= length:
                linear_masses.append(peptide_mass_all_amino_acids(subpeptide))
            cyclo_masses.append(peptide_mass_all_amino_acids(subpeptide))
    linear_masses.sort()
    cyclo_masses.sort()
    return linear_masses, cyclo_masses



def linear_and_cyclo_peptide_spectrum(peptide):
    linear_masses = [0, peptide_mass(peptide)]
    cyclo_masses = [0, peptide_mass(peptide)]
    length = len(peptide)
    peptide = peptide + peptide  # makes it easier to get subpeptides that wrap around
    for i in range(length):
        for j in range(1, length):
            subpeptide = peptide[i:i + j]
            if i + j <= length:
                linear_masses.append(peptide_mass(subpeptide))
            cyclo_masses.append(peptide_mass(subpeptide))
    linear_masses.sort()
    cyclo_masses.sort()
    return linear_masses, cyclo_masses



def peptide_mass(subpeptide):
    mass = 0
    for amino_acid in subpeptide:
        mass += amino_acid_mass[amino_acid]
    return mass



def peptide_mass_all_amino_acids(subpeptide):
    total_mass = 0
    for mass in subpeptide:
        total_mass += mass
    return total_mass



def brute_force_cyclopeptide_sequencer(spectrum):
    spectrum = [int(mass) for mass in spectrum]
    all_peptides = all_peptides_for_mass(spectrum[-1])
    for peptide in all_peptides:
        if spectrum == cyclopeptide_to_spectrum(peptide):
            return peptide
    return None



def all_peptides_for_mass(mass):
    all_peptides = []
    max_length = int(mass / 57)
    for i in range(1, max_length+1):
        all_peptides_length_i = list(itertools.permutations(amino_acid_mass.keys(), i))
        for peptide in all_peptides_length_i:
            peptide_string = ''.join(peptide)
            if peptide_mass(peptide_string) == mass:
                all_peptides.append(peptide_string)
    return all_peptides



def cyclopeptide_sequencing_ideal_spectrum(spectrum):
    parent_mass = spectrum[-1]
    candidate_peptides = [('',-1,-1)]
    final_peptides = []
    while (len(candidate_peptides) > 0):
        candidate_peptides = expand_peptides(candidate_peptides)
        for peptide_tuple in candidate_peptides:
            if peptide_mass(peptide_tuple[0]) == parent_mass:
                if (cyclopeptide_to_spectrum(peptide_tuple[0]) == spectrum) and (peptide_tuple[0] not in final_peptides):
                    final_peptides.append(peptide_tuple[0])
                candidate_peptides.remove(peptide_tuple)
            elif not is_consistent_with_ideal_spectrum(peptide_tuple[0], spectrum):
                candidate_peptides.remove(peptide_tuple)
    return final_peptides



def expand_peptides(peptides):
    new_peptides = []
    for peptide_tuple in peptides:
        for amino_acid in amino_acid_mass.keys():
            new_peptides.append((peptide_tuple[0] + amino_acid, -1, -1)) # adds a 2-tuple with a peptide and its score
    return new_peptides



def expand_peptides_all_amino_acids(peptides, amino_acid_alphabet):
    new_peptides = []
    for peptide_tuple in peptides:
        for amino_acid in amino_acid_alphabet:
            new_peptide_sequence = peptide_tuple[0][:]
            new_peptide_sequence.append(amino_acid)
            new_peptides.append((new_peptide_sequence, -1, -1))  # adds a 2-tuple with a peptide and its score
    return new_peptides



def is_consistent_with_ideal_spectrum(peptide, spectrum):
    experimental_spectrum = copy.copy(spectrum)
    theoretical_spectrum = cyclopeptide_to_spectrum(peptide)
    for mass in theoretical_spectrum:
        try:
            experimental_spectrum.remove(mass)
        except ValueError:
            return False
    return True



def linear_and_cyclo_peptide_score(peptide, spectrum):
    experimental_spectrum = copy.copy(spectrum)
    theoretical_linear_spectrum, theoretical_cyclo_spectrum = linear_and_cyclo_peptide_spectrum_all_amino_acids(peptide)
    linear_score = 0
    cyclo_score = 0
    for mass in theoretical_cyclo_spectrum:
        if mass in experimental_spectrum:
            experimental_spectrum.remove(mass)
            cyclo_score += 1
            if mass in theoretical_linear_spectrum:
                linear_score += 1
    return linear_score, cyclo_score



def leaderboard_cyclopeptide_sequencer(n, spectrum):
    n = int(n)
    spectrum = [int(item) for item in spectrum]
    global max_score_peptide
    parent_mass = spectrum[-1]
    leaderboard = [([0],1,1)]
    leader_peptide = ([0],1,1)
    amino_acid_alphabet = list(set(amino_acid_mass.values()))

    while len(leaderboard) > 0:
        leaderboard = expand_peptides_all_amino_acids(leaderboard, amino_acid_alphabet) # keep old peptides and add new expanded peptides
        leaderboard_copy = copy.copy(leaderboard)
        for peptide_tuple in leaderboard_copy: # filter out peptides that are more massive than the spectrum allows
            if peptide_mass_all_amino_acids(peptide_tuple[0]) > parent_mass:
                leaderboard.remove(peptide_tuple)

        leaderboard = trim_leaderboard(leaderboard, spectrum, n)
        if len(leaderboard) > 0:
            max_score_peptide = max(leaderboard, key=lambda item: item[2])
        if max_score_peptide[2] > leader_peptide[2]:
            leader_peptide = max_score_peptide
    return leader_peptide[0][1:]



def spectral_convolution(spectrum):
    convolution = []
    for i in range(len(spectrum)):
        for j in range(i+1, len(spectrum)):
            difference = abs(spectrum[i] - spectrum[j])
            convolution.append(difference)
    return convolution



def filter_convolution(convolution, m):
    multiplicities = {}
    for mass in convolution:
        if 57 <= mass <= 200:
            if mass not in multiplicities:
                multiplicities[mass] = 1
            else:
                multiplicities[mass] += 1
    multiplicities = {k: v for k, v in sorted(multiplicities.items(), key=lambda item: item[1], reverse=True)}
    multiplicities = list(multiplicities)
    return multiplicities[:m]



def convolution_cyclopeptide_sequencing(m, n, spectrum):
    m = int(m)
    n = int(n)
    spectrum = [int(item) for item in spectrum]
    global max_score_peptide
    parent_mass = spectrum[-1]
    leaderboard = [([0],1,1)]
    leader_peptide = ([0],1,1)
    convolution = spectral_convolution(spectrum)
    amino_acid_alphabet = filter_convolution(convolution, m)

    while len(leaderboard) > 0:
        leaderboard = expand_peptides_all_amino_acids(leaderboard, amino_acid_alphabet) # keep old peptides and add new expanded peptides
        leaderboard_copy = copy.copy(leaderboard)
        for peptide_tuple in leaderboard_copy: # filter out peptides that are more massive than the spectrum allows
            if peptide_mass_all_amino_acids(peptide_tuple[0]) > parent_mass:
                leaderboard.remove(peptide_tuple)

        leaderboard = trim_leaderboard(leaderboard, spectrum, n)
        if len(leaderboard) > 0:
            max_score_peptide = max(leaderboard, key=lambda item: item[2])
        if max_score_peptide[2] > leader_peptide[2]:
            leader_peptide = max_score_peptide
    return leader_peptide[0][1:]




def trim_leaderboard(leaderboard, spectrum, n):
    new_leaderboard = []
    for peptide_tuple in leaderboard:
        if peptide_tuple[1] < 0:

            linear_score, cyclo_score = linear_and_cyclo_peptide_score(peptide_tuple[0], spectrum)
            new_leaderboard.append((peptide_tuple[0], linear_score, cyclo_score))
        else:
            new_leaderboard.append(peptide_tuple)

    new_leaderboard = sorted(new_leaderboard, key=lambda item: item[1], reverse=True)

    # for lists shorter than n
    if len(new_leaderboard) <= n:
        return new_leaderboard

    # includes ties with the nth place
    lowest_score = new_leaderboard[n-1][1]
    i = n
    while i < len(new_leaderboard) and lowest_score == new_leaderboard[i][1]:
        i += 1

    return new_leaderboard[:i]
















'''tyrocidine_b1 = 'VKLFPWFNQY'

spectrum1 = cyclopeptide_to_spectrum(tyrocidine_b1)
spectrum2 = '0 97 99 114 128 147 147 163 186 227 241 242 244 260 261 262 283 291 333 340 357 385 389 390 390 405 430 430 447 485 487 503 504 518 543 544 552 575 577 584 632 650 651 671 672 690 691 738 745 747 770 778 779 804 818 819 820 835 837 875 892 917 932 932 933 934 965 982 989 1030 1039 1060 1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322'
spectrum2 = spectrum2.split()
spectrum2 = [int(mass) for mass in spectrum2]
spectrum3 = '0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 333 340 347 357 385 388 389 390 390 405 430 430 435 447 485 487 503 504 518 543 544 552 575 577 584 599 608 631 632 650 651 653 671 672 690 691 717 738 745 747 770 778 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1031 1039 1060 1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322'
spectrum3 = spectrum3.split()
spectrum3 = [int(mass) for mass in spectrum3]

print(spectrum2)
print(spectral_convolution(cyclopeptide_to_spectrum('NQEL')))
print(convolution_cyclopeptide_sequencing(spectrum3, 12, 1000))
print(leaderboard_cyclopeptide_sequencer(spectrum3, 1000))'''