"""
Created on Tue Jan 30:
@author: Nasserallah Manai, S5700698, University of Groningen, The Netherlands
@email: n.manai@student.rug.nl

This python script is used to obtain statistical information such as: the number of unique haplotypes / the number of
segregation sites and their positions from a dataset in fasta format.  In our case, we will use it for a dataset named
bp_seq.fasta which contains genetic sequences of green sea turtles.  we create multiple functions
to facilitate our calculations and information extraction.
    :return txt file containing a summary table of everything calculated in our script + a list of all of our unique
     haplotypes with unique IDs/ number of occurrences in the main dataset/ Sequence


"""

# imports
from Bio import SeqIO
from collections import Counter


# functions


def nucleotides_diff(sequence_one, sequence_two: str):
    """
    This function takes two sequences of same length (395) and compares each nucleotide while returning a list
    containing the positions where differences were observed.
    :param sequence_one: First genetic sequence of length 395
    :param sequence_two: Second genetic sequence of length 395
    :return: a list of positions where differences were observed
    """
    seq_positions = list()
    for position in range(395):
        if sequence_one[position] != sequence_two[position]:
            seq_positions.append(position + 1)
    return seq_positions


def comparer(hap_list: list):
    """
    This function takes our list of haplotypes (unique DNA sequences observed) and compares them one by one to each
    other by calling on the previously created nucleotide_diff function.  The final result is a list of lists of
    different positions for each individual comparison made.
    :param hap_list: a list containing the DNA sequences of our unique haplotypes
    :return: A list of lists containing the positions where nucleotide differences where observed for each individual
    comparison
    """
    holder = list()

    for i in range(len(hap_id) - 1):  # I ran out of names for counters and am using typical names I used to use in java
        for j in range(i + 1, len(hap_list)):  # and C++ :')
            holder.append(nucleotides_diff(hap_list[i][3], hap_list[j][3]))

    return holder


def clean(pos_list: list):
    """
    This function takes the output from the comparer function and cleans it by first putting every element of all the
    sub-lists in one list, then removing duplicates and sorta in ascending order. The output is a final list with all
    the positions where differences were observed. its length is the number of segregating sites

    :param pos_list: the list returned from our previous function "comparer".
    :return: a final list with all the positions where differences were observed. its length is the number of
    segregation sites
    """
    return_list = list()
    for sublist in pos_list:
        for position in sublist:
            return_list.append(position)

    return_list = list(dict.fromkeys(return_list))
    return_list.sort()

    return return_list


def pi_calc(compared: list):
    """
    This function takes the lists "seq_pos" from comparer and appends the length of each sublist to a new list. this
    list will then be used in function nuc_diversity to calculate the nucleotide diversity.

    :param compared: the list returned from our previous function "comparer".
    :return: a list of the number of differences observed in each comparison made.
    """

    product = list()
    for element in compared:
        product.append(len(element) / 395)

    return product


def nuc_diversity(list1, list2: list):
    """
    This function uses two lists as input to calculate the nucleotide diversity of our data. it does so by extracting
    haplotype frequency numbers from a list1 and plugging them into the equation (can be found here):
    https://en.wikipedia.org/wiki/Nucleotide_diversity
    :param list1: This list is one containing the haplotype frequencies (hap_id in our case), as these frequencies are
    necessary to calculate the nucleotide diversity.
    :param list2: list containing all the Pi values created by the pi_calc function
    :return: The nucleotide diversity
    """
    holder2 = list()
    final = list()

    for i in range(len(list1) - 1):
        for j in range(i + 1, len(list1)):
            holder2.append(list1[i][2] * list1[j][2])

    for k in range(len(list2)):
        final.append(holder2[k] * list2[k])

    tot = sum(final)

    diversity1 = (395 / 394) * tot

    return diversity1


def hap_diversity(lista: list):
    """
    This function calculates the haplotype diversity of our data. it does so by extracting haplotype frequency numbers
    from our inputted list (lista) and plugging them into the equation (can be found here):
    https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13643 (section 2.1 (3))
    :param lista: This list is one containing the haplotype frequencies (hap_id in our case), as these frequencies are
    necessary to calculate the haplotype diversity.
    :return: Haplotype diversity
    """
    temp_list = list()
    for counter in range(len(lista)):
        temp_list.append(lista[counter][2] ** 2)

    # Uncomment to check
    # print(temp_list)
    summation = sum(temp_list)

    diversity2 = (395 / 394) * (1 - summation)

    return diversity2


# main script


if __name__ == "__main__":

    initial_dict = {}

    with open('bp_seqs.fas', "r") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            initial_dict[seq_record.description] = str(seq_record.seq)

    """
    The block of code from 147-168 shows how I used the data that I extracted from the .fas file into the dictionary
    initial_dict, to create a list containing all the unique haplotypes (i.e their DNA sequences), an ID I gave them, 
    the number of times each one was observed in the dataset, and the frequency of each.
    """

    valuesList = []
    for key, values in initial_dict.items():
        valuesList.append(values)
        freq = Counter(valuesList)
        uniqueValues = list(freq.keys())
        uniqueValues.sort()

    total_entries = len(initial_dict)
    print(f'Total number of entries in the fasta file are: {total_entries}')

    temp_1 = dict(freq)
    temp_2 = sorted(temp_1.items(), key=lambda t: (-t[1], t[0]))
    haplos = dict(temp_2)

    hap_id = list()
    id = 1

    for key, values in haplos.items():
        hap_id.append(('HAP' + str(id), values, (values / total_entries), key))
        id += 1

    print(hap_id)

    # The next section is where I call on my created functions to do all the work

    hap_num = len(hap_id)
    print(f'The number of unique sequences in this data set is: {hap_num}')

    seq_pos = comparer(hap_id)

    pi_list = pi_calc(seq_pos)

    # next four commented print statements are to check (uncomment to check)
    # print(seq_pos)
    # print(len(seq_pos))

    # print(pi_list)
    # print(len(pi_list))

    nucleotide_diversity_temp = nuc_diversity(hap_id, pi_list)
    nucleotide_diversity = round(nucleotide_diversity_temp, 4)
    print(f'Our nucleotide diversity is : {nucleotide_diversity}')

    haplotype_diversity_temp = hap_diversity(hap_id)
    haplotype_diversity = round(haplotype_diversity_temp, 4)
    print(f'Our haplotype diversity is: {haplotype_diversity}')

    cleaned_list = clean(seq_pos)
    seg_num = len(cleaned_list)

    # Uncomment print statements to check
    # print(cleaned_list)

    print(f'The number of segregating sites is: {seg_num}')

    # This section prints out our final output (the table) to a .txt file named final_table.txt

    f = open("final_table.txt", "a")

    f.write("Number_of _Unique_Haplotypes\tNumber_of_Segregating_Sites\tNucleotide_Diversity\tHaplotype_Diversity\n")
    f.write(f'{hap_num}\t{len(cleaned_list)}\t{nucleotide_diversity}\t{haplotype_diversity}\n\n')

    f.write(f'ID\t#_of_Occurrences\tFrequency(fractional)\tSequence\n')
    for hap in range(hap_num):
        f.write(f'{hap_id[hap][0]}\t{hap_id[hap][1]}\t{round(hap_id[hap][2], 4)}\t{hap_id[hap][3]}\n')
    f.close()

    print("All Done!")

# A messy and not fully functioning alternative to find segregation sites just for fun! This part gave me the correct
# number of segregating sites,but the positions of the differences were off and wonky, even tho I was using the index
# method/function.


"""
haplo_seqs = list(haplos.keys())


nucleotide_comp = list()

for seq in haplo_seqs:
    for letter in seq:
         nucleotide_comp.append(letter)


#print (len(nucleotide_comp))

seg_site = [0] * 395

for counter in range(395):
    for count in range(1, 24):
        if (nucleotide_comp[counter] == nucleotide_comp[counter + (395 * count)]):
            seg_site[counter] += 0
        else:
            seg_site[counter] += 1


print(seg_site)

counter = 0
seg_index = list()

for nucleotide in seg_site:
    if nucleotide > 0:
        counter += 1
        seg_index.append(seg_site.index(nucleotide) + 1)

print(counter)
print(seg_index)
"""