import collections
from Bio import SeqIO
import sys, time, math
import subprocess as sp

list1, list2, list3, list3a, final_list, new_list, filename = [], [], [], [], [], [], ''
# list1 contains the sequences names, list2 has the sequences, list3 has all string of letters at positions
# final_list contains the amino acid sequence with the most common letters that appear throughout the database

def getting filename:
    # Reads in the name of the fasta file in the command line if mentioned, otherwise 'KenyaHA2016-2017.fasta' is the default
    if len(sys.argv) > 1:
        filename = open(sys.argv[1])
    else:
        print('Did not receive a file name, try again!')

# Method is defined that reads through the FASTA file using the Biopython module and extracts the sequence and the name
# of the sequence into 2 separate lists
def getting_aa(name, list_one, list_two, list_three):
    i = 0
    fasta = list(SeqIO.parse(name, "fasta"))
    for each in fasta:    # The loop stores the sequence and the sequence id into 2 lists
        ids = fasta[i].id   # Takes the ID of the sequence
        seq = fasta[i].seq  # Takes the sequence in
        description = fasta[i].description

        # Taking in just the name of the Strain and its isoltaion year for future reference
        start = description.find('Strain Name:') + 12
        end = description.find('Protein Name:') - 1
        name = description[start:end]

        list_one.append(name)
        list_two.append(seq)
        list_three.append(ids[:ids.index('|')])
        i = i+1

# Method to extract individual letters of the sequence and them make them into a string so that they can be compared later
def finding_letters(list2, list3):
    j = 0
    while j < len(list2[0]):    # the while loop to read through each sequence
        k = 0
        letter = ''
        while k < len(list2):   # the while loop to read through each of letter of the sequence
            letter = list2[k][j]+letter # a string that contains letters of all sequences at a certain position
            k += 1
        j += 1
        list3.append(letter)    # Adds the collection of letters into a list

# Method to find out the most common letter at each position of the string in the list
def finding_frequency(list3, final_list):
    m = 0
    while m < len(list3):   # loop that reads through the list
        string = list3[m]
        most_common = (collections.Counter(string).most_common(1)[0])   # counter to find out the most frequent letter
        final_list.append(most_common[0])     # extracts and adds the most common letter to the final_list
        m += 1

# Method that compares all of the database sequences to the predicted sequence
def finding_sequence(list_of_sequences, list_two, list_one, list_three, list_four):
    count = 0
    counter = 0
    for sequences in list_of_sequences:     # reads through all the sequences one by one
        i = 0
        while i < len(sequences):   # While loop reads through each letter of the sequence
            if sequences[i] == list_two[i]:     # A count increases when there is a match of sequences
                counter += 1
            count += 1
            i += 1
        result = (float(counter)/float(count))*100  # Equation calculates the percentage of amino acids that are matches
        list_three.append(result)   # The percentages are added to a list for future reference

    for items in list_three:    # For loop to find out which sequence alignment has the best match
        max = items
        if items >= max:
            max = items     # if loop constantly updates the max found within the list

    index = list_three.index(max)
    new_num = ("%.2f" % round(list_three[index],2))     # rounds the decimal number to 2 decimal places
    print('The highest sequence identity match was found for %s (%s)' % (list_one[index], list_four[index]))
    print('The sequence alignment was found to be at %s (%s/%s)' % (new_num, counter, count))

    print("\n{:<22}|{:>10}".format("Influenza", " Identity Match")) # column widths: 12, 11 & 6
    newList = sorted(list_three, reverse=True)
    time.sleep(0.3)
    z = 0
    for sequences in list_one:
        index2 = list_three.index(newList[z])
        print("{:<22}|  {:<10}".format(list_one[index2], ("%.2f" % round(list_three[index2],2))))
        z += 1

# Method that saves the newly created sequence into a FASTA file so that it can be used on blast
# *** CURRENTLY INCOMPLETE
def save_sequence(string):
    output = open('sample.fasta', 'w')#opening a file called output_text to be able to write to it
    # output.write(">gb:00000000|ncbiId:NA|UniProtKB:NA|Organism:Influenza A virus (Predicted)|Strain Name:NA|Protein Name:HA Hemagglutinin|Gene Symbol:HA|Segment:NA|Subtype:H1N1|Host:Human \n")
    z = 1
    newString = ''
    while z < len(string)+1:
        if (z%70 != 0):
            newString = newString+string[z-1]
        else:
            newString = newString+string[z-1]+'\n'
        z += 1
    print('\n'+newString+'\n')
    output.write(newString)

def run_blast(subject, query, datatype, outputfile):
	sp.Popen(['/usr/local/ncbi-blast-2.7.1+/bin/blastp', '-db', subject, '-query', 'sample.fasta', '-out', outputfile])

while i < 5:
# Calling of all the methods mentioned
getting_aa(filename, list1, list2, list3a)
finding_letters(list2, list3)
time.sleep(0.3)
finding_frequency(list3, final_list)

# The final_list is converted into a string, so that it can be printed as one
produced_sequence = ''.join(final_list)
# print ('\n%s\n'% produced_sequence)     # prints the created sequence to be seen

save_sequence(produced_sequence)
time.sleep(0.3)
finding_sequence(list2, final_list, list1, new_list, list3a)

# time.sleep(1)
# run_blast(filename, '-query', 'prot', 'comparison.txt')
