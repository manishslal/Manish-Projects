__author__ = 'NR'
'''
Perform permutation test on AML/ALL gene expression data
'''

import csv
from math import sqrt
from random import shuffle

########################
# Function Definitions #
########################

def get_diagnosis_columns(diag, header):
    '''get list of indices corresponding to given diagnosis in header line
    '''
    return [index for index, value in enumerate(header) if diag == value]

def get_labels_from_file(infile):
    '''Get diagnoses from header line in csv file
    '''
    with open(infile) as f:
        header_line = next(f)  # Get the first line in the file
        header_line = header_line.strip()  # line contains white characters that need to be stripped first
        diagnoses = header_line.split(',')  # contains diagnosis labels.
        return diagnoses

def get_expression_from_csv(infile):
    '''get expression data from csv. **Uses csv packages**. Returns a 2D list'''
    with open(infile) as f:
        next(f)  # skip the first line
        lines = [line for line in csv.reader(f, quoting=csv.QUOTE_NONNUMERIC,
                                             quotechar='"')]  # read in data using csv package
    return lines

def get_diagnosis_columns(diag, header):
    '''get list of indices corresponding to given diagnosis in header line'''
    return [index for index, value in enumerate(header) if diag == value]

def mean(values):
    ''' Calculate the mean of a sequence of numbers'''
    return float(sum(values)) / len(values)

def sample_variance(values):
    ''' Calculate the sample variance'''
    sample_mean = mean(values)
    return sum([(value - sample_mean) ** 2 for value in values]) / (len(values) - 1)

def get_test_statistic(x, sample1, sample2):
    '''Calculates a two sample test statistic. Assumes unequal sample variances.
    Sample1 and Sample2 are the indices that correspond each sample'''
    x1 = [x[i] for i in sample1]
    x2 = [x[j] for j in sample2]
    difference_in_means = abs(mean(x1) - mean(x2))
    denominator = sqrt((sample_variance(x1) / len(x1)) + (sample_variance(x2) / len(x2)))
    return difference_in_means / denominator

def shuffle_diagnosis_labels(diagnosis_list):
    '''shuffle a list of diagnosis lists. Returns a list of the same
    length as the input list'''
    # NOTE: shuffle works in place and returns None. Because of this,
    # you have to duplicate the diagnosis list then shuffle.
    new_list = diagnosis_list[:]  # duplicate diagnosis list. Why do you think there's a slice operator in the expression?
    shuffle(new_list)  # shuffle the new list
    return new_list

################
# MAIN PROGRAM #
################

# read in data and get diagnosis indices
true_labels =  get_labels_from_file('leukemia_big.csv')
genes = get_expression_from_csv('leukemia_big.csv')
ALL = get_diagnosis_columns('ALL', true_labels)
AML = get_diagnosis_columns('AML', true_labels)

x = genes[0]
print get_test_statistic(genes[0], ALL, AML)
i = 0
while i < 1000:
    shuffled_labels = shuffle_diagnosis_labels(true_labels)
    print(get_test_statistic(genes[0], shuffled_labels, true_labels))
    #print(get_test_statistic)
    i += 0
print true_labels
print shuffled_labels
