import operator as op
import argparse
import os
import sys
import numpy as np
import pandas as pd
import math as mt
import matplotlib.pyplot as plt
import csv

# Hier is waar ek dit n command line tool maak :)
"""
parser = argparse.ArgumentParser(prog='pHMMdi', usage='%(prog)s [options] path',
                                 description='pHMMdi is a profile Hidden Markov Modeller that provides the user with '
                                             'a diagramatic representation of the results generated from the Forward '
                                             'and Viterbi algorithms',
                                 epilog='If we knew what it was we were doing,it would not be called research, would it? - Albert Einstein')
parser.add_argument('Path', metavar='path', type=str, help='the path to program file')
parser.add_argument('MSA', metavar='Fasta_file', type=str, action='store',
                    help='Make sure the file is in the correct format otherwise it really wont work...')
parser.add_argument('Datatype', metavar='datatype', type=str, action='store', help='Either P or N')
parser.add_argument('Query', metavar='Query_sequence', type=str, action='store',
                    help='Make sure the length of this sequence matches the sequence length in the fasta files')
parser.add_argument('-o', '--output', action='store', help='output file name containing generated data in table format (.csv file)')
parser.add_argument('-p', '--pseudocounts', action='store', help='pseudocounts method used (N, 0.1 or Pascal)')
parser.add_argument('-t', '--threshold', action='store', help='threshold (to calculate conserved regions from)')
parser.add_argument('-f', '--format', action='store', help='specify a format for the diagram output (svg or png)')
args = parser.parse_args()
user_specified_path = args.Path
user_specified_file = args.MSA
user_specified_datatype = args.Datatype
user_specified_query = args.Query
user_specified_settings = vars(args)
#print(user_specified_settings)

fasta_file = user_specified_settings['MSA']
sample_type = user_specified_settings['Datatype']
if sample_type == "N":
    pass
elif sample_type == "P":
    pass
else:
    sys.exit("This wasn't one of the available options, please try again")
query_sequence = user_specified_settings['Query']  # full disclaimer, this program only works with query sequences of 
                                                   # length equal or shorter than sequence length of sequences in fasta 
if user_specified_settings['pseudocounts'] != None:
    pseudo = user_specified_settings['pseudocounts']
else:
    pseudo = 'Pascal'
if user_specified_settings['threshold']:
    thres = user_specified_settings['threshold']
    if int(thres) > 100:
        sys.exit("Please pick a value between 0 and 100")
    if int(thres) <= 0:
        sys.exit("Please pick a number between 0 and 100")
else:
    thres = 50
if user_specified_settings['output'] != None:
    file_name = user_specified_settings['output']
else:
    file_name = 'pHMM_tables.csv'
if user_specified_settings['format'] != None:
    format = user_specified_settings['format']
else:
    format = 'svg'
"""
fasta_file = 'C:/Users/anika/Documents/Bioinformatics Modules and Assignments/HMM_projek/nmsa2.fasta'
sample_type = 'N'
query_sequence = '-A-ATGA'
pseudo = 'N'
thres = 50
file_name = 'pHMM_tables.csv'

class ReadAlignment:

    def __init__(self):
        pass

    def read_file(self, filepath):
        f = open(filepath,'r')
        name = []
        sequences = []
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                name.append(line[1:].rstrip('\n'))
                if len(current_sequence) > 0:
                    sequences.append(current_sequence)
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        sequences.append(current_sequence)
        amount_of_seq = len(sequences)
        if amount_of_seq > 500000:
            sys.exit("This might take a while, are you sure you need to check so many sequences...")
        len1 = len(sequences[0])
        if len1 <= 1:
            sys.exit("Maybe check your fasta file... It seems it doesn't contain sequences...")
        if len1 > 1000000:
            print('This might take some time... Are you sure you need an alignment this long?')
            print('You can restart with shorter sequences at any time...')
        counter = 0
        for seq in sequences:
            if (len(seq) == len1):
                counter += 1
            else:
                counter = counter
        if (counter == amount_of_seq):
            print("Thanks, all the sequences are the same length. Step one is completed correctly :) You're doing great!")
            return (name, sequences, amount_of_seq)
        else:
            sys.exit("This is not gonna work... Let's try again. Maybe check your sequence length :)")

    def checker(self):
        if sample_type == 'N':
            if 'J' in sequences or 'O' in sequences or 'U' in sequences:
                sys.exit("Check sample type")
            else:
                pass
        elif sample_type == 'P':
            if 'U' in sequences:
                sys.exit("Check sample type")
            else:
                pass
        else:
            pass

class EmissionsTable:

    def __init__(self):
        pass

    def split(self):
        # Function 1: Splitting the string of sequences into a shaped array
        sequences2 = np.asarray(sequences)
        tot_len = len(sequences2)
        seq_len = len(sequences2[0])
        counter = 0
        nuc = []
        for item in sequences2:
            list_1 = list(sequences2[counter])
            counter += 1
            nuc += list_1
            list_1 = []
        splitted = np.array(nuc, dtype='str').reshape(tot_len, seq_len)
        # Appending the array of splitted string to a list
        new = []
        counter2 = 0
        lys = []
        for array in range(len(splitted[counter2])):
            for j in range(len(splitted)):
                seq = splitted[j]
                new.append(seq[array])
            lys.append(new)
            new = []
        # for loop to split into columns
        length1 = len(lys)
        str1 = []
        counter = 0
        for item in range(len(lys)):
            split = lys[counter]
            counter += 1
            for letter in split:
                str1.append(letter)
        str2 = np.array_split(str1, length1)  # split the separated and sorted characters into groups
        str2 = np.array(str2)  # converting the list to an array to see whether the count function works now
        # print("str2:")
        # print(str2)
        return str2

    def check_type(self, str2):
        # Function 2: Determining whether we are working with DNA or RNA or protein
        origin = []
        sub = []
        for x in str2:
            if sample_type == 'P':
                a = op.countOf(x, "A")
                sub.append(a)
                b = op.countOf(x, "B")
                sub.append(b)
                c = op.countOf(x, "C")
                sub.append(c)
                d = op.countOf(x, "D")
                sub.append(d)
                e = op.countOf(x, "E")
                sub.append(e)
                f = op.countOf(x, "F")
                sub.append(f)
                g = op.countOf(x, "G")
                sub.append(g)
                h = op.countOf(x, "H")
                sub.append(h)
                i = op.countOf(x, "I")
                sub.append(i)
                k = op.countOf(x, "K")
                sub.append(k)
                l = op.countOf(x, "L")
                sub.append(l)
                m = op.countOf(x, "M")
                sub.append(m)
                n = op.countOf(x, "N")
                sub.append(n)
                p = op.countOf(x, "P")
                sub.append(p)
                q = op.countOf(x, "Q")
                sub.append(q)
                r = op.countOf(x, "R")
                sub.append(r)
                s = op.countOf(x, "S")
                sub.append(s)
                t = op.countOf(x, "T")
                sub.append(t)
                v = op.countOf(x, "V")
                sub.append(v)
                w = op.countOf(x, "W")
                sub.append(w)
                y = op.countOf(x, "Y")
                sub.append(y)
                z = op.countOf(x, "Z")
                sub.append(z)
            elif sample_type == 'N':
                element_exist = "U" in x
                if element_exist == True:  # for RNA
                    a = op.countOf(x, "A")
                    sub.append(a)
                    u = op.countOf(x, "U")
                    sub.append(u)
                    g = op.countOf(x, "G")
                    sub.append(g)
                    c = op.countOf(x, "C")
                    sub.append(c)
                else:  # for DNA # aware when AA's are added this could become problematic
                    a = op.countOf(x, "A")
                    sub.append(a)
                    t = op.countOf(x, "T")
                    sub.append(t)
                    g = op.countOf(x, "G")
                    sub.append(g)
                    c = op.countOf(x, "C")
                    sub.append(c)
            else:
                print("This wasn't one of the available options, please try again :)")
            origin.append(sub)
            sub = []  # clearing before continuing
        return origin

    def conserved(self, origin):
        # Function 3: To determine conserved regions
        conserved = []
        index = 0
        for o in origin:
            con = sum(origin[index])
            conserved.append(con)
            con = 0
            index += 1
        index = 0
        match = []
        insertion = []
        threshold = (int(thres) / 100) * len(sequences)
        collen = len(origin[0])
        for j in conserved:
            if j >= threshold:
                match.append(index)
            else:
                insertion.append(index)
            index += 1
        #print("conserved e:", match, "inserts:", insertion)
        index = 0
        return (insertion, threshold)

    def binary_table(self,str2):
        mini_tab = []
        table = []
        for element in str2:
            for letter in element:
                if letter == '-':
                    mini_tab.append(0)
                else:
                    mini_tab.append(1)
            table.append(mini_tab)
            mini_tab = []
        bin_table = np.asarray(table)
        df_table = pd.DataFrame(bin_table)
        # print("df_table")
        # print(df_table)
        return df_table

    def ranges(self, insertion):
        # function to determines whether the insertion column should become an insertion region
        insertion = sorted(set(insertion))
        gaps = [[s, e] for s, e in zip(insertion, insertion[1:]) if s+1 < e]
        edges = iter(insertion[:1] + sum(gaps, []) + insertion[-1:])
        insert_reg = list(zip(edges, edges))
        return insert_reg

    def counts(self, origin, insert_reg):

        # Function 4: Creating a dataframe (counts table)
        origin2 = np.transpose(origin)
        if sample_type == 'N':
            rows = ['A', 'T', 'G', 'C']
        else:  # sample is protein
            rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Z']
        df = pd.DataFrame(origin2, index=rows)
        # c_i header list original (not the combined one)
        c_i_list = []
        # get list of numerics for insert regions
        insert = []
        for tuple in insert_reg:
            for t in range(tuple[0], tuple[1]+1):
                insert.append(t)
        for i in range(len(origin)):
            if i in insert:
                c_i_list.append("i")
            else:
                c_i_list.append("c")
        # df.columns = c_i_list
        # Creating a second df that is not a counts table but the emissions table
        df2 = df[:]
        df.columns = c_i_list
        c_l = []
        for k in insert_reg:
            col_list = []
            for i in range(k[0], k[1] + 1):  # right so difference between values for emission and transition is
                # that for emission I need an emission probs for the insert regions, that will mean that emissions
                # table will have more columns than transition table, so please don't freak out when you see the output
                col_list.append(i)
            c_l.append(col_list)
        for l in range(len(c_l)):
            sub = c_l[l]
            df2[sub[0]] = df2[c_l[l]].sum(axis=1)
            # loop through remaining insertion columns to drop columns
        for l in range(len(c_l)):
            subb = c_l[l]
            for i in range(subb[0] + 1, subb[len(subb) - 1] + 1):
                # print(i, end="index to be dropped")
                del df2[i]
        # write df and df2 to file
        with open('dfs_Original.csv', 'w') as f:
                f.write('df')
                f.write("\n")
                df.to_csv(f)
        with open('dfs_Concatenated.csv', 'w') as f:
                f.write('df2')
                f.write("\n")
                df2.to_csv(f)
        return(df2)

    def emission(self, df2):

        # Function 5: Determining the nucleotide frequencies
        totfreq1 = df2.sum(axis=1)
        totfreq2 = totfreq1.sum()
        totfreq3 = []
        for i in range(len(totfreq1)):
            totfreq3.append(totfreq1[i]/totfreq2)
        # Setting up a frequencies table
        index = 0
        for k in df2.columns:
            coltot = sum(df2[k])
            # coltot = df[index].sum()
            df2[k] /= coltot
            index += 1
        freq_a = np.array(df2)
        pssm_pseudo = []
        pssm_log = []
        pssm_pseudo_tot = []
        log_tot = []
        number_of_seq = len(freq_a[0])
        for row in freq_a:
            if pseudo == 'N':
                index1 = 0
                ps1 = np.sqrt(number_of_seq)
                counts = np.sum(freq_a, axis=0)
                for x in range(len(freq_a[index1])):
                    index2 = 0
                    numerator = row[x] + (ps1*totfreq3[index2])
                    denominator = (counts[index2] + ps1)
                    new_val = numerator/denominator
                    pssm_pseudo.append(new_val)
                    log = -(mt.log10(new_val))
                    pssm_log.append(log)
                    index2 += 1
                pssm_pseudo_tot.append(pssm_pseudo)
                log_tot.append(pssm_log)
                pssm_log = []
                pssm_pseudo = []
            elif pseudo == '0.1':
                index1 = 0
                ps1 = 0.1
                counts = np.sum(freq_a, axis=0)
                for x in range(len(freq_a[index1])):
                    index2 = 0
                    #print(freq_a[index1])
                    #print(counts[index2])
                    numerator = row[x] + (ps1*totfreq3[index2])
                    denominator = (counts[index2] + ps1)
                    new_val = numerator/denominator
                    pssm_pseudo.append(new_val)
                    log = -(mt.log10(new_val))
                    pssm_log.append(log)
                    index2 += 1
                pssm_pseudo_tot.append(pssm_pseudo)
                log_tot.append(pssm_log)
                pssm_log = []
                pssm_pseudo = []
            elif pseudo == 'Pascal':
                index1 = 0
                ps1 = 1
                ps2 = number_of_seq
                counts = np.sum(freq_a, axis=0)
                for x in range(len(freq_a[index1])):
                    index2 = 0
                    numerator = row[x] + (ps1)
                    denominator = (counts[index2] + ps2)
                    new_val = numerator / denominator
                    pssm_pseudo.append(new_val)
                    log = -(mt.log10(new_val))
                    pssm_log.append(log)
                    index2 += 1
                pssm_pseudo_tot.append(pssm_pseudo)
                log_tot.append(pssm_log)
                pssm_log = []
                pssm_pseudo = []
            else:
                sys.exit('Unacceptable, check the pseudocount you selected')
            index1 += 1
        print("Emissions probabilities table created :)")
        return(pssm_pseudo_tot, log_tot)


class TransitionTable:

    def __init__(self):
        pass

    def bin_table(self):
        a = EmissionsTable.split(sequences)  # importing from emissions_table.py to get my sequences
        # creating a binary table to work with (-) are (0's) and letters are (1)
        mini_tab = []
        table = []
        for element in a:
            for letter in element:
                if letter == '-':
                    mini_tab.append(0)
                else:
                    mini_tab.append(1)
            table.append(mini_tab)
            mini_tab = []
        binary_table = np.asarray(table)
        df_table = pd.DataFrame(binary_table)
        df_table = df_table.transpose()
        return(binary_table, df_table)
        # Okayyy so this turns my data into a dataframe consisting of 0's and 1's to represent dashes and letters

# So now I neeeeed to determine conserved regions
    def concerved_region(self, binary_table):
        conserved = []
        index = 0
        for o in binary_table:
            con = sum(binary_table[index])
            conserved.append(con)
            con = 0
            index += 1
        index = 0
        match = []
        insertion = []
        collen = len(binary_table[0])
        threshold = (int(thres) / 100) * len(sequences)
        for j in conserved:
            if j >= threshold:
                match.append(index)
            else:
                insertion.append(index)
            index += 1
        #print("conserved t:", match, "inserts", insertion)
        return(match, insertion)

# So this creates an array containing the current states of the alignment
    def current_states(self, match, df_table):
        array_row = []
        new_array = []
        i = 0
        j = 0
        # print(df_table)
        # print(match)
        for col in df_table:
            # print(df_table.columns[i])
            # print(match[j])
            if df_table.columns[i] == match[j]:
                for nuc in df_table[i]:
                    if nuc == 1:
                        array_row.append("M")
                    else:
                        array_row.append("D")
                i += 1
                j += 1
                new_array.append(array_row)
                array_row = []
            else:
                for nuc in df_table[i]:
                    if nuc == 1:
                        array_row.append("I")
                    else:
                        array_row.append('-')
                i += 1
                new_array.append(array_row)
                array_row = []
        #print("new_array aka current states:", new_array)
        return(new_array)   #Array containing current states

    def transition_probs(self, new_array):
        # Okay so what I want to do now is to determine transition probabilities using the array
        # But first lets create some states
        mm = 0
        md = 0
        mi = 0
        dd = 0
        dm = 0
        di = 0
        ii = 0
        im = 0
        id = 0
        nuc = 0  # counter to loop through nucleotides in alignment
        temp_next = []  # temp list to keep values when running through insertion transitions
        total_iterations = []  # list to append the iterations to
        for s in range(len(new_array)-1):  # loop through all the alignments one by one
            previous = temp_next
            temp_next = []
            current = new_array[s]  # stores current alignment
            next = new_array[s+1]  # stores next alignment
            if 'I' in previous:  # meaning there was an insertion state and a temp list was created
                current = previous
                if '-' in next:  # if it goes to another insertion
                    for a in current:
                        if current[nuc] == 'M' and 'I' in next[nuc]:
                            mi += 1
                            temp_next.append('I')
                        elif current[nuc] == 'M' and next[nuc] == '-':
                            temp_next.append('M')
                        elif current[nuc] == 'D' and next[nuc] == 'I':
                            di += 1
                            temp_next.append('I')
                        elif current[nuc] == 'D' and next[nuc] == '-':
                            temp_next.append('D')
                        elif current[nuc] == 'I' and next[nuc] == 'I':
                            ii += 1
                            temp_next.append('I')
                        else:  # if current == 'I' AND next == '-'
                            temp_next.append('I')
                        nuc += 1
                else:  # if it goes to a conserved region
                    for a in current:
                        if current[nuc] == 'M' and next[nuc] == 'M':
                            mm += 1
                        elif current[nuc] == 'M' and next[nuc] == 'D':
                            md += 1
                        elif current[nuc] == 'D' and next[nuc] == 'D':
                            dd += 1
                        elif current[nuc] == 'D' and next[nuc] == 'M':
                            dm += 1
                        elif current[nuc] == 'I' and next[nuc] == 'M':
                            im += 1
                        else:  # if current == 'I' AND next == 'D'
                            id += 1
                        nuc += 1
            else:  # hopefully meaning if previous is empty
                if 'M' in current:  # checking whether current state is conserved
                    if 'M' in next:  # checking whether next state is conserved
                        for a in current:  # running through all the nucleotides in current and next to compare and determine transition
                            if current[nuc] == next[nuc] == 'M':
                                mm += 1
                            elif current[nuc] == next[nuc] == 'D':
                                dd += 1
                            elif current[nuc] == 'M' and next[nuc] == 'D':
                                md += 1
                            else:
                                dm += 1
                            nuc += 1  # increment nucleotide counter to move on and compare the next two nucleotides
                    else:  # if next state is an insert
                        previous = []
                        for a in current:  # running through all the nucleotides in current and next to compare and determine transition
                            if current[nuc] == 'M' and 'I' in next[nuc]:
                                mi += 1
                                temp_next.append('I')  # appending to a temp_next list as complete transitions can't be calculated yet
                            elif current[nuc] == 'D' and 'I' in next[nuc]:
                                di += 1
                                temp_next.append('I')
                            elif current[nuc] == 'M' and next[nuc] == '-':
                                temp_next.append('M')
                            else:
                                temp_next.append('D')
                            nuc += 1
                else:  # meaning we are starting with an insert state
                    if 'M' in next:  # checking whether next state is conserved
                        for a in current:  # running through all the nucleotides in current and next to compare and determine transition
                            if current[nuc] == 'I' and next[nuc] == 'M':
                                im += 1
                            elif current[nuc] == 'I' and next[nuc] == 'D':
                                id += 1
                            elif current[nuc] == '-' and next[nuc] == 'D':  # if there is nothing in state i0 we revert back to m0
                                md += 1
                            else:
                                mm += 1
                            nuc += 1
                    else:  # if next state is an insert -- this should repeat until next state = 'M'
                        previous = []
                        for a in current:  # running through all the nucleotides in current and next to compare and determine transition
                            if current[nuc] == 'I' and 'I' in next[nuc]:
                                ii += 1
                                temp_next.append('I')  # appending to a temp_next list as complete transitions can't be calculated yet
                            elif current[nuc] == '-' and 'I' in next[nuc]:
                                mi += 1
                                temp_next.append('I')
                            elif current[nuc] == '-' and next[nuc] == '-':
                                temp_next.append('M')
                            else:
                                temp_next.append('D')
                            nuc += 1
            nuc = 0
            total_iterations.append([mm,md,mi,dd,dm,di,ii,im,id])
            mm = 0
            md = 0
            mi = 0
            dd = 0
            dm = 0
            di = 0
            ii = 0
            im = 0
            id = 0
        return(total_iterations)

    def concatenate_it(self, total_iterations):
        concatenated = []
        skip = False
        for sub in range(len(total_iterations)):
            subs = total_iterations[sub]
            subs_next = []
            if skip is not True:
                if sub < len(total_iterations)-1:
                    subs_next = total_iterations[sub+1]
                    if (subs_next[2] > 0 or subs_next[5] > 0 or subs_next[6] > 0 or subs_next[7] > 0 or subs_next[8] > 0) and sub > 0:
                        new = [a + b for a, b in zip(subs, subs_next)]
                        concatenated.append(new)
                        sub = sub + 1
                        skip = True
                    else:
                        concatenated.append(subs)
                else:
                    concatenated.append(subs)
            else:
                skip = False
        return(concatenated)

    def logg_odd_scores(self, concatenated):
        # nou gou om pseudocounts te add en uit te sorteer
        transition_array = []
        temp_transition_array = []
        ps1 = 1
        ps2 = 3  # this is what they used in literature (http://www.lcqb.upmc.fr/julianab/teaching/SBAS/2019/cours8.pdf)
        # This is hopefully Pascal's pseudocount method
        # print("state", concatenated)
        for state in concatenated:  # loop through all states
            index1 = 0  # position of transition in state
            index2 = 0  # position of state in array
            for transition in state:  # loop through every state in the transition maw mm until id
                current_transition = transition
                if index1 <= 2:  # for states mm, md, mi
                    new_tran = (current_transition+ ps1) / (concatenated[index2][0] + concatenated[index2][1] + concatenated[index2][2] + ps2)
                    index1 += 1
                    temp_transition_array.append(new_tran)
                    new_tran = []
                elif index1 > 2 and index1 <= 5:  # for states dd, dm, di
                    new_tran = (current_transition + ps1) / (concatenated[index2][3] + concatenated[index2][4] + concatenated[index2][5] + ps2)
                    index1 += 1
                    temp_transition_array.append(new_tran)
                    new_tran = []
                else:  # for states ii, id, im
                    new_tran = (current_transition + ps1) / (concatenated[index2][6] + concatenated[index2][7] + concatenated[index2][8] + ps2)
                    index1 += 1
                    temp_transition_array.append(new_tran)
                    new_tran = []
            index2 += 1
            transition_array.append(temp_transition_array)
            temp_transition_array = []
        print("Transition probabilities table created :)")
        return(transition_array)


class Algorithms:

    def __init__(self):
        self.transition = t_df3.to_numpy()
        self.emissions = pssm_p.to_numpy()
        # print("transition table")
        # print(self.transition)
        # print("emissions table")
        # print(self.emissions)
        pass

    def querysequence(self, insert_reg):
        seq_len = len(self.emissions)
        if sample_type == 'P':
            query = query_sequence
            query = query.upper()
            #print(query)
            if 'J' in query or 'O' in query or 'U' in query or 'X' in query:
                #print(query)
                sys.exit("Contains invalid character... maybe check your file?")  # elif check whether any special char except '-' is present OR any numbers
            else:
                print("No invalid characters identified in query sequence")
        elif sample_type == 'N':
            query = query_sequence
            query = query.upper()
            #print("query", query)
            if 'B' in query or 'D' in query or 'E' in query or 'F' in query or 'H' in query or 'I' in query \
                        or 'J' in query or 'K' in query or 'L' in query or 'M' in query or 'N' in query or 'O' in query \
                        or 'P' in query or 'Q' in query or 'R' in query or 'S' in query or 'V' in query or 'W' in query \
                        or 'X' in query or 'Y' in query or 'Z' in query:
                sys.exit("Contains invalid character... maybe check your file?")  # elif check whether any special char except '-' is present OR any numbers
            else:
                print("No invalid characters identified in query sequence")
        else:
            sys.exit("Uhmmm... just retype the query sequence...")
        s = 0
        i = 0
        region = 0
        bi_o = []
        for emission in query:
            if region < len(insert_reg):
                if i == insert_reg[region][0] or i >= insert_reg[region][0] and i <= insert_reg[region][1]:
                    if i > insert_reg[region][0]:
                        s -= 1
                if i == insert_reg[region][1]:
                    region += 1
            emission_o = query[i]  # the letter being omitted at a specific time
            i += 1
            if sample_type == 'N':
                if emission_o == 'A':
                    bi_o.append(self.emissions[0][s])
                elif emission_o == 'T':
                    bi_o.append(self.emissions[1][s])
                elif emission_o == 'U':
                    bi_o.append(self.emissions[1][s])
                elif emission_o == 'G':
                    bi_o.append(self.emissions[2][s])
                elif emission_o == 'C':
                    bi_o.append(self.emissions[3][s])
                else:  # if - is emitted
                    bi_o.append(1)  # not 0 otherwise my calculations won't work
            else:  # type is p
                if emission_o == 'A':
                    bi_o.append(self.emissions[s][0])
                elif emission_o == 'B':
                    bi_o.append(self.emissions[s][1])
                elif emission_o == 'C':
                    bi_o.append(self.emissions[s][2])
                elif emission_o == 'D':
                    bi_o.append(self.emissions[s][3])
                elif emission_o == 'E':
                    bi_o.append(self.emissions[s][4])
                elif emission_o == 'F':
                    bi_o.append(self.emissions[s][5])
                elif emission_o == 'G':
                    bi_o.append(self.emissions[s][6])
                elif emission_o == 'H':
                    bi_o.append(self.emissions[s][7])
                elif emission_o == 'I':
                    bi_o.append(self.emissions[s][8])
                elif emission_o == 'K':
                    bi_o.append(self.emissions[s][9])
                elif emission_o == 'L':
                    bi_o.append(self.emissions[s][10])
                elif emission_o == 'M':
                    bi_o.append(self.emissions[s][11])
                elif emission_o == 'N':
                    bi_o.append(self.emissions[s][12])
                elif emission_o == 'P':
                    bi_o.append(self.emissions[s][13])
                elif emission_o == 'Q':
                    bi_o.append(self.emissions[s][14])
                elif emission_o == 'R':
                    bi_o.append(self.emissions[s][15])
                elif emission_o == 'S':
                    bi_o.append(self.emissions[s][16])
                elif emission_o == 'T':
                    bi_o.append(self.emissions[s][17])
                elif emission_o == 'V':
                    bi_o.append(self.emissions[s][18])
                elif emission_o == 'W':
                    bi_o.append(self.emissions[s][19])
                elif emission_o == 'Y':
                    bi_o.append(self.emissions[s][20])
                elif emission_o == 'Z':
                    bi_o.append(self.emissions[s][21])
                else:  # emission_o == '-'
                    bi_o.append(1)  # to allow the calculation to continue
            s += 1
        #print("list with the emission probabilities for the specific sequence:", bi_o)
        # print("length of list with emission probabilities for the query sequence", len(bi_o))
        return (query, bi_o)

    def viterbi(self, query, bi_o):
        most_probable_path = []
        # print("t_df3")
        # print(self.transition)
        # initialization
        # probability_dist = 1
        state = 0
        em_prob = 0
        aMj_1Mj = self.transition[state][0]
        aMj_1Dj = self.transition[state][1]
        aMjIj = self.transition[state][2]
        aDj_1Dj = self.transition[state][3]
        aDj_1Mj = self.transition[state][4]
        aDjIj = self.transition[state][5]
        aIjIj = self.transition[state][6]
        aIj_1Mj = self.transition[state][7]
        aIj_1Dj = self.transition[state][8]
        log_values = []
        # initialization for real now
        VjM_i = 0
        log_values.append(VjM_i)
        VjI_i = 0
        VjD_i = 0
        most_probable_path.append(log_values)
        log_values = []
        max_state = []
        max_state.append((0,0))
        #print(len(self.transition))
        l = len(self.transition)
        # for the recursion part I want to loop through all the states except the first and last
        for st in range(0, l-1):
            # where the em_prob = em_count + ps/em_freq + ps
            # so just following my logic here:
            m = (-mt.log10(bi_o[em_prob])) + max((VjM_i + (-(mt.log10(aMj_1Mj)))),
                                                      (VjI_i + (-(mt.log10(aIj_1Mj)))),
                                                      (VjD_i + (-(mt.log10((aDj_1Mj))))))
            i = (-mt.log10(bi_o[em_prob])) + max((VjM_i + (-(mt.log10(aMjIj)))), (VjI_i + (-(mt.log10(aIjIj)))),
                                                      (VjD_i + (-(mt.log10((aDjIj))))))
            d = max((VjM_i + (-mt.log10(aMj_1Dj))), (VjI_i + (-(mt.log10(aIj_1Dj)))),
                        (VjD_i + (-(mt.log10((aDj_1Dj))))))
            log_values.append(d)
            log_values.append(m)
            log_values.append(i)
            VjM_i = m
            VjI_i = i
            VjD_i = d
            # find max for best path
            max_st = max(m, i, d)
            index_array = np.array([m, i, d])
            max_st_index = np.argmax(index_array)
            #print(index_array)
            max_state.append((max_st_index, max_st))
            em_prob += 1
            most_probable_path.append(log_values)
            log_values = []
            state += 1
            aMj_1Mj = self.transition[state][0]
            aMj_1Dj = self.transition[state][1]
            aMjIj = self.transition[state][2]
            aDj_1Dj = self.transition[state][3]
            aDj_1Mj = self.transition[state][4]
            aDjIj = self.transition[state][5]
            aIjIj = self.transition[state][6]
            aIj_1Mj = self.transition[state][7]
            aIj_1Dj = self.transition[state][8]
        # termination
        V_end = max((VjM_i + (-(mt.log10(aMj_1Mj)))), (VjI_i + (-(mt.log10(aIj_1Mj)))), (VjD_i + (-(mt.log10((aDj_1Mj))))))
        most_probable_path.append(V_end)
        max_state.append((0, V_end))
        pad = []
        for state in max_state:
            if state[0] == 0:
                pad.append("match")
            elif state[0] == 1:
                pad.append("insert")
            else:
                pad.append("deletion")
        print(pad)
        # print("most probable path:", most_probable_path)
        # print("length:", len(most_probable_path))
        # print("best path:", str(max_state))
        return (max_state, pad)

    def forward(self, query, bi_o):
        all_possibilities = []
        # initialization
        probability_dist = 1
        # recursion
        st = 0
        em_prob = 0
        tr_prob = 0
        tr_mm = self.transition[tr_prob][0]
        tr_md = self.transition[tr_prob][1]
        tr_mi = self.transition[tr_prob][2]
        tr_dd = self.transition[tr_prob][3]
        tr_dm = self.transition[tr_prob][4]
        tr_di = self.transition[tr_prob][5]
        tr_ii = self.transition[tr_prob][6]
        tr_im = self.transition[tr_prob][7]
        tr_id = self.transition[tr_prob][8]
        log_values = []
        # initialization for real now
        FjM_i = 0
        log_values.append(FjM_i)
        FjI_i = 0
        FjD_i = 0
        all_possibilities.append(log_values)
        log_values = []
        for state in range(1, (len(query) - 1)):  # for the recursion part I want to loop through all the states except the first and last
            st += 1
            # where the em_prob = em_count + ps/em_freq + ps
            m = -(mt.log10(bi_o[em_prob]) + -(mt.log10(
                    (tr_mm * mt.exp(FjM_i)) + (tr_im * mt.exp(FjI_i)) + (tr_dm * mt.exp(FjD_i)))))
            log_values.append(m)
            i = -(mt.log10(bi_o[em_prob]) + -(mt.log10(
                    (tr_mi * mt.exp(FjM_i)) + (tr_ii * mt.exp(FjI_i)) + (tr_di * mt.exp(FjD_i)))))
            log_values.append(i)
            d = -(mt.log10((bi_o[em_prob])) + -(mt.log10(
                    (tr_md * mt.exp(FjM_i)) + (tr_id * mt.exp(FjI_i)) + (tr_dd * mt.exp(FjD_i)))))
            log_values.append(d)
            FjI_i = i
            FjM_i = m
            FjD_i = d
            em_prob += 1
            all_possibilities.append(log_values)
            log_values = []
        # termination
        F_end = [-(mt.log10(bi_o[em_prob]) + (-mt.log10(
                    (tr_mm * mt.exp(FjM_i)) + (tr_im * mt.exp(FjI_i)) + (tr_dm * mt.exp(FjD_i)))))]
        all_possibilities.append(F_end)
        # print("all possibilities:", all_possibilities)
        # print("length:", len(all_possibilities))
        print("Both algorithms are completed, the csv file should be available in the output directory.")
        return (all_possibilities)


class Model:

    def __init__(self):
        print("Please be patient while the models are created.")
        self.scaler = len(match)/0.6
        self.y = [0, 2, 5]
        extend_y = len(match)/0.6
        self.y.append(extend_y)
        self.match2 = [-5]
        extention = len(match)*1.5
        for x in match:
            self.match2.append(x*3)
        extend_x = self.match2[-1] + extention
        self.match2.append(extend_x)
        pass

    def transition_visualizer(self, transition):
        colour = "black"
        figure, axes = plt.subplots(1, sharex=True, sharey=True)
        #gs = figure.add_gridspec(3, hspace=0)
        # initial state for m0 and i0
        width_mm = transition[0][0]
        width_mi = transition[0][2]
        width_md = transition[0][1]
        width_dd = transition[0][3]
        width_dm = transition[0][4]
        width_di = transition[0][5]
        width_ii = 10 + transition[0][6]
        width_im = transition[0][7]
        width_id = transition[0][8]
        t_size = 'x-small'
        if len(match) > 12:
            t_size = 'xx-small'
        rectangle = plt.Rectangle((2 - 1, 3), width=2, height=1, fc='None', ec="black", linewidth=1.5)
        state_m0 = plt.text(1.5, 3.2, 'm0', size='x-small')
        diamond = plt.Rectangle((2, 7), width=1.25, height=1.25, angle=45, facecolor="None", edgecolor="black", linewidth=1.25)
        state_i = plt.text(1.75, 7.5, 'i0', size='x-small')
        plt.arrow(2, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_im)  # im
        plt.arrow(3,3.5,1.5,0,head_width=0.5,head_length=0.5, fc =colour, ec=colour, linewidth=width_mm)  # mm
        plt.arrow(2, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_mi)  # mi
        plt.arrow(2, 8.7, 3, 2.2, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_id)  # id
        plt.arrow(2, 4, 3.5, 6.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_md)  # md
        axes.set_xlim(1,5)
        axes.set_ylim(1,5)
        axes.plot([1],[7.3],marker=r'$\circlearrowright$',ms=10, color=colour)
        plt.arrow(2, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_im)  # im
        colour = 'black'
        axes.add_artist(rectangle)
        axes.add_artist(diamond)
        # then for the middle part
        x = 6
        counter = 0
        print("amount of match states predicted:", len(match))
        #("forward:", len(transition))
        j = 1
        for i in range(1, len(match)-1):
            if j < len(transition):
                width_mm = transition[i][0]*3
                width_mi = transition[i][2]*3
                width_md = transition[i][1]*3
                width_dd = transition[i][3]*3
                width_dm = transition[i][4]*3
                width_di = transition[i][5]*3
                width_ii = 10 + transition[i][6]*3
                width_im = transition[i][7]*3
                width_id = transition[i][8]*3
            t_size = 'x-small'
            if len(match) > 12:
                t_size = 'xx-small'
            rectangle = plt.Rectangle((x-1, 3), width=2.2, height=1.1, fc='None', ec="black", linewidth=1.25)
            state_m = plt.text(x-0.5, 3.2, s=('m' + str(i)), size=t_size)
            diamond = plt.Rectangle((x, 7), width=1.35, height=1.35, angle=45, facecolor="None", edgecolor="black", linewidth=1.25)
            state_i = plt.text(x - 0.3, 7.5, s=('i' + str(i)), size=t_size)
            circle = plt.Circle((x, 12), radius=1.1, fc='None', ec='black', linewidth=1.25)
            state_d = plt.text(x-0.5, 11.75, s=('d' + str(i)), size=t_size)
            if i < (len(match)):
                counter += 1
                plt.arrow(x + 1, 3.5, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_mm)  # mm
                plt.arrow(x + 1, 12, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_dd)  #dd
                plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_mi)  # mi
                plt.arrow(x, 11, 0, -1.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_di)  # di
                plt.arrow(x, 8.7, 3, 2.2, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_id)  # id
                plt.arrow(x, 4, 3.5, 6.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_md)  # md
                plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_im)  # im
                plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_mi)  # mi
                plt.arrow(x, 11, 3.6, -6.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_dm)  # dm
                # axes.set_xlim(1, 135)
                # axes.set_ylim(1, 35)
                axes.plot([x - 1], [7.3], marker=r'$\circlearrowright$', ms=width_ii, color=colour)
                x += 4
                j += 1
            axes.set_xlim(0.5, len(self.match2))
            axes.set_ylim(0, self.scaler/2)
            axes.plot([x - 1], [7.3], marker=r'$\circlearrowright$', ms=width_ii, color='black')
            axes.set_aspect(1)
            axes.add_artist(rectangle)
            axes.add_artist(diamond)
            axes.add_artist(circle)
        t_size = 'x-small'
        if len(match) > 12:
            t_size = 'xx-small'
        rectangle = plt.Rectangle((x-1, 3), width=2.2, height=1.1, fc='None', ec="black", linewidth=1.25)
        state_m = plt.text(x-0.5, 3.2, s=('m' + str(i+1)), size=t_size)
        diamond = plt.Rectangle((x, 7), width=1.35, height=1.35, angle=45, facecolor="None", edgecolor="black", linewidth=1.25)
        state_i = plt.text(x - 0.3, 7.5, s=('i' + str(i+1)), size=t_size)
        circle = plt.Circle((x, 12), radius=1.1, fc='None', ec='black', linewidth=1.25)
        state_d = plt.text(x-0.5, 11.75, s=('d' + str(i+1)), size=t_size)
        axes.set_aspect(1)
        axes.add_artist(rectangle)
        axes.add_artist(diamond)
        axes.add_artist(circle)
        # to end it all off
        width_mm = transition[-2][0]*3
        width_mi = transition[-2][2]*3
        width_md = transition[-2][1]*3
        width_dd = transition[-2][3]*3
        width_dm = transition[-2][4]*3
        width_di = transition[-2][5]*3
        width_ii = 10 + transition[-2][6]*3
        width_im = transition[-2][7]*3
        width_id = transition[-2][8]*3
        plt.arrow(x + 1, 3.5, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_mm)  # mm
        plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_mi)  # mi lyne
        plt.arrow(x, 11, 0, -1.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_di)  # di
        plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_im)  # im
        plt.arrow(x, 11, 3.6, -6.5, head_width=0.5, head_length=0.5, fc=colour, ec=colour, linewidth=width_dm)  # dm
        x += 4
        rectangle = plt.Rectangle((x-1, 3), width=2, height=1, fc='None', ec="black", linewidth=1.5)
        state_m0 = plt.text(x-0.6, 3.2, 'END', size='x-small')
        axes.add_artist(rectangle)
        plt.xticks(self.match2)
        plt.yticks(self.y)
        if len(match) < 10:
            plt.xlim(0, 25)
            plt.ylim(1, 15)
        elif len(match) < 30:
            plt.xlim(0, 50)
            plt.ylim(1, 20)
        elif len(match) < 65:
            plt.xlim(0, 75)
            plt.ylim(1, 30)
        else:
            plt.xlim(0, 100)
            plt.ylim(1, 30)
        axes.set_title( "Profile HMM representing MSA" )
        figure.tight_layout()
        plt.axis('off')
        if format == 'png':
            plt.savefig("transitions_diagram.png")
        else:
            plt.savefig("transitions_diagram.svg")
        diagram = plt.show()
        return diagram

    def viterbi_vizualizer(self, pad):
        pad2 = []
        for state in pad:
            #print("state", state)
            if state[0] == 0:
                pad2.append('match')
            elif state[0] == 2:
                pad2.append('deletion')
            else:  # if insert
                pad2.append('match')
                pad2.append('insert')
        #print("pad2:", pad2)
        figure, axes = plt.subplots(1, sharex=True, sharey=True)
        # shapes
        current = 0
        next = 1
        x = 2
        for state in range(0, (len(pad2) - 1)):
            t_size = 'x-small'
            if len(match) > 12:
                t_size = 'xx-small'
            rectangle = plt.Rectangle((x - 1, 3), width=2, height=1, fc='None', ec="black", linewidth=1.5)
            state_m0 = plt.text(x - 0.5, 3.2, s=('m' + str(current)), size=t_size)
            axes.add_artist(rectangle)
            diamond = plt.Rectangle((x, 7), width=1.25, height=1.25, angle=45, facecolor="None", edgecolor="black",
                                    linewidth=1.25)
            state_i = plt.text(x - 0.3, 7.5, s=('i' + str(current)), size=t_size)
            axes.add_artist(diamond)
            if state == 0:
                pass
            else:
                circle = plt.Circle((x, 12), radius=1.1, fc='None', ec='black', linewidth=1.25)
                state_d = plt.text(x - 0.5, 11.75, s=('d' + str(current)), size=t_size)
                axes.add_artist(circle)
            x += 4
            next += 1
            current += 1
        # lines
        current = 0
        next = 1
        previous = -1
        x = 2
        for state in range(0, (len(pad2)-1)):
            width_mm = 1
            colour_mm = 'black'
            width_mi = 1
            colour_mi = "black"
            width_md = 1
            colour_md = "black"
            width_dd = 1
            colour_dd = "black"
            width_dm = 1
            colour_dm = "black"
            width_di = 1
            colour_di = "black"
            size = 10
            colour_ii = "black"
            width_im = 1
            colour_im = "black"
            width_id = 1
            colour_id = "black"
            if state == (len(pad2)-2):
                plt.arrow(x + 1, 3.5, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour_mm, ec=colour_mm,
                          linewidth=width_mm)  # mm
                plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour_mi, ec=colour_mi,
                          linewidth=width_mi)  # mi
            else:
                plt.arrow(x + 1, 3.5, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour_mm, ec=colour_mm,
                          linewidth=width_mm)  # mm
                plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour_mi, ec=colour_mi,
                          linewidth=width_mi)  # mi
                plt.arrow(x, 4, 3.5, 6.5, head_width=0.5, head_length=0.5, fc=colour_md, ec=colour_md,
                          linewidth=width_md)  # md
            if state == (len(pad2)-1):
                pass
            else:
                if state == (len(pad2) - 2):
                    plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour_im, ec=colour_im, linewidth=width_im)  # im
                    axes.plot([x - 1], [7.3], marker=r'$\circlearrowright$', ms=size, color=colour_ii)
                else:
                    plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour_im, ec=colour_im, linewidth=width_im)  # im
                    plt.arrow(x, 8.7, 3, 2.2, head_width=0.5, head_length=0.5, fc=colour_id, ec=colour_id, linewidth=width_id)  # id
                    axes.plot([x-1], [7.3], marker=r'$\circlearrowright$', ms=size, color=colour_ii)
            if state == 0 or state == (len(pad2)-1):
                pass
            else:
                if state == (len(pad2) - 2):
                    plt.arrow(x, 11, 3.6, -6.5, head_width=0.5, head_length=0.5, fc=colour_dm, ec=colour_dm,
                              linewidth=width_dm)  # dm
                    plt.arrow(x, 11, 0, -1.5, head_width=0.5, head_length=0.5, fc=colour_di, ec=colour_di,
                              linewidth=width_di)  # di
                else:
                    plt.arrow(x + 1, 12, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour_dd, ec=colour_dd,
                              linewidth=width_dd)  # dd
                    plt.arrow(x, 11, 0, -1.5, head_width=0.5, head_length=0.5, fc=colour_di, ec=colour_di,
                              linewidth=width_di)  # di
                    plt.arrow(x, 11, 3.6, -6.5, head_width=0.5, head_length=0.5, fc=colour_dm, ec=colour_dm,
                              linewidth=width_dm)  # dm
            if pad2[current] == 'match' and pad2[next] == 'match':
                width_mm = 5
                colour_mm = "red"
                plt.arrow(x + 1, 3.5, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour_mm, ec=colour_mm,
                          linewidth=width_mm)  # mm
            if pad2[current] == 'match' and pad2[next] == 'insert':
                width_mi = 5
                colour_mi = "red"
                plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour_mi, ec=colour_mi,
                          linewidth=width_mi)  # mi
                width_im = 5
                colour_im = "red"
                plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour_im, ec=colour_im,
                          linewidth=width_im)  # im
            if pad2[current] == 'match' and pad2[next] == 'deletion':
                width_md = 5
                colour_md = "red"
                plt.arrow(x, 4, 3.5, 6.5, head_width=0.5, head_length=0.5, fc=colour_md, ec=colour_md,
                          linewidth=width_md)  # md
            if pad2[current] == 'deletion' and pad2[next] == 'deletion':
                width_dd = 5
                colour_dd = "red"
                plt.arrow(x + 1, 12, 1.5, 0, head_width=0.5, head_length=0.5, fc=colour_dd, ec=colour_dd,
                          linewidth=width_dd)  # dd
            if pad2[current] == 'deletion' and pad2[next] == 'match':
                width_dm = 5
                colour_dm = "red"
                plt.arrow(x, 11, 3.6, -6.5, head_width=0.5, head_length=0.5, fc=colour_dm, ec=colour_dm,
                          linewidth=width_dm)  # dm
            if pad2[current] == 'deletion' and pad2[next] == 'insert':
                width_di = 5
                colour_di = "red"
                plt.arrow(x, 11, 0, -1.5, head_width=0.5, head_length=0.5, fc=colour_di, ec=colour_di,
                          linewidth=width_di)  # di
            if pad2[current] == 'insert':
                if state == 0:
                    if pad2[current] == 'insert' and pad2[next] == 'deletion':
                        width_id = 5
                        colour_id = "red"
                        plt.arrow(x, 8.7, 3, 2.2, head_width=0.5, head_length=0.5, fc=colour_id, ec=colour_id,
                                  linewidth=width_id)  # id
                    if pad2[current] == 'insert' and pad2[next] == 'match':
                        width_im = 5
                        colour_im = "red"
                        plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour_im, ec=colour_im,
                                  linewidth=width_im)  # im
                    if pad2[current] == 'insert' and pad2[next] == 'insert':
                        size = 12
                        colour_ii = "red"
                        axes.plot([x - 1], [7.3], marker=r'$\circlearrowright$', ms=size, color=colour_ii)
                else:
                    if pad2[current] == 'insert' and pad2[next] == 'deletion':
                        width_id = 5
                        colour_id = "red"
                        plt.arrow(x, 8.7, 3, 2.2, head_width=0.5, head_length=0.5, fc=colour_id, ec=colour_id,
                                  linewidth=width_id)  # id
                    if pad2[current] == 'insert' and pad2[next] == 'match':
                        width_im = 5
                        colour_im = "red"
                        plt.arrow(x, 7, 3, -2.5, head_width=0.5, head_length=0.5, fc=colour_im, ec=colour_im,
                                  linewidth=width_im)  # im
                        width_mi = 5
                        colour_mi = "red"
                        plt.arrow(x, 4, 0, 2.5, head_width=0.5, head_length=0.5, fc=colour_mi, ec=colour_mi,
                                  linewidth=width_mi)  # mi
                    if pad2[current] == 'insert' and pad2[next] == 'insert':
                        width_id = 5
                        colour_id = "red"
                        axes.plot([x - 1], [7.3], marker=r'$\circlearrowright$', ms=size, color=colour_ii)
            x += 4
            next += 1
            current += 1
            previous += 1
        rectangle = plt.Rectangle((x - 1, 3), width=2, height=1, fc='None', ec="black", linewidth=1.5)
        state_m0 = plt.text(x - 0.5, 3.2, s=('END'), size=t_size)
        axes.add_artist(rectangle)
        plt.xticks(self.match2)
        plt.yticks(self.y)
        if len(match) < 10:
            plt.xlim(0, 25)
            plt.ylim(1, 15)
        elif len(match) < 30:
            plt.xlim(0, 50)
            plt.ylim(1, 20)
        elif len(match) < 65:
            plt.xlim(0, 75)
            plt.ylim(1, 30)
        else:
            plt.xlim(0, 100)
            plt.ylim(1, 30)
        figure.suptitle("Profile HMM representing MSA")
        figure.tight_layout()
        plt.axis('off')
        if format == 'png':
            plt.savefig("viterbi_diagram.png")
        else:
            plt.savefig("viterbi_diagram.svg")
        diagram2 = plt.show()
        return diagram2



if __name__=='__main__':
    read_alignment = ReadAlignment()
    r = read_alignment.read_file(fasta_file)
    sequence_names = r[0]
    sequences = r[1]
    amount = r[2]
    # print("imported sequences from fasta file:")
    # print(sequences)
    read_alignment.checker()
    t_and_e_object = EmissionsTable()
    str2 = t_and_e_object.split()
    origin = t_and_e_object.check_type(str2)
    binary_table = t_and_e_object.binary_table(str2)
    insertion, threshold = t_and_e_object.conserved(origin)
    insert_reg = t_and_e_object.ranges(insertion)
    e_freq_counts = t_and_e_object.counts(origin, insert_reg)
    (tot_pseudo, tot_log) = t_and_e_object.emission(e_freq_counts)
    pssm_p = pd.DataFrame(tot_pseudo)
    transition_object = TransitionTable()
    binary_table, df_table = transition_object.bin_table()
    match, insertion = transition_object.concerved_region(binary_table)
    new_array = transition_object.current_states(match, df_table)
    total_iterations = transition_object.transition_probs(new_array)
    concatenated = transition_object.concatenate_it(total_iterations)
    transition_array = transition_object.logg_odd_scores(concatenated)
    # creating a df so that I can actually see what's going on
    t_df1 = pd.DataFrame(total_iterations, columns=['mm', 'md', 'mi', 'dd', 'dm', 'di', 'ii', 'im', 'id'])
    t_df2 = pd.DataFrame.swapaxes(t_df1, axis1="index", axis2="columns")
    # print("df2")
    # print(df2)  # this dataframe is the values before preforming calculations and before merging states
    t_df3 = pd.DataFrame(transition_array, columns=['mm', 'md', 'mi', 'dd', 'dm', 'di', 'ii', 'im', 'id'])
    t_df4 = pd.DataFrame.swapaxes(t_df3, axis1="index", axis2="columns")
    # print(df4)  # this is what happens when I merge stuff
    arrayy = t_df4.to_numpy()
    concatenated = pd.DataFrame(concatenated)
    algorithms = Algorithms()
    query_sequence, bi_o = algorithms.querysequence(insert_reg)
    all_pos_list = algorithms.forward(query_sequence, bi_o)
    best_path_list, states = algorithms.viterbi(query_sequence, bi_o)
    print("The likeliness that the query sequences matches the MSA is:", all_pos_list[-1])
    all_pos = pd.DataFrame(all_pos_list, columns=['m', 'i', 'd'])
    most_probable = pd.DataFrame(best_path_list)
    model = Model()
    # print(transition_array)
    model.transition_visualizer(transition_array)
    model.viterbi_vizualizer(best_path_list)
    list_dfs = [e_freq_counts, pssm_p, t_df2, t_df4, all_pos, most_probable]
    headings = ["Table 1: Emission Frequency counts (without pseudo counts)",
                "Table 2: Emissions probabilities (with pseudo counts)",
                "Table 3: Transition frequencies counts", "Table 4: Transition probabilities",
                "Table 5: Forward probabilities", "Table 6: Viterbi probabilities for best path"]
    with open(file_name, 'w') as f:
        for i in range(len(list_dfs)):
            f.write(headings[i])
            f.write("\n")
            (list_dfs[i]).to_csv(f)
            f.write("\n")
    print("Check the output directory for your tables")
