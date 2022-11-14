# pHMMgraph: Profile Hidden Markov Model graphs

## About the Tool
pHMMgraph takes a Multiple Sequence Alignment (MSA) as input in the form of a fasta file and produces a visually informative diagram 
that shows the most probable alignment of a sequence to a model as well as the weighted transition probabilities. 

## How to use this tool
* Download the pHMMgraph.py file from the github repository
* Install the dependancies using the requirements.txt file
* Run the command below within the project directory: `pip install -r requirements.txt`
* Obtain a MSA using a tool like muscle
* Prepare a query sequence (of the same length than the sequences in the MSA)
* Open the terminal (or CMD)
* Go to the directory where the downloaded program file was stored
* Run the program using the command `python3 pHMMgraph.py` followed by the _positional arguments_
* Example format: `python3 pHMMgraph.py filepath_to_programfile filepath_to_fasta_file sample_type query_sequence`
* Example: 
`python pHMMgraph.py C:\PycharmProjects\pHMMgraph.py "C:\PycharmProjects\demonstration_protein.fasta" P FK-R-AK`  
* For help use `python3 pHMMgraph.py -h` or `python3 pHMMgraph.py --help`

## Positional and Optional Arguments
```
pHMMgraph is a profile Hidden Markov Modeller that provides the user with a diagramatic representation of 
the results generated from the Forward and Viterbi algorithms

positional arguments:
  path                  the path to program file
  Fasta_file            Make sure the file is in the correct format otherwise it really wont work...
  datatype              Either P or N (protein or nucleotide)
  Query_sequence        Make sure the length of this sequence matches the sequence length in the fasta files

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output file name containing generated data in table format (.csv file)
  -p PSEUDOCOUNTS, --pseudocounts PSEUDOCOUNTS
                        pseudocounts method used (N, 0.1 or Pascal)
  -t THRESHOLD, --threshold THRESHOLD
                        threshold (to calculate conserved regions from)
  -f FORMAT, --format FORMAT
                        specify a format for the diagram output

If we knew what it was we were doing,it would not be called research, would it? - Albert Einstein

```
## Implementation
### Getting started
pHMMgraph can be downloaded from the GitHub repository. It runs on any platform with Python. It takes positional and optional arguments. The positional arguments, starts with the file path to the multiple sequence alignment file, in fastA format, the sample type, which is either N (nucleotide) or P (protein), and the query symbol sequence. The optional arguments have to be provided after ‘-‘, and the letter of the option that the user would like to specify. The options are -h (--help) to see the help message which contains information about the positional and optional arguments; -o (--output) to specify the file name for the output file, written in csv format; -p (--pseudo counts) to specify that transition and emission frequencies should be corrected with pseudo counts; -t (--threshold) to specify the threshold used to calculate the conserved regions, and lastly -f (--format), to specify the output format of the diagrams. 
The default setting for the optional arguments are: -o pHMM_tables.csv; -p ‘Pascal’; -t 50 and -f svg.  
### Classes and functions
The program consists of five classes: ReadAlignment, EmissionsTable, TransitionTable, Algorithms and Model. These classes were named according to the function they are performing in the program and encode functions that performs specific tasks within the class. 
#### ReadAlignment
The first class, ReadAlignment consists of two functions (excluding __init__()): read_file and checker. Read_file takes the file path to the fastA file, reads the file into memory and ensures that all of the sequences are of the same length, stores the amount of sequences and confirms that the fastA file is not empty. Checker does not require any input and performs a quick initial scan to determine whether the specified sample type is valid. 
#### EmissionsTable
EmissionTable’s first class is split, and it splits the sequences into a comma seperated string array, whereafter check_type does a second check to determine that the sequences do not contain any invalid characters. Thereafter, the conserved function determines the conserved and insert regions by comparing the column-wise percentage of the number of aligned symbols relative to the number of sequences, to the specified threshold. This is determined using a binary table (from function bin_table) where letters are replaced with 1's and -'s with 0's. The letters are then counted after which the total column count is divided by the amount of sequences and then multiplied by the specified threshold. This is followed by ranges, which is used to determine the length of an insertion region. When the insertion regions are determined, consecutive insert columns are concatenated into a combined insert region column as isert regions share an emission probability. An emission frequency counts table is created by the counts function. The count function counts the occurance of a specific letter in a column of the alignment. It counts all of the expected letters for either DNA, RNA or proteins. Lastly, the function emission then converts the count frequencies table into the emissions probability table. The emission probabilities are calculated using a combination of the emission frequency count and the preferred pseudo count multiplied by the emitted nucleotide frequency divided by the sum of the total possible emissions and the pseudo count.
#### TranstionTable
The next class is TransitionTable.  The first function, bin_table, creates a new binary table specific to the state transitions. To calculate the transition probabilities, the current_states function uses the lists with the column numbers of match and insertion states to create a new list containing the letters ‘m’, ‘i’ and ‘d’ showing the underlying states of the MSA. Transition_probs creates another frequency counts table, but this time for the state transitions. It is created using the combination of the current and the next letter in the aforementioned list, adding a count to the transition state corresponding to the combination of the current and next letter. Concatenate_it, includes the insertion counts into the preceding match state as the insert region following a match state falls into the same state on the pHMM. Lastly, logg_odd_scores is called to create a transition probability table. The transition probabilities is calculated per state by adding the pseudo count to the transition frequency count and then dividing it by the sum of the possible state transitions (which serves as the pseudo count- ensuring it accounts for at least one possible transition per state) and the actual amount of transitions for the state. There is an option to add log odd scores at this step, but it is not used for the algorithm by default, since the log odd scores are already incorporated into the equation used in the algorithm for efficiency.
#### Algorithms
The next class, Algorithms, starts by calling querysequence, which matches the conserved and insert regions for the query sequence to the MSA’s. A new list is created containing the emission probabilities for the symbols emitted by the query sequence.  Thereafter, Viterbi is called to determine the most probable state path through the model that could produce the query sequence. This is done by initializing the algorithm as M0 = 0. This is followed by the forward algorithm, that predicts the probability that the query sequence was produced by the model. 
#### Model
The last class is called Model and starts off with a function called transition_visualizer that creates a diagrammatic representation of the pHMM with the line width adjusted given the transition probabilities between states. Lastly, viterbi_visualizer is called to create a second model that shows the most probable path through the alignment using red arrows instead of the normal black arrows. 
After the user called the program in the terminal and provided input, the program will run and provide the user with feedback as it finishes tasks. Messages appear when a task is complete to make it easier to find the step where something went wrong if the program did not run to completion.


## How to interpret the output
The csv file consists of 6 dataframes. 
The first dataframe represent the emission frequencies counts table.
The second, the emission probabilities with pseudocounts added.
The third, the transition counts table.
The fourth, the transition probabilities table.
Then the probabilities to be at each state at a given time as calculated using the Forward Algorithm.
The final value here, shows the probability that the query sequence was produced by the model.
Lastly, the most probable path through the alignment calculated using the Viterbi algorithm.

Then the diagrams should hopefully speak for themselves as the purpose of this project was to make them visually informative.
But shortly, the red arrows shows the path that most likely produced the query sequence given the pHMM.
And the thickness of the lines shows the transition probabilities and how likely it is to be at a specific state at a given time. 

Disclaimers: 
* The CSV file to which you write the tables out cannot be open while running the program otherwise it will give an error. 
* If you do not specify a new output file for a new run, it will overwrite the csv from the previous run as it writes to the same file. 
* The diagram will get unreadably small if the sequence length of the sequences in the MSA is longer than 100. 

## About the Creator
Name: Anika du Plessis 

Current Occupation: Prospective BSc (Hons) Bioinformatics and Computational Biology

Institution: Stellenbosch University

Previous Tertiary Education: BSc Biochemistry and Microbiology (NWU)

Please contact me for more information on the inner workings of this tool or if any bugs crawl out

26703173@sun.ac.za
