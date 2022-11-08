# pHMMgraph: Profile Hidden Markov Model graphs

## About the Tool
pHMMgraph takes a Multiple Sequence Alignment (MSA) as input in the form of a fasta file and produces a visually informative diagram 
that shows the most probable alignment of a sequence to a model as well as the weighted transition probabilities. 

## How to use this tool
* Download the pHMMdi.py file from the github repository
* Obtain a MSA using a tool like muscle
* Prepare a query sequence (of the same length than the sequences in the MSA
* Open the terminal (or CMD)
* Go to the directory where the downloaded program file was stored
* Run the program using the command `python3 pHMMdi.py` followed by the _positional arguments_
* For help use `python3 pHMMdi.py -h` or `python3 pHMMdi.py --help`

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
