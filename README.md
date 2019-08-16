# BA-Thesis
My BA Thesis concerning bioinformatical analysis of IDP (intrinsically disordered proteins).

IDP are proteins that lack a fixed or ordered three-dimensional structure. 
My goal was to check if such a proteins differ from normal (that is structurized) proteins in amino acid composition. If yes, how? 


Because there is not much experimentally confirmed IDP proteins (less than 3000 sequences in disprot database), results would be statisticaly not significant. So I've decided to use program ([IUPred2a](https://iupred2a.elte.hu)) which predicts if protein or its part is IDP or not. Dataset of sequences which I have used is UniRef50. It is nonredundant and clustered database of all known sequences.
After predicting "IDPcity" of approximately 171 000 sequences I have splited obtained sequences into two groups: IDP and nonIDP using **afterIUPred.py**. In the next step I extracted from each sequence n-meres and counted them (I've received almost all possible permutations of 5-meres) by **divider.py**. 
Because of the fact that in some organisms aminoacids have different usage I've calculated z-scores of every n-mere using the size of permutation group for which it belong. This was done using **someStats.py** script.
Hydrophobicity of n-meres was calculated by **hydroph.R** with [Peptides](https://github.com/cran/Peptides) package. 
Obtained data were next comapred (**comparator.py**).


# Programms used:

## afterIUPred.py
Script takes the data from iupred2a and divide it to two files. In one there are fragments of proteins which were classified as IDP-like, and in the other there are fragments classified as nonIDP. One protein sequence can be splited many time, in particular into very small fragments (for example of length 1). So the script takes advantage of running mean to make the predictive function smooth and get more longer sequences.

```bash
$ ./afterIUPred.py -i iupred_20-30k 
$ head nonIDP_10-20k.fasta 
>UniRef50_A8PWH1 Structure-specific endonuclease subunit SLX1 n=1 Tax=Malassezia globosa (strain ATCC MYA-4612 / CBS 7966) TaxID=425265 RepID=SLX1_MALGO
MRISSVVEHRTPRVYVCYCLRSLSRPNQTYIGS
QHNGLVKQGAFYTRMARPWTMDVVVYGFPSKLAALQFEWSWQKPHASRHLR
ATAAVYAGRSSKPLFPATRS
RVRMRSSTVPEYKFLVLRALLASEPFCFWHLHVGFYSEYAYGVWQFMDRANPTRYSVSRITRRPLPPSYPPVACDF
SYPVLPEATTETLGLTWEQLEHAP
LRSHTSFLEALDQDASALEAHLVHASRKDMSSSICGLCGGHINRHVPLSYTHCPHACDAVFHLTCLARYSLEQETRAHARTFCLPTSAWCPMCQRPPVPWPEIVRRVFRRAELKVSM

>UniRef50_A9BER7 GMP synthase [glutamine-hydrolyzing] n=547 RepID=GUAA_PETMO
MEKILVIDYGSQYTQLLAKRIRDLGVFSEVIQYDDNISLSNVKGIILSGGPDSVYNIDAPDISDEILNAELPILGICYGMQLIAKKLGGKVEQRGIAEYGKTKINITDQSLLFKKIPSTFNVWMSHKDMVTKVPEKFKITSLTSNNIISSFENESENIYCIQFHPEVRHTEFGIDILKNFIQGICGLKGSWTLMDFVENKIKEIKDTIGDKKAIIALSGGVDSSVAAVLTHRAIGNNLKAIFVNHGFLRMNEVEEVESTFRDYMGLNLTTVDAQERFLSKLKGVTDPEQKRKIIGEEFIRVFEQEAKKEEGCEYLIQGTIYSDVIESAKSGKKTFKIKSHHNVGGLPEDIDLKIVEPLKELFKDEVRSVGEILGLPREILYRHPFPGPGLAIRIMGEINDEKLTILKKVDNIFINTLKETGWYDKVWQAFAVLIPVKTVGITGDKRSYGYVAALRSVDSVEGMTADWSKVPFEILDLVSSRITNEVEEITRVVYDISSKPPATIEWE
```


## divider.py
This script returns n-meres sequences of proteins which were provided. 
It takes as an input .fasta file with many protein sequences and n which is the length of meres.
Additionaly this program counts the number of occurances of all meres.

As an example:
```bash
./divider.py -i test.fasta -n 5 
```
returns a file 
```bash
$ cat test_5-meres
AAAAA	5490
LLLLL	3391
HTGEK	3059
TGEKP	2937
GEKPY	2593
SSSSS	2589
QQQQQ	2407
CGKAF	2075
IHTGE	1731
GGGGG	1363
ECGKA	1358
...
```


## someStats.py
This program is for anylysis of short n-meres-peptides in proteins. It returns z-scores of n-meres-peptides based on the size of its permutation group. 
Input file should contain nmeres with its counts in dataset.

```bash
$ ./someStats.py -i test_5-meres -n 5
```
Output file is all possible n-mers with its z-score.
```bash
$ cat test_5-meres_ZScores
HTGEK	496.95
TGEKP	415.91
GEKPY	407.37
CGKAF	395.1
IHTGE	312.81
ECGKA	269.92
RIHTG	252.2
KPYEC	223.42
...
```

## hydroph.R
Script for assigning a hydrophobicity for a column of n-mer sequences
```bash
$ ./hydroph.R test_5-meres_ZScores > test_5-meres_Z_hydr

HTGEK 	496.95 	0.426 
TGEKP 	415.91 	0.672 
GEKPY 	407.37 	0.616 
CGKAF 	395.1 	-0.342 
IHTGE 	312.81 	-0.08 
ECGKA 	269.92 	0.248 
RIHTG 	252.2 	0.136 
KPYEC 	223.42 	0.266 
KPYKC 	215.15 	0.38 
HQRIH 	183.05 	0.146 

```



## comparator.py
This program compares two files with n-meres and their z-scores (and optionally hydrophobicity).
The result is a plot comparing two sets of proteins or file with merged columns.
```bash
$ ./comparator.py -i nonIDP_Z_hydr -j IDP_Z_hydr -p plot -b
```





