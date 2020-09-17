# Sequences in R

## Introduction

DNA sequences are comprised of **nucleotides**. Each nucleotide contains a phosphate group, a sugar group and a nitrogen base. The four types of nitrogen bases are adenine (**A**), thymine (**T**), guanine (**G**) and cytosine (**C**). Find out more about nucleotides here: https://www.genome.gov/genetics-glossary/Nucleotide. In the DNA code, these are organized to form various "genomic elements" - introns, exons, promoters, etc. Some of the elements are regulatory, and some are **transcribed** into an intermediate RNA format (mRNA) that can be further **translated** into proteins that are made up of Amino Acids. Each Amino Acid consists of a codon - a tri-nucleotide sequence that encodes one of 20 amino acids:

| A - C                   | Q - I                   | L - P                   | S - V                |
| ----------------------- | ----------------------- | ----------------------- | -------------------- |
| alanine - ala - A       | glutamine - gln - Q     | leucine - leu - L       | serine - ser - S     |
| arginine - arg - R      | glutamic acid - glu - E | lysine - lys - K        | threonine - thr - T  |
| asparagine - asn - N    | glycine - gly - G       | methionine - met - M    | tryptophan - trp - W |
| aspartic acid - asp - D | histidine - his - H     | phenylalanine - phe - F | tyrosine - tyr - Y   |
| cysteine - cys - C      | isoleucine - ile - I    | proline - pro - P       | valine - val - V     |

You can learn more about the transcription and translation process in our Introduction to Transcriptomics Course.

## Load a sequence as a DNA string

There are several ways you can store a sequence of letters in R, for example I can just make a string and call it something, like:`

```R
string1 <- "Gene sequence"
string2 <- "ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT"
```

once something is a string, I can perform several operations on this string like calculate it's length, find the number of characters of a specific type, or take a subset of a string at a given position. Some of these functions are easy in base R, but additional libraries have been developed like dplyr from tidyverse to deal with strings efficiently as well as specialized libraries to deal with strings of genomic data and genomic files like FASTA. In base R, you can write `nchar(string1)` to calculate the length, or use the tidyverse library to do the same: `str_length(string1)`

You can learn more about sequences in R here: https://www.r-bloggers.com/how-to-work-with-strings-in-base-r-an-overview-of-20-methods-for-daily-use/

## Convert DNA code to RNA

Let's assume you have 2 gene sequences and you would like to know what they would look like if they were trascribed into RNA. Remeber that RNA replaces the thymine nucleotide (T) with uracil (U). In R, you can use the function `chartr` to convert an "old" character to a "new character" in a given string. For example:

```R
string2 <- "ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT"

RNAstring2 <- chartr(old="T", new="U", string2)

print(RNAstring2)
```

The output will be:

"AUUUUAAGUAGUUAAGCCAGUGCCCGAUGCAAAGCGGUCCAUGAAUGCAAUGCCCUUGUAUAUAUAUGUGAU"

## Convert RNA to AA

In regular strings, you would have to manually go through the code and convert each try-nucleotide combination into an associated Amino Acid Letter, but there is a library in R called "[Biostrings]( https://kasperdanielhansen.github.io/genbioconductor/html/Biostrings.html)" that automates this type of functionality.

```R
library(Biostrings)

dna1 <- DNAString(string2)

translate(dna1)
```

The result for this sequence is:
seq: ILSS*ASARCKAVHECNALVYICD

you can read more about translation, codon sequences and working with open readng frames here: https://rdrr.io/bioc/Biostrings/man/translate.html



# Aligning and Comparing Sequences

Often times, sequences are hard to analyze because they can consist of 10s, 100s or even thousands of nucleotides. In many cases, we don't even know exactly what the sequences do and if we find some pattern what that might do to the product of that gene. One useful way to understand sequence variation is by performing sequence alignment.



## Pairwise Alignment

To illustrate this idea, let's align one sequence to another and see how similar they are. Let's imagine 3 sequnces in this case, 2 that are fairly similar (length is 72 and 77 bp) and another one that is much shorter (20 bp), so it has a similar segment, but overall the sequence is very different. First, let's align the similar sequences :

```R
library(DECIPHER)
library(Biostrings)

#load data as a DNA string
Patient1 <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT")
Patient2 <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC")
Patient3 <- DNAString("CAAAGCGGTCCATGAGATGC")

myAlign <- pairwiseAlignment(Patient1, Patient2)
print(myAlign)
```

The result of this code will look like this:

```bash
> Global PairwiseAlignmentsSingleSubject (1 of 1)
> pattern: ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGA-ATGCAATGCCCTTGTATATATATGTGAT----
> subject: ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC
> score: 102.6865 
```

What if we compare the long and the short sequences?

```r
myAlign <- pairwiseAlignment(Patient1, Patient3)
print(myAlign)
```

Now we will get the following:

```
> print(myAlign)
> Global PairwiseAlignmentsSingleSubject (1 of 1)
> pattern: ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT
> subject: -----------------------------CAAAGCGGTCCATGAG----ATGC-------------------
> score: -206.2459 
```

You can see that by default the alignment is global - one that tries to align the whole sequnece as mucha s possible by introducing gaps inbetween and at the ends of the whol sequence. There are additional alignment types, including global, local and overlap. Local alignment will seek to keep only the parts that algign, for example:

```r
localAlign <- pairwiseAlignment(Patient1, Patient3, type = "local")
print(localAlign)
```

the output would be:

```
> print(localAlign)
Local PairwiseAlignmentsSingleSubject (1 of 1)
pattern: [30] CAAAGCGGTCCATGA
subject:  [1] CAAAGCGGTCCATGA
score: 29.72634 
```

For overlap, the alignment would be:

```R
overlapAlign <- pairwiseAlignment(Patient1, Patient3, type = "overlap")
print(overlapAlign)
```

the output would be:

```
> print(overlapAlign)
Overlap PairwiseAlignmentsSingleSubject (1 of 1)
pattern: [30] CAAAGCGGTCCATGA-ATGC
subject:  [1] CAAAGCGGTCCATGAGATGC
score: 23.65337 
```

## Multiple Sequence Alignment:

What if we wanted to align multiple sequences? IN this case, it would be cumbersome to try to align all pairs, so we turn to a method called "multiple sequnece alignment". This method can also be used to compare multiple long and short sequences to find matching portions, introduce gaps and find a consensus sequence. In order to perform multiple sequence alignment in R without runnin g the traditional algorithms like MUSCLE, we can use the DECIPHER package from bioconductor.

```R
library(DECIPHER)

# multiple sequence alignment
Patient <- list()
DNASeq[1] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT")
DNASeq[2] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC")
DNASeq[3] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGATAA")

seqs <- DNAStringSet(unlist(DNASeq))

# perform the alignment
aligned <- AlignSeqs(seqs)

# print out the aligned sequences
print(aligned)
```

as a result, you will get:

```
> aligned
  A DNAStringSet instance of length 5
    width seq
[1]    78 ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGA--ATGCAATGCCCTTGTATATATATGTGAT----
[2]    78 ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAG-ATGCAATGCCCTTGTATATATATGTGATTTAC
[3]    78 ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGA--ATGCAATGCCCTTGTATATATATGTGATAA--
[4]    78 ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGTCTGCAATGCCCTTGTATATATATGTGAT----
[5]    78 ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGA--ATGCAATGCCCTTGTATATATATGTGAT----
```

You can also save the output as a FASTA file using the function `writeXStringSe()` and visualize the alignment using `BrowseSeqs()`

```R
# view the alignment in a browser (optional)
myFASTA <- "my.fasta"
writeXStringSet(aligned, myFASTA)

# view or save the alignment in an HTML file for a browser view
TF <- tempfile("plot___", fileext = ".html")
TF <- BrowseSeqs(aligned, highlight=0, htmlFile = TF)
file.copy(TF, './')
```

Read more about the DECIPHER package for multiple sequnce alignment here: https://rdrr.io/bioc/DECIPHER/f/inst/doc/ArtOfAlignmentInR.pdf

Try to run this full code:

```R
library(DECIPHER)

# multiple sequence alignment
Patient <- list()
DNASeq[1] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT")
DNASeq[2] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC")
DNASeq[3] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGATAA")

seqs <- DNAStringSet(unlist(DNASeq))

# perform the alignment
aligned <- AlignSeqs(seqs)

# print out the aligned sequences
print(aligned)

# view the alignment in a browser (optional)
myFASTA <- "my.fasta"
writeXStringSet(aligned, myFASTA)

# view or save the alignment in an HTML file for a browser view
TF <- tempfile("plot___", fileext = ".html")
TF <- BrowseSeqs(aligned, highlight=0, htmlFile = TF)
file.copy(TF, './')
```

# Assignment:

Now it is time to test ywhat you have learned. In this activity you will modify a python script to detect mutant patients and determine the recommended medicine dosage for mutant and wild type patients
In this exercise, we will solve the following case study based on the short lecture on pharmacogenomics and with the help of a python sequence alignment script.
Below are the partial gene sequences of CYP3A5 gene from 5 different individuals -

<u>Patient 1:</u> ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT
<u>Patient 2:</u> ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC
<u>Patient 3:</u> ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGATAA
<u>Patient 4:</u> ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGTCTGCAATGCCCTTGTATATATATGTGAT
<u>Patient 5:</u> ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT

Here is the wild type sequence for this gene - ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT

The second A in the sequence GAATGC is the position 6986. This mutation/SNP is denoted as 6986 A>G.

Consider that all 5 of these individuals are kidney transplant patients and their doctor needs to decide the dosage for immunosuppressant based on their genetic makeup. Here are the dosage guidelines from CPIC (Stanford Medical School):

![img](https://lh4.googleusercontent.com/ZJJ_eBCuXCE5LdLbZWnMZXchLr2XzAnW3I7kn5LKF6q7AReu__IpFN0nNBTp4KsEDYdIJDyS9UtbOG8oAM5UBX90iLPWypBEdqwQFkRrFogagnL7gwo1fSJO5udg62wucmwk__zi)![Figure](https://user-images.githubusercontent.com/71350730/93373096-7cae9d80-f872-11ea-9a75-35b8787eb24d.png)  

\## Table : Dosage recommendation for Tacrolimus based on CYP3A5phenotype by CPIC (Stanford Medical School) [Sources link: https://cpicpgx.org/content/guideline/publication/tacrolimus/2015/25801146.pdf ].

Submit your answer below:

```R
#Submit your assignment here




```

