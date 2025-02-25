# BentonAdvancedPython_Portfolio
This is my portfolio from my advanced python work.

# Advanced Python Final Project

## Sequence Objectives

We begin our Advance Python work with learning Sequence Objectives. After loading the data in terminal, you can load the project packs here to begin our work.
```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```



We can give numbers to the letters in the sequence.
```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G




We can also print the length of each sequence.
```python
print(len(my_seq))
```

    5




You can use the numbers to see where the letter is in the sequence.
```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T




We can do a count code, which give us how many times the letter/letters you're looking for appears in the sequence.
```python
Seq("AAAA").count("AA")
```




    2





Random sequence used for example.
```python
my_seq = Seq("GAATCCAGTACTAGTGGAATTACGTAGTAAGGTCCCATCCA")
```



We can count the length of our sequence.
```python
len(my_seq)
```




    41





We can count how many of each letter in the sequence.
```python
my_seq.count("G")
```




    9





We can measure the "G" "C" count which is important for designing primers.
```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    43.90243902439025





Now we need to import the gc fraction to begin our next example.
```python
from Bio.SeqUtils import gc_fraction
```



We made up a sequence again for this example.
```python
my_seq = Seq("GAATCCAGTAGTGGAATTACGTAGTAAGGTCCCATCCA")
```


We can get the percentage of the "G" and "C" in the sequence.
```python
gc_fraction(my_seq)
```




    0.4473684210526316




Now we're going to slice sequences into parts.
```python
my_seq[4:12]
```




    Seq('CCAGTAGT')





We can also slice out cut outs of our sequences.
```python
my_seq[0::3]
```




    Seq('GTAAGAATTGCAC')




```python
my_seq[1::3]
```




    Seq('ACGGGTCAAGCTA')




You can locate a specific gene and chose our start position.
```python
my_seq[2:3]
```




    Seq('A')




We can use negatives to start from the other end.
```python
my_seq[::-1]
```




    Seq('ACCTACCCTGGAATGATGCATTAAGGTGATGACCTAAG')




If you have a seq object and want to turn it back into a string, you can use this code.
```python
str(my_seq)
```




    'GAATCCAGTAGTGGAATTACGTAGTAAGGTCCCATCCA'





Now we will start working with fasta files.
```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


You can print out the information of the file using the print code.
```python
print(fasta_format_string)
```

    >Name
    GAATCCAGTAGTGGAATTACGTAGTAAGGTCCCATCCA
    



You can add two scripts together.
```python
seq1 = Seq("ACGT")
seq2 = Seq("AGTAAG")
```


Adding the two together
```python
seq1 + seq2
```




    Seq('ACGTAGTAAG')




They will add together in the order you type them in, so this is the reverse of the first example.
```python
seq2 + seq1
```




    Seq('AGTAAGACGT')





Starting contig
```python
contigs = [Seq("ATG"), Seq("ATCCGA"), Seq("TTAGCA")]
```


```python
spacer = Seq("N" *10)
```



Now we can use the joint function to take the spacer option we made and we compute it to join with the contigs.
```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCGANNNNNNNNNNTTAGCA')





Sometimes you have case sensitive stuff, so there is a code to help with that so that you won't have that problem.
```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




Can make everything uppercase.
```python
dna_seq.upper()
```




    Seq('ACGTACGT')




Can make everything lowercase. These are used to easily manipulate the strings, so you don't have to go through the sequence individually.
```python
dna_seq.lower()
```




    Seq('acgtacgt')





Since the sequence is in all uppercase then when typed in lowercase it will be false.
```python
"gtac" in dna_seq
```




    False



 You can change the lowercase to uppercase by putting in the code for upper and saving it.
```python
dna_seq = dna_seq.upper()
```

Now the sequence will be true for uppercase.
```python
"GTAC" in dna_seq
```




    True





Loading our example sequence.
```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAATCGC")
```


Getting the complementing base pairs for our sequence.
```python
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTAGCG')




Can get the reverse complementing base pairs. It will do the reverse and use the -1 together.
```python
my_seq.reverse_complement()
```




    Seq('GCGATTTCGATCCTATATAGGCCCATCGATC')




You can also get the reverse compliment of your protein sequence. It's pulling the compliment from ambiguity codes.
```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




Create a random coding DNA sequence.
```python
coding_dna = Seq("ATGCCATGTAGTCTGATCCATGTGTAACCATTAAACTGTGTA")
```


Checking to make sure our code is there.
```python
coding_dna
```




    Seq('ATGCCATGTAGTCTGATCCATGTGTAACCATTAAACTGTGTA')




Create a template DNA sequence.
```python
template_dna = coding_dna.reverse_complement()
```


Now we can view our reverse complement which is how transciption is done.
```python
template_dna
```




    Seq('TACACAGTTTAATGGTTACACATGGATCAGACTACATGGCAT')



Now we will transcribe it into a protein.
```python
coding_dna
```




    Seq('ATGCCATGTAGTCTGATCCATGTGTAACCATTAAACTGTGTA')




Can transcibe the T base pairs to U using this code to transribe it.
```python
messenger_rna = coding_dna.transcribe()
```


Now we can view our transcribed data.
```python
messenger_rna
```




    Seq('AUGCCAUGUAGUCUGAUCCAUGUGUAACCAUUAAACUGUGUA')




We can get the same results using this code instead.
```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGCCAUGUAGUCUGAUCCAUGUGUAACCAUUAAACUGUGUA')




You can also do reverse transcription.
```python
messenger_rna.back_transcribe()
```




    Seq('ATGCCATGTAGTCTGATCCATGTGTAACCATTAAACTGTGTA')




```python
messenger_rna
```




    Seq('AUGCCAUGUAGUCUGAUCCAUGUGUAACCAUUAAACUGUGUA')




We can also translate the sequence. RNA is turning into protein. The "*" are stop codons.
```python
messenger_rna.translate()
```




    Seq('MPCSLIHV*PLNCV')




Doing translation of a mitochondria genome.
```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MPCSLIHV*PLNCV')



You can use this code to give you the same results.
```python
coding_dna.translate(table = 2)
```




    Seq('MPCSLIHV*PLNCV')



Translate the neclotons to the first stop codon.
```python
coding_dna.translate(to_stop = True)
```




    Seq('MPCSLIHV')




```python
coding_dna.translate(table =2, to_stop=True)
```




    Seq('MPCSLIHV')




You can change the stop symbol to whatever you want.
```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MPCSLIHV!PLNCV')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTGCACTTTCCCTGGTTCTGGTCGCTCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


You can get the Seq using this code.
```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLHFPWFWSLHGSTGCGNYVSPVSKITDRRS**SWLLLGWRSLARPR...SPL')



We can pull from before the first stop codon.
```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLHFPWFWSLHGSTGCGNYVSPVSKITDRRS')





Tells python to start with methylanine since this is a complete gene. If the sequence is not valid, then python will spit out an error. Since mine is not valid, am error occurs.
```python
# gene.translate(table = "Bacterial", cds = True)
```


Let's look at codon tables.
```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```

We can visual the tables by printing them.
```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



We can also see our mitochondrial table.
```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



We can ask what the stop codons are in our mitochondrial table.
```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




We can also see the start codons.
```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




We can use sequence objects to help navigate though our code.
```python
seq = Seq("ACGT")
```


Testing to see if the codons equal the seq object.
```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




Sequences can be represented by creating a sequence object for none followed by the sequence. We create a well defined length.
```python
unknown_seq = Seq(None, 10)
```


We can get the information from the code.
```python
unknown_seq
```




    Seq(None, length=10)




You can also get the length using this code.
```python
len(unknown_seq)
```




    10




We can pull information from seq objectives.
```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


The seq will only give you information based off what you have. Since we don't have more than 20 base pairs, then the seq can't give us what we need from out input.
```python
seq[1000:1020]
```




    Seq(None, length=20)



You have to code for what's defined in your sequence.
```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




This shows that you can have partial information for a chromosome combined with placeholders for length.
```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




You can also created these partial unknown sequenced by appending.
```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


We can also produce mutations in our seq. We have to download the program needed first.
```python
from Bio.Seq import MutableSeq
```


You have to import mutable Seq.
```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')



This is how we can cause mutations in our sequence.
```python
mutable_seq[5] = "C"
```

Now we see that we have changed the fifth codon to "C".
```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




We can also remove the codons from our sequence if needed using this code.
```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')



You can get the reverse of your codon.
```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')



If you want to change your code and keep protection from making it mutable.
```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')



If you don't want to used Seq object, you can use this instead.
```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate 
```


You have to create a string.
```python
my_string = "GCTGTTATGGGTCGTTGGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCCAACGACCCATAACAGC'



If you don't want to download Biopython, then you can transcribe your data.
```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGGAAGGGTGGTCGTGCTGCTGGTTAG'



 How you get your stop codon using this method.
```python
translate(my_string)
```




    'AVMGRWEGWSCCWL'

This concludes our sequencing objectives.

## Sequence Annotations
Now we will begin looking at our work from sequence annotations.


