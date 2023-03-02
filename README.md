# SSR Detector
Short script written in Perl language using regular expression techniques to detect SSR in the sequences provided in FASTA file format and save resulted output in tab delimited text file.

> Sample command line call: `perl ssr.pl fasta.txt > results.txt`

This code has been developed for the study of [Development of 1000 Microsatellite Markers across the Date Palm (Phoenix dactylifera L.) Genome](https://www.ishs.org/ishs-article/882_29) which won the prize of Khalifa International Date Palm Award as the most Distinguished Research/Studies in 2012 ([Winner Book 2012-2013-2014](https://www.kidpa.ae/sites/default/files/2012-2013-2014%20english.pdf#page=26)).

## Procedure Parameters:
- **$min_base and $max_base:** Minimum and maximum SSR base length (2 and 6 for example, so we are looking after di-, tri-, tetra-, penta-, and hexa- nucleotide repeats).
- **$min_repeat:** Minimum number of repeats to tag it as SSR (5 for example).
- **$max_split:** Maximum number of stretches may included in one compound microsatellite (3 for example).
- **$max_comp:** Maximum length of nucleotides between stretches of microsatellite consists of more than one type of repeat unit (4 for example).
- **$max_error_ratio:** The maximum percentage of accepted non-repeated nucleotides within the microsatellite motifs (0.20 for example, stands for 20%).

## Output Format:
- **ID** (string): Sequence ID as represented in the FASTA file.
- **Category** [Simple | Compound]: A microsatellite is referred as “simple”, if a single type of repeat unit repeats several times; a “compound‟ microsatellite consists of stretches of more than one type of repeat unit.
- **Base** (string): Nucleotide repeats.
- **Position** (integer): Start position.
- **Length** (integer): Total length.
- **Type** [Perfect | Imperfect]: A "perfect" microsatellite does not contain mutations or interruptions; an “imperfect‟ microsatellite contains mutations or interruptions.
- **Error** (%percentage): The percentage of accepted non-repeated nucleotides within the microsatellite motifs.
- **SSR** (string): Full detected SSR string.
