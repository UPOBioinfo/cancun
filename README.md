# Scripts used to analyze *Acinetobacter baumannii* CRISPR-Cas systems

## Perl scripts
* bindConting2ab.pl - bind contigs of phage-plasmids from a Blast report
* bindConting2ab_phageplasmid.pl - bind contigs from a Blast report
* convertffn2referencesFromRoary.pl - search for CDS sequences correspondig to pangenome references
* countVirusVariantsBySpacers.pl - count the number of viral variants targeted by the spacers
* createFileForPeaks1.pl - create a distribution plot of protospacers
* createFileForPeaks1_genes.pl - create a distribution plot of protospacers (by genes)
* createGroupsOfRedundantSequences.pl - cluster redundant sequences
* extractCancunSequences.pl - extract genomic islands
* extractSpacerTypes.pl - classify spacers by CRISPR-Cas type
* extractSpacerTypes_withDirections.pl - classify spacers by CRISPR-Cas type (oriented based on the repeat sequences)
* getPamSequences.pl - extract PAM sequences
* maskCRISPRinViruses.pl - mask CRISPR arrays with Ns
* typeSpacersFromFasta.pl
* typeSpacersFromFasta_frequency.pl - extract and classify spacers (adding frequencies)

## R scripts (figures)
* fig2.R
* fig3.R
* fig4.R
* fig5.R
* fig6.R
* fig7.R
* supplfigs2.R
