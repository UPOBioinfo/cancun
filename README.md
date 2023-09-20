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

## R scripts (figure script: input files)
* fig2.R: all.gggenes
* fig3.R: islands.xlsx, p1virus.ab, p1virus_incomplete.ab, ../types/ifa.ab, ../types/ifb.ab
* fig4.R: metadata_ab_def_caudovirus.tsv, tree.tree, metadata_st79.tsv, ppcrispr/data/CP087336/polymorphisms_move.phy_phyml_tree.txt, metadata_ppcrispr.tsv, pp_completos/data/CP087336/polymorphisms_move.phy_phyml_tree.txt, metadata_pp.tsv
* fig5.R: nspacers3.tsv, nspacers_types.tsv, rp_ifa1.seq, rp_ifa2.seq, rp_ifb.seq
* fig6.R: variants_spacers2.tsv
* fig7.R: pam_p1virus_ifa1_12.seq, pam_p1virus_ifa2_12.seq, pam_p1virus_ifb_12.seq, spacers_vs_ab00002_80_80_2.tsv, spacers_vs_p1virusffn_ditribution.tsv, ../gggenes/p1virus.gggenes, pam_freq_p1virus_cds.tsv
* supplfigs2.R: ../../spacers_new/virus_100_phage.id, ../../strains.ab, phage_matrix_p1virus.tsv, metadata_ab_phigaro.tsv

# Additional files
* ref_phigaro_ab.fasta: prophages (references)
* all_ab_phigaro_out_sup8kb.fasta: prophages (all)
* ab_phigaro.asoc: prophages (prophage clusters)
* spacers_ab_uniq.fasta: non-redundant spacers
* spacers_ab_uniq_types2.tsv: spacer types
