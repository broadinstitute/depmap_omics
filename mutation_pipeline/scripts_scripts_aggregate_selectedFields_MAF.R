## reads selected fields from maf files and aggregates them 
args<-commandArgs(TRUE)
print(args)

outfn = args[1];
inFNs = args[2];

SelectFields = as.character(args[[3]])  ## use SelectFields = 'Chromosome,Start_position,t_alt_count,t_ref_count,Variant_Type' for PoN filter
selFields=unlist(strsplit(SelectFields,','))


#inFNs = '/xchip/cle/analysis/mghandi/CCLE/R/ExonUsage/fns.txt'

fns = as.character(read.csv(inFNs, header=FALSE)[,])
fe = sapply(fns, file.exists)
samples = fns[which(fe)]

nsamples = length(samples)

readsel = function(fn, selFields){
  f = read.delim(fn, sep='\t', comment.char = '#', header = TRUE, quote='')
  ii = match(selFields,colnames(f))
  if(length(which(is.na(ii)))>0){
    cat('\n Warning: some fields not found: ', fn,' ', selFields[which(is.na(ii))])
  }
  res = f[,ii]
}

for(i in 1:nsamples){
  ri = readsel(samples[i],selFields)
  write.table(ri,  file=outfn,  sep='\t', row.names = FALSE, quote = FALSE, append = (i>1), col.names = !(i>1))
}
#"Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Chromosome","Start_position","End_position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS","dbSNP_Val_Status","Genome_Change","Annotation_Transcript","Tumor_Sample_Barcode","Protein_Change","t_alt_count","t_ref_count","realign_judgment","keep_min","keep_mean","keep_max","n_tie","reject_min","reject_mean","reject_max","keep_scores","reject_scores","Passed_Filters","CGC_Cancer_Germline_Mut","HGNC_RefSeq_supplied_by_NCBI_","HGNC_Primary_IDs","reads_mapped_away","mate_not_mapped_away","Verification_Status","Score","n_ref_count","HGNC_EnsemblGeneID","HGVS_protein_change","dbNSFP_Ancestral_allele","UniProt_Site","ESP_Chromosome","Codon_Change","dbNSFP_MutationAssessor_pred","dbNSFP_phastCons46way_placental_rankscore","ClinVar_TYPE","dbNSFP_GERP_RS_rankscore","dbNSFP_Interpro_domain","CCLE_ONCOMAP_overlapping_mutations","dbNSFP_Uniprot_id","RC","covered","observed_in_normals_count","Transcript_Exon","n_alt_count","failure_reasons","UniProt_Natural_Variations","ESP_CG","CGC_CancerSyndrome","gnomADg_AF","gnomADg_AC","cDNA_Change","tumor_f","UniProt_Experimental_Info","total_reads","COSMIC_overlapping_mutation_descriptions","UniProt_Region","CGC_Tumor_Types_Germline","Refseq_mRNA_Id","GO_Cellular_Component","GO_Molecular_Function","gnomADg_AF_FIN","Ensembl_so_term","COSMIC_overlapping_mutations","SwissProt_acc_Id","COSMIC_total_alterations_in_gene","COSMIC_overlapping_primary_sites","CGC_CancerGermlineMut","gnomADg_AF_NFE","DrugBank","COSMIC_n_overlapping_mutations","CGC_Cancer_Molecular_Genetics","CGC_Cancer_Syndrome"



