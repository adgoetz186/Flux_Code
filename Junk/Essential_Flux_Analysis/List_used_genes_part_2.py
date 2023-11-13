import json
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
import os
import numpy as np
import collections
os.chdir("../..")

# Genes which are in the metabolic network but not in the scRNAseq data
non_matching_list = ['ADSS1', 'ADSS2', 'GPAT4', 'ADH1A', 'ADH1B', 'PTGR3', 'ARG1', 'CDO1', 'ALDOB', 'G6PC1', 'GLUD2', 'MMUT', 'NME2', 'PDHA2', 'PCK1', 'PLPP1', 'PLPP3', 'PLPP2', 'PRPS1L1', 'TKTL2', 'GGT2P', 'GGT6', 'AGXT', 'NAPRT', 'PYCR3', 'ATP5MGL', 'ATP5MG', 'DMAC2L', 'ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5F1E', 'ATP5PB', 'ATP5MC2', 'ATP5PD', 'ATP5ME', 'ATP5PF', 'ATP5MF', 'ATP5PO', 'ATP5MC1', 'ATP5MC3', 'COX6A2', 'CYTB', 'SLC5A3', 'AQP4', 'SLC16A8', 'SLC7A3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'COX7A1', 'COX1', 'COX2', 'COX3', 'SLC22A7', 'AQP10']

# Genes which are in the model, obtained from part 1
genes_in_model_symbol = ['AOC2', 'AOC3', 'AOC1', 'SDS', 'PDHX', 'DLD', 'DLST', 'OGDH', 'HPD', 'KDSR', 'HAAO', 'GOT1', 'GOT2', 'AACS', 'AADAT', 'ABAT', 'ALDH9A1', 'HADHB', 'ACAA2', 'ACAT1', 'ACAT2', 'HADHA', 'ACACA', 'ACACB', 'ACLY', 'ACADM', 'ACADSB', 'IVD', 'ACAD8', 'ACOX1', 'FASN', 'ACO1', 'IREB2', 'ACO2', 'ACSS2', 'CYP2E1', 'ADA', 'AK7', 'AK1', 'AK2', 'AK5', 'AK4', 'AK3', 'AMD1', 'ADK', 'APRT', 'ADSL', 'ADSS1', 'ADSS2', 'GPAT4', 'AGPAT4', 'AGPAT5', 'MBOAT2', 'AGPAT1', 'AGPAT2', 'AGPAT3', 'ETNPPL', 'AGXT2', 'AHCY', 'AHCYL1', 'ATIC', 'AKR1C4', 'GPT2', 'GPT', 'ADHFE1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'PTGR3', 'ALDH3A2', 'AKR1A1', 'AKR1B1', 'AKR7A2', 'AMPD1', 'AMPD2', 'AMPD3', 'ARG1', 'ARG2', 'ASL', 'ASS1', 'ASRGL1', 'ASNS', 'GAD1', 'CAD', 'ALDH1A2', 'ALDH1A1', 'ALDH7A1', 'UPB1', 'TM7SF2', 'NSDHL', 'HSD17B4', 'MSMO1', 'CAT', 'CPS1', 'CDIPT', 'CDS1', 'CDS2', 'CHPT1', 'CEPT1', 'PHOSPHO1', 'PCYT1A', 'PCYT1B', 'CHKA', 'CHKB', 'CRLS1', 'CS', 'CRAT', 'CROT', 'CTPS1', 'CTPS2', 'CDO1', 'CTH', 'CBS', 'CMPK1', 'CMPK2', 'DGKA', 'DGKB', 'DGKD', 'DGKE', 'DGKG', 'DGKH', 'DGKI', 'DGKQ', 'DGKZ', 'DCTD', 'CDA', 'AICDA', 'SCD', 'GUK1', 'DGUOK', 'DHCR24', 'CYP7A1', 'DHCR7', 'DEGS1', 'DEGS2', 'DHFR', 'DPYS', 'QDPR', 'GGPS1', 'FDPS', 'MVD', 'DERA', 'DTYMK', 'DPYD', 'TK1', 'TK2', 'UPP2', 'PNP', 'DUT', 'ITPA', 'EBP', 'ECHS1', 'AUH', 'ENO1', 'ENO3', 'ENO2', 'ETFA', 'ETFB', 'ETFDH', 'ETNK1', 'ETNK2', 'OLAH', 'ACSBG2', 'ACSL1', 'ELOVL2', 'ELOVL6', 'ELOVL5', 'ALDOA', 'ALDOB', 'ALDOC', 'ALDH1L1', 'AFMID', 'FTCD', 'ALDH1L2', 'MTHFD1', 'MTHFD1L', 'FH', 'FAH', 'GPD1', 'ALDH18A1', 'GNPDA2', 'GNPDA1', 'H6PD', 'G6PD', 'G6PC3', 'G6PC2', 'G6PC1', 'UGP2', 'GAPDHS', 'GAPDH', 'GART', 'ALDH3A1', 'ALDH1A3', 'ALDH3B1', 'ALDH3B2', 'GCSH', 'GLDC', 'AMT', 'GFPT1', 'GFPT2', 'SHMT1', 'SHMT2', 'GAD2', 'GLUD1', 'GLUD2', 'GLS2', 'GLS', 'PPAT', 'GCDH', 'GCAT', 'HAO1', 'HAO2', 'GK', 'GK2', 'GMPS', 'PGD', 'GNMT', 'GPAM', 'GSR', 'HPRT1', 'CA1', 'CA12', 'CA14', 'CA2', 'CA3', 'CA4', 'CA6', 'CA7', 'CA9', 'CA13', 'CA5A', 'CA5B', 'HADH', 'HSD17B10', 'EHHADH', 'HKDC1', 'GCK', 'HK1', 'HK2', 'HK3', 'HGD', 'HIBADH', 'HAL', 'KYNU', 'HMGCS1', 'HMGCS2', 'HMGCLL1', 'HMGCL', 'GRHPR', 'HSD11B1', 'HSD11B2', 'HSD17B8', 'HSD17B1', 'HSD17B7', 'HSD17B2', 'HSD3B2', 'HSD3B7', 'IDH3A', 'IDH3B', 'IDH3G', 'IDH1', 'IDH2', 'BCAT1', 'BCAT2', 'IMPDH1', 'IMPDH2', 'MIOX', 'IDI2', 'IDI1', 'SUCLG1', 'SUCLG2', 'SUCLA2', 'AMDHD1', 'KHK', 'KMO', 'LDHD', 'LDHAL6B', 'LDHAL6A', 'UEVLD', 'LDHA', 'LDHB', 'LDHC', 'GLO1', 'LSS', 'CEL', 'LPL', 'MGLL', 'LIPC', 'SC5D', 'GSTZ1', 'MCCC1', 'MCCC2', 'MLYCD', 'MDH1B', 'MDH1', 'MDH2', 'ME2', 'ME1', 'ME3', 'MAT1A', 'MAT2A', 'MAT2B', 'MTR', 'MVK', 'IMPA1', 'IMPA2', 'MTMR1', 'MTMR2', 'ISYNA1', 'MCEE', 'MMUT', 'ALDH6A1', 'MTAP', 'MTHFD2', 'MTHFR', 'MTHFD2L', 'ACY3', 'ASPA', 'NME7', 'NME6', 'NME1', 'NME2', 'NME3', 'NME4', 'NT5C', 'NT5C3A', 'NT5C1B', 'NT5E', 'NT5C2', 'GUCA1A', 'OTC', 'OXCT2', 'OXCT1', 'DBT', 'BCKDHB', 'BCKDHA', 'UMPS', 'ODC1', 'OAT', 'PYCR1', 'PLD1', 'PLD2', 'ACMSD', 'PC', 'DLAT', 'PDHA1', 'PDHB', 'PDHA2', 'PCK1', 'PCYT2', 'PEMT', 'PFKP', 'PFKL', 'PFKM', 'PHGDH', 'GPI', 'PGK1', 'PGK2', 'PGLS', 'BPGM', 'PGAM1', 'PGAM2', 'PGM1', 'PGM2', 'PTPMT1', 'PGS1', 'PAH', 'PLCB1', 'PLCE1', 'PLCZ1', 'PLCXD2', 'PLCH2', 'PLCH1', 'PLCB2', 'PLCB3', 'PLCB4', 'PLCD1', 'PLCD3', 'PLCD4', 'PLCL1', 'PLCG1', 'PLCG2', 'PMVK', 'LHPP', 'PPA1', 'PPA2', 'PLPP1', 'PLPP3', 'PLPP2', 'PCCA', 'PCCB', 'ACAD9', 'ACAD10', 'ACAD11', 'ACADS', 'PFAS', 'PRODH2', 'PRODH', 'PRPS1', 'PRPS1L1', 'PRPS2', 'PISD', 'PSAT1', 'PSPH', 'PTDSS1', 'PTDSS2', 'PKM', 'PKLR', 'UPP1', 'RDH16', 'RDH5', 'RDH8', 'RDH11', 'RDH10', 'RDH12', 'RDH13', 'RDH14', 'SDR16C5', 'RRM1', 'RRM2', 'RRM2B', 'RPE', 'RPIA', 'SARDH', 'SORD', 'SPTLC1', 'SPTLC2', 'SPTLC3', 'SGMS1', 'SRM', 'SQLE', 'FDFT1', 'SDHA', 'SDHB', 'SDHC', 'SDHD', 'TALDO1', 'PCBD1', 'SDSL', 'TKT', 'TKTL1', 'TKTL2', 'TYMP', 'TYMS', 'TPI1', 'TXNRD1', 'TDO2', 'IDO2', 'IDO1', 'TAT', 'UGDH', 'UROC1', 'UCK2', 'UCK1', 'UCKL1', 'XYLB', 'DCXR', 'ALDH4A1', 'NIT2', 'ACSS1', 'ASPG', 'GGT1', 'GGT2P', 'GGT5', 'GGT6', 'GGT7', 'GGTLC1', 'GGTLC2', 'GGTLC3', 'SUOX', 'AGXT', 'ALDH5A1', 'TMEM91', 'NAPRT', 'QPRT', 'ALDH2', 'ALDH1B1', 'XDH', 'CPT1A', 'CPT2', 'TST', 'MPST', 'HIBCH', 'PYCR2', 'CYP51A1', 'CA8', 'ACOX3', 'CYP2B6', 'POR', 'CYP4F3', 'CYP4F22', 'CPT1C', 'CPT1B', 'HMGCR', 'PYCR3', 'SLC25A21', 'SLC4A4', 'SLC4A5', 'SLC25A12', 'SLC25A13', 'SLC19A1', 'SLC16A1', 'SLC16A7', 'SLC16A3', 'SLC25A20', 'SLC25A10', 'SLC25A11', 'SLC13A3', 'SLC36A1', 'SLC36A2', 'SLC38A3', 'SLC38A5', 'SLC25A26', 'SLC3A2', 'SLC7A6', 'SLC25A15', 'SLC25A2', 'ATP5MGL', 'ATP5MG', 'DMAC2L', 'NSF', 'ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5F1E', 'ATP5PB', 'ATP5MC2', 'ATP5PD', 'ATP5ME', 'ATP5PF', 'ATP5MF', 'ATP5PO', 'ATP5MC1', 'ATP5MC3', 'SLC25A4', 'SLC25A5', 'SLC25A6', 'SLC15A2', 'STARD3', 'SLC25A1', 'SLC26A6', 'SLC25A16', 'UQCRB', 'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRH', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6A2', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A', 'COX7B2', 'COX8C', 'CYC1', 'UQCRQ', 'UQCR11', 'UQCR10', 'CYTB', 'SLC43A2', 'SLC7A11', 'DHODH', 'SLC25A19', 'GPD2', 'SLC37A4', 'SLC5A1', 'SLC5A2', 'SLC5A3', 'SLC5A9', 'SLC5A11', 'SLC5A10', 'SLC5A12', 'SLC25A18', 'SLC25A22', 'SLC17A6', 'SLC17A7', 'SLC17A8', 'SLC6A9', 'AQP1', 'AQP2', 'AQP4', 'AQP5', 'AQP8', 'MIP', 'SLC4A1', 'SLC4A2', 'SLC4A3', 'SLC4A9', 'SLC12A4', 'SLC12A6', 'SLC12A7', 'CYB5D1', 'SLC16A8', 'SLC5A8', 'AQP9', 'SLC22A5', 'SLC6A14', 'SLC7A1', 'SLC7A2', 'SLC7A3', 'SLC38A4', 'NDUFA13', 'NDUFA11', 'NDUFB11', 'NDUFA12', 'TUSC3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'NDUFA1', 'NDUFA10', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFC1', 'NDUFC2', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFV3', 'SLC4A8', 'SLC4A10', 'SLC5A5', 'SLC9A1', 'SLC9A2', 'SLC9A3', 'SLC9A4', 'SLC9A5', 'SLC24A1', 'SLC24A2', 'RHAG', 'RHBG', 'SLC12A1', 'SLC12A2', 'ATP10A', 'SLC25A3', 'SLC17A4', 'SLC17A3', 'SLC17A2', 'SLC17A1', 'SLC34A1', 'SLC34A2', 'ATP8A1', 'SCP2', 'SLC7A10', 'SLC3A1', 'SLC7A7', 'SLC26A1', 'SLC26A2', 'SLC26A7', 'SLC26A8', 'SLC26A9', 'SLC26A11', 'SLC26A3', 'SLC13A4', 'NNT', 'SLC29A1', 'SLC14A1', 'SLC14A2', 'CP', 'SLC7A5', 'SLC7A9', 'SLC43A1', 'SLC6A5', 'SLC4A7', 'SLC20A1', 'SLC20A2', 'SLC34A3', 'SLC6A7', 'COX4I2', 'COX7A1', 'COX7A2L', 'COX6B2', 'COX1', 'COX2', 'COX3', 'ABCD1', 'SLC22A7', 'CFTR', 'SLC22A6', 'SLC22A8', 'SLC22A11', 'SLCO1A2', 'SLCO1B1', 'AQP10', 'AQP3', 'AQP7', 'SLC36A4', 'PAICS', 'ACADL', 'SFXN1']

# Gives a list of the locations of the non matching genes in the total gene list
index_of_non_matching_list = [genes_in_model_symbol.index(i) for i in non_matching_list]

# Gets scRNAseq gene list and data matrix
scRNAseq_gene_list = list(lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["Features"])
scRNAseq_mat = lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["data"].todense()

# Takes the data for t = 0 hours, this is often not needed as scRNAseq data typically has no temporal dependence.
scRNAseq_mat_0h = scRNAseq_mat[:,:702]

# Opens saved results from entrez database query
with open("Data/Input/Symbol_Handling/Model_scRNA_mismatch_dictionary/A549_essential_recon2_2.json","r") as erf:
	annotations = json.load(erf)["DocumentSummarySet"]["DocumentSummary"]

# Generates list of lists of aliases for the non matching genes
# alias_lists index corresponds to index of non_matching_list
alias_lists = [annotations[i]["OtherAliases"].split(", ") for i in range(len(annotations))]

# for each entry of alias lists scans through each alias and checks to see if that alias is in the scRNAseq gene names
# if it is it replaces the symbol from the original with the correct alias
for aliase_list_index in range(len(alias_lists)):
	for i in alias_lists[aliase_list_index]:
		if i in scRNAseq_gene_list:
			loi = index_of_non_matching_list[aliase_list_index]
			genes_in_model_symbol[loi] = i

# Creates new matrix for just the metabolic gene expression
# For each entry in genes_in_model_symbol it tries to find the right row in the scRNAseq data
# If it finds the right row it puts it at the row number corresponding to the entry of the gene in genes_in_model_symbol
# If it doesnt, this means the gene symbol along with its aliases arnt in the scRNAseq list, we assume that
# to mean there are no reads for that data
metabolic_scRNAseq_mat_0h = np.zeros((len(genes_in_model_symbol),np.shape(scRNAseq_mat_0h)[1]))
for i in range(len(genes_in_model_symbol)):
	try:
		metabolic_scRNAseq_mat_0h[i] = scRNAseq_mat_0h[scRNAseq_gene_list.index(genes_in_model_symbol[i])]
	except ValueError:
		continue

name_array = np.reshape(np.array(genes_in_model_symbol),(-1,1))
array_to_save = np.hstack((name_array,metabolic_scRNAseq_mat_0h))
np.savetxt("Data/Intermediate/scRNAseq_for_model/A549_recon_2_2.txt",array_to_save,delimiter=",", fmt="%s")
