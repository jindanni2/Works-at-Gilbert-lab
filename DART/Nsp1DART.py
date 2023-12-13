#!/usr/bin/env python3
import pandas as pd
import re
import numpy as np

def pileup(sample_files):
# Extracting "Plus_reads" counts from bbmap coverage.txt files, for every sample
# ID indices:
#		Variant pool: 0 - 11403 (excluded from pile-up)
#		Short and long: 11404 - 51594
#		SARS-CoV-2 UTRs: 51595 - 51598
#		300-nt spike-in control: 51599 - 51698	(excluded from pile-up)
	path = "read_counts_demultiplexed/"
	first_file = True
	df = pd.DataFrame()
	
	for sample in sample_files:
		data = pd.read_csv(f"{path}{sample_files[sample]}",sep = '\t')
		if first_file:
			df["#ID"] = data["#ID"]
			first_file = False
		df[sample] = data["Plus_reads"]
	df_sub = df.iloc[11404:51599,]
	df_sub.to_csv("Raw_reads_novar.tsv",sep = "\t", index = False)
	#df.to_csv("Raw_reads_noind.tsv",sep = "\t", index = False)

def total_mreads_ctrl(sample_fileheads):
# Extracting uniquelly mapped read numbers (convert to million reads) for each sample from STAR Log files. Also extract total spike-in ctrl reads and calculate CPM for the ctrl of each sample. Returns two values: a dictionary of uniquelly mapped reads numbers in million reads; a dictionary of ctrl CPM for each sample.

	UMRN_dict = {}
	ctrlCPM_dict = {}

	for sample in sample_fileheads:
		
		# Extract UMRN and convert to million reads
		UMRN_dict[sample] = 0
		for bar in range(1,6):
			STAR_file = "STARlog/"+sample_fileheads[sample]+f"__barcode_Bar{bar}_1.fastq.gzLog.final.out"
			fh = open(STAR_file,"r")
			for line in fh:
				line = line.lstrip()
				if line.startswith("Uniquely mapped reads number"):
					line = line.rstrip()
					UMRN_dict[sample] += float(line.split("\t")[-1])/1000000
					break
			fh.close()

		# Extract and calculate spike-in control CPM
		ctrlCPM_dict[sample] = 0
		if sample.startswith("IVT"):
			ctrlCPM_dict[sample] = np.nan
		else:
			ctrl_file = "read_counts_demultiplexed/"+sample_fileheads[sample]+"__barcode_Bar3_1.fastq.gz_coverage.txt"
			raw_data = pd.read_csv(ctrl_file, sep = "\t")
			ctrl_data = raw_data[raw_data["#ID"].str.startswith("300nt_Standard")]
			total_ctrl_reads = ctrl_data["Plus_reads"].sum()
			total_ctrl_CPM = total_ctrl_reads/UMRN_dict[sample]
			ctrlCPM_dict[sample] = total_ctrl_CPM

	return(UMRN_dict,ctrlCPM_dict)


def raw_to_CPM(raw_file,UMRN_dict):
# Take raw reads from piled-up tsv file and a dictionary containing uniquelly mapped read number (million reads) for all samples. Write CPM values to a new tsv file.

	raw_df = pd.read_csv(raw_file,sep = "\t")
	CPM_df = pd.DataFrame()
	#print(raw_df)
	CPM_df["#ID"] = raw_df["#ID"]
	
	for sample in UMRN_dict:
		CPM_df[sample] = raw_df[sample]/UMRN_dict[sample]
	
	#print(CPM_df)
	CPM_df.to_csv("CPM.tsv", sep = "\t", index = False)

def Norm_ctrl(ctrl_dict):
# Takes a dictionary containing ctrl spike-in CPM values for each sample library. Normalize values to the first sample.

	norm_ctrl_dict = {}
	first_sample = True

	for sample in ctrl_dict:
		if first_sample:
			norm_ctrl_dict[sample] = 1.0
			denominator = ctrl_dict[sample]
			first_sample = False
		else:
			norm_ctrl_dict[sample] = ctrl_dict[sample]/denominator

	return(norm_ctrl_dict)

def RRS_Retention(CPM_file,norm_dict):
# Takes a datasheet containing CPM values for all monosome, input, and IVT samples, and a dictionary of normalization factor from spike-in control of each sample. Write RRS, log2RRS, Retention, log2Retention to new tsv files.

	# Initiate new data frames
	RRS_df = pd.DataFrame()
	Retention_df = pd.DataFrame()
	log2RRS_df = pd.DataFrame()
	log2Retention_df = pd.DataFrame()
	
	# Read CPM files, adding ID columns to new data frames	
	CPM_df = pd.read_csv(CPM_file, sep = "\t")
	RRS_df["#ID"] = CPM_df["#ID"]
	Retention_df["#ID"] = CPM_df["#ID"]
	log2RRS_df["#ID"] = CPM_df["#ID"]
	log2Retention_df["#ID"] = CPM_df["#ID"]
	
	# Pair sample with correct IVT
	corresponding_IVT = {}
	for sample in norm_dict:
		if sample.endswith("Rep-1") or sample.endswith("Rep-2") or sample.endswith("Rep-3"):
			corresponding_IVT[sample] = "IVT_bar1"
		elif sample.endswith("Rep-4") or sample.endswith("Rep-5") or sample.endswith("Rep-6"):
			corresponding_IVT[sample] = "IVT_bar2"
	
	#print(corresponding_IVT)
	# Calculate RRS and Retention，and their logarithms
	for sample in norm_dict:
		if sample.startswith("monosome"):
			RRS_df[sample[9:]] = CPM_df[sample]/norm_dict[sample]/CPM_df[corresponding_IVT[sample]]
			log2RRS_df[sample[9:]] = np.log2(RRS_df[sample[9:]])
		elif sample.startswith("input"):
			Retention_df[sample[6:]] = CPM_df[sample]/norm_dict[sample]/CPM_df[corresponding_IVT[sample]]
			log2Retention_df[sample[6:]] = np.log2(Retention_df[sample[6:]])

	# Write tsv files
	RRS_df.to_csv("RRS_norm_1.tsv", sep = "\t", index = False)
	log2RRS_df.to_csv("log2RRS_norm_1.tsv", sep = "\t", index = False)
	Retention_df.to_csv("Retention_norm_1.tsv", sep = "\t", index = False)
	log2Retention_df.to_csv("log2Retention_norm_1.tsv", sep = "\t", index = False)


def RRS_prime():
# Calculates RRS_prime for each UTR. RRS_prime is defined as (norm CPM in monosome)/(norm CPM in translation input), which is equavalent to RRS/Retention.

	# Read RRS and Retention files, initiate new data frame, fill in UTR IDs
	RRS_df = pd.read_csv("RRS_norm_1.tsv", sep = "\t")
	Retention_df = pd.read_csv("Retention_norm_1.tsv", sep = "\t")

	RRSprime_df = pd.DataFrame()
	RRSprime_df["#ID"] = RRS_df["#ID"]

	# Extract sample names from RRS dataframe
	sample_list = []
	for col in RRS_df.columns:
		if col != "#ID":
			sample_list.append(col)

	# Calculate RRSprime for each sample
	for sample in sample_list:
		RRSprime_df[sample] = RRS_df[sample]/Retention_df[sample]
	
	# Write tsv file
	RRSprime_df.to_csv("RRSprime.tsv", sep = "\t", index = False)

####### YOU ARE HERE ###########

def Read_filter():
	pass

def main():
	sample_files = {"monosome-Rep-2":"80S-Rep-2-5__barcode_Bar1_1.fastq.gz_coverage.txt",
"monosome-Rep-3":"80S-Rep-3-6__barcode_Bar1_1.fastq.gz_coverage.txt",
"monosome-Rep-5":"80S-Rep-2-5__barcode_Bar2_1.fastq.gz_coverage.txt",
"monosome-Rep-6":"80S-Rep-3-6__barcode_Bar2_1.fastq.gz_coverage.txt",
"monosome_Nsp1-Rep-1":"80SNsp1-Rep-1-4__barcode_Bar1_1.fastq.gz_coverage.txt",
"monosome_Nsp1-Rep-2":"80SNsp1-Rep-2-5__barcode_Bar1_1.fastq.gz_coverage.txt",
"monosome_Nsp1-Rep-3":"80SNsp1-Rep-3-6__barcode_Bar1_1.fastq.gz_coverage.txt",
"monosome_Nsp1-Rep-4":"80SNsp1-Rep-1-4__barcode_Bar2_1.fastq.gz_coverage.txt",
"monosome_Nsp1-Rep-5":"80SNsp1-Rep-2-5__barcode_Bar2_1.fastq.gz_coverage.txt",
"monosome_Nsp1-Rep-6":"80SNsp1-Rep-3-6__barcode_Bar2_1.fastq.gz_coverage.txt",
"input-Rep-1":"Input-Rep1__barcode_Bar1_1.fastq.gz_coverage.txt",
"input-Rep-2":"Input-Rep2__barcode_Bar1_1.fastq.gz_coverage.txt",
"input-Rep-3":"Input-Rep3__barcode_Bar1_1.fastq.gz_coverage.txt",
"input-Rep-4":"Input-Rep4__barcode_Bar2_1.fastq.gz_coverage.txt",
"input-Rep-5":"Input-Rep5__barcode_Bar2_1.fastq.gz_coverage.txt",
"input-Rep-6":"Input-Rep6__barcode_Bar2_1.fastq.gz_coverage.txt",
"input_Nsp1-Rep-1":"Input-Nsp1-Rep1__barcode_Bar1_1.fastq.gz_coverage.txt",
"input_Nsp1-Rep-2":"Input-Nsp1-Rep2__barcode_Bar1_1.fastq.gz_coverage.txt",
"input_Nsp1-Rep-3":"Input-Nsp1-Rep3__barcode_Bar1_1.fastq.gz_coverage.txt",
"input_Nsp1-Rep-4":"Input-Nsp1-Rep4__barcode_Bar2_1.fastq.gz_coverage.txt",
"input_Nsp1-Rep-5":"Input-Nsp1-Rep5__barcode_Bar2_1.fastq.gz_coverage.txt",
"input_Nsp1-Rep-6":"Input-Nsp1-Rep6__barcode_Bar2_1.fastq.gz_coverage.txt",
"IVT_bar1":"IVTbar-1-4__barcode_Bar1_1.fastq.gz_coverage.txt",
"IVT_bar2":"IVTbar-2-5__barcode_Bar2_1.fastq.gz_coverage.txt"}

	#ctrl_files = {}
	#for sample in sample_files:
	#	if not sample.startswith("IVT"):
	#		ctrl_files[sample] = re.sub(r"Bar\d_1","Bar3_1",sample_files[sample])
	#sample_list = list(sample_files.keys())
	
#	pileup(sample_files)
	
	sample_fileheads = {"monosome-Rep-2":"80S-Rep-2-5",
"monosome-Rep-3":"80S-Rep-3-6",
"monosome-Rep-5":"80S-Rep-2-5",
"monosome-Rep-6":"80S-Rep-3-6",
"monosome_Nsp1-Rep-1":"80SNsp1-Rep-1-4",
"monosome_Nsp1-Rep-2":"80SNsp1-Rep-2-5",
"monosome_Nsp1-Rep-3":"80SNsp1-Rep-3-6",
"monosome_Nsp1-Rep-4":"80SNsp1-Rep-1-4",
"monosome_Nsp1-Rep-5":"80SNsp1-Rep-2-5",
"monosome_Nsp1-Rep-6":"80SNsp1-Rep-3-6",
"input-Rep-1":"Input-Rep1",
"input-Rep-2":"Input-Rep2",
"input-Rep-3":"Input-Rep3",
"input-Rep-4":"Input-Rep4",
"input-Rep-5":"Input-Rep5",
"input-Rep-6":"Input-Rep6",
"input_Nsp1-Rep-1":"Input-Nsp1-Rep1",
"input_Nsp1-Rep-2":"Input-Nsp1-Rep2",
"input_Nsp1-Rep-3":"Input-Nsp1-Rep3",
"input_Nsp1-Rep-4":"Input-Nsp1-Rep4",
"input_Nsp1-Rep-5":"Input-Nsp1-Rep5",
"input_Nsp1-Rep-6":"Input-Nsp1-Rep6",
"IVT_bar1":"IVTbar-1-4",
"IVT_bar2":"IVTbar-2-5"}
	UMRN_dict, ctrlCPM_dict = total_mreads_ctrl(sample_fileheads)
	#print(ctrlCPM_dict)

	#raw_to_CPM("Raw_reads_novar.tsv", UMRN_dict)
	norm_ctrl_dict = Norm_ctrl(ctrlCPM_dict)
	RRS_Retention("CPM.tsv", norm_ctrl_dict)

#main()

RRS_prime()		
