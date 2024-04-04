#!/usr/bin/env python
# coding: utf-8

# # Process Fusion Calls from LongGF, JAFFA, FusionSeeker
# # All Samples

# In[ ]:


### STRANDNESS ### 

import pandas as pd
# read in the gtf annotation file into a dataframe named: reference_df

print("reading the reference gtf file for known strandness look up")


ref_gtf_annotation = sys.argv[1] # path/to/gencode.v44.annotation.bed

reference_df = pd.read_csv(ref_gtf_annotation, sep='\t')

# Select the desired columns
reference_df = reference_df.iloc[:, [0, 1, 2, 5]]

# Rename the columns
reference_df.columns = ['Chromosome', 'Start', 'End', 'Strand']

# Display the resulting DataFrame
print(reference_df)


# # Input and GF Criteria Output Paths

# In[ ]:


# PARAMETERS

import sys
import os

sample_number = sys.argv[2] 
print("Gene Fusion Criteria for ", sample_number, "\n")


# GENE FUSION CRITERIA OUTPUT FOLDER
gf_criteria_path = sys.argv[3] # path/to/06_GF_criteria/ include the "/" at the end 
new_folder = gf_criteria_path + sample_number + "/"

# Create the new folder
os.makedirs(new_folder)

print("NEW PATH", new_folder)

artifact_filter_path = new_folder


# LONGGF
longgf_parameters = ['100-10-100', '100-25-100', '25-25-25', '50-25-50', '75-25-75', '10-10-10', '25-10-25', '50-10-50', '75-10-75', '10-25-10']
input_path = sys.argv[4] # path to 05_LongGF/ folder (include the slash at the end)


artifact_filter_longgf_output = artifact_filter_path + sample_number + "/" + sample_number + "_longgf.tsv"



# JAFFA
jaffa_input_path = sys.argv[5] # path to 05_JAFFA folder as input (include the slash at the end) 

input_jaffa = jaffa_input_path + sample_number + "/jaffa_results.csv"


# FUSION SEEKER 
fusionseek_input_path = sys.argv[6] # path to 05_FusionSeeker folder (include the slash at the end) 
input_fusion_seek = fusionseek_input_path + sample_number + "/fusionseeker_out/confident_genefusion.txt"


# input BAM file - sorted by position 
bam_path = sys.argv[7]



# Genes of Interest File
genes_of_interest = sys.argv[8] 



# number of supporting reads cutoff 
num_supporting_reads = sys.argv[9] 


# strandness check output file 
strandness_check_output_file_path = artifact_filter_path + "all_fusions_genes_of_interest_strandness.tsv"

gf_summary_file = artifact_filter_path + "all_fusions_criteria_summary_no_readIDs.tsv"

gf_summary_file_readIDS = artifact_filter_path + "all_fusions_criteria_summary.tsv"



# Literature Cross Reference 
# Specify the file path for the list of known fusion names
mitelman_db_gfs = sys.argv[10] 

chimer_db_gfs = sys.argv[11]

cosmic_db_gfs = sys.argv[12] 



print("input files are saved")


# # Getting Known Strandness

# In[ ]:


def get_known_strandness(consolidated_df, reference_df, artifact_filter_path, name_of_output):
    # Check if consolidated_df is not None and is a pandas DataFrame
    if consolidated_df is None or not isinstance(consolidated_df, pd.DataFrame):
        print("Input DataFrame is empty or not valid. Creating an empty file with the header.")
        
        # Create an empty DataFrame with the header
        empty_df = pd.DataFrame(columns=[
            'Fusion Name', 'Gene 1 Symbol', 'Gene 1 Breakpoint', 'Gene 1 Strand',
            'Gene 2 Symbol', 'Gene 2 Breakpoint', 'Gene 2 Strand',
            'Supporting Reads', 'Read IDs', 'Gene 1 Known Strandness', 'Gene 2 Known Strandness'
        ])

        # Save the empty DataFrame to the specified output path
        empty_file_path = artifact_filter_path + name_of_output + "_strandness.tsv"
        empty_df.to_csv(empty_file_path, sep='\t', index=False)
        
        print(f"Empty DataFrame saved to: {empty_file_path}")
        
        return empty_file_path  # Or any value that indicates a skip


    print("Look up the known strandness from the gtf annotation reference and store in the summary dataframe")

    # Create empty lists to store the strandness values for Gene 1 and Gene 2
    gene1_strandness_values = []
    gene2_strandness_values = []

    print("Wait for loop to finish...")

    # Loop through the values in the "Gene 1 Breakpoint" and "Gene 2 Breakpoint" columns
    for gene1_column, gene2_column in zip(consolidated_df["Gene 1 Breakpoint"], consolidated_df["Gene 2 Breakpoint"]):
        # Gene 1 lookup
        gene1_parse = gene1_column.split(":")
        gene1_chromosome = gene1_parse[0]
        gene1_position = int(gene1_parse[1])

        gene1_filtered_df = reference_df[(reference_df['Chromosome'] == gene1_chromosome) & 
                                          (reference_df['Start'] <= gene1_position) & 
                                          (reference_df['End'] >= gene1_position)]

        if not gene1_filtered_df.empty:
            gene1_strandness = gene1_filtered_df['Strand'].values[0]
        else:
            gene1_strandness = "No Match"

        # Gene 2 lookup
        gene2_parse = gene2_column.split(":")
        gene2_chromosome = gene2_parse[0]
        gene2_position = int(gene2_parse[1])

        gene2_filtered_df = reference_df[(reference_df['Chromosome'] == gene2_chromosome) & 
                                          (reference_df['Start'] <= gene2_position) & 
                                          (reference_df['End'] >= gene2_position)]

        if not gene2_filtered_df.empty:
            gene2_strandness = gene2_filtered_df['Strand'].values[0]
        else:
            gene2_strandness = "No Match"

        # Append the strandness values for both genes
        gene1_strandness_values.append(gene1_strandness)
        gene2_strandness_values.append(gene2_strandness)

    # Add the strandness values as new columns in the input DataFrame
    consolidated_df.insert(6, "Gene 2 Known Strandness", gene2_strandness_values, True)
    consolidated_df.insert(3, "Gene 1 Known Strandness", gene1_strandness_values, True)

    print("Known Strandness of Fusion Genes")
    # print(consolidated_df)

    
    # Save the updated DataFrame to a new tab-delimited file
    consolidated_output_file_with_strandness = artifact_filter_path + name_of_output + "_unique_strandness.tsv"
    consolidated_df.to_csv(consolidated_output_file_with_strandness, sep='\t', index=False)

    print('Updated summary dataframe with known strandness saved to:', consolidated_output_file_with_strandness)
    
    return consolidated_output_file_with_strandness

    
print("function stored")    


# # LongGF 

# In[ ]:


# LongGF

import pandas as pd

def process_longgf_output(longgf_output, summary_file):
    fusion_name_lst = []
    gene1_symbol_lst = [] 
    gene2_symbol_lst = []
    support_reads_lst = []
    gene1_break_lst = []
    gene2_break_lst = []

    read_id_process_lst = []
    read_id_lst = []

    gene1_strand_process_lst = []
    gene1_strand_lst = []

    gene2_strand_process_lst = []
    gene2_strand_lst = []
    
    gene1_position_process_lst = []
    gene1_position_lst = []
    
    gene2_position_process_lst = []
    gene2_position_lst = []
    

    with open(longgf_output, 'r') as file:
        process_lines = False  # Initialize a flag to False
        for line in file:
            if line.startswith("Potential_fusion"):
                process_lines = True
                continue  

            if process_lines:   
                if line.startswith("SumGF"):
                    fusion_name = line.split()[1]
                    fusion_name_lst.append(fusion_name)

                    gene1_symbol = line.split()[1].split(":")[0]
                    gene1_symbol_lst.append(gene1_symbol)

                    gene2_symbol = line.split()[1].split(":")[1]
                    gene2_symbol_lst.append(gene2_symbol)

                    support_reads = line.split()[2]
                    support_reads_lst.append(support_reads)

                    gene1_break = line.split()[3]
                    gene1_break_lst.append(gene1_break)

                    gene2_break = line.split()[4]
                    gene2_break_lst.append(gene2_break)

                # Process the read ids 
                if not line.startswith("GF") and not line.startswith("SumGF"):
                    
                    read_id_name = line.split('/')[1].split(':')[0]
                    read_id_process_lst.append(read_id_name)

                    gene1_strand = line.split("(")[1][0]
                    gene1_strand_process_lst.append(gene1_strand)
                    
                    
                    gene1_position = line.split(")")[0].split(":")[-1]
                    gene1_position_process_lst.append(gene1_position)
                    
                    gene2_position = line.split(")")[1].split(":")[-1].split("/")[1]
                    gene2_position_process_lst.append(gene2_position)
                    
                    gene1_strand = line.split("(")[2][0]
                    gene2_strand_process_lst.append(gene1_strand)

    # Append Read IDs, separate by comma, and store in lists 
    start = 0
    for num in support_reads_lst:
        group_ids = read_id_process_lst[start:start + int(num)]
        read_id_lst.append(",".join(group_ids))

        group_gene1_strand = gene1_strand_process_lst[start:start + int(num)]
        gene1_strand_lst.append(",".join(group_gene1_strand))

        group_gene2_strand = gene2_strand_process_lst[start:start + int(num)]
        gene2_strand_lst.append(",".join(group_gene2_strand)) 
        
        group_gene1_position = gene1_position_process_lst[start:start + int(num)]
        gene1_position_lst.append(",".join(group_gene1_position))
        
        group_gene2_position = gene2_position_process_lst[start:start + int(num)]
        gene2_position_lst.append(",".join(group_gene2_position))

        start += int(num)

    # Combine everything into a DataFrame 
    longgf_df = pd.DataFrame({'Fusion Name': fusion_name_lst, 'Gene 1 Symbol': gene1_symbol_lst, 'Gene 1 Breakpoint': gene1_break_lst, 'Gene 1 Strand': gene1_strand_lst, 
                              'Gene 1 Input Position': gene1_position_lst, 'Gene 2 Symbol': gene2_symbol_lst, 'Gene 2 Breakpoint': gene2_break_lst, 'Gene 2 Strand': gene2_strand_lst, 
                              'Gene 2 Input Position': gene2_position_lst,'Supporting Reads': support_reads_lst, 'Read IDs': read_id_lst})

    # Save as a tab-delimited file 
    longgf_df['GF Program'] = 'LongGF'
    longgf_df.to_csv(summary_file, sep='\t', index=False)


def combine_longgf_outputs(longgf_parameters, input_path, output_file, all_output_filenames):
    for num in longgf_parameters:
        input_file = input_path + "longgf." + num + ".txt"
        process_longgf_output(input_file, output_file)
        all_output_filenames.append(output_file)
        print("Finished processing LongGF output with parameters: ", num)
       
        
def consolidate_outputs(all_output_filenames, sample_number, artifact_filter_path):
    data_frames = [pd.read_csv(file_path, sep='\t') for file_path in all_output_filenames]
    
    # Concatenate the DataFrames
    consolidated_df = pd.concat(data_frames, ignore_index=True)
    
    num_rows_LGF = consolidated_df.shape[0]
    print("Number of fusions in LongGF:", num_rows_LGF)
    
    # Save the consolidated DataFrame to a new tab-delimited file
    summary_output = artifact_filter_path + sample_number + "_all_longgf.tsv"
    consolidated_df.to_csv(summary_output, sep='\t', index=False)


    # Group by 'Fusion Name' and keep the row with the highest supporting reads
    
    # Convert 'Supporting Reads' to numeric if it's not already
    consolidated_df['Supporting Reads'] = consolidated_df['Supporting Reads'].astype(int)

    # Group by multiple columns and keep the row with the highest supporting reads
    grouped_df = consolidated_df.sort_values('Supporting Reads', ascending=False).groupby(['Fusion Name', 'Gene 1 Breakpoint', 'Gene 2 Breakpoint']).first().reset_index()

    
    #previous way 
    #consolidated_df['Supporting Reads'] = consolidated_df['Supporting Reads'].astype(int)
    #consolidated_df = consolidated_df.sort_values('Supporting Reads', ascending=False).drop_duplicates(subset='Fusion Name', keep='first')

    # Save the consolidated DataFrame to a new tab-delimited file
    summary_output = artifact_filter_path + sample_number + "_unique_longgf.tsv"
    grouped_df.to_csv(summary_output, sep='\t', index=False, header=True)

    print()
    print("Finished combining to a summary ")

    # print number of fusions in total from the different parameters, no repeats, filtered on the number of supporting reads
    num_entries = grouped_df.shape[0]
    print("Number of Fusions:", num_entries)
    
    return grouped_df  # Add this line to explicitly return the DataFrame

    
print("functions stored")


# # JAFFA 

# In[ ]:


import os
import pandas as pd

def process_jaffa_output(input_jaffa, output_path):
    # Check if the input file exists
    if not os.path.isfile(input_jaffa):
        # Create an empty DataFrame with the columns
        empty_df = pd.DataFrame(columns=[
            'Fusion Name', 'Gene 1 Symbol', 'Gene 1 Breakpoint', 'Gene 1 Strand',
            'Gene 2 Symbol', 'Gene 2 Breakpoint', 'Gene 2 Strand',
            'Supporting Reads', 'Read IDs'
        ])

        # Save the empty DataFrame to the specified output path
        empty_df.to_csv(artifact_filter_path + sample_number + "_all_JAFFA.tsv", sep='\t', index=False, header=True)
        empty_df.to_csv(output_path, sep='\t', index=False, header=True)

        print()
        print("No input file found. Created an empty file with the header.")
        return

    # Create empty lists for each column
    fusion_name_lst = []
    gene1_symbol_lst = []
    gene1_break_lst = []
    gene1_strand_lst = []
    gene2_symbol_lst = []
    gene2_break_lst = []
    gene2_strand_lst = []
    support_reads_lst = []
    read_id_lst = []

    # Read the confident fusion input, tab-delimited file
    with open(input_jaffa, 'r') as file:
        next(file)  # Skip header

        for line in file:
            parts = line.strip().split(',')

            # Process the data and append to respective lists
            fusion_name_lst.append(parts[1])
            gene1_symbol_lst.append(parts[1].split(":")[0])
            gene1_break_lst.append(parts[2] + ":" + parts[3])
            gene1_strand_lst.append(parts[4])
            gene2_symbol_lst.append(parts[1].split(":")[1])
            gene2_break_lst.append(parts[5] + ":" + parts[6])
            gene2_strand_lst.append(parts[7])
            support_reads_lst.append(parts[10])
            read_id_lst.append(parts[14])

    # Create the DataFrame without specifying column names
    data = {
        'Fusion Name': fusion_name_lst,
        'Gene 1 Symbol': gene1_symbol_lst,
        'Gene 1 Breakpoint': gene1_break_lst,
        'Gene 1 Strand': gene1_strand_lst,
        'Gene 2 Symbol': gene2_symbol_lst,
        'Gene 2 Breakpoint': gene2_break_lst,
        'Gene 2 Strand': gene2_strand_lst,
        'Supporting Reads': support_reads_lst,
        'Read IDs': read_id_lst
    }

    df = pd.DataFrame(data)
    
####    
    num_rows_JAFFA = df.shape[0]
    print("Number of fusions in the JAFFA Output:", num_rows_JAFFA)
    
    # Save the consolidated DataFrame to a new tab-delimited file
    summary_output_JAFFA = artifact_filter_path + sample_number + "_all_JAFFA.tsv"
    df.to_csv(summary_output_JAFFA, sep='\t', index=False, header=True)


    # Group by 'Fusion Name' and keep the row with the highest supporting reads
    
    # Convert 'Supporting Reads' to numeric if it's not already
    df['Supporting Reads'] = df['Supporting Reads'].astype(int)

    # Group by multiple columns and keep the row with the highest supporting reads
    grouped_df = df.sort_values('Supporting Reads', ascending=False).groupby(['Fusion Name', 'Gene 1 Breakpoint', 'Gene 2 Breakpoint']).first().reset_index()

 ####   
    
    # Group by 'Fusion Name' and keep the row with the highest supporting reads
    #df['Supporting Reads'] = df['Supporting Reads'].astype(int)
    #df = df.sort_values('Supporting Reads', ascending=False).drop_duplicates(subset='Fusion Name', keep='first')

    # Save the grouped_df DataFrame to a new tab-delimited file
    grouped_df['GF Program'] = 'JAFFA'
    grouped_df.to_csv(output_path, sep='\t', index=False, header=True)

    print(grouped_df)
    
    print()
    print("Finished combining Jaffa data to a summary")
    
    return grouped_df 

print("function stored")


# # FUSION SEEKER

# In[ ]:


import pandas as pd
import pysam
import os

def process_fusionseeker_output(input_fusion_seek, output_path, bam_path):
    # Check if the input file exists
    print(input_fusion_seek)
    
    if not os.path.isfile(input_fusion_seek):
        # Create an empty DataFrame with the columns
        empty_df = pd.DataFrame(columns=[
            'Fusion Name', 'Gene 1 Symbol', 'Gene 1 Breakpoint', 'Gene 1 Strand',
            'Gene 2 Symbol', 'Gene 2 Breakpoint', 'Gene 2 Strand',
            'Supporting Reads', 'Read IDs'
        ])

        # Save the empty DataFrame to the specified output path
        empty_df.to_csv(output_path, sep='\t', index=False)

        print()
        print("No input file found. Created an empty file with the header.")
        return

    # Read the confident fusion input, tab-delimited file
    with open(input_fusion_seek, 'r') as file:
        # Check if the file has only the header
        first_line = file.readline().strip()
        print(first_line)
        
        second_line = file.readline().strip() if not file.readline().strip() == '' else None

        if second_line is None:
            if first_line.lower().startswith('fusion_name'):
                # Only header is present, create an empty DataFrame with the columns
                empty_df = pd.DataFrame(columns=[
                    'Fusion Name', 'Gene 1 Symbol', 'Gene 1 Breakpoint', 'Gene 1 Strand',
                    'Gene 2 Symbol', 'Gene 2 Breakpoint', 'Gene 2 Strand',
                    'Supporting Reads', 'Read IDs'
                ])

                # Save the empty DataFrame to the specified output path
                empty_df.to_csv(output_path, sep='\t', index=False)

                print()
                print("Input file has only the header. Created an empty file with the header.")
                return



    # Create empty lists for each column
    fusion_name_lst = []
    gene1_symbol_lst = []
    gene1_break_lst = []
    gene1_strand_lst = []
    gene2_symbol_lst = []
    gene2_break_lst = []
    gene2_strand_lst = []
    support_reads_lst = []
    read_id_lst = []

    # Read the confident fusion input, tab-delimited file
    with open(input_fusion_seek, 'r') as file:
        next(file)  # Skip header

        for line in file:
            
            parts = line.strip().split('\t')
            #print(parts)

            # Process the data and append to respective lists
            fusion_name_lst.append(parts[1] + ":" + parts[2])
            gene1_symbol_lst.append(parts[1])
            
            gene1_break_lst.append(parts[4] + ":" + parts[5])
            
            gene1_strand_lst.append('N/A')  # Set to 'N/A' because strand is not available
            
            gene2_symbol_lst.append(parts[2])
            
            gene2_break_lst.append(parts[6] + ":" + parts[7])
            
            gene2_strand_lst.append('N/A')  # Set to 'N/A' because strand is not available
            
            support_reads_lst.append(parts[3])
            
            read_id_lst.append(parts[8])

    # Create the DataFrame without specifying column names
    data = {
        'Fusion Name': fusion_name_lst,
        'Gene 1 Symbol': gene1_symbol_lst,
        'Gene 1 Breakpoint': gene1_break_lst,
        'Gene 1 Strand': gene1_strand_lst,
        'Gene 2 Symbol': gene2_symbol_lst,
        'Gene 2 Breakpoint': gene2_break_lst,
        'Gene 2 Strand': gene2_strand_lst,
        'Supporting Reads': support_reads_lst,
        'Read IDs': read_id_lst
    }

    df = pd.DataFrame(data)
    
#previous
    # Group by 'Fusion Name' and keep the row with the highest supporting reads
    #df['Supporting Reads'] = df['Supporting Reads'].astype(int)
    #df = df.sort_values('Supporting Reads', ascending=False).drop_duplicates(subset='Fusion Name', keep='first')

    # Print the number of fusions in total from the different parameters, no repeats,
    # filtered on the number of supporting reads

    num_entries = df.shape[0]
    print("Number of Fusions in FusionSeeker data:", num_entries)

    # Find the strandness of the genes from the BAM file, based on gene positions and readIDs
    bam_file = pysam.AlignmentFile(bam_path, 'rb')

    # Create empty lists to store strand information
    gene1_strands = []
    gene2_strands = []

    for index, row in df.iterrows():
        gene1_bp_chrm = row['Gene 1 Breakpoint'].split(":")[0]        
        gene1_bp_location = int(row['Gene 1 Breakpoint'].split(":")[1])
        
        gene2_bp_chrm = row['Gene 2 Breakpoint'].split(":")[0]
        gene2_bp_location = int(row['Gene 2 Breakpoint'].split(":")[1])

        
        
        # Loop over the read IDs for "Gene 1"
        gene1_readIDs_loop_over = set(row['Read IDs'].split(","))  # Use a set to ensure unique read IDs for "Gene 1"

        # Loop over the read IDs for "Gene 2"
        gene2_readIDs_loop_over = set(row['Read IDs'].split(","))  # Use a set to ensure unique read IDs for "Gene 2"

        # Create empty lists to store strand information for the current row
        gene1_strands_row = []
        gene2_strands_row = []

        gene1_readIDs_processed = set()  # To keep track of processed read IDs for "Gene 1"
        gene2_readIDs_processed = set()  # To keep track of processed read IDs for "Gene 2"

        for read in bam_file.fetch(gene1_bp_chrm, start=gene1_bp_location-100, end=gene1_bp_location+100):
            if read.qname in gene1_readIDs_loop_over and read.qname not in gene1_readIDs_processed:
                strand = '+' if read.is_forward else '-'
                gene1_strands_row.append(strand)
                gene1_readIDs_processed.add(read.qname)  # Mark this read ID as processed for "Gene 1"

        for read in bam_file.fetch(gene2_bp_chrm, start=gene2_bp_location-100, end=gene2_bp_location+100):
            if read.qname in gene2_readIDs_loop_over and read.qname not in gene2_readIDs_processed:
                strand = '+' if read.is_forward else '-'
                gene2_strands_row.append(strand)
                gene2_readIDs_processed.add(read.qname)  # Mark this read ID as processed for "Gene 2"

        # Join the strand information for the current row and store in the respective columns
        gene1_strands.append(','.join(gene1_strands_row))
        gene2_strands.append(','.join(gene2_strands_row))

    # Store the strand information in the DataFrame
    df['Gene 1 Strand'] = gene1_strands
    df['Gene 2 Strand'] = gene2_strands

    df['GF Program'] = 'FusionSeeker' 
    

    # Save the updated DataFrame to a file
    #df.to_csv(output_path, sep='\t', index=False)
    
    
    
####    
    num_rows_FS = df.shape[0]
    print("Number of fusions in the FusionSeeker Output:", num_rows_FS)
    
    # Save the consolidated DataFrame to a new tab-delimited file
    summary_output_FS = artifact_filter_path + sample_number + "_all_FusionSeeker.tsv"
    df.to_csv(summary_output_FS, sep='\t', index=False)


    # Group by 'Fusion Name' and keep the row with the highest supporting reads
    
    # Convert 'Supporting Reads' to numeric if it's not already
    df['Supporting Reads'] = df['Supporting Reads'].astype(int)

    # Group by multiple columns and keep the row with the highest supporting reads
    grouped_df = df.sort_values('Supporting Reads', ascending=False).groupby(['Fusion Name', 'Gene 1 Breakpoint', 'Gene 2 Breakpoint']).first().reset_index()
    grouped_df.to_csv(output_path, sep='\t', index=False, header=True)

###


    
    
    print("done")

    return grouped_df

print("function stored")



# # LONGGF PROCESSING

# In[ ]:


# LONGGF PROCESSING
name_of_output = "longgf"

output_file = input_path + sample_number + "/" + "longgf_summary.tsv"

print(output_file)
all_output_filenames = []  

print(input_path)

input_path_sub = input_path + sample_number + "/"

# Summary file
combine_longgf_outputs(longgf_parameters, input_path_sub, output_file, all_output_filenames)

# Consolidate
consolidated_df = consolidate_outputs(all_output_filenames, sample_number, artifact_filter_path)


# Call get_known_strandness with the consolidated DataFrame
artifact_filter_longgf_output = get_known_strandness(consolidated_df, reference_df, artifact_filter_path, name_of_output)

# Check if the function returned None, indicating a skip
if artifact_filter_longgf_output is None:
    print("Process skipped.")
else:
    print("Process completed.")
    


# # Strandness Check 

# In[ ]:


# strandness check: create two files: correct and incorrect

# conditions: 

# Gene 1 Gene 2 Gene 1 Ref Gene 2 Ref Status
# + + + + Yes
# - - - - Yes
# + + - - Flip
# - - + + Flip
# + - + - Yes
# - + - + Yes
# + - - + Flip
# - + + - Flip
# + + - + No
# + + + - No
# - - - + No
# - - + - No
# + - + + No
# - + + + No
# + - - - No
# - + - - No



# In[ ]:


# Function to determine status
def determine_status(gene1_strands, gene2_strands, gene1_ref, gene2_ref):
    status_list = []
    
    for gene1, gene2 in zip(gene1_strands, gene2_strands):
        if gene1 == '+' and gene2 == '+' and gene1_ref == '+' and gene2_ref == '+':
            status_list.append('Yes')
        elif gene1 == '-' and gene2 == '-' and gene1_ref == '-' and gene2_ref == '-':
            status_list.append('Yes')
            
        elif gene1 == '+' and gene2 == '+' and gene1_ref == '-' and gene2_ref == '-':
            status_list.append('Flip')
        elif gene1 == '-' and gene2 == '-' and gene1_ref == '+' and gene2_ref == '+':
            status_list.append('Flip')
            
        elif gene1 == '+' and gene2 == '-' and gene1_ref == '+' and gene2_ref == '-':
            status_list.append('Yes')
        elif gene1 == '-' and gene2 == '+' and gene1_ref == '-' and gene2_ref == '+':
            status_list.append('Yes')
            
        elif gene1 == '+' and gene2 == '-' and gene1_ref == '-' and gene2_ref == '+':
            status_list.append('Flip')
        elif gene1 == '-' and gene2 == '+' and gene1_ref == '+' and gene2_ref == '-':
            status_list.append('Flip')
            
        else:
            status_list.append('Invalid combination')
    

    return status_list


print("function stored")


# In[ ]:


# Strandness Processing 

# Specify the path to the LongGF summary strandness output file
longgf_strandness_output = artifact_filter_path + "longgf_unique_strandness.tsv"

# Read the file into a DataFrame
longgf_strandness_df = pd.read_csv(longgf_strandness_output, sep='\t')

# Display the DataFrame
print(longgf_strandness_df)




# Add new columns to the DataFrame
longgf_strandness_df['Gene_Orientation'] = ''
longgf_strandness_df['Invalid_Read_Count'] = ''



# Iterate over the rows of the DataFrame
for index, row in longgf_strandness_df.iterrows():


    gene_order = []
    
    fusion_name = row['Fusion Name']
    print(f"\nProcessing Fusion: {fusion_name}")

    # Reference strands for Gene 1 and Gene 2
    reference_strand_gene1 = row['Gene 1 Known Strandness']
    #print("Reference gene 1: ", reference_strand_gene1)
    reference_strand_gene2 = row['Gene 2 Known Strandness']
    #print("Reference gene 2: ", reference_strand_gene2)

    gene1_strands = row['Gene 1 Strand'].split(',')
    #print(gene1_strands)
    gene2_strands = row['Gene 2 Strand'].split(',')
    #print(gene2_strands)
    
    
    
    
    # Handle NaN values in 'Gene 1 Input Position'
    gene1_pos = row['Gene 1 Input Position']
    print(gene1_pos)
    
    
    if pd.notna(gene1_pos) and isinstance(gene1_pos, str):
        gene1_pos = gene1_pos.split(',')
    else:
        # Set default values when position information is missing
        gene1_pos = ['0-0']

    # Handle NaN values in 'Gene 2 Input Position'
    gene2_pos = row['Gene 2 Input Position']
    if pd.notna(gene2_pos) and isinstance(gene2_pos, str):
        gene2_pos = gene2_pos.split(',')
    else:
        # Set default values when position information is missing
        gene2_pos = ['0-0']
        
    # Initialize 'Gene_Orientation' column for the fusion
    #gene_orientation = None

    
    
    
    # Iterate over the strands and apply determine_status function
    status_list = determine_status(gene1_strands, gene2_strands, reference_strand_gene1, reference_strand_gene2)
    print("STATUS..........", status_list)
    
    # Count occurrences of 'Invalid Combination' in status_list
    invalid_read_count = status_list.count('Invalid combination')
    print("I N V A L I D: ", invalid_read_count)
    
    # Update the 'Invalid_Read_Count' column in the DataFrame
    longgf_strandness_df.at[index, 'Invalid_Read_Count'] = invalid_read_count

    
    #print(set(status_list))
    
    

    
    # Check if 'Invalid combination' is in the set(status_list)
    if 'Invalid combination' in set(status_list):
        
        status_set = set(status_list)

        # Create a new list without 'Invalid combination'
        filtered_status_list = [status for status in status_list if status != 'Invalid combination']
        print()
        print()
        print()

    
    
        filtered_status_set = set(filtered_status_list)


        # If the set is not empty after removal, continue processing
        if filtered_status_set:
            print('removed invalid, continuing with:', filtered_status_set)

            
            
            if filtered_status_set == {'Yes'}:
                # Check positions to determine if gene1 is located at the beginning or end 
                print('only yes')
                
                
                if len(gene1_pos) == 1:
                    pos1_start, pos1_end = map(int, gene1_pos[0].split('-'))
                    pos2_start, pos2_end = map(int, gene2_pos[0].split('-'))
                    
                    if pos1_start == 0 and pos1_end == 0 and pos2_start == 0 and pos2_end == 0:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation

                    
                    elif pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation


                    # Check positions to determine if gene1 is located at the beginning or end 
                    elif pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation   

                        
                        
                elif len(gene1_pos) > 1:
                    
                    #need to find the index of yes on the original status_list 
                    yes_index_num = status_list.index("Yes")
                    print(yes_index_num)
                
                    pos1_start, pos1_end = map(int, gene1_pos[yes_index_num].split('-'))
                    #print("position gene 1 start", pos1_start)

                    pos2_start, pos2_end = map(int, gene2_pos[yes_index_num].split('-'))
                    #print("position gene 2 start", pos2_start)

                    # Check positions to determine if gene1 is located at the beginning or end 
                    if pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation


                    # Check positions to determine if gene1 is located at the beginning or end 
                    if pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation   
                    

                    
                    
                    
                    
            elif filtered_status_set == {'Flip'}:
                # Check positions to determine if gene1 is located at the beginning or end 
                print('only flip')
                
                if len(gene1_pos) == 1:
                    pos1_start, pos1_end = map(int, gene1_pos[0].split('-'))
                    pos2_start, pos2_end = map(int, gene2_pos[0].split('-'))
                    
                    if pos1_start == 0 and pos1_end == 0 and pos2_start == 0 and pos2_end == 0:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation   
                        
                    
                    elif pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation


                    # Check positions to determine if gene1 is located at the beginning or end 
                    elif pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation   


                elif len(gene1_pos) > 1:
                    # Need to find the index of 'Flip' on the original status_list 
                    flip_index_num = status_list.index("Flip")
                    print(flip_index_num)

                    pos1_start, pos1_end = map(int, gene1_pos[flip_index_num].split('-'))
                    pos2_start, pos2_end = map(int, gene2_pos[flip_index_num].split('-'))

                    # Check positions to determine if gene1 is located at the beginning or end 
                    if pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation

                    elif pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation
            
                  

                
                
                    
            elif filtered_status_set == {'Yes', 'Flip'}:
                # Check positions to determine if gene1 is located at the beginning or end 
                print('only yes and flip')
                
                if len(gene1_pos) == 1:
                
                    pos1_start, pos1_end = map(int, gene1_pos[0].split('-'))
                    pos2_start, pos2_end = map(int, gene2_pos[0].split('-'))
                    
                    
                    if pos1_start == 0 and pos1_end == 0 and pos2_start == 0 and pos2_end == 0:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation

                    
                    elif pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation


                    # Check positions to determine if gene1 is located at the beginning or end 
                    elif pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation   

                        
                        
                elif len(gene1_pos) > 1:
                    
                    #need to find the index of yes on the original status_list 
                    yes_index_num = status_list.index("Yes")
                    print(yes_index_num)
                
                    pos1_start, pos1_end = map(int, gene1_pos[yes_index_num].split('-'))
                    #print("position gene 1 start", pos1_start)

                    pos2_start, pos2_end = map(int, gene2_pos[yes_index_num].split('-'))
                    #print("position gene 2 start", pos2_start)

                    # Check positions to determine if gene1 is located at the beginning or end 
                    if pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation


                    # Check positions to determine if gene1 is located at the beginning or end 
                    elif pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation              

            
            

        else:
            # If the set is empty after removal, set the gene_orientation to 'Invalid'
            gene_orientation = 'Invalid'
            longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation
            continue


            
            
            
            
            
            
            
            
 ### 
            
    elif set(status_list) == {'Yes'}:  
        # Check positions to determine if gene1 is located at the beginning or end 
        pos1_start, pos1_end = map(int, gene1_pos[0].split('-'))
        print("position gene 1 start", pos1_start)

        pos2_start, pos2_end = map(int, gene2_pos[0].split('-'))
        print("position gene 2 start", pos2_start)

        # Check positions to determine if gene1 is located at the beginning or end 
        if pos1_start < pos2_start:
            gene_orientation = 'Yes'
            print('Checking orientation', gene_orientation)
            longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation
            
            
            

    elif set(status_list) == {'Flip'}:
        # Check positions to determine if gene1 is located at the beginning or end 
        pos1_start, pos1_end = map(int, gene1_pos[0].split('-'))
        print("position gene 1 start", pos1_start)

        pos2_start, pos2_end = map(int, gene2_pos[0].split('-'))
        print("position gene 2 start", pos2_start)

        # Check if all position values are zero
        if pos1_start == 0 and pos1_end == 0 and pos2_start == 0 and pos2_end == 0:
            gene_orientation = 'Flip'
            print('All position values are zero, setting gene orientation to Yes')
        else:
            # Check positions to determine if gene1 is located at the beginning or end 
            if pos1_start < pos2_start:
                gene_orientation = 'Yes'
                print('Checking orientation', gene_orientation)
            elif pos1_start > pos2_start:
                gene_orientation = 'Flip'
                print('Checking orientation', gene_orientation)

        # Set the gene_orientation in the DataFrame
        print(gene_orientation)
        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation

        
        
        
    elif set(status_list) == {'Yes', 'Flip'}:
        print('only yes and flip HERE')
        
        
        #need to find the index of yes on the original status_list 
        yes_index_num = status_list.index("Yes")
        print(yes_index_num)

                    
        print("THIS", status_list[yes_index_num])

                

        if status_list[yes_index_num] == "Yes":
            print("gene1_pos:", gene1_pos)  # Print the entire list for debugging
            if 0 <= yes_index_num < len(gene1_pos):
                print("gene1_pos[yes_index_num]:", gene1_pos[yes_index_num])  # Print the specific element
                pos1_start, pos1_end = map(int, gene1_pos[yes_index_num].split('-'))
                print("pos1_start:", pos1_start)

                print("gene2_pos:", gene2_pos)  # Add similar checks for gene2_pos
                if 0 <= yes_index_num < len(gene2_pos):
                    pos2_start, pos2_end = map(int, gene2_pos[yes_index_num].split('-'))
                    print("pos2_start:", pos2_start)

                    if pos1_start == 0 and pos1_end == 0 and pos2_start == 0 and pos2_end == 0:
                        gene_orientation = 'Yes'
                    elif pos1_start < pos2_start:
                        gene_orientation = 'Yes'
                        print('Checking orientation', gene_orientation)
                    elif pos1_start > pos2_start:
                        gene_orientation = 'Flip'
                        print('Checking orientation', gene_orientation)
                else:
                    print("yes_index_num is out of range for gene2_pos")
            else:
                print("yes_index_num is out of range for gene1_pos") 



        # Set the gene_orientation in the DataFrame
        print(gene_orientation)
        longgf_strandness_df.at[index, 'Gene_Orientation'] = gene_orientation



    
# Print the filtered DataFrame
print(longgf_strandness_df)

# # Save the filtered DataFrame to a new tab-delimited file
# filtered_orientation_output_file = artifact_filter_path + "longgf_orientation_check.tsv"
# longgf_strandness_df.to_csv(filtered_orientation_output_file, sep='\t', index=False)



# # Make Orientation Flips if needed from LongGF outputs

# In[ ]:


# Condition: If Gene_Orientation is 'Flip', swap Gene 1 and Gene 2 information
flip_condition = longgf_strandness_df['Gene_Orientation'] == 'Flip'

# Swap Gene 1 Symbol with Gene 2 Symbol
longgf_strandness_df.loc[flip_condition, 'Gene 1 Symbol'], longgf_strandness_df.loc[flip_condition, 'Gene 2 Symbol'] = (
    longgf_strandness_df.loc[flip_condition, 'Gene 2 Symbol'], longgf_strandness_df.loc[flip_condition, 'Gene 1 Symbol']
)

# Swap Gene 1 Breakpoint with Gene 2 Breakpoint
longgf_strandness_df.loc[flip_condition, 'Gene 1 Breakpoint'], longgf_strandness_df.loc[flip_condition, 'Gene 2 Breakpoint'] = (
    longgf_strandness_df.loc[flip_condition, 'Gene 2 Breakpoint'], longgf_strandness_df.loc[flip_condition, 'Gene 1 Breakpoint']
)

# Swap Gene 1 Known Strandness with Gene 2 Known Strandness 
longgf_strandness_df.loc[flip_condition, 'Gene 1 Known Strandness'], longgf_strandness_df.loc[flip_condition, 'Gene 2 Known Strandness'] = (
    longgf_strandness_df.loc[flip_condition, 'Gene 2 Known Strandness'], longgf_strandness_df.loc[flip_condition, 'Gene 1 Known Strandness']
)

# Swap Gene 1 Strand with Gene 2 Strand
longgf_strandness_df.loc[flip_condition, 'Gene 1 Strand'], longgf_strandness_df.loc[flip_condition, 'Gene 2 Strand'] = (
    longgf_strandness_df.loc[flip_condition, 'Gene 2 Strand'], longgf_strandness_df.loc[flip_condition, 'Gene 1 Strand']
)


# Swap Gene 1 Input Position with Gene 2 Input Position
longgf_strandness_df.loc[flip_condition, 'Gene 1 Input Position'], longgf_strandness_df.loc[flip_condition, 'Gene 2 Input Position'] = (
    longgf_strandness_df.loc[flip_condition, 'Gene 2 Input Position'], longgf_strandness_df.loc[flip_condition, 'Gene 1 Input Position']
)

# Concatenate Gene 1 Symbol with ":" and Gene 2 Symbol in 'Fusion Name' column
longgf_strandness_df['Fusion Name'] = longgf_strandness_df['Gene 1 Symbol'] + ':' + longgf_strandness_df['Gene 2 Symbol']


# Replace 'Flip' with 'Yes' in 'Gene_Orientation' column since we made the change already 
longgf_strandness_df['Gene_Orientation'] = longgf_strandness_df['Gene_Orientation'].replace('Flip', 'Yes', inplace=False)



print(longgf_strandness_df)


# Save the filtered DataFrame to a new tab-delimited file
filtered_orientation_output_file = artifact_filter_path + "longgf_unique_orientation_check.tsv"
longgf_strandness_df.to_csv(filtered_orientation_output_file, sep='\t', index=False)


print("LongGF Orientation Check Complete")


# # JAFFA PROCESSING

# In[ ]:


# JAFFA PROCESSING
name_of_output = "jaffa"

output_file_J = artifact_filter_path + sample_number + "_jaffa.tsv"

# call function to process jaffa output
consolidated_df_jaffa = process_jaffa_output(input_jaffa, output_file_J)

# Check if the function returned None, indicating a skip
if consolidated_df_jaffa is None:
    print("Process skipped.")
else:
    print("Process completed.")


# # FUSION SEEKER PROCESSING

# In[ ]:


#FUSION SEEKER PROCESSING
name_of_output = "fusionseeker"

output_file_FS = artifact_filter_path + sample_number + "_fusionseeker.tsv"

# call function to process FusionSeeker output
consolidated_df_FS = process_fusionseeker_output(input_fusion_seek, output_file_FS, bam_path)

# Check if the function returned None, indicating a skip
if consolidated_df_FS is None:
    print("Process skipped.")
else:
    print("Process completed.")


# # Combine Fusion Program Summaries into One File 

# In[ ]:


# Concatenate the DataFrames into a complete summary file 
consolidated_LGF_J_FS_df = pd.concat([longgf_strandness_df, consolidated_df_jaffa, consolidated_df_FS], ignore_index=True)

print('Concatenated into a summary dataframe')

# Count the total number of fusions
num_fusions = consolidated_LGF_J_FS_df.shape[0]
print('Number of Fusions in the Consolidated Summary Report:', num_fusions)

# Save the consolidated DataFrame to a new tab-delimited file
consolidated_output_file = artifact_filter_path + "all_fusions_LGF_J_FS.tsv"
consolidated_LGF_J_FS_df.to_csv(consolidated_output_file, sep='\t', index=False)


consolidated_LGF_J_FS_df


# # Check if at least one gene is in the specified file of Genes of Interest

# In[ ]:


# check if these are in your genes_of_interest file

# Read fusion genes from the text file
with open(genes_of_interest, 'r') as file:
    fusion_genes_list = [line.strip() for line in file]

# Add a column 'Genes of Interest' and set values based on whether either gene symbol is in fusion_genes_list
consolidated_LGF_J_FS_df['Genes of Interest'] = consolidated_LGF_J_FS_df.apply(
    lambda row: 'Yes' if row['Gene 1 Symbol'] in fusion_genes_list or row['Gene 2 Symbol'] in fusion_genes_list else 'No',
    axis=1
)
   
print(consolidated_LGF_J_FS_df)

# Save the consolidated DataFrame to a new tab-delimited file
consolidated_output_file = artifact_filter_path + "all_fusions_genes_of_interest.tsv"
consolidated_LGF_J_FS_df.to_csv(consolidated_output_file, sep='\t', index=False)


# Filter DataFrame to keep only rows where 'Genes of Interest' is 'Yes'
gene_of_interest_filtered_df = consolidated_LGF_J_FS_df[consolidated_LGF_J_FS_df['Genes of Interest'] == 'Yes']

# Print the filtered DataFrame
print(gene_of_interest_filtered_df)

# Save the filtered DataFrame to a new tab-delimited file
filtered_output_file = artifact_filter_path + "all_fusions_genes_of_interest_filtered.tsv"
gene_of_interest_filtered_df.to_csv(filtered_output_file, sep='\t', index=False)



# # Filter based on Number of Supporting Reads

# In[ ]:


# Assuming "Supporting Reads" is a numeric column
filtered_df_supporting_reads = gene_of_interest_filtered_df[gene_of_interest_filtered_df['Supporting Reads'] >= num_supporting_reads]

filtered_df_supporting_reads_sorted = filtered_df_supporting_reads.sort_values(by='Supporting Reads', ascending=False)


# Save the filtered DataFrame to a new tab-delimited file
filtered_output_file = artifact_filter_path + "all_fusions_num_support_reads_filtered.tsv"
filtered_df_supporting_reads_sorted.to_csv(filtered_output_file, sep='\t', index=False)

print('filtered by num supporting reads: ', num_supporting_reads)



# In[ ]:


filtered_df_supporting_reads


# # Cross Reference to Literature and Final Summary Outputs

# In[ ]:


## Load dataframe and reference files

reference_file1 = pd.read_csv(mitelman_db_gfs, delimiter='\t')  
reference_file2 = pd.read_csv(chimer_db_gfs, delimiter='\t')  
reference_file3 = pd.read_csv(cosmic_db_gfs, delimiter='\t')  


# Create a function to check if fusion is present in reference files
def check_fusion_in_references(row):
    fusion1 = f"{row['Gene 1 Symbol']}:{row['Gene 2 Symbol']}"
    fusion2 = f"{row['Gene 2 Symbol']}:{row['Gene 1 Symbol']}"

    found_in_file1 = fusion1 in reference_file1.values or fusion2 in reference_file1.values
    found_in_file2 = fusion1 in reference_file2.values or fusion2 in reference_file2.values
    found_in_file3 = fusion1 in reference_file3.values or fusion2 in reference_file3.values

    result = []

    if found_in_file1:
        result.append(fusion1 if fusion1 in reference_file1.values else fusion2)
    if found_in_file2:
        result.append(fusion1 if fusion1 in reference_file2.values else fusion2)
    if found_in_file3:
        result.append(fusion1 if fusion1 in reference_file3.values else fusion2)

    return ', '.join(result) if result else ''


filtered_df_supporting_reads_sorted['Literature Cross Reference'] = filtered_df_supporting_reads_sorted.apply(check_fusion_in_references, axis=1)

# Filter rows where 'Literature Cross Reference' is not empty
literature_cross_reference_df_sorted = filtered_df_supporting_reads_sorted[filtered_df_supporting_reads_sorted['Literature Cross Reference'] != '']

# Sort by Number of Supporting Reads 
literature_cross_reference_df_sorted = literature_cross_reference_df_sorted.sort_values(by='Supporting Reads', ascending=False)

# Save the filtered DataFrame to a new tab-delimited file
output_file = artifact_filter_path + "literature_cross_referenced_output.tsv"
literature_cross_reference_df_sorted.to_csv(output_file, sep='\t', index=False)

# Display the resulting DataFrame
print(literature_cross_reference_df_sorted)









# # Prioritize by known in literature and then another by number of supporting reads

# In[ ]:


# Summary File #1

# Select the desired columns for the final output
selected_columns = ["Fusion Name", "Gene 1 Symbol", "Gene 1 Breakpoint", "Gene 2 Symbol", "Gene 2 Breakpoint", "Supporting Reads", "Read IDs", "Literature Cross Reference", 'GF Program']

# Concatenate the selected columns from both DataFrames
final_df_with_readIDs = filtered_df_supporting_reads_sorted[selected_columns]

# Save the final DataFrame to a new tab-delimited file
final_output_file = artifact_filter_path + "final_output.tsv"
final_df_with_readIDs.to_csv(final_output_file, sep='\t', index=False)

# Display the final DataFrame
print(final_df_with_readIDs)


# In[ ]:


# Without the ReadIDs for Analysis of the results

# Summary File #2 - no read ids

# Select the desired columns for the final output
selected_columns_no_readIDs = ["Fusion Name", "Gene 1 Symbol", "Gene 1 Breakpoint", "Gene 2 Symbol", "Gene 2 Breakpoint", "Supporting Reads", "Literature Cross Reference", 'GF Program']

# Concatenate the selected columns from both DataFrames
final_df_NO_readIDs = filtered_df_supporting_reads_sorted[selected_columns_no_readIDs]

# Save the final DataFrame to a new tab-delimited file
final_output_file_no_readIDs = artifact_filter_path + "final_output_no_readIDs.tsv"
final_df_NO_readIDs.to_csv(final_output_file_no_readIDs, sep='\t', index=False)

# Display the final DataFrame
print(final_df_NO_readIDs)



print("Gene Fusion Criteria Completed")


