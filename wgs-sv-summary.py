import vcfpy
import matplotlib.pyplot as plt
import pandas as pd
import re
import os
import shutil

input_vcf='/Users/ash/Downloads/p3l.vcf'
output_dir='/Users/ash/Downloads'
directory_path = os.path.join(output_dir, 'SV_summary')

# Check if the directory exists
if os.path.exists(directory_path):
    # Force remove the directory and its contents
    shutil.rmtree(directory_path)
    print(f"Removed existing directory: {directory_path}")
else:
    print(f"Directory does not exist: {directory_path}")

# Create the directory
os.makedirs(directory_path)
print(f"Created directory: {directory_path}")

#Read the sniffles output VCF file
reader = vcfpy.Reader.from_path(input_vcf)

# Step 2: Create a list to hold the records
records = []

# Step 3: Iterate over each record in the VCF and collect relevant data
for record in reader:
    records.append({
        'CHROM': record.CHROM,
        'POS': record.POS,
        'ID': record.ID,
        'REF': record.REF,
        'ALT': record.ALT,
        'QUAL': record.QUAL,
        'FILTER': record.FILTER,
        'INFO': record.INFO
    })
#Convert the list of records to a DataFrame
raw_sv_df = pd.DataFrame(records)

###Extracting info
def extract_info(info):
    info_dict = {}
    for item in info.split(','):
        key_value = item.split('=')
        if len(key_value) == 2:
            info_dict[key_value[0]] = key_value[1]
    return info_dict


raw_sv_df['SVTYPE'] = raw_sv_df['INFO'].apply(lambda x: x.get('SVTYPE'))
raw_sv_df['SVLEN'] = raw_sv_df['INFO'].apply(lambda x: int(x.get('SVLEN', 0)))

#Removing unknown chromosomes
raw_sv_df = raw_sv_df[~raw_sv_df['CHROM'].str.startswith('chrUn')]
##filtering out random chromsoomes
raw_sv_df['CHROM']=raw_sv_df['CHROM'].str.replace('chr', '', regex=False)
raw_sv_df=raw_sv_df[raw_sv_df['CHROM'].isin(map(str, range(1, 23)))]

###creating separate df for SV and chr rearrnagements
sv_df = raw_sv_df[raw_sv_df['SVTYPE'].isin(['INS', 'DEL', 'DUP','INV'])]
chrom_rearrangs = raw_sv_df[raw_sv_df['SVTYPE'] == 'BND']

#function for categorizing into bins
def categorize_svlens(svlen):
    if 1000 <= svlen < 2500:
        return '1000-2500'
    elif 2500 <= svlen < 5000:
        return '2500-5000'
    elif svlen >= 5000:
        return '>5000'
    else:
        return '<1000'



################ PLotting raw data

# Raw data: Apply categorization
sv_df['Len_Category'] = sv_df['SVLEN'].abs().apply(categorize_svlens)

# Count occurrences of each SVTYPE for each CHROM
sv_df['CHROM']=sv_df['CHROM'].astype(int)

# Raw data: Group by SVTYPE and Category to get counts
raw_counts_df = sv_df.groupby(['SVTYPE', 'Len_Category']).size().reset_index(name='Count')


############## PLOTTING RAW DATA ##################

# Pivot the DataFrame to get counts for each SVTYPE in separate columns
pivot_df = raw_counts_df.pivot(index='Len_Category', columns='SVTYPE', values='Count').fillna(0)

# Create a figure with subplots
fig, axs = plt.subplots(2, 2, figsize=(20, 20))  # Create a 2x2 grid of subplots

# Plot each SVTYPE in a separate subplot
axs[0, 0].bar(pivot_df.index, pivot_df['DEL'], color='steelblue')
axs[0, 0].set_title('DEL Counts')
axs[0, 0].set_xlabel('Length Category')
axs[0, 0].set_ylabel('Count')
axs[0, 0].tick_params(axis='x', rotation=45)

    
axs[0, 1].bar(pivot_df.index, pivot_df['DUP'], color='orange')
axs[0, 1].set_title('DUP Counts')
axs[0, 1].set_xlabel('Length Category')
axs[0, 1].set_ylabel('Count')
axs[0, 1].tick_params(axis='x', rotation=45)

axs[1, 0].bar(pivot_df.index, pivot_df['INS'], color='green')
axs[1, 0].set_title('INS Counts')
axs[1, 0].set_xlabel('Length Category')
axs[1, 0].set_ylabel('Count')
axs[1, 0].tick_params(axis='x', rotation=45)

axs[1, 1].bar(pivot_df.index, pivot_df['INV'], color='red')
axs[1, 1].set_title('INV Counts')
axs[1, 1].set_xlabel('Length Category')
axs[1, 1].set_ylabel('Count')
axs[1, 1].tick_params(axis='x', rotation=45)

for i in ['0,0','0,1','1,0','1,1']:
    a=int(i.split(',')[0])
    b=int(i.split(',')[1])
    for bar in axs[a,b].patches:
        yval = bar.get_height()
        axs[a,b].text(bar.get_x() + bar.get_width() / 2, yval, int(yval), ha='center')
 
#saving image in output dir
plt.savefig(output_dir+'/SV_summary/raw_sv_dist.png')
        
############## STACKED bar chart for SVTYPES

raw_chrom_count_df = sv_df.groupby(['CHROM', 'SVTYPE']).size().unstack(fill_value=0)

# Plotting
raw_chrom_count_df.plot(kind='bar', stacked=True, figsize=(10, 6))


plt.title('Structural Variant Types by Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Count of SVTYPE')
plt.xticks(rotation=0)  # Rotate x-axis labels for better readability
plt.legend(title='SVTYPE')
plt.tight_layout()  # Adjust layout to prevent overlap

#saving image in output dir
plt.savefig(output_dir+'/SV_summary/raw_chrom_dist.png')

####################################################################################
#######################    Filtration and summarizing                 ###################
####################################################################################
#####################################################

#Filtering: Filter for 'PASS' in FILTER column and INFO starting with 'PRECISE'
filtered_df = sv_df[
    (sv_df['FILTER'].apply(lambda x: 'PASS' in x)) & 
    (sv_df['INFO'].apply(lambda x: any(k.startswith('PRECISE') for k in x.keys())))
]

# Step 5: Group by SVTYPE and Category to get counts
filter_counts_df = filtered_df.groupby(['SVTYPE', 'Len_Category']).size().reset_index(name='Count')

############################################################
############ PLOTTING filtered data

# Pivot the DataFrame to get counts for each SVTYPE in separate columns
pivot_filter_df = filter_counts_df.pivot(index='Len_Category', columns='SVTYPE', values='Count').fillna(0).astype(int)

# Create a figure with subplots
fig, axs = plt.subplots(2, 2, figsize=(20, 20))  # Create a 2x2 grid of subplots

# Plot each SVTYPE in a separate subplot
axs[0, 0].bar(pivot_filter_df.index, pivot_filter_df['DEL'], color='steelblue')
axs[0, 0].set_title('DEL Counts')
axs[0, 0].set_xlabel('Length Category')
axs[0, 0].set_ylabel('Count')
axs[0, 0].tick_params(axis='x', rotation=45)

    
axs[0, 1].bar(pivot_filter_df.index, pivot_filter_df['DUP'], color='orange')
axs[0, 1].set_title('DUP Counts')
axs[0, 1].set_xlabel('Length Category')
axs[0, 1].set_ylabel('Count')
axs[0, 1].tick_params(axis='x', rotation=45)

axs[1, 0].bar(pivot_filter_df.index, pivot_filter_df['INS'], color='green')
axs[1, 0].set_title('INS Counts')
axs[1, 0].set_xlabel('Length Category')
axs[1, 0].set_ylabel('Count')
axs[1, 0].tick_params(axis='x', rotation=45)

axs[1, 1].bar(pivot_filter_df.index,pivot_filter_df['INV'], color='red')
axs[1, 1].set_title('INV Counts')
axs[1, 1].set_xlabel('Length Category')
axs[1, 1].set_ylabel('Count')
axs[1, 1].tick_params(axis='x', rotation=45)

for i in ['0,0','0,1','1,0','1,1']:
    a=int(i.split(',')[0])
    b=int(i.split(',')[1])
    for bar in axs[a,b].patches:
        yval = bar.get_height()
        axs[a,b].text(bar.get_x() + bar.get_width() / 2, yval, int(yval), ha='center')
      
#saving image in output dir
plt.savefig(output_dir+'/SV_summary/filter_sv_dist.png')

############## STACKED bar chart for SVTYPES

# Count occurrences of each SVTYPE for each CHROM

filtered_df['CHROM']=filtered_df['CHROM'].astype(int)
filter_chrom_count_df = filtered_df.groupby(['CHROM', 'SVTYPE']).size().unstack(fill_value=0)

# Plotting
filter_chrom_count_df.plot(kind='bar', stacked=True, figsize=(10, 6))


plt.title('Structural Variant Types by Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Count of SVTYPE')
plt.xticks(rotation=0)  # Rotate x-axis labels for better readability
plt.legend(title='SVTYPE')
plt.tight_layout()  # Adjust layout to prevent overlap

#saving image in output dir
plt.savefig(output_dir+'/SV_summary/filter_chrom_dist.png')


########################################################################################################
########## Getting a sheet for translocations ###########
########################################################################################################

##### making a df for plotting
raw_chr_rearrngs_granges=pd.DataFrame()

raw_chr_rearrngs_granges['chromosome1']=chrom_rearrangs['CHROM']
raw_chr_rearrngs_granges['chromosome2']=chrom_rearrangs[chrom_rearrangs['SVTYPE'] == 'BND']['ALT'].apply(lambda x: re.search(r'chr(\d+)', x[0].mate_chrom).group(1))
raw_chr_rearrngs_granges['bp1']=chrom_rearrangs[chrom_rearrangs['SVTYPE'] == 'BND']['POS']
raw_chr_rearrngs_granges['bp2']=chrom_rearrangs[chrom_rearrangs['SVTYPE'] == 'BND']['ALT'].apply(lambda x: x[0].mate_pos)
raw_chr_rearrngs_granges['rearrangements'] = raw_chr_rearrngs_granges.apply(
    lambda row: 'interchromosomal' if row['chromosome1'] != row['chromosome2'] else 'intrachromosomal',
    axis=1
)



############# Filtering the chr_rearrangs


#Filtering: Filter for 'PASS' in FILTER column and INFO starting with 'PRECISE'
filter_chrom_rearrangs = chrom_rearrangs[
    (chrom_rearrangs['FILTER'].apply(lambda x: 'PASS' in x)) & 
    (chrom_rearrangs['INFO'].apply(lambda x: any(k.startswith('PRECISE') for k in x.keys())))
]

##### making a df for plotting
filter_chr_rearrngs_granges=pd.DataFrame()

filter_chr_rearrngs_granges['chromosome1']=filter_chrom_rearrangs['CHROM']
filter_chr_rearrngs_granges['chromosome2']=filter_chrom_rearrangs[filter_chrom_rearrangs['SVTYPE'] == 'BND']['ALT'].apply(lambda x: re.search(r'chr(\d+)', x[0].mate_chrom).group(1))
filter_chr_rearrngs_granges['bp1']=filter_chrom_rearrangs[filter_chrom_rearrangs['SVTYPE'] == 'BND']['POS']
filter_chr_rearrngs_granges['bp2']=filter_chrom_rearrangs[filter_chrom_rearrangs['SVTYPE'] == 'BND']['ALT'].apply(lambda x: x[0].mate_pos)
filter_chr_rearrngs_granges['rearrangements'] = filter_chr_rearrngs_granges.apply(
    lambda row: 'interchromosomal' if row['chromosome1'] != row['chromosome2'] else 'intrachromosomal',
    axis=1
)


################################################################
#### SV - deletions, duplications

raw_sv_plot_data=pd.DataFrame()

raw_sv_plot_data['chrom']= sv_df[sv_df['SVTYPE'] != 'BND']['CHROM']
raw_sv_plot_data['SVTYPE']=sv_df[sv_df['SVTYPE'] != 'BND']['SVTYPE']
raw_sv_plot_data['SVLEN']=sv_df[sv_df['SVTYPE'] != 'BND']['SVLEN'].abs()
raw_sv_plot_data['POS']=sv_df[sv_df['SVTYPE'] != 'BND']['POS']
raw_sv_plot_data['END']=[info.get('END') for info in sv_df[sv_df['SVTYPE'] != 'BND']['INFO']]
raw_sv_plot_data['Len_Category'] =sv_df[sv_df['SVTYPE'] != 'BND']['SVLEN'].abs().apply(categorize_svlens)


############### Filtered SV data
filter_sv_plot_data=pd.DataFrame()

filter_sv_plot_data['chrom']= filtered_df[filtered_df['SVTYPE'] != 'BND']['CHROM']
filter_sv_plot_data['SVTYPE']=filtered_df[filtered_df['SVTYPE'] != 'BND']['SVTYPE']
filter_sv_plot_data['SVLEN']=filtered_df[filtered_df['SVTYPE'] != 'BND']['SVLEN'].abs()
filter_sv_plot_data['POS']=filtered_df[filtered_df['SVTYPE'] != 'BND']['POS']
filter_sv_plot_data['END']=[info.get('END') for info in filtered_df[filtered_df['SVTYPE'] != 'BND']['INFO']]
filter_sv_plot_data['Len_Category'] =filtered_df[filtered_df['SVTYPE'] != 'BND']['SVLEN'].abs().apply(categorize_svlens)


########### Summarizing it in excel sheet ##########
sheet_name=output_dir+'/SV_summary/SV-summary.xlsx'

#filter_sv_plot_data.to_csv(output_dir+'/SV_summary/filter_sv.csv', index=False)
#raw_chr_rearrngs_granges.to_csv(output_dir+'/SV_summary/raw_chr_rearrangs.csv', index=False)
#filter_chr_rearrngs_granges.to_csv(output_dir+'/SV_summary/filter_chr_rearrangs.csv', index=False)
#aw_sv_plot_data.to_csv(output_dir+'/SV_summary/raw_sv.csv', index=False)


# Create an Excel writer object and save DataFrames to sheets
with pd.ExcelWriter(sheet_name, engine='openpyxl') as writer:
    raw_counts_df.to_excel(writer, sheet_name='raw_SV_counts', index=False)
    raw_chrom_count_df.to_excel(writer, sheet_name='raw_chrom_counts', index=True)
    raw_sv_plot_data.to_excel(writer, sheet_name='raw_sv', index=False)  
    filter_sv_plot_data.to_excel(writer, sheet_name='filtered_sv', index=False)
    raw_chr_rearrngs_granges.to_excel(writer, sheet_name='raw_chr_rearrangs', index=False)
    filter_chr_rearrngs_granges.to_excel(writer, sheet_name='filter_chr_rearrangs', index=False)
    
    