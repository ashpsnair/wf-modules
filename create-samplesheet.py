import pandas as pd
import re
import argparse
import os

# Function to create the sample sheet
def create_sample_sheet(input_file):
    # Read file names from the provided input file
    with open(input_file, 'r') as f:
        file_names = [line.strip() for line in f if line.strip()]  # Read lines and strip whitespace

    # Initialize an empty list to hold the data for the DataFrame
    data = []

    # Sex mapping based on patient information
    sex_mapping = {
        "T01": "XY", "N01": "XY",
        "T02": "XY", "N02": "XY",
        "T03": "XX", "N03": "XX",
        "T04": "XY", "N04": "XY",
        "T05": "XY", "N05": "XY",
        "T06": "XX", "N06": "XX",
        "T07": "XY", "N07": "XY",
        "T08": "XY", "N08": "XY",
        "T09": "XY", "N09": "XY",
        "T10": "XY", "N10": "XY",
        "T11": "XY", "N11": "XY",
        "T12": "XY", "N12": "XY",
        "T13": "XY", "N13": "XY",
        "T14": "XY", "N14": "XY"
    }

    # Process each file name
    for file_name in file_names:
        # Split the file name by '_'
        parts = file_name.split('_')
        
        # Extract relevant information
        sample = parts[0]  # e.g., “N01” or “T02”
        patient = int(re.search(r'\d+', sample).group())
        
        # Determine patient and status
        if sample.startswith('N'):
            status = 0
        elif sample.startswith('T'):
            status = 1
        
        lane = parts[-2]  # Get the lane number after 'L'

        fastq_1 = file_name.split(lane)[0] + lane + '_1.fq.gz'
        fastq_2 = file_name.split(lane)[0] + lane + '_2.fq.gz'
        
        # Assign sex based on the mapping; default to None if not found
        sex = sex_mapping.get(sample, None)

        # Append data to the list
        data.append({
            'patient': patient,
            'sex': sex,  # Use the mapped sex value
            'status': status,
            'sample': sample,
            'lane': lane,
            'fastq_1': fastq_1,
            'fastq_2': fastq_2,
        })

    # Create a DataFrame from the data list
    df = pd.DataFrame(data)

    # Remove duplicate rows from the DataFrame
    df.drop_duplicates(inplace=True)  # This will modify df in place

    # Write the DataFrame to a CSV file in the same directory as the input file
    output_csv_file = os.path.join(os.path.dirname(input_file), 'samplesheet.csv')
    df.to_csv(output_csv_file, index=False)  # Write to CSV without row indices

    print(f'DataFrame has been written to {output_csv_file}')

# Main function to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description='Create a sample sheet CSV from a list of file names.')
    parser.add_argument('-input', required=True, help='Input text file containing list of file names')
    
    args = parser.parse_args()
    
    create_sample_sheet(args.input)

if __name__ == '__main__':
    main()