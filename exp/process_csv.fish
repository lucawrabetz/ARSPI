#!/usr/bin/env fish

# Check if the user provided the directory containing CSV files as an argument
if test (count $argv) -ne 1
    echo "Usage: $argv[0] <directory>"
    exit 1
end

set csv_dir $argv[1]

# Check if the provided directory exists
if not test -d $csv_dir
    echo "Directory '$csv_dir' does not exist."
    exit 1
end

# Loop through all CSV files in the provided directory
for file in $csv_dir/*.csv
    # Check if the file exists
    if test -f $file
        # Run the command python sheet_cat.py <filename> for each CSV file
        python (dirname (status -f))"/sheet_cat.py" $file
    end
end

