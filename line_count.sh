#!/bin/bash
# Set the search criteria
search_dir=$1
extension=".noseq.tbl"
# Find all files in the search directory and its subdirectories
files=$(find $search_dir -type f -name "*$extension")
# Loop through each file and print the line count
for file in $files
do
  line_count=$(wc -l "$file" | awk '{print $1}')
  echo "$file: $line_count lines"
done
