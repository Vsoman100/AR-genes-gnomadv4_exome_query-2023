
# Set the filenames
output_file="${1}_concatenated_GCF.tsv"
header_file="header.txt"
temp_file="temp_concatenated.txt"

# Create the temporary file
> "$temp_file"

# Initialize a variable to store the header
header=""

# Iterate over the files
for file in *_GCF.tsv; do
  # Remove the header from each file and append the content to the temporary file
  tail -n +2 "$file" >> "$temp_file"

  # Extract the header from the first file
  if [ -z "$header" ]; then
    header=$(head -n 1 "$file")
  fi
done

# Concatenate the header and the temporary file into the final concatenated file
{ echo "$header"; cat "$temp_file"; } > "$output_file"

# Remove the temporary file
rm "$temp_file"


