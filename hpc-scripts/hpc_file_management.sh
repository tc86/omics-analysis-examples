#!/bin/bash

## This is a script to organize the specific downstream files from running a WES pipeline so these can be utilized to create CLIA style pdf reports later
## Author: Tamrin Chowdhury
## Date: 12/19/2024


# Paths and directories

LIST_PATH="~/raw_data/batch1/sample_list.txt"

SOURCE_ROOT="~/raw_data/batch1"

DEST_ROOT="~/output"

LOG_FILE="~/output/log/sample_move_log.txt"



# Create log file or clear if it exists

: > "$LOG_FILE"



# Function for logging

log_message() {

    echo "$(date +"%Y-%m-%d %T") - $1" | tee -a "$LOG_FILE"

}



# Function to move file if it does not already exist in the destination

move_file() {

    local src="$1"

    local dest="$2"

    local logfile="$3"



    if [[ -e "$dest/$(basename "$src")" ]]; then

        log_message "Skipped $(basename "$src"): already exists in $dest" | tee -a "$logfile"

    else

        mv "$src" "$dest" && log_message "Moved $(basename "$src") to $dest" || log_message "Failed to move $(basename "$src") to $dest"

    fi

}



# Read the list file line by line

while IFS=$'\t' read -r case_id tumor_sample normal_sample pair_name; do

    # Strip carriage return characters if any

    case_id=$(echo "$case_id" | tr -d '\r')

    tumor_sample=$(echo "$tumor_sample" | tr -d '\r')

    normal_sample=$(echo "$normal_sample" | tr -d '\r')

    pair_name=$(echo "$pair_name" | tr -d '\r')



    if [[ -z "$case_id" || -z "$tumor_sample" || -z "$normal_sample" || -z "$pair_name" ]]; then

        log_message "Skipping incomplete line: $case_id, $tumor_sample, $normal_sample, $pair_name"

        continue

    fi



    # Define source directories for each sample and pair

    tumor_dir="$SOURCE_ROOT/$tumor_sample"

    normal_dir="$SOURCE_ROOT/$normal_sample"

    pair_dir="$SOURCE_ROOT/$pair_name"



    # Define destination directories

    ss_output_dir="$DEST_ROOT/ss_output/$case_id"

    cnv_dir="$DEST_ROOT/cnv/$case_id"

    loh_dir="$DEST_ROOT/loh/$case_id"

    report_dir="$DEST_ROOT/report/$case_id"

    cnv_plot_dir="$DEST_ROOT/cnv_plot/$case_id"

    metric_dir="$DEST_ROOT/metric/$case_id"

    paired_output_dir="$DEST_ROOT/paired_output/$case_id"



    # Create destination directories if they don't exist

    mkdir -p "$ss_output_dir" "$cnv_dir" "$loh_dir" "$report_dir" "$cnv_plot_dir" "$metric_dir" "$paired_output_dir"



    # Move tumor and normal sample directories

    if [[ -d "$ss_output_dir/$tumor_sample" ]]; then

        log_message "Skipped $tumor_sample directory: already exists in $ss_output_dir" | tee -a "$LOG_FILE"

    else

        mv "$tumor_dir" "$ss_output_dir/" && log_message "Moved $tumor_dir to $ss_output_dir" || log_message "Failed to move $tumor_dir to $ss_output_dir"

    fi

    

    if [[ -d "$ss_output_dir/$normal_sample" ]]; then

        log_message "Skipped $normal_sample directory: already exists in $ss_output_dir" | tee -a "$LOG_FILE"

    else

        mv "$normal_dir" "$ss_output_dir/" && log_message "Moved $normal_dir to $ss_output_dir" || log_message "Failed to move $normal_dir to $ss_output_dir"

    fi



    # Move and organize files from the pair directory

    move_file "$pair_dir/${pair_name}.cnvedit.txt" "$cnv_dir" "$LOG_FILE"

    move_file "$pair_dir/${pair_name}.lohedit.txt" "$loh_dir" "$LOG_FILE"

    move_file "$pair_dir/${pair_name}.report.xlsx" "$report_dir" "$LOG_FILE"

    move_file "$pair_dir/${pair_name}_bafplot.png" "$cnv_plot_dir" "$LOG_FILE"

    move_file "$pair_dir/${pair_name}_edits.txt" "$metric_dir" "$LOG_FILE"



    # Move and rename working directory

    if [[ -d "$paired_output_dir/$tumor_sample" ]]; then

        log_message "Skipped working directory: already exists as $paired_output_dir/$tumor_sample" | tee -a "$LOG_FILE"

    else

        mv "$pair_dir/working" "$paired_output_dir/$tumor_sample" && log_message "Moved and renamed working directory to $paired_output_dir/$tumor_sample" || log_message "Failed to move and rename working directory to $paired_output_dir/$tumor_sample"

    fi



done < "$LIST_PATH"



log_message "Script execution completed."
