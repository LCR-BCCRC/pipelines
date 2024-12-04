#!/usr/bin/env python

import subprocess
import sys
import os

## Logging files
global stdout
global stderr
global stdout_f
global stderr_f

stdout = snakemake.log["stdout"]
stderr = snakemake.log["stderr"]

def increaseSleepInterval(batch_file):
    """
    Increase sleep interval between batch commands by 5 seconds
    """
    os.system(f'sleep=$(grep "Sleep" {batch_file} | cut -d " " -f2) && new_sleep=$(($sleep + 5000)) && sed -i "s/setSleepInterval $sleep/setSleepInterval $new_sleep/g" {batch_file}')

def runIGV(batch_file, igv, status, message, attempt):
    os.system(f'echo "Snapshot may be {status}. {message}... Rerunning IGV... Attempt {str(attempt)}:" >> {stdout} 2>> {stderr}')
    os.system(f'maxtime=$(($(wc -l < {batch_file}) * 60 + 15)) && timeout --foreground $maxtime xvfb-run -s "-screen 0 1980x1020x24" {igv["server_num"]} {igv["server_args"]} {igv["igv"]} -b {batch_file} >> {stdout} 2>> {stderr}')

def getImageQualities(snapshot, batch_file, igv, summary_file, attempts = 0):
    height = None
    width = None
    kurtosis = None
    skewness = None
    corrupt = True
    while attempts < 4 and corrupt:
        attempts += 1
        try:
            height = str(subprocess.check_output(f"identify -format '%h' {snapshot}", shell=True)).split("'")[1].split("\\n")[0]
            width = str(subprocess.check_output(f"identify -format '%w' {snapshot}", shell=True)).split("'")[1].split("\\n")[0]
            kurtosis, skewness = [float(value.split(": ")[1]) for value in str(subprocess.check_output(f"identify -verbose {snapshot} | grep -E 'kurtosis|skewness' | tail -n 2", shell=True)).replace("\\n'","").split("\\n      ")]
            corrupt = False
        except:
            status = "corrupt"
            message = ""
            if attempts < 3:
                runIGV(batch_file, igv, status, message, attempts)
    if attempts == 4 and corrupt:
        # Quit because snapshot is corrupt and can't run quality control
        qc_status = "failed"
        logFailedSnapshots(snapshot, summary_file, qc_status)
        sys.exit()

    quality_dict = {
        "height": height,
        "width": width,
        "kurtosis": kurtosis,
        "skewness": skewness
    }

    return quality_dict

def handleIncorrectDimensions(snapshot, img_values, thresholds, batch_file, igv, failed_summary):
    status = "in incorrect dimensions"
    attempts = 0

    # Increase sleep interval 
    increaseSleepInterval(batch_file)

    while float(img_values["width"]) == 640 and attempts < 5:
        attempts += 1
        previous_server_arg = igv["server_num"]

        if igv["server_num"] == "--auto-servernum" or int(igv["server_num"].replace("-n ","")) >= 99:
            new_server_arg = "-n 1"
        else:
            new_server_arg = f'-n {str(float(igv["server_num"]) + 1)}'

        messsage = f'Current snapshot width is {img_values["width"]}, but at least 1020px is expected. This might occurr if xvfb-run is unable to connect to current server ({previous_server_arg}) due to a server lock. Switching server numbers... Attempting new server argument: {new_server_arg}'

        igv["server_num"] = new_server_arg

        # Rerun IGV
        runIGV(batch_file, igv, status, message, attempts)

        # Update image values
        img_values = getImageQualities(snapshot, batch_file, igv, failed_summary)

    return attempts, img_values

def handleTruncated(snapshot, img_values, thresholds, batch_file, igv, dim_attempts, failed_summary):
    status = "truncated"
    attempts = 0

    # Increase sleep interval
    if dim_attempts == 0:
        increaseSleepInterval(batch_file)

    # Rerun IGV until height is no longer truncated
    while float(img_values["height"]) in thresholds["truncated"] and attempts < 3:
        attempts += 1
        message = f'Current snapshot height is {img_values["height"]}.'

        # Rerun IGV
        runIGV(batch_file, igv, status, message, attempts)
        
        # Update image values
        img_values = getImageQualities(snapshot, batch_file, igv, failed_summary)

    return attempts, img_values

def handleBlank(snapshot, img_values, thresholds, batch_file, igv, truncated_attempts, dim_attempts, failed_summary):
    status = "blank"
    blank = True
    attempts = 0

    if dim_attempts == 0 and truncated_attempts == 0:
        increaseSleepInterval(batch_file)
    
    while blank and attempts < 3:
        attempts += 1
        message = f'Current snapshot values are: {img_values["height"]} height, {img_values["kurtosis"]} kurtosis, and {img_values["skewness"]} skewness. Snapshots with these values may be blank. Blank snapshots may be due to errors reading BAM file headers, Java address bind errors, or other errors that occur during the IGV run. Rerunning with increased sleep interval.'

        # Rerun IGV
        runIGV(batch_file, igv, status, message, attempts)

        # Get updated values
        img_values = getImageQualities(snapshot, batch_file, igv, failed_summary)

        # Check if still blank
        blank = is_blank(img_values, thresholds)
    
    return attempts, img_values

def is_blank(img_values, thresholds):
    blank_check = any(
        (
            (
                float(img_values["height"]) == float(height_threshold) or
                ("<" in height_threshold and float(img_values["height"]) < float(height_threshold.replace("<",""))) or
                (">" in height_threshold and float(img_values["height"]) > float(height_threshold.replace(">","")))
            ) and
            (float(img_values["kurtosis"]) > float(thresholds[height_threshold]["kurtosis"])) and
            (float(img_values["skewness"]) < float(thresholds[height_threshold]["skewness"]))
        )
        for height_threshold in list(thresholds)
    )

    return blank_check

def qualityControl(snapshot, batch_file, igv, img_values, thresholds, failed_summary, attempts=0):
    # Set default values for attempts so sleep value is not perpetually increased
    dimension_attempts = 0
    truncated_attempts = 0
    blank_attempts = 0

    # Set default qc status
    qc_status = "pass"

    # Check width
    if float(img_values["width"]) == 640:
        dimension_attempts, img_values = handleIncorrectDimensions(snapshot, img_values, thresholds, batch_file, igv, failed_summary)

    # Handle truncated attempts
    if float(img_values["height"]) in thresholds["truncated"]:
        truncated_attempts, img_values = handleTruncated(snapshot, img_values, thresholds, batch_file, igv, dimension_attempts, failed_summary)

    # Check if blank
    blank_thresholds = thresholds["blank"]

    blank = is_blank(img_values, blank_thresholds)

    if blank:
        blank_attempts, img_values = handleBlank(snapshot, img_values, blank_thresholds, batch_file, igv, truncated_attempts, dimension_attempts, failed_summary)

    # Check final values and log failed/suspicious
    if float(img_values["width"]) == 640:
        os.system(f'echo "Snapshot width is {img_values["width"]}. Improper dimensions should be fixed and rerun. Check snapshot {snapshot}" >> {stdout}')
        qc_status = "fail"

    if float(img_values["height"]) in thresholds["failed"]:
        os.system(f'echo "Snapshot height is {img_values["height"]} and may still be truncated or improperly loaded. Check snapshot {snapshot}" >> {stdout}')
        qc_status = "fail"

    if qc_status == "pass" and (dimension_attempts >= 5 or truncated_attempts >= 3 or blank_attempts >= 3):
        qc_status = "suspicious"

    return qc_status

def logFailedSnapshots(snapshot, summary_file, qc_status):
    outline = "\t".join([snakemake.wildcards["sample_id"], 
                        snakemake.wildcards["seq_type"],
                        snakemake.wildcards["genome_build"],
                        snakemake.wildcards["gene"],
                        snakemake.wildcards["chromosome"],
                        snakemake.wildcards["start_position"],
                        snakemake.wildcards["preset"],
                        qc_status,
                        snapshot])
    with open(summary_file, "a") as handle:
        handle.write(outline + "\n")

def main():
    stdout_f = open(snakemake.log["stdout"], "a")
    stderr_f = open(snakemake.log["stderr"], "a")

    # Output file
    outfile = snakemake.output["snapshot_qc"]

    ## Quality control variables
    snapshot = snakemake.input["snapshot"]
    qc_thresholds = snakemake.params["thresholds"]

    ## Batch scripts

    batch_script = snakemake.params["batch_script"]
    merged_batch = snakemake.params["merged_batch"]
    batch_temp = snakemake.params["batch_temp"]
    # Set up the temporary batch script file
    os.system(f'cat {batch_script} > {batch_temp} && echo "exit" >> {batch_temp}')

    ## Variables for running IGV 
    igv_exec = {
        "igv": snakemake.params["igv"],
        "server_num": snakemake.params["server_number"],
        "server_args": snakemake.params["server_args"]
    }

    ## Summary file to append failed snapshots to
    f_summary = snakemake.input["failed_summary"]

    # Get image qualities
    img_values = getImageQualities(snapshot, batch_temp, igv_exec, f_summary)

    description_line = f'Initial image values are:\nHeight: {img_values["height"]}\nWidth: {img_values["width"]}\nKurtosis: {img_values["kurtosis"]}\nSkewness:{img_values["skewness"]}\n'
    stdout_f.write(description_line)

    # Control the qualities
    results = qualityControl(snapshot, batch_temp, igv_exec, img_values, qc_thresholds, f_summary)

    if results != "pass":
        logFailedSnapshots(snapshot, f_summary, results)

    if results != "fail":
        os.system(f'touch {outfile}')

    # Cleanup
    os.remove(batch_temp)
    stdout_f.close()
    stderr_f.close()

if __name__ == "__main__":
    main()
