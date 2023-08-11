# eplink-spred2bids
Snakemake workflow for downloading and converting EpLink EEG data to BIDS

## Prerequisites:

You need to have BrainCODE SPRED (XNAT) access to the Eplink EPL31 project. 
Put your username and password in environment variables named `SPRED_USER` and `SPRED_PASS` respectively.

## Instructions:

1. Clone this repository and install with `pip install <path_to_cloned_repo>`

2. Update the config file to change the `tmp_download` folder to a local disk with large enough space. 

3. Run snakemake with a dry-run first:
```
snakemake -np
```

4. If everything looks fine, run with the specified number of parallel cores, e.g. 4 as below:
```
snakemake --cores 4
```


## Instructions for bids EEG and iEEG:

### Step 1: Download zip files, or symlink to an existing `zips` folder

```
snakemake all_download_ephys_zips
```
Note: this is 4TB and may take a couple days.


### Step 2: Extract the zip files 

```
snakemake all_extract_ephys_zips
```

### Step 3: Create the consolidated `bids_eeg` or `bids_ieeg` datasets


```
snakemake all_bids_eeg
snakemake all_bids_ieeg
```

Note: all zip files must be extracted before running this part, since it uses the files to infer the wildcards (ie to tell between eeg and ieeg) 

