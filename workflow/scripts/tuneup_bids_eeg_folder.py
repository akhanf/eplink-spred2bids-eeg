from snakemake.shell import shell
from glob import glob
from os.path import getsize
import pandas as pd
import numpy as np

shell("mkdir -p {snakemake.output.eeg_dir}")

#make symlinks for edf files instead of copying
for edf in glob(f'{snakemake.input.raw_dir}/*.edf'):
    shell("ln -srv {edf} {snakemake.output.eeg_dir}")

#for all other files (json), copy the actual files so we can edit them if needed..
tsv_json_files = glob(f'{snakemake.input.raw_dir}/*.json')+glob(f'{snakemake.input.raw_dir}/*.tsv')
shell("cp -RLv {tsv_json_files} {snakemake.output.eeg_dir}")

#convert all files from dos to unix (symlinks are skipped automatically):
shell("dos2unix {snakemake.output.eeg_dir}/*")

#rename annotations.tsv to events.tsv:
shell("rename annotations.tsv events.tsv {snakemake.output.eeg_dir}/*")

#check for events.tsv files with file size <2 bytes (near-empty files containing single character)
for tsv in glob(f'{snakemake.output.eeg_dir}/*events.tsv'):
    if getsize(tsv) < 2:
        shell("cp {snakemake.input.events_hdr_only} {tsv}")

#for tsv files, we will read in with pandas, then write out, replacing '\n' with '; ' 
# since this seems to occur frequently in the event description, but causes the validator to fail (since 
# it seems the parsing of the tsv files is not as smart as pandas)...

# also need to replace empty strings with n/a (validator doesn't like empty strings)

for tsv in glob(f'{snakemake.output.eeg_dir}/*events.tsv'):
    df = pd.read_csv(tsv,sep='\t')
    df.event = df.event.str.replace('\n','; ')
    df.event = df.event.replace(np.nan,'n/a') #replace nan (from empty string) with 'n/a'
    df.to_csv(tsv,sep='\t',index=False)


