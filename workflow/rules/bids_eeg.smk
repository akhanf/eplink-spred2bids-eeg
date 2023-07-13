from snakebids.utils.snakemake_io import glob_wildcards
import pandas as pd

rule all_bids_eeg:
    input:
        'bids_eeg/dataset_description.json'

rule all_bids_ieeg:
    input:
        'bids_ieeg/dataset_description.json'


#load zip lists from tsv files - these must first be created 
zip_list={}
for suffix in ['eeg','ieeg']:
    zip_list[suffix] = pd.read_table(f'resources/subject_sessions_{suffix}.tsv',dtype={'subject':str})[['subject','session','site','run','task']].to_dict(orient='list')

ieeg_dirs = expand('bids_ieeg/sub-EPL31{site}{subject}/ses-{session}/ieeg',
                        zip,
                        **zip_list['ieeg'],
                        allow_missing=True)

ieeg_dirs = expand('bids_ieeg/sub-EPL31{site}{subject}/ses-{session}/ieeg',
                        zip,
                        **zip_list['ieeg'],
                        allow_missing=True)


print(ieeg_dirs)



rule create_eeg_ieeg_scans_tsv:
    """ this parses the contents of the downloaded and extracted zip files to create a tsv file for the zip lists"""
    output:
        tsv='resources/subject_sessions_{datatype}.tsv'
    params:
        datatype='{datatype}'
    script:
        '../scripts/create_tsv_from_downloaded_eeg_ieeg.py'

rule create_bids_eeg_folder:
    input:
        raw_dir = 'raw/site-{site}/sub-{subject}/EEG/EPL31_{site}_{subject}_{visit}_SE{sesnum}_EEG'
    output:
        eeg_dir = directory('bids_{eeg_type}/sub-EPL31{site}{subject}/ses-V{visit}SE{sesnum}/{eeg_type}')
    shell: 
        "mkdir -p {output.eeg_dir} && "
        "for f in $(find {input.raw_dir} -type f); "
        "do "
        " ln -srv  $f {output.eeg_dir};"
        "done && "
        "rename annotations.tsv events.tsv {output.eeg_dir}/*" #rename annotations.tsv to events.tsv


def get_bids_eeg_dirs(wildcards):
    eeg_dirs = expand(f'bids_{wildcards.eeg_type}'+'/sub-EPL31{site}{subject}/ses-{session}/'+f'{wildcards.eeg_type}',
                    zip,
                    **zip_list[wildcards.eeg_type])
    print(eeg_dirs)
    return eeg_dirs

rule create_dataset_json_eeg:
    input:
        eeg_dirs = get_bids_eeg_dirs,
        json = 'resources/dataset_description_template.json',
        events_json = 'resources/events.json'
    output:
        json = 'bids_{eeg_type}/dataset_description.json',
        events_json = 'bids_{eeg_type}/events.json'
    shell: 'cp {input.json} {output.json} && cp {input.events_json} {output.events_json}'



