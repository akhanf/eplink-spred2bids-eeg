from snakebids.utils.snakemake_io import glob_wildcards

def get_raw_dirs_(filetype):
    raw = []
    for site in config['sites']:
        for subject in get_subjects(site):
            
            exp_list = f'resources/site-{site}_sub-{subject}_experiments.txt'
        
            with open(exp_list) as fd:
                for exp in fd.read().splitlines():
                    if exp in config['ignore_experiments']:
                        continue
                    if exp.split('_')[-1] == filetype:
                        rawdir = f'raw/site-{site}/sub-{subject}/{filetype}/{exp}'
                        raw.append(rawdir) 
    return raw                      




rule all_bids_eeg:
    input:
        'bids_eeg/dataset_description.json'

rule all_bids_ieeg:
    input:
        'bids_ieeg/dataset_description.json'



def get_wildcards_from_downloaded_ieeg():
    wildcards = glob_wildcards('raw/site-{site}/sub-{subject}/EEG/EPL31_{site}_{subject}_{visit}_SE{sesnum}_EEG/sub-EPL31{site}{subject}_ses-{session}_task-{task}_run-{run}_ieeg.edf')
    
    zip_list = dict()
    zip_list['subject'] = wildcards.subject
    zip_list['site'] = wildcards.site
    zip_list['session'] = [f'V{visit}SE{sesnum}' for visit,sesnum in zip(wildcards.visit,wildcards.sesnum) ]
    zip_list['run'] = wildcards.run
    zip_list['task'] = wildcards.task

    return zip_list

def get_wildcards_from_downloaded_eeg():
    wildcards = glob_wildcards('raw/site-{site}/sub-{subject}/EEG/EPL31_{site}_{subject}_{visit}_SE{sesnum}_EEG/sub-EPL31{site}{subject}_ses-V{visit}SE{sesnum}_task-{task}_run-{run}_eeg.edf')

    zip_list = dict()
    zip_list['subject'] = wildcards.subject
    zip_list['site'] = wildcards.site
    zip_list['session'] = [f'V{visit}SE{sesnum}' for visit,sesnum in zip(wildcards.visit,wildcards.sesnum) ]
    zip_list['run'] = wildcards.run
    zip_list['task'] = wildcards.task

    return zip_list


zip_list=dict()
zip_list['ieeg'] = get_wildcards_from_downloaded_ieeg()
zip_list['eeg'] = get_wildcards_from_downloaded_eeg()



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

rule create_dataset_json_eeg:
    input:
        eeg_dirs = lambda wildcards: expand(expand('bids_{eeg_type}/sub-EPL31{site}{subject}/ses-{session}/{eeg_type}',
                        zip,
                        **zip_list[wildcards.eeg_type], 
                        allow_missing=True),
                        eeg_type=wildcards.eeg_type),
        json = 'resources/dataset_description_template.json',
        events_json = 'resources/events.json'
    output:
        json = 'bids_{eeg_type}/dataset_description.json',
        events_json = 'bids_{eeg_type}/events.json'
    shell: 'cp {input.json} {output.json} && cp {input.events_json} {output.events_json}'



