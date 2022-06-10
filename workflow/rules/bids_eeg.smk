from snakebids.utils.snakemake_io import glob_wildcards


def get_wildcards_from_downloaded_ieeg():
    wildcards = glob_wildcards('raw/site-{site}/sub-{subject}/EEG/sub-EPL31{site}{subject}_ses-{session}_task-{task}_run-{run}_ieeg.edf')
    
    zip_list = dict()
    zip_list['subject'] = wildcards.subject
    zip_list['site'] = wildcards.site
    zip_list['session'] = wildcards.session
    zip_list['run'] = wildcards.run
    zip_list['task'] = wildcards.task

    return zip_list

def get_wildcards_from_downloaded_eeg():
    wildcards = glob_wildcards('raw/site-{site}/sub-{subject}/EEG/sub-EPL31{site}{subject}_ses-{session}_task-{task}_run-{run}_eeg.edf')
    zip_list = dict()
    zip_list['subject'] = wildcards.subject
    zip_list['site'] = wildcards.site
    zip_list['session'] = wildcards.session
    zip_list['run'] = wildcards.run
    zip_list['task'] = wildcards.task

    return zip_list




zip_list=dict()
zip_list['ieeg'] = get_wildcards_from_downloaded_ieeg()
zip_list['eeg'] = get_wildcards_from_downloaded_eeg()




rule create_bids_eeg_folder:
    input:
        raw_dir = 'raw/site-{site}/sub-{subject}/EEG'
    params:
        glob = 'sub-EPL31{site}{subject}_ses-{session}_*'
    output:
        eeg_dir = directory('bids_{eeg_type}/sub-EPL31{site}{subject}/ses-{session}/{eeg_type}')
    shell: 
        'mkdir -p {output.eeg_dir} && '
        'cp -v {input.raw_dir}/{params.glob} {output.eeg_dir}'

rule create_dataset_json_eeg:
    input:
        eeg_dirs = lambda wildcards: expand(expand('bids_{eeg_type}/sub-EPL31{site}{subject}/ses-{session}/{eeg_type}',
                        zip,
                        **zip_list[wildcards.eeg_type], 
                        allow_missing=True),
                        eeg_type=wildcards.eeg_type),
        json = 'resources/dataset_description_template.json'
    output:
        json = 'bids_{eeg_type}/dataset_description.json'
    shell: 'cp {input.json} {output.json}'


rule :
    input: 
#        raw_dir = lambda wildcards: expand('raw/site-{site}/sub-{subject}/EEG',site=wildcards.site,subject=get_subjects(wildcards.site))
#rule grab_ieeg_prebids:
#    input:
#        'raw/site-{site}/sub-{subject}/EEG/sub-EPL31LHS0003_ses-V02SE16_task-full_run-01_ieeg.edf'
#        get_inputs
#    output:
#        'bids-combined/sub-{subject}/ses-{session}/ieeg/sub-{subject}_ses-{session}_task-{task}_run-{run}_ieeg.edf'


