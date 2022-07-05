from snakebids.utils.snakemake_io import glob_wildcards

rule extract_zip_eeg:
    """extract zipfile, that optionally contains another zipfile that needs 7za to extract"""
    input:
        zipfile = ancient('zips/site-{site}/sub-{subject}/{filetype}/EPL31_{site}_{subject}_{visit}_SE{sesnum}_{filetype}.zip'),
        unzip_exec = 'ext_bin/7za' #7zip binary required for the enclosed zipfile, compile from https://github.com/jinfeihan57/p7zip
    output:
        raw_dir = directory('raw/site-{site}/sub-{subject}/{filetype,EEG}/EPL31_{site}_{subject}_{visit}_SE{sesnum}_{filetype}')
    shadow: 'minimal'
    group: 'unzip'
    shell:
        "mkdir -p {output.raw_dir} temp_unzipped && "
        "unzip -j -d temp_unzipped {input.zipfile} && "
        "for zip in $(find temp_unzipped -type f -name '*.zip'); "
        "do"
        " {input.unzip_exec} x -bb3 -y -otemp_unzipped ${{zip}}; "
        "done && "
        "mv $(find temp_unzipped -type f ! -name '*.zip') {output.raw_dir}"




def get_wildcards_from_downloaded_ieeg():
    wildcards = glob_wildcards('raw/site-{site}/sub-{subject}/EEG/EPL31_{site}_{subject}_{visit}_SE{sesnum}_EEG/sub-EPL31{site}{subject}_ses-{session}_task-{task}_run-{run}_ieeg.edf')
    
    zip_list = dict()
    zip_list['subject'] = wildcards.subject
    zip_list['site'] = wildcards.site
    zip_list['session'] = wildcards.session
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


print(zip_list)

rule create_bids_eeg_folder:
    input:
        raw_dir = 'raw/site-{site}/sub-{subject}/EEG/EPL31_{site}_{subject}_{visit}_SE{sesnum}_EEG'
    output:
        eeg_dir = directory('bids_{eeg_type}/sub-EPL31{site}{subject}/ses-V{visit}SE{sesnum}/{eeg_type}')
    shell: 
        'mkdir -p {output.eeg_dir} && '
        'for f in $(find {input.raw_dir} -type f ); '
        'do '
        ' ln -srv  $f {output.eeg_dir};'
        'done'

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


