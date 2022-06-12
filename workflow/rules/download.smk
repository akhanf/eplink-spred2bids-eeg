import xnat

envvars:
    'SPRED_USER',
    'SPRED_PASS'


localrules: get_subject_list

rule get_subject_list:
    params:
        site_id = 'EPL31_{site}'
    output:
        subj_list = 'resources/subjects_{site}.txt'
    run:
        session = xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
        subjects = [row[0] for row in session.projects[params.site_id].subjects.tabulate(columns=['label'])]
        with open(output.subj_list, "w") as out:
            for s in subjects:
                out.write(s.split('_')[2]+'\n') 
        session.disconnect()

rule get_imaging_sessions:
    params:
        site_id = 'EPL31_{site}',
        subj_id = 'EPL31_{site}_{subject}'
    output:
        exp_list = 'resources/site-{site}_sub-{subject}_experiments.txt'
    run:
        session = xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
        experiments = [row[0] for row in session.projects[params.site_id].subjects[params.subj_id].experiments.tabulate(columns=['label'])]

        with open(output.exp_list, "w") as out:
            for exp in experiments:
                out.write(exp+'\n') 
        session.disconnect()




rule get_indiv_zip:
    input:
        exp_list = 'resources/site-{site}_sub-{subject}_experiments.txt'
    params:
        remote_path = '/data/projects/EPL31_{site}/experiments/EPL31_{site}_{subject}_{visit}_SE{sesnum}_{filetype}',
        zipdir = 'zips/site-{site}/sub-{subject}/{filetype}',
        debug_args = {'loglevel':'DEBUG','debug':True} if config.get('debug_download',False) else {},
    output:
        zipfile = protected('zips/site-{site}/sub-{subject}/{filetype}/EPL31_{site}_{subject}_{visit}_SE{sesnum}_{filetype}.zip')
    run:
        shell('mkdir -p {params.zipdir}')
        with xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'],**params.debug_args) as session:
            experiment = session.create_object(params.remote_path)
            experiment.download(output.zipfile) 
    




rule extract_zip_eeg_lhs:
    """LHS has zip in zip file; extract zip file, then extract the zip file in that; creates a flat bids tree"""
    input:
        zipfile = 'zips/site-{site}/sub-{subject}/{filetype}/EPL31_{site}_{subject}_{visit}_SE{sesnum}_{filetype}.zip'
    output:
        raw_dir = directory('raw/site-{site,LHS}/sub-{subject}/ses-V{visit}SE{sesnum}_{filetype,EEG}')
    shadow: 'minimal'
    shell:
        'mkdir -p {output.raw_dir} temp_zip && '
        'unzip -j -d temp_zips {input.zipfile} && '
        'for zip in $(ls temp_zips/*.zip); '
        'do'
        ' unzip -j -d {output.raw_dir} ${{zip}}; '
        'done'

rule extract_zip_eeg_twh_hsc:
    """TWH and HSC just have a zip file with the bids data"""
    input:
        zipfile = 'zips/site-{site}/sub-{subject}/{filetype}/EPL31_{site}_{subject}_{visit}_SE{sesnum}_{filetype}.zip'
    output:
        raw_dir = directory('raw/site-{site,TWH|HSC}/sub-{subject}/ses-V{visit}SE{sesnum}_{filetype,EEG}')
    shadow: 'minimal'
    shell:
        'mkdir -p {output.raw_dir} && '
        'unzip -j -d temp_zips {input.zipfile}'





#TODO: update this rule to expand over the indiv MR zips (usually only one, but could be two?)
rule make_dicom_tar:
    input:
        zip_dir = 'zips/site-{site}/sub-{subject}/{filetype}'
    params:
        file_match = '*/scans/*/resources/DICOM/files/*',
        temp_dir = os.path.join(config['tmp_download'],'raw/site-{site}/sub-{subject}/{filetype}')
    output:
        tar = 'raw/site-{site}/sub-{subject}/{filetype,MR}/sub-{subject}.tar'
    group: 'dl'
    shell:
        'mkdir -p {params.temp_dir} && '
        'for zip in $(ls {input.zip_dir}/*.zip); '
        'do '
        ' unzip -d {params.temp_dir} ${{zip}} {params.file_match}; '
        'done && '
        'tar -cvf {output.tar} {params.temp_dir} && '
        'rm -rf {params.temp_dir}'

       
rule tar_to_bids:
    input:
        tar = 'raw/site-{site}/sub-{subject}/MR/sub-{subject}.tar',
        heuristic = lambda wildcards: config['tar2bids'][wildcards.site],
        container = 'resources/singularity/tar2bids.sif'

    params:
        temp_bids_dir = 'raw/site-{site}/sub-{subject}/MR/temp_bids',
        heudiconv_tmpdir = os.path.join(config['tmp_download'],'{site}','{subject}')
    output:
        dir = directory('bids/site-{site}/sub-{subject}')
    group: 'dl'
    shell: 
        "mkdir -p {params.heudiconv_tmpdir} && "
        "singularity run -e {input.container} -o {params.temp_bids_dir} "
        " -w {params.heudiconv_tmpdir} -h {input.heuristic} -N 1 -T 'sub-{{subject}}' {input.tar} || true && " #force clean exit for tar2bids
        "mkdir -p {output.dir} && rsync -av {params.temp_bids_dir}/sub-{wildcards.subject}/ {output.dir} && "
        "rm -rf {params.heudiconv_tmpdir}"


rule create_dataset_json:
    input:
        dir = lambda wildcards: expand('bids/site-{site}/sub-{subject}',site=wildcards.site,subject=get_subjects(wildcards.site)),
        json = 'resources/dataset_description_template.json'
    output:
        json = 'bids/site-{site}/dataset_description.json'
    shell: 'cp {input.json} {output.json}'



