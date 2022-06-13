

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



