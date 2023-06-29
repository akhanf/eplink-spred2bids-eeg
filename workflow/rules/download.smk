import xnat

envvars:
    'SPRED_USER',
    'SPRED_PASS'

def get_subjects(site):
    if os.path.exists(f'resources/subjects_{site}.txt'):
        with open(f'resources/subjects_{site}.txt','r') as f:
            return [s.replace('\n','') for s in f.readlines()]

def get_all_experiment_lists():
    exps = []
    for site in config['sites']:
        for subject in get_subjects(site):
            exps.append(f'resources/site-{site}_sub-{subject}_experiments.txt')
    return exps


def get_zips(filetype):
    zips = []
    for site in config['sites']:
        for subject in get_subjects(site):
            
            exp_list = f'resources/site-{site}_sub-{subject}_experiments.txt'
        
            with open(exp_list) as fd:
                for exp in fd.read().splitlines():
                    if exp in config['ignore_experiments']:
                        continue
                    if exp.split('_')[-1] == filetype:
                        zipfile = f'zips/site-{site}/sub-{subject}/{filetype}/{exp}.zip'
                        zips.append(zipfile) 
    return zips                      


rule all_download_ephys_zips:
    input: 
        get_zips('EEG')



rule get_indiv_zip:
    input:
        exp_list = ancient('resources/site-{site}_sub-{subject}_experiments.txt') 
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
    


def get_raw_dirs(filetype):
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
    print(raw)
    return raw                      



rule all_extract_ephys_zips:
    input:
        get_raw_dirs('EEG')



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





