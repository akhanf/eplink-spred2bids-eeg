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
    





