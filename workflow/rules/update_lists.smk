import xnat

envvars:
    'SPRED_USER',
    'SPRED_PASS'


localrules: get_subject_list


rule all_update_subject_lists:
    input:
        expand('resources/subjects_{site}.txt',site=config['sites'])


rule get_subject_list:
    params:
        site_id = 'EPL31_{site}',
        exp_list = lambda wildcards: 'resources/site-{site}_sub-{subject}_experiments.txt'
    output:
        subj_list = 'resources/subjects_{site}.txt'
    run:
        session = xnat.connect(config['spred_url'],user=os.environ['SPRED_USER'],password=os.environ['SPRED_PASS'])
        subjects = sorted([row[0] for row in session.projects[params.site_id].subjects.tabulate(columns=['label'])])
        with open(output.subj_list, "w") as out_subjlist:
            for s in subjects:
                out_subjlist.write(s.split('_')[2]+'\n') 

                experiments = [row[0] for row in session.projects[params.site_id].subjects[s].experiments.tabulate(columns=['label'])]
                with open(params.exp_list.format(subject=s.split('_')[2],site=wildcards.site), "w") as out_explist:
                    for exp in experiments:
                        out_explist.write(exp+'\n') 


        session.disconnect()


