import pandas as pd
from snakebids.utils.snakemake_io import glob_wildcards

datatype = snakemake.params.datatype

wildcards = glob_wildcards('raw/site-{site}/sub-{subject}/EEG/EPL31_{site}_{subject}_{visit}_SE{sesnum}_EEG/sub-EPL31{site}{subject}_ses-{session}_task-{task}_run-{run}'+f'_{datatype}.edf')

zip_list = dict()
zip_list['subject'] = wildcards.subject
zip_list['site'] = wildcards.site
zip_list['session'] = [f'V{visit}SE{sesnum}' for visit,sesnum in zip(wildcards.visit,wildcards.sesnum) ]
zip_list['run'] = wildcards.run
zip_list['task'] = wildcards.task

df = pd.DataFrame(zip_list)
df.sort_values(by=['site','subject','session','task','run']).to_csv(snakemake.output.tsv,index=False,sep='\t')


