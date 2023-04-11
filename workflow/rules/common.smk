import pandas as pd

def get_mem_mb(wildcards, threads):
    return threads * 6000

def r1_from_sample(wildcards):
    return samplesheet.loc[wildcards.sample]['read1']

def r2_from_sample(wildcards):
    return samplesheet.loc[wildcards.sample]['read2']

def get_control(wildcards):
    rep_index = samplesheet.loc[wildcards.sample]['group']
    controls = samplesheet[samplesheet['condition'] == "Control"]
    sample = controls.loc[controls['group'] == rep_index].index.tolist()[0]
    return "bedalignments/"+ sample + ".bedgraph"
