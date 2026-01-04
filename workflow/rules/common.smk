
def r1_from_sample(wildcards):
    return samplesheet.loc[wildcards.sample]['read1']

def r2_from_sample(wildcards):
    return samplesheet.loc[wildcards.sample]['read2']

def get_bedg_control(wildcards):
    rep_index = samplesheet.loc[wildcards.sample]['group']
    controls = samplesheet[samplesheet['condition'] == "Control"]
    sample = controls.loc[controls['group'] == rep_index].index.tolist()[0]
    return "bedalignments/"+ sample + ".bedgraph"

def get_bam_control(wildcards):
    rep_index = samplesheet.loc[wildcards.sample]['group']
    controls = samplesheet[samplesheet['condition'] == "Control"]
    sample = controls.loc[controls['group'] == rep_index].index.tolist()[0]

    # Respect the dedup configuration - use appropriate path
    if config["dedup"]:
        return "dedupalignments/" + sample + ".sorted.qflt.bam"
    else:
        return "alignments/" + sample + ".sorted.qflt.bam"

