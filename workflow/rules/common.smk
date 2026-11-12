
def r1_from_sample(wildcards):
    return samplesheet.loc[wildcards.sample]['read1'].strip()

def r2_from_sample(wildcards):
    return samplesheet.loc[wildcards.sample]['read2'].strip()

def get_bedg_control(wildcards):
    rep_index = samplesheet.loc[wildcards.sample]['group']
    controls = samplesheet[samplesheet['condition'] == "Control"]
    sample = controls.loc[controls['group'] == rep_index].index.tolist()[0].strip()
    return "bedalignments/{}.bedgraph".format(sample)

def get_bam_control(wildcards):
    rep_index = samplesheet.loc[wildcards.sample]['group']
    controls = samplesheet[samplesheet['condition'] == "Control"]
    sample = controls.loc[controls['group'] == rep_index].index.tolist()[0].strip()
    # quality-filtered BAM (deduplicated if dedup applies to the control)
    return "dedupalignments/{}.sorted.qflt.bam".format(sample)

