
def get_mem_mb(wildcards, threads):
    return threads * 8000

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
    # Match the dedup setting used for experimental samples and use quality-filtered BAMs
    if config["dedup"]:
        return "dedupalignments/{}.sorted.qflt.bam".format(sample)
    else:
        return "alignments/{}.sorted.qflt.bam".format(sample)

