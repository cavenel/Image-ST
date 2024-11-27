PaSTa is a nextflow-based end-to-end image analysis pipeline for decoding image-based spatial transcriptomics data. It performs imaging cycle registration, cell segmentation and transcripts peak decoding. It is currently supports analysis of three types of ST technology:

- in-situ sequencing-like encoding
- MERFISH-like encoding
- RNAScope-like labelling

Demo run
---

Check this HackMD from I2K2024 workshop: https://hackmd.io/w4DeWEDxRlKwIPTDCc77XA

FAQ
---

1. My HOME dir is full when running Singularity image conversion on HPC.

A quick solution is to manually specify singularity dir by setting:


```
singularity cache clean
export SINGULARITY_CACHEDIR=./singularity_image_dir
export NXF_SINGULARITY_CACHEDIR=./singularity_image_dir
```

2. How do I modify parameters to specific process/step?

By following nf-core standard, it is possible to add any parameters to the main script using ext.args=”--[key] [value]” in the run.config file.

An example is

withName: POSTCODE {
    ext.args = "--channel_names 'DAPI,Cy5,AF488,Cy3,AF750'"
}

3. Cannot download pretrained model for the deep-learning tools (Spotiflow/CellPose)

> Exception: URL fetch failure on https://drive.switch.ch/index.php/s/6AoTEgpIAeQMRvX/download: None -- [Errno -3] Temporary failure in name resolution
Or CellPose
urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>

Mostly likely you've reached max download (?), wait a bit and try later OR manually download those models and update the configuration file.