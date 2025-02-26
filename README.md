PaSTa is a nextflow-based end-to-end image analysis pipeline for decoding image-based spatial transcriptomics data. It performs imaging cycle registration, cell segmentation and transcripts peak decoding. It is currently supports analysis of three types of ST technology:

- in-situ sequencing-like encoding
- MERFISH-like encoding
- RNAScope-like labelling

Prerequisites:
1. Nextflow. Installation guide: https://www.nextflow.io/docs/latest/getstarted.html
2. Docker or Singularity. Installation guide: https://docs.docker.com/get-docker/ or https://sylabs.io/guides/3.7/user-guide/quick_start.html

Demo run with GitPod
---

Check this HackMD from I2K2024 workshop: https://hackmd.io/w4DeWEDxRlKwIPTDCc77XA

Basic run
---
1. Clone the repository
```
git clone -b v0.3.2 https://github.com/cellgeni/Image-ST.git
```
2. Prepare the run.config file *
```
process {
        withName: CELLPOSE {
                maxForks = 1
                ext.args = "--channels [0,0]"
                storeDir = "./output/naive_cellpose_segmentation/"
        }

        withName: POSTCODE {
                memory = {160Gb * task.attempt}
                ext.args = "--channel_names 'DAPI,Cy5,AF488,Cy3,AF750'"
                storeDir = "./output/PoSTcode_decoding_output"
        }

        withName: TO_SPATIALDATA {
                memory = {20.Gb * task.attempt}
                ext.args = "--feature_col 'Name' --expansion_in_pixels 30 --save_label_img False"
                queue = "teramem"
        }

        withName: MERGE_OUTLINES {
                storeDir = "./output/cellpose_segmentation_merged_wkt/"
        }

        withName: BIOINFOTONGLI_MICROALIGNER {
                memory = {50.Gb * task.attempt}
                storeDir = "./output/registered_stacks"
        }

        withName: BIOINFOTONGLI_TILEDSPOTIFLOW {
                maxForks = 1
                memory = {100.Gb * task.attempt}
                storeDir = "./output/spotiflow_peaks/"
        }

        withName: BIOINFOTONGLI_MERGEPEAKS {
                memory = {100.Gb * task.attempt}
                storeDir = "./output/spotiflow_peaks/"
        }

        withName: BIOINFOTONGLI_CONCATENATEWKTS {
                memory = {100.Gb * task.attempt}
                storeDir = "./output/spotiflow_peaks/"
        }

        withName: EXTRACT_PEAK_PROFILE {
                memory = {100.Gb * task.attempt}
                storeDir = "./output/peak_profiles/"
        }

        withLabel: gpu {
                containerOptions = "--nv"
        }
}
```
3. Prepare the parameters file (e.g. iss.yaml)
```
images:
   - ['id': "ISS_exp9",
       [
         "cycle1.ome.tiff",
         "cycle2.ome.tiff",
         "cycle3.ome.tiff",
         "cycle4.ome.tiff",
         "cycle5.ome.tiff",
         "cycle6.ome.tiff",
       ]
     ]
cell_diameters: [15]
chs_to_call_peaks: [1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29]
codebook:
  - ['id': "ISS_exp9", "./codebook.csv", "./dummy.txt"]
segmentation_method: "CELLPOSE"

out_dir: "./output"
```
4. Run the pipeline
```
nextflow run ./Image-ST/main.nf -profile lsf,singularity -c run.config -params-file iss.yaml -entry RUN_DECODING -resume
```
5. Check the output in the specified storeDir.

Spin up Napari with napari-spatialdata plugin installed (https://spatialdata.scverse.org/projects/napari/en/latest/notebooks/spatialdata.html)

Then use the following command to visualize the output
```
from napari_spatialdata import Interactive
import spatialdata as spd

data = spd.read_zarr([path-to-.sdata-folder])
Interactive(data)
```

*: You may leave the process block empty if you want to use the default parameters.

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