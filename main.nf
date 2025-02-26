include { DECODE } from './workflows/main' 
include { RNASCOPE } from './workflows/main' 

params.images = [
    [['id': "run_id"], [
        "cycle1",
        "cycle2",
    ]],
]

params.chs_to_call_peaks = [1, 2]

params.codebook = [["id":'' ], "codebook.csv"]
params.segmentation_method = "CELLPOSE"


workflow RUN_RNASCOPE {
    RNASCOPE()
}

workflow RUN_DECODING {
    DECODE()
}