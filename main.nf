include { DECODE } from './workflows/main' 
include { RNASCOPE } from './workflows/main' 

params.images = [
    [['id': "run_id"], [
        "cycle1",
        "cycle2",
    ]],
]

workflow RUN_RNASCOPE {
    RNASCOPE()
}

workflow RUN_DECODING {
    DECODE()
}