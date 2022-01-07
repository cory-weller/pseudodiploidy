#!/usr/bin/env bash

# script to retrieve raw data set


getData () {
    filePath=${1}
    fileURL=${2}
    fileSHA512=${3}

    if [[ ! -f ${filePath} ]]; then
        wget -O ${filePath} ${fileURL}
    else
        echo "${filePath} already downloaded"
    fi
    checksum=$(sha512sum ${filePath} | cut -d " " -f 1)
    if [[ ${checksum} != "${fileSHA512}" ]]; then
        echo "WARNING: ${filePath} checksum failed, contents may differ from those originally used"
    else
        echo "${filePath} checksum passed"
    fi
}

getData \
    "data/external/SGD-features.tab" \
    "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab" \
    "cc3c54abb0106993b4d09227f4748e46460b4a7614518b6513863d7c175d9e9bdde773a676729437f7c2cba622cc40c98dbc83bf38b84149e05aedf056639373"

getData \
    "src/R.sif" \
    "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119121&authkey=ALJdmNxRUqQ10B4"  \
    "43007b40347d8d521afd2cac2f99a8fcb39b74a1f0cc73c02c18e2da3fef17b9bf11b20e46b7b53ac05943e1e30a57fb5d58f18b1c3efc5500ce7e9bd892c096"

