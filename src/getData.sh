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
    "data/external/SGD_features.tab" \
    "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab" \
    "cc3c54abb0106993b4d09227f4748e46460b4a7614518b6513863d7c175d9e9bdde773a676729437f7c2cba622cc40c98dbc83bf38b84149e05aedf056639373"

getData \
    "src/R.sif" \
    "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119114&authkey=AH63WVNezJbyWqM"  \
    "5e1acf95474bdfa14e9a25f46002186f4b724ae7457e7f7c18b02c702a0d4ef3338fb1eac5ce1c78d25d622379ce80f97320cd94c10697ba831034bd66eeab4f"
