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
    "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119113&authkey=AP5vqj2WsnZf1mw" \
    "dec245404d4e3b86371a4495d31fc9cf9c5a8299421a7763c4e75efde72ed7812cfef05f67df3da4905a7332fe44c5dae5ca16bc5f361c824acf40ba77711273"
