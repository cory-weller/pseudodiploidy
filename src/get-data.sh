#!/usr/bin/env bash

# script to retrieve raw data


getFile () {
    destinationPath=${1}
    fileURL=${2}
    fileSHA512=${3}

    if [[ ! -f ${destinationPath} ]]; then
        wget -O ${destinationPath} ${fileURL}
        if [ $? -ne 0 ]; then
            echo "Failed to download $(basename ${destinationPath})"
        fi
    else
        echo "${destinationPath} already downloaded"
    fi

    if [[ -f ${destinationPath} ]]; then
        checksum=$(sha512sum ${destinationPath} | cut -d " " -f 1)
        if [[ ${checksum} != "${fileSHA512}" ]]; then
            echo "WARNING: ${destinationPath} checksum failed, contents may differ from those originally used."
            echo "File may be currupted outdated. Consider removing the file and running this script again."
        else
            echo "${destinationPath} checksum passed"
        fi
    fi
}


getFile \
    "data/external/SGD-features.tab" \
    "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab" \
    "cc3c54abb0106993b4d09227f4748e46460b4a7614518b6513863d7c175d9e9bdde773a676729437f7c2cba622cc40c98dbc83bf38b84149e05aedf056639373"
    # tab-delimited features file for S288C

getFile \
    "src/pseudodiploidy.sif" \
    "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119161&authkey=ALZOoR1P0_-ih_M"  \
    "af79d106d832f1aeffe2cb5bd19fc03228e745c9549bb1ce67e8ca3ee6a52111181d9f9a0c9f694d180b4d81e1e7119d020a4e0e73dd88b26c0d925699d77b99"
    # singularity image file containing tabix plus R tools for differential expression analysis

getFile \
    "data/external/S288C.vcf.zip" \
    "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119158&authkey=AEMV-omDAGOoe7U"  \
    "b3abea664cf2835b2451dc9649739aa3b392723e9cdb5188ff1edcd31c399a785254dfc9a3d6098ceb221ce3cd6ba0bc7863109be8111df54d914e823dbb2606"
    # Zipfile containing S288C vcfs, bgzipped and indexed by tabix
