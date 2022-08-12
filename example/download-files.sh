#!/bin/bash
set -e
TEMP_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t 'holoaverage-files-temp')
echo "Temporarily downloading files into: $TEMP_DIR"
curl --output "$TEMP_DIR/GaN_holographic_focal_series.zip" https://depositonce.tu-berlin.de/bitstream/11303/7448/2/GaN_holographic_focal_series.zip
unzip -d . "$TEMP_DIR/GaN_holographic_focal_series.zip" > /dev/null
if [ ! -f ./GaN_holographic_focal_series/README.txt ]; then
    echo "Something went wrong, expecting GaN_holographic_focal_series/README.txt file."
    exit 1
fi
rm -rf "$TEMP_DIR"
