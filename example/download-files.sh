#!/bin/bash
set -e
TEMP_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t 'holoaverage-files-temp')
echo "Temporarily downloading files into: $TEMP_DIR"
curl -L --max-redirs 5 --output "$TEMP_DIR/GaN_holographic_focal_series.zip" https://depositonce.tu-berlin.de/bitstreams/e8a777a0-110b-4232-9ed1-0959418724a6/download
unzip -d . "$TEMP_DIR/GaN_holographic_focal_series.zip" > /dev/null
if [ ! -f ./GaN_holographic_focal_series/README.txt ]; then
    echo "Something went wrong, expecting GaN_holographic_focal_series/README.txt file."
    exit 1
fi
rm -rf "$TEMP_DIR"
