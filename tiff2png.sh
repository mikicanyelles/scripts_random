#!/bin/bash
#
# DESCRIPTION
# Convert all images in a folder from .tiff to .png using the 'convert' utility
#   available on macOS and on Linux through the installation of imagemagick


for i in *.tiff; do
    convert $i ${i%.*}.png
    echo "$i done"
done
