#!/bin/tcsh
./cube_to_target \
--no_ridges \
--smoothing_over_ocean \
--grid_descriptor_file='J05_topo_grid.nc' \
--intermediate_cs_name='gmted2010_bedmachine-ncube0540-220518.nc' \
--output_grid='J05' \
--smoothing_scale=240.0 \
-u 'Nicholas Kevlahan kevlahan@mcmaster.ca' \
-q '.'
