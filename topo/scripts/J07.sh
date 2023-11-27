#!/bin/tcsh
./cube_to_target \
--no_ridges \
--smoothing_over_ocean \
--grid_descriptor_file='J07_topo_grid.nc' \
--source_data_identifier='gmted2010_modis_bedmachine-ncube0540-220518.nc' \
--intermediate_cs_name='gmted2010_bedmachine-ncube0540-220518.nc' \
--output_grid='J07' \
--smoothing_scale=60.0 \
-u 'Nicholas Kevlahan kevlahan@mcmaster.ca' \
-q '.'
