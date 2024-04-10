# Programmable filter
# Select data sets to perform computations on using "command+click"
#
# Data
all_time     = inputs[0].CellData['mass_average']
partial_time = inputs[1].CellData['mass_average']

# Compute difference 
output.CellData.append(abs(all_time - partial_time), 'difference')

# Max difference
max_diff = max(output.CellData['difference'])
max_diff = scientific_notation="{:.1e}".format(max_diff)

# RMS difference
N = sum(val > 0 for val in output.CellData['difference'])

rms_all = sqrt(sum(inputs[0].CellData['mass_average']**2)/(N-1))

rms_diff = sqrt(sum(output.CellData['difference']**2)/(N-1)) / rms_all * 100 # in percent
rms_diff = scientific_notation="{:.1e}".format(rms_diff)

print(rms_diff, max_diff)  
