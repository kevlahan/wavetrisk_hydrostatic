# Highlight data using CMD+click
approx = inputs[0].CellData['OMEGA_average']
exact  = inputs[1].CellData['mass_average']

# Compute difference 
output.CellData.append(approx - exact, 'difference')
diff = output.CellData['difference']

# RMS difference
N = sum(val > 0 for val in output.CellData['difference'])

rms_exact = sqrt(sum(exact**2)/(N-1))
rms_diff  = sqrt(sum(diff**2)/(N-1)) / rms_exact * 100 # in percent

output = scientific_notation="{:.1f}".format(rms_diff)
print("rms error = "+output+"%")   
