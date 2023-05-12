# script test_readColumnData.py
# -*- coding: utf-8 -*-

import numpy as np
import PhyPraKit as ppk

# -----example Code illustrating usage --------------------
if __name__ == "__main__":
  import sys, numpy as np, matplotlib.pyplot as plt
  from PhyPraKit import odFit, labxParser, readColumnData
  from scipy import interpolate
  from scipy import signal
  import kafe
  from kafe import ASCII, LaTeX, FitFunction
  from kafe.function_library import linear_2par
  


I=(0.01,0.04,0.16,0.49)
xdata=I


beta=(0.0164864152,0.0474360799,0.1678556441,0.4980079681)

bkor=np.zeros(4)
for i in range(4):
    bkor[i]=(beta[i]-0.0041797283)
ydata=bkor




yer=0.00000001

def generate_datasets(output_file_path1):
    '''The following block generates the Datasets and writes a file for
    each of them.'''

    import numpy as np  # need some functions from numpy

    my_datasets = []


    my_datasets.append(kafe.Dataset(data=(xdata, ydata)))
    my_datasets[-1].add_error_source('y', 'simple', yer)
    
    my_datasets[0].write_formatted(output_file_path1)
    

############
# Workflow #
############

# Generate the Dataseta and store them in files

#generate_datasets('Aufg2_bkor.dat')

# Initialize the Datasets
my_datasets = [kafe.Dataset(title=r'bkor')]

# Load the Datasets from files
my_datasets[0].read_from_file(input_file='Aufg2_bkor.dat')
# Create the Fits
my_fits = [kafe.Fit(dataset,
                    linear_2par,
                    fit_label="Linear regression " )
           for dataset in my_datasets]

# Do the Fits
for fit in my_fits:
    fit.do_fit()

# Create the plots
my_plot = kafe.Plot(my_fits[0])

mn  = my_fits[0].get_parameter_values(rounding=False)



#mn=np.around(mn, decimals=3)
print mn
#my_plot.axes.annotate(r'$\pm_{eich} $' + str(mnUB), xy=(3.82, 4.229), size=10, ha='left')

my_plot.axis_labels = ['$Strom \ I^2 \  [A^2]$', r'$ \beta_{kor} \  [1/s]$']


# Draw the plots
my_plot.plot_all()


###############
# Plot output #
###############

# Save the plots
my_plot.save('Aufgabe2_bkor.pdf')




# Show the plots
#my_plot.show()
