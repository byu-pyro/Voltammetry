# modified_python_script.py
print('py file opened')
import os
import sys
#import glob
from NW_functions import *
print('py file 1')
# Define the parent directory as the current directory
#parent_dir = os.getcwd()

# Print the current working directory
#print(f"Current working directory: {parent_dir}")

print('py file 2')
# Get the file name from the command line arguments
file_name = sys.argv[1]
print(file_name)
print('py file 3')

file = open(f'{file_name}', 'rb')
sim = pickle.load(file)

print('py file 4')

# Redirect stdout to a file
#sys.stdout = open(sim['file_name'] + '.out', 'w')
    
print('py file 5')


if sim['status'] == 'not started':
    print(f"running {file_name}")
    run(sim, False)
    print('py file 6a')
    
else:
    print(f"not running {file_name}")
    print(f"sim status is{sim['status']}")
    print('py file 6b')
    
print('py file 7')
t = sim['data']['t'][-1]
E = sim['data']['E'][-1]
i = sim['data']['i'][-1]
C = sim['data']['C'][-1, :, :]
N = sim['data']['N'][-1]
capture_rate = sim['data']['capture_rate']
compute_time = sim['data']['compute_time'][-1]

x = sim['parameters']['x']
T = sim['parameters']['T']
v = sim['parameters']['v']
n = sim['parameters']['n']
Eo = sim['parameters']['Eo']
D = get_D_ions(T, sim['parameters']['diffusivity_parameters_ions'])
z = sim['parameters']['z_ions']
stoic_rxn = sim['parameters']['stoic_rxn']
j_analyte = sim['parameters']['j_analyte']
E_start = sim['parameters']['E_start']
E_star = sim['parameters']['E_star']

ideal_solid_activity = sim['assumptions']['ideal_solid_activity']
ideal_ion_activity = sim['assumptions']['ideal_ion_activity']
no_migration = sim['assumptions']['no_migration']
no_ohmic_losses = sim['assumptions']['no_ohmic_losses']
iR_compensation = sim['assumptions']['iR_compensation']
No = sim['assumptions']['No']
dEo = sim['assumptions']['dEo']
L_RE = sim['assumptions']['L_RE']
R_ohmic = sim['assumptions']['R_ohmic']
print('py file 8')
fig, ax = plt.subplots(1, 1)
ax.set_xlabel(r'Potential/V vs. $E^*$')
ax.set_ylabel(r'Current Density/mA $cm^2$')
ax.set_title(sim['file_name'])
line, = ax.plot([], [])
line.set_data(sim['data']['E'] - E_star, sim['data']['i'])
ax.relim()
ax.autoscale_view(True, True, True)
fig.canvas.draw()
fig.canvas.flush_events()
fig.tight_layout()
# Set the limits of x and y axes based on your data
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
# Save the figure before showing it
fig.savefig(sim['file_name'] + '.png')
fig.show()
print('py file 9')

sys.stdout.close()
