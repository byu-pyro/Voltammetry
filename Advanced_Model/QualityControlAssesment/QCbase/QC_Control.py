##Great idea, I would recommend the following simulations

 
#Run each of the following at 
#* high concentration (10 mol% UCl3)
#* low concentration (0.1 wt% UCl3) 
#* fast scan rate (1 V/s) 
#* low scan rate (10 mV/s)
 

#- Diffusion only (no migration, ideal activites, no ohmic losses)  ###0000
#- Test solid activity only                                         ###0001
#- Test ion activity only                                           ###0010
#- Test migration only                                              ###0100
#- Test migration with complexation                                 ###0100b
#- Test ohmic losses                                                ###1000
#- Test ohmic losses with 85% IR compensation                       ###1000b
#- Test solid activity, ion activity, migration with complexation, and ohmic losses with IR compensation. ###1111
 
#Then lets say you add radial, you would add a radial test on its own and then add radial to the number 8.

from functions import *


##### High Concentartion Fast Scan ######

# salts used
salts = np.array(["LiCl", "KCl", "UCl3"])

# molecular weight of each salt
# salt mass fractions     ??? This is from Michaels Raw code, keping the Li/Kcl ratio from that raw code y_salts = np.array([0.287, 0.703, 0.01])

MW_salts = np.array([42.394, 74.551, 344.388])  # g/mol
OGy_salts = np.array([0.287, 0.703, 0.01])

mol0 = OGy_salts[0]/MW_salts[0]
mol1= OGy_salts[1]/MW_salts[1]
mol2= OGy_salts[2]/MW_salts[2]

totalmol = mol0+mol1+mol2
print(totalmol)

molpecent_salts= np.array([mol0, mol1, mol2])
print(molpecent_salts)

mol_balanced0 = mol0/(mol1+mol0)

molper_UCl3 = 0.1
molper_other_salt = 1-molper_UCl3
molper_Li = molper_other_salt*mol_balanced0 
molper_K = molper_other_salt*(1-mol_balanced0)
print(molper_Li+molper_K+molper_UCl3)


mass_UCl3 = molper_UCl3*MW_salts[2]
mass_Li = molper_Li*MW_salts[0]
mass_K = molper_K*MW_salts[1]
new_totalmass = mass_UCl3+mass_Li +mass_K
new_massper_UCl3 =mass_UCl3/new_totalmass
new_massper_Li = mass_Li/new_totalmass
new_massper_K = mass_K/new_totalmass

y_salts = np.array([new_massper_Li, new_massper_K, new_massper_UCl3])
print(y_salts)


# density parameters of the salt based on a linear function (density at 0K, density change per degree K)
density_parameters_salts = np.array([
    [1.88, -4.32E-04],  # LiCl
    [2.14, -5.83E-04],  # KCl
    [13.66, -79.4E-04, ],  # UCl3
])

# ions in the salts
ions = np.array(["Li+", "K-", "[UCl6]3-", "Cl-"])

# charge of the ions
z_ions = np.array([1, 1, -3, -1])

# stoichiometry of the ions in the salt
stoic_ions = np.array([
    # Li+,      K+,     [UCl6]3-,   Cl+
    [1,         0,      0,          1],  # LiCl
    [0,         1,      0,          1],  # KCl
    [0,         0,      1,          -3],  # UCl3, this has -3 for chlorine because it forms the complex
])

# diffusivity parameters of the salt based on the arrheneus equation(pre-exponential factor, activation energy)
diffusivity_parameters_ions = np.array([
    [1.65E-3, -2.516E4],  # Li+
    [1.93E-3, -2.812E4],  # K+
    [1.37E-3, -2.420E4],  # Cl-
    [5.92E-4, -2.044E4],  # [UCl6]3-
])

# Temperature
T_C = 500  # C
T = T_C + 273.15  # K

# The scan rate
v = 1.00 #V/s

# number of electrons transferred
n = 3

# standard potential
Eo = 0.0

# the index of the analyte
j_analyte = 2

# stoichiometry of the reduction  reaction
stoic_rxn = np.array([0, 0, -1, 6])

# Grid Space
x = get_gridspace(L=0.5, n_x=70, alpha=1.1)

# how often do you want a data point (in V)
capture_rate = 0.0001

def simcreation(IC,ID):
    if ID == f'0000':
        ideal_solid_activity = True #assume ideal solid activity?
        ideal_ion_activity = True
        no_migration = True #assume no migration?
        no_ohmic_losses = True
        
        # percent of estimated ohmic resistance to account for
        iR_compensation = 0 #0.85 #make sure no_ohmic losses is False if not 0

    elif ID == f'0001':
        ideal_solid_activity = False #assume ideal solid activity?
        ideal_ion_activity = True
        no_migration = True #assume no migration?
        no_ohmic_losses = True
                # percent of estimated ohmic resistance to account for
        iR_compensation = 0 #0.85 #make sure no_ohmic losses is False if not 0
        
    elif ID == f'0010':
        ideal_solid_activity = True #assume ideal solid activity?
        ideal_ion_activity = False
        no_migration = True #assume no migration?
        no_ohmic_losses = True
                # percent of estimated ohmic resistance to account for
        iR_compensation = 0 #0.85 #make sure no_ohmic losses is False if not 0
        
    elif ID == f'0100':
        ideal_solid_activity = True #assume ideal solid activity?
        ideal_ion_activity = True
        no_migration = False #assume no migration?
        no_ohmic_losses = True
                # percent of estimated ohmic resistance to account for
        iR_compensation = 0 #0.85 #make sure no_ohmic losses is False if not 0
        
    elif ID == f'0100b':
        ideal_solid_activity = True #assume ideal solid activity?
        ideal_ion_activity = True
        no_migration = False #assume no migration?
        no_ohmic_losses = True
                # percent of estimated ohmic resistance to account for
        iR_compensation = 0 #0.85 #make sure no_ohmic losses is False if not 0
        
        #Insert Complexation Argument Once Michael Responds
        
    elif ID == f'1000':
        ideal_solid_activity = True #assume ideal solid activity?
        ideal_ion_activity = True
        no_migration = True #assume no migration?
        no_ohmic_losses = False
                # percent of estimated ohmic resistance to account for
        iR_compensation = 1 #0.85 #make sure no_ohmic losses is False if not 0
        
    elif ID == f'1000b':
        ideal_solid_activity = True #assume ideal solid activity?
        ideal_ion_activity = True
        no_migration = True #assume no migration?
        no_ohmic_losses = False
                # percent of estimated ohmic resistance to account for
        iR_compensation = 0.85 #make sure no_ohmic losses is False if not 0

    elif ID == f'1111':
        ideal_solid_activity = False #assume ideal solid activity?
        ideal_ion_activity = False
        no_migration = False #assume no migration?
        no_ohmic_losses = False
                # percent of estimated ohmic resistance to account for
        iR_compensation = 1 #0.85 #make sure no_ohmic losses is False if not 0
        
    else:
        print('ID not recognized')
        pass


    # percent of estimated ohmic resistance to account for
    iR_compensation = 0 #0.85 #make sure no_ohmic losses is False if not 0

    # amount of deposit to have full coverage of the electrode
    No = 4e-9  # mol/cm^2

    # difference between the pure and dilute reference state
    dEo = -0.113  # V

    # distance between the working electrode and reference electrode
    L_RE = 0.25 # cm



    # Calculatin Input
    # --------------------------------------------------------------------------------------------------------


    # initial concentration of all the ions
    C_ions = get_C_ions(T, y_salts, MW_salts, density_parameters_salts, stoic_ions)

    # resistance between the working electrode and the reference electrode
    R_ohmic = get_R_ohmic(L_RE, C_ions, T, z_ions, diffusivity_parameters_ions)

    # the reference potential for the plots
    E_star = get_E_star(C_ions[j_analyte], Eo, n, T, ideal_ion_activity, dEo)

    # The starting potential for the scans
    E_start = E_star + 0.3

    # get the amount of deposit at the surface
    N = np.linspace(No / 1e8, No, 10000)  # an array of guess values
    # E at all of these N
    N_guess = N[np.argmin(np.abs(Nernst(N, C_ions[j_analyte], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)
                                 - E_start))]  # closest guess
    N_start = get_N(N_guess, E_start, C_ions[j_analyte], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)

    #Initial concentraiton profile
    C_start = np.ones((len(ions), len(x))) * C_ions[:, None]


    # Initialize Simulation File
    # --------------------------------------------------------------------------------------------------------------------
    sim = {
        'status': 'not started',
        'file_name': f'QCbase{IC}{ID}',
        'path': './',
        'parameters': {
            'T': T,
            'salts': salts,
            'y_salts': y_salts,
            'MW_salts': MW_salts,
            'density_parameters_salts': density_parameters_salts,
            'ions': ions,
            'z_ions': z_ions,
            'stoic_ions': stoic_ions,
            'diffusivity_parameters_ions': diffusivity_parameters_ions,
            'j_analyte': j_analyte,
            'E_start': E_start,
            'v': v,
            'n': n,
            'Eo': Eo,
            'stoic_rxn': stoic_rxn,
            'x': x,
            'E_star': E_star},
        'assumptions': {
            'ideal_solid_activity': ideal_solid_activity,
            'ideal_ion_activity': ideal_ion_activity,
            'no_migration': no_migration,
            'no_ohmic_losses': no_ohmic_losses,
            'iR_compensation': iR_compensation,
            'No': No,
            'dEo': dEo,
            'L_RE': L_RE,
            'R_ohmic': R_ohmic},
        'data': {
            'capture_rate': capture_rate,
            't': np.array([0.0]),  # 1D array for time
            'E': np.array([E_start]),  # 1D array from potential
            'i': np.array([0.0]),  # 1D array for current
            'C': np.array([C_start]),  # 3D array of concentration (time, ions, distance)
            'N': np.array([N_start]),
            'compute_time': np.array([0.0])}}  # 1D array from amount of deposit

    file = open(f'{sim["path"]}{sim["file_name"]}', 'wb')
    pickle.dump(sim, file)

ID_List = np.array([f'0000',f'0001', f'0010', f'0100', f'1000', f'0100b', f'1000b', f'1111'])

for i in range(len(ID_List)):
    IC = f'_HF_'
    ID = ID_List[i]
    simcreation(IC,ID)
    


## -- Slow Scan -- ##

# The scan rate
v = 0.010 #V/s

for i in range(len(ID_List)):
    IC = f'_HS_'
    ID = ID_List[i]
    simcreation(IC,ID)

##### Low Concentration Slow Scan #####

OGy_salts = np.array([0.287, 0.703, 0.01])
massfper = OGy_salts[0]/(OGy_salts[0]+OGy_salts[1])
massfU = 0.001
massfrest = 1-massfU 
massfLi = massfper*massfrest
massfK = (1-massfper)*massfrest
y_salts = np.array([massfLi, massfK, massfU])

for i in range(len(ID_List)):
    IC = f'_LS_'
    ID = ID_List[i]
    simcreation(IC,ID)

## -- Fast Scan -- ##

v = 1.0 #V/s
for i in range(len(ID_List)):
    IC = f'_LF_'
    ID = ID_List[i]
    simcreation(IC,ID)