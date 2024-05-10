import os
import sys
import numpy as np
from NW_functions import *
import matplotlib.pyplot as plt

# Redirect standard output to a file
sys.stdout = open('QA_Output.txt', 'w')

print( " ----- Quality Assment Post Processing Has Begun ----- ")

ID_List =  np.array([f'0000',f'0001', f'0010', f'0100', f'1000', f'0100b', f'1000b', f'1111']) #np.array('1111', '0o0001')
IC_List =   np.array([f'HF_' , f'HS_' , f'LF_' , f'LS_'])     #np.array('HS_', "LS_")

for j in range(len(IC_List)):
    IC = IC_List[j]
    for i in range(len(ID_List)):
        ID = ID_List[i]
        test_file = f"NW_QCtest_{IC}{ID}"
        ref_file = f"QCbase_{IC}{ID}"
        print("===============================================================================")
        print("-------------------------------------------------------------------------------")
        print("===============================================================================")
        print(f"------------------- Preparing file: {test_file} --------------------------")
        print("===============================================================================")
        print("-------------------------------------------------------------------------------")
        print("===============================================================================")
        print("   ")
        print("   ")
        print("   ")
        print("file ", test_file, " is being processed...")
        
        test_file_path = f"./{test_file}"
        ref_file_path = f"../QCbase/{ref_file}"
        try:
            # Get the size of the file
            test_file_size = os.path.getsize(test_file_path)
            ref_file_size = os.path.getsize(ref_file_path)
            
            print(f"The size of the test file is {test_file_size} bytes.")
            print(f"The size of the test file is {ref_file_size} bytes.")
            
            file1 = open(test_file_path, 'rb')
            print('file2 opened')
            file2 = open(ref_file_path, 'rb')
            print('file2 opened')
            
            sim = pickle.load(file1)
            print('file1 loaded')
            ref = pickle.load(file2)
            print('file2 loaded')
            
            ## Data from New Functions module
            t = sim['data']['t'][-1]
            E = sim['data']['E'][-1]
            i = sim['data']['i'][-1]
            C = sim['data']['C'][-1, :, :]
            N = sim['data']['N'][-1]
            capture_rate = sim['data']['capture_rate']
            compute_time = sim['data']['compute_time'][-1]

            print('sim variables saved')

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

            print('sim parameters saved')

            ideal_solid_activity = sim['assumptions']['ideal_solid_activity']
            ideal_ion_activity = sim['assumptions']['ideal_ion_activity']
            no_migration = sim['assumptions']['no_migration']
            no_ohmic_losses = sim['assumptions']['no_ohmic_losses']
            iR_compensation = sim['assumptions']['iR_compensation']
            No = sim['assumptions']['No']
            dEo = sim['assumptions']['dEo']
            L_RE = sim['assumptions']['L_RE']
            R_ohmic = sim['assumptions']['R_ohmic']

            print('sim assumptions saved')
            
            sim_variables = [t, E, i, C, N, capture_rate, compute_time]
            print('yes1')
            sim_parameters = [x, T, v, n, Eo, D, z, stoic_rxn, j_analyte, E_start, E_star]
            print('yes2')
            sim_assumptions = [ideal_solid_activity, ideal_ion_activity, no_migration, no_ohmic_losses, iR_compensation, No, dEo, L_RE, R_ohmic]

            print('sim arrays created')
            
            ## Data from Controls 
            rt = ref['data']['t'][-1]
            rE = ref['data']['E'][-1]
            ri = ref['data']['i'][-1]
            rC = ref['data']['C'][-1, :, :]
            rN = ref['data']['N'][-1]
            rcapture_rate = ref['data']['capture_rate']
            rcompute_time = ref['data']['compute_time'][-1]

            print('ref variables saved')

            rx = ref['parameters']['x']
            rT = ref['parameters']['T']
            rv = ref['parameters']['v']
            rn = ref['parameters']['n']
            rEo = ref['parameters']['Eo']
            rD = get_D_ions(T, ref['parameters']['diffusivity_parameters_ions'])
            rz = ref['parameters']['z_ions']
            rstoic_rxn = ref['parameters']['stoic_rxn']
            rj_analyte = ref['parameters']['j_analyte']
            rE_start = ref['parameters']['E_start']
            rE_star = ref['parameters']['E_star']
            print('ref parameters saved')

            rideal_solid_activity = ref['assumptions']['ideal_solid_activity']
            rideal_ion_activity = ref['assumptions']['ideal_ion_activity']
            rno_migration = ref['assumptions']['no_migration']
            rno_ohmic_losses = ref['assumptions']['no_ohmic_losses']
            riR_compensation = ref['assumptions']['iR_compensation']
            rNo = ref['assumptions']['No']
            rdEo = ref['assumptions']['dEo']
            rL_RE = ref['assumptions']['L_RE']
            rR_ohmic = ref['assumptions']['R_ohmic']
            print('ref assumptions saved')
            
            ref_variables = [rt, rE, ri, rC, rN, rcapture_rate, rcompute_time]
            ref_parameters = [rx, rT, rv, rn, rEo, rD, rz, rstoic_rxn, rj_analyte, rE_start, rE_star]
            ref_assumptions = [rideal_solid_activity, rideal_ion_activity, rno_migration, rno_ohmic_losses, riR_compensation, rNo, rdEo, rL_RE, rR_ohmic]
            print('ref arrays created')
            print("  ")
            print("---------------------------")
            print("Variable IDs are as follows:")
            print("0 = t")
            print("1 = E")
            print("2 = i")
            print("3 = C")
            print("4 = N")
            print("5 = capture_rate")
            print("6 = compute_time")
            print("-----------------------------")
            


            
            print("   ")
            print("-----------------------------------------------------------------------------")
            print("====================  Parameters and Assumptions  ======== [ ", f"{test_file}", " ] =====")
            print("-----------------------------------------------------------------------------")
            print("      ")
            
            print("Temprature: ", f"Sim - {T} K , ", f"Ref - {rT} K ") 
            print("Scan Rate: ", f"Sim - {v} V/s , ", f"Ref - {rv} V/s ") 
            print("Number of Electrons Transfered: ", f"Sim - {n} , ", f"Ref - {rn} ")
            print("E_start: ", f"Sim - {E_start} , ", f"Ref - {rE_start} ")
            print("E_star: ", f"Sim - {E_star} , ", f"Ref - {rE_star} ")
            
            print("Difusion Parameters:")
            print("----Sim----")
            print(D)
            print("----Ref----")
            print(rD)
            print("  ")
            
            print("Z ions:")
            print("----Sim----")
            print(z)
            print("----Ref----")
            print(rz)
            print("  ")
            
            print("Stochiometry:")
            print("----Sim----")
            print(stoic_rxn)
            print("Analyte: ", j_analyte )
            print("----Ref----")
            print(rstoic_rxn)
            print("Analyte: ", rj_analyte )
            print("  ")
            

            print("Grid Spacing")
            print("----Sim----")
            print(x)
            print("----Ref----")
            print(rx)
            print("  ")
                
            print("Assumption IDs are as follows: ")
            print("0 = ideal_solid_activity")
            print("1 = ideal_ion_activity")
            print("2 = no_migration")
            print("3 = no_ohmic_losses")
            print("4 = No")
            print("5 = dEo")
            print("6 = L_RE")
            print("7 = R_ohmic")
            print("-----------------------------")
            
            print("Assumptions:")
            print("----Sim----")
            print(sim_assumptions)
            print("----Ref----")
            print(ref_assumptions)
            print("  ")

            
            print("   ")
            print("-----------------------------------------------------------------------------")
            print("====================  Variable Comparison Begins  ======== [ ", f"{test_file}", " ] =====")
            print("-----------------------------------------------------------------------------")
            print("      ")
            print("      ")
            
            print("Variable IDs are as follows:")
            print("0 = t")
            print("1 = E")
            print("2 = i")
            print("3 = C")
            print("4 = N")
            print("5 = capture_rate")
            print("6 = compute_time")
            print("-----------------------------")
            
            
            # Create an empty list to store the SSEs
            SSEs = []

            for k in range(len(sim_variables)):
                print("Variable ID: ", "[", k, "]")
                Err = sim_variables[k] - ref_variables[k]
                SE = Err**2
                SSE = np.sum(SE)
                # Append the SSE to the list
                SSEs.append(SSE)
                print(f"For Variable ID {k} the SSE is {SSE}")

            # Now, SSEs is a list of all SSEs
            print(f"All SSEs: {SSEs}")
            print(test_file)
            print("Generating Bar Graph...")
            plt.figure()  # Create a new figure
            plt.bar(range(len(SSEs)), SSEs)
            plt.xlabel('Variable ID')
            plt.ylabel('SSE')
            plt.title(f'Bar Graph of all the SSEs for: {test_file}')
            plt.savefig(f'SSE_{test_file}.png')
            plt.show()
            
            print("Generating Potential V Current Plot...")
            plt.figure()  # Create a new figure
            fig, ax = plt.subplots(1, 1)
            ax.set_xlabel(r'Potential/V vs. $E^*$')
            ax.set_ylabel(r'Current Density/mA $cm^2$')
            ax.set_title(sim['file_name'])

            # Plot sim data
            line_sim, = ax.plot(sim['data']['E'] - E_star, sim['data']['i'], 'r:', label='sim')
            # Plot ref data
            line_ref, = ax.plot(ref['data']['E'] - rE_star, ref['data']['i'], 'b--', label='ref')

            ax.relim()
            ax.autoscale_view(True, True, True)
            fig.canvas.draw()
            fig.canvas.flush_events()
            fig.tight_layout()

            # Add a legend
            ax.legend()

            # Save the figure before showing it
            fig.savefig(f'Comp_{test_file}.png')
            fig.show()
            
            print("file end...")
            print("   ")
            print("   ")
            print("   ")
            print("   ")
            print("   ")

    
            
    
        except:
            print(f"Error: The file {test_file} or {ref_file} does not exist. Moving to the next file. Or Some Other Error Occurred.")
            #text = f"Error: {e}"
            #print(text)
            print("   ")
            print("   ")
            print("   ")
            print("   ")
            print("   ")


