##########################################################################################
##########################################################################################
######################################  LIBRARIES  #######################################
##########################################################################################
##########################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import simpson
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

##########################################################################################
##########################################################################################
######################################  INPUT DATA  ######################################
##########################################################################################
##########################################################################################

#subfolder_input_data_name = "examples/Plummer_example/Plummer_Model_np_(1000)_a_(1000.0 pc)_M_(1.0e+08 Msun)"
images_output_folder_name = "images"

##########################################################################################
##########################################################################################
####################################  USER PARAMTERS  ####################################
##########################################################################################
##########################################################################################


#The output units are kpc, Msun, years

conversion_velocity = 977775188.069539330899715424 # kpc/year to km/s; x kpc/year = x * conversion_velocity km/s
conversion_energy = 1.901572149306448e+54 * 10**-46# 1 [J] = 1.901572149306448e+54 [Msun kpc^2/year^2]
conversion_time = 1.0e-6 # from year to MegaYear

conversion_au = 2.063e+8 # 1kpc = 2.063e+8 au

conversion_distance = 1.0 # 1kpc 
conversion_to_vel_cm = 1.0 * conversion_velocity / (1.0e8) # * km/s

##########################################################################################
##########################################################################################
####################################  READING PARAMETERS  ################################
##########################################################################################
##########################################################################################

#QUANTIZER

times_example = np.array([])

kinetic_energy_error = np.array([])
potential_energy_error = np.array([])
total_energy_error = np.array([])
virial_energy_error = np.array([])


energies_error = np.array([])
module_VCM_error = np.array([])

rcut90 = np.array([])
rcut75 = np.array([])
rcut50 = np.array([])
rcut30 = np.array([])
rcut10 = np.array([])



	
# Reading data
name_output_folder = "output/Simulation_Maxdt[My](1000.0)_lv(5+12)_Np(1000)"
parameter_name = name_output_folder + '/Parameters.dat'
data_name = name_output_folder + '/Data/Data_'

f = open(parameter_name,"r")
lines = f.readlines()
f.close()

parameters = []
values = []
for i in range(len(lines)):
	#if (len(lines[i].split())==2):
	if (i < 38):
		parameters = parameters + [lines[i].split()[0]]
		if(float(lines[i].split()[1])==int(float(lines[i].split()[1]))):
			values = values + [int(float(lines[i].split()[1]))]
		else:
			values = values + [float(lines[i].split()[1])]

Nparticles = values[parameters.index('Nparticles_final')]
Iterations = values[parameters.index('Number_outputs')] + 1
vtype = values[parameters.index('vtype')]
datatype = np.dtype("f" + str(2**(vtype+1)))

Error_energy_Total = np.zeros(Iterations)
times = np.zeros(Iterations)

Kinetic = np.zeros(Iterations)
Potential = np.zeros(Iterations)
Total_energy = np.zeros(Iterations)
Virial_energy = np.zeros(Iterations)

Momentum_x = np.zeros(Iterations)
Momentum_y = np.zeros(Iterations)
Momentum_z = np.zeros(Iterations)

Error_velCM_x = np.zeros(Iterations)
Error_velCM_y = np.zeros(Iterations)
Error_velCM_z = np.zeros(Iterations)

Momentum_x_aux = np.zeros(Iterations)
Momentum_y_aux = np.zeros(Iterations)
Momentum_z_aux = np.zeros(Iterations)

cut90 = int((Nparticles+1)*90/100-1); 
cut75 = int((Nparticles+1)*75/100-1); 
cut50 = int((Nparticles+1)*50/100-1); 
cut30 = int((Nparticles+1)*30/100-1); 
cut10 = int((Nparticles+1)*10/100-1); 
		
rcut90 = np.zeros(Iterations)
rcut75 = np.zeros(Iterations)
rcut50 = np.zeros(Iterations)
rcut30 = np.zeros(Iterations)
rcut10 = np.zeros(Iterations)


#Reading the QUANTIZER data
for j in range(Iterations):
	data = np.fromfile(data_name + str(j) + '.bin', dtype = datatype)


	posX = data[1*Nparticles:2*Nparticles] * conversion_distance	#First
	posY = data[2*Nparticles:3*Nparticles] * conversion_distance	#Second
	posZ = data[3*Nparticles:4*Nparticles] * conversion_distance	#Third
	# velX = data[4*Nparticles:5*Nparticles] * conversion_velocity
	# velY = data[5*Nparticles:6*Nparticles] * conversion_velocity
	# velZ = data[6*Nparticles:7*Nparticles] * conversion_velocity


	r = np.sqrt(posX**2 + posY**2 + posZ**2)
	rsort = sorted(r)
	rcut90[j] = rsort[cut90]
	rcut75[j] = rsort[cut75]
	rcut50[j] = rsort[cut50]	
	rcut30[j] = rsort[cut30]
	rcut10[j] = rsort[cut10]

	Kinetic[j] = data[7*Nparticles] * conversion_energy	#Forth
	Potential[j] = data[7*Nparticles+1] * conversion_energy	#Fifth

	Virial_energy[j] = (2*data[7*Nparticles] + data[7*Nparticles+1])*conversion_energy

	Total_energy[j] = data[7*Nparticles+2] * conversion_energy	#Sixth

	Error_energy_Total[j] =  np.abs((Total_energy[j] - Total_energy[0])/Total_energy[0])

	#Momentum
	Momentum_x[j] = data[7*Nparticles+3]	#Seventh
	Momentum_y[j] = data[7*Nparticles+4]	#eigth
	Momentum_z[j] = data[7*Nparticles+5]	#ninth

	# Error in the velocity of the CM
	Error_velCM_x[j] =  (Momentum_x[j] - Momentum_x[0]) 
	Error_velCM_y[j] =  (Momentum_y[j] - Momentum_y[0]) 
	Error_velCM_z[j] =  (Momentum_z[j] - Momentum_z[0]) 

	times[j] = data[7*Nparticles+6] * conversion_time	#Tenth

	#Error_velCM_module = np.sqrt(Error_velCM_x**2 + Error_velCM_y**2 + Error_velCM_z**2) * conversion_to_vel_cm


##########################################################################################
##########################################################################################
###########################################  PLOTS  ######################################
##########################################################################################
##########################################################################################

##########################################################################################
##################################### Figure 1



# CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
#                   '#f781bf', '#a65628', '#984ea3',
#                   '#999999', '#e41a1c', '#dede00']

# legend_elements = [Line2D([0], [0], color='k', lw=4, label='p-p solution'),
#                    Line2D([0], [0],linestyle="None", marker='.', color='k', label='QUANTIZER', markersize=18)]


fontsize_legend = 25
fontsize_labels = 30
font_tick_size = 10
size_panel_title = 35
# title_text_position_x = 0.01
# title_text_position_y = 1.18

# title_text_position_x_sub4 = 0.01
# title_text_position_y_sub4 = 1.22

#legend_label_size = 1.5
#legend_label_text_sep = 0.2
labelspacing = 0.1

# legend_position_1 = (0.01,0.62)
# legend_position_2 = (0.01,0.62)
# legend_position_3 = (0.01,0.8)
# legend_position_4 = (0.01,0.76)



# second_legend_pos_1 = (0.68,0.62)
# second_legend_pos_2 = (0.68,0.62)

#vertical_space_between_subplots = 0
horizontal_space_between_subplots = 0.5

#tick_width = 2
tick_size = 15

#columnspacing_fig=0.6
borderpad_fig=0.05



fig = plt.gcf()
fig.set_size_inches(16, 10)

#Lagrangian Radii
ax = fig.add_subplot(221)
plt.title("Lagrangian radius",fontsize=size_panel_title)

plt.plot(times,rcut90,label='r90')
plt.plot(times,rcut75,label='r75')
plt.plot(times,rcut50,label='r50')
plt.plot(times,rcut30,label='r30')
plt.plot(times,rcut10,label='r10')

ax.legend(fontsize=fontsize_legend,labelspacing=labelspacing, borderpad=borderpad_fig)
plt.xlabel('Time [Myr]',fontsize=fontsize_labels)
plt.ylabel('[kpc]' ,fontsize=fontsize_labels)

plt.yticks(fontsize=font_tick_size)
plt.xticks(fontsize=font_tick_size)

#Energy
ax = fig.add_subplot(222)
plt.title("Energies",fontsize=size_panel_title)

plt.plot(times,Kinetic,label='Kinetic')
plt.plot(times,Potential,label='Potential')
plt.plot(times,Total_energy,label='Total_energy')
plt.plot(times,Virial_energy,label='Virial energy')

ax.legend(fontsize=fontsize_legend,labelspacing=labelspacing, borderpad=borderpad_fig)
plt.xlabel('Time [Myr]',fontsize=fontsize_labels)
plt.ylabel('[$10^{46}$ J]',fontsize=fontsize_labels)

plt.yticks(fontsize=font_tick_size)
plt.xticks(fontsize=font_tick_size)

#Energy Error
ax = fig.add_subplot(223)
plt.title("Energy error",fontsize=size_panel_title)

plt.plot(times,Error_energy_Total,label='Error Total energy')

ax.legend(fontsize=fontsize_legend,labelspacing=labelspacing, borderpad=borderpad_fig)
plt.xlabel('Time [Myr]',fontsize=fontsize_labels)

plt.yticks(fontsize=font_tick_size)
plt.xticks(fontsize=font_tick_size)

#System Momentum
ax = fig.add_subplot(224)
plt.title("System Momentum",fontsize=size_panel_title)

plt.plot(times,Momentum_x,label='Momentum x')
plt.plot(times,Momentum_y,label='Momentum y')
plt.plot(times,Momentum_z,label='Momentum z')

ax.legend(fontsize=fontsize_legend,labelspacing=labelspacing, borderpad=borderpad_fig)
plt.xlabel('Time [Myr]',fontsize=fontsize_labels)
plt.ylabel('[Msun kpc/year]',fontsize=fontsize_labels)

plt.yticks(fontsize=font_tick_size)
plt.xticks(fontsize=font_tick_size)


plt.subplots_adjust(hspace=horizontal_space_between_subplots)
plt.savefig(images_output_folder_name + "/example.png",bbox_inches='tight',dpi=300)
plt.show()	#Activite here to show the plots directly from python