#!/usr/bin/env python3
# import time
import numpy as np
import matplotlib.pyplot as plt
import yaml
from yaml.loader import SafeLoader
import time
import imageio
import glob
import os

#from Landlab
from landlab import RasterModelGrid, imshow_grid, imshowhs_grid
from landlab.io import read_esri_ascii
from landlab.io.netcdf import write_raster_netcdf

#Hillslope geomorphology
from landlab.components import ExponentialWeatherer
from landlab.components import DepthDependentTaylorDiffuser
from landlab.components import DepthDependentDiffuser

#Fluvial Geomorphology and Flow routing
from landlab.components import FlowDirectorMFD #trying the FlowDirectorMFD
from landlab.components import FlowAccumulator, Space, FastscapeEroder, PriorityFloodFlowRouter
from landlab.components.space import SpaceLargeScaleEroder

from ss_fault_function import ss_fault

#reading the file with parameters
config = yaml.safe_load(open('parameters.yaml','r'))#,Loader=yaml.FullLoader)

model_name= config['saving']['model_name']
alt_name= config['comments']['alt_name']
save_location = 'output_model_run/%s'%model_name # location to save output

try:
    os.chdir(save_location)
except:
    os.mkdir('output_model_run/%s'%model_name)
    os.chdir(save_location)

print('saving into: %s' %(os.getcwd()), file=open('out_%s.txt' %model_name, 'w')) # print to a file

#grid
ymax=config['shape']['ymax']
xmax=config['shape']['xmax']
dxy=config['shape']['dxy']

nrows = int(ymax/dxy)
ncols = int(xmax/dxy)

#TECTONICS AND TIME PARAMETERS
total_slip_1= config['tectonics']['total_slip_1']
total_slip_2= config['tectonics']['total_slip_2']
total_model_time= config['time']['total_model_time']
dt= config['time']['dt']
iterations= np.arange(0,total_model_time,dt)
print(iterations)
desired_slip_per_event_1=(total_slip_1/total_model_time)*dt
desired_slip_per_event_2=(total_slip_2/total_model_time)*dt
shrink = 0.5


# Geomorphic parameters

# uplift
uplift_rate= config['geomorphology']['uplift_rate']

#Hillsope Geomorphology for DDTD component
Sc=config['geomorphology']['Sc']
Hstar= config['geomorphology']['Hstar'] # characteristic transport depth, m
V0= config['geomorphology']['V0'] #transport velocity coefficient changed this
D= Hstar*V0#V0 *Hstar  #effective(maximum) diffusivity
P0=config['geomorphology']['P0']
run_off=config['geomorphology']['run_off']

#Fluvial Erosion for SPACE Large Scale Eroder
K_sed=config['geomorphology']['K_sed'] #sediment erodibility
K_br= config['geomorphology']['K_br'] #bedrock erodibility
F_f=config['geomorphology']['F_f']#fraction of fine sediment
phi= config['geomorphology']['phi'] #sediment porosity
H_star=config['geomorphology']['H_star'] #sediment entrainment lenght scale
Vs= config['geomorphology']['Vs'] #velocity of sediment
m_sp= config['geomorphology']['m_sp'] #exponent ondrainage area stream power
n_sp= config['geomorphology']['n_sp'] #exponent on channel slope in the stream power framework
sp_crit_sed= config['geomorphology']['sp_crit_sed'] #sediment erosion threshold
sp_crit_br= config['geomorphology']['sp_crit_br'] #bedrock erosion threshold

#Initialize video maker
# writer = imageio.get_writer(f'/Users/taranguiz/Research/Lateral_advection/output_model_run/{model_name}/{alt_name}.mp4', fps=20)

def add_file_to_writer(filename):
    image=imageio.imread(filename)
    writer.append_data(image)

def get_file_sequence(time_step):
    num_digits = len(str(total_model_time))
    file_sequence_num = time_step * 10 ** -num_digits
    float_convert = f'{{:.{num_digits}f}}'
    file_sequence_num = float_convert.format(file_sequence_num).split('.')[1]
    return file_sequence_num

#READING STEADY STATE TOPO
(mg,z)=read_esri_ascii('/Users/taranguiz/Research/NZ_Tararua_SSF/steady_state_1/finaltopo_topographic__elevation.asc', name="topographic__elevation")
(mg1,soil_0)=read_esri_ascii('/Users/taranguiz/Research/NZ_Tararua_SSF/steady_state_1/finaltopo_soil__depth.asc', name='soil__depth')
(mg2,bed_0)=read_esri_ascii('/Users/taranguiz/Research/NZ_Tararua_SSF/steady_state_1/finaltopo_bedrock__elevation.asc', name='bedrock__elevation')

mg.add_field("soil__depth", soil_0, at='node')
mg.add_field("bedrock__elevation", bed_0, at='node')

mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=False, left_is_closed=True,
                                       right_is_closed=True, top_is_closed=True)
shrink = 0.5
imshow_grid(mg,z, cmap='viridis', shrink=shrink)
plt.title('Initial Topography')
plt.show()
# write_raster_netcdf(f'{model_name}0000000.nc', mg, format="NETCDF4")
# initial_topo_img = f'/Users/taranguiz/Research/Lateral_advection/output_model_run/{model_name}/{model_name}0000000.png'
# plt.savefig(initial_topo_img, dpi=300, facecolor='white')
# add_file_to_writer(initial_topo_img)

soil=mg.at_node['soil__depth']
bed= mg.at_node['bedrock__elevation']

# print (mg.at_node['topographic__elevation'][mg.core_nodes])
# print (mg.at_node['bedrock__elevation'][mg.core_nodes])
# print (mg.at_node['soil__depth'][mg.core_nodes])

fault_loc_y_1=int(mg.number_of_node_rows / 4.)
fault_nodes_1 = np.where(mg.node_y==(fault_loc_y_1*10))[0]
print(fault_nodes_1)

fault_loc_y_2=int(mg.number_of_node_rows / 1.5)
fault_nodes_2 = np.where(mg.node_y==(fault_loc_y_2*10))[0]
print(fault_nodes_2)

#instantiate components
expweath=ExponentialWeatherer(mg, soil_production__maximum_rate=P0, soil_production__decay_depth=Hstar)

# Hillslope with Taylor Diffuser
ddtd=DepthDependentTaylorDiffuser(mg,slope_crit=Sc,
                                  soil_transport_velocity=V0,
                                  soil_transport_decay_depth=Hstar,
                                  nterms=2,
                                  dynamic_dt=True,
                                  if_unstable='warn', courant_factor=0.1)
#Flow Router
fr=PriorityFloodFlowRouter(mg, flow_metric='D8', suppress_out=True, runoff_rate=run_off)
#SPACE Large Scale
space= SpaceLargeScaleEroder(mg,
                             K_sed=K_sed,
                             K_br=K_br,
                            F_f=F_f,
                            phi=phi,
                            H_star=H_star,
                            v_s=Vs,
                            m_sp=m_sp,
                            n_sp=n_sp,
                            sp_crit_sed=sp_crit_sed,
                             sp_crit_br=sp_crit_br)


# #fluvial array (look up table)
# fluvial_freq=config['climate']['fluvial_freq'] #how often the humid period occurs
# fluvial_len=config['climate']['fluvial_len'] #how long the humid period last
# fluvial_0=np.arange(fluvial_freq,total_model_time, fluvial_freq)
# fluvial_n=fluvial_0+fluvial_len
# fluvial_times=np.vstack((fluvial_0,fluvial_n)).T
# print(fluvial_times)


#### print parameters to file ####
print('name of the model is: %s' %alt_name, file=open('out_%s.txt' %model_name, 'a'))
print('uplift_rate: %s' %uplift_rate, file=open('out_%s.txt' %model_name, 'a'))
print('Sc: %s' %Sc, file=open('out_%s.txt' %model_name, 'a'))
print('Hstar: %s' %Hstar, file=open('out_%s.txt' %model_name, 'a'))
print('V0: %s' %V0, file=open('out_%s.txt' %model_name, 'a'))
print('P0: %s' %P0, file=open('out_%s.txt' %model_name, 'a'))
print('run_off: %s' %run_off, file=open('out_%s.txt' %model_name, 'a'))
print('K_sed: %s' %K_sed, file=open('out_%s.txt' %model_name, 'a'))
print('K_br: %s' %K_br, file=open('out_%s.txt' %model_name, 'a'))
print('F_f: %s' %F_f, file=open('out_%s.txt' %model_name, 'a'))
print('phi: %s' %phi, file=open('out_%s.txt' %model_name, 'a'))
print('H_star: %s' %H_star, file=open('out_%s.txt' %model_name, 'a'))
print('Vs: %s' %Vs, file=open('out_%s.txt' %model_name, 'a'))
print('m_sp: %s' %m_sp, file=open('out_%s.txt' %model_name, 'a'))
print('n_sp: %s' %n_sp, file=open('out_%s.txt' %model_name, 'a'))
print('sp_crit_sed: %s' %sp_crit_sed, file=open('out_%s.txt' %model_name, 'a'))
print('sp_crit_br: %s' %sp_crit_br, file=open('out_%s.txt' %model_name, 'a'))
# print('total_slip: %s' %total_slip, file=open('out_%s.txt' %model_name, 'a'))
# print('method: %s' %method, file=open('out_%s.txt' %model_name, 'a'))
print('total_model_time: %s' %total_model_time, file=open('out_%s.txt' %model_name, 'a'))
print('dt: %s' %dt, file=open('out_%s.txt' %model_name, 'a'))
# print('fluvial_freq: %s' %fluvial_freq, file=open('out_%s.txt' %model_name, 'a'))
# print('fluvial_len: %s' %fluvial_len, file=open('out_%s.txt' %model_name, 'a'))
print('',file=open('out_%s.txt' %model_name, 'a'))


#things for the loop
time=0 #time counter
f=0 #index counter for fluvial
h=0 #index counter for hillslope
accumulate_1=0
accumulate_2=0

while time < total_model_time:

      z[mg.core_nodes]+= (uplift_rate*dt) #do uplift all the time
      bed[mg.core_nodes]+= (uplift_rate*dt)
      expweath.calc_soil_prod_rate()
      ddtd.run_one_step(dt)
      fr.run_one_step()
      space.run_one_step(dt)

      accumulate_1 += desired_slip_per_event_1
      accumulate_2 += desired_slip_per_event_2

      print('is accumulating')

      if accumulate_1>= mg.dx:
         ss_fault(grid=mg, fault_loc_y=fault_loc_y_1, total_slip=total_slip_1,
                  total_time=total_model_time, accumulate=accumulate_1)
         accumulate_1 = accumulate_1% mg.dx

         print('one slip fault 1')
      if accumulate_2 >= mg.dx:
          ss_fault(grid=mg, fault_loc_y=fault_loc_y_2, total_slip=total_slip_2,
                   total_time=total_model_time, accumulate=accumulate_2)
          accumulate_2 = accumulate_2 % mg.dx

          print('one slip fault 2')




      if time>0 and time%200 == 0:
          imshow_grid(mg, z, cmap='viridis', shrink=shrink)
          plt.title('Topography after ' + str(time) + ' years')
          plt.show()
          # print(mg.at_node.keys())
          # write_raster_netcdf(f'{model_name}{get_file_sequence(time)}.nc', mg, format="NETCDF4",
          #                     names=['surface_water__discharge', 'drainage_area',
          #                            'bedrock__elevation', 'soil__depth',
          #                            'sediment__flux', 'soil_production__rate', 'topographic__steepest_slope'])
          # loop_topo_img  = f'/Users/taranguiz/Research/Lateral_advection/output_model_run/{model_name}/{model_name}{get_file_sequence(time)}.png'
          # plt.savefig(
          #     loop_topo_img,
          #     dpi=300, facecolor='white'
          # )
          # add_file_to_writer(loop_topo_img)

          plt.clf()

      time = time + dt
      print(time)

print(time)
imshow_grid(mg,z, cmap='viridis', shrink=shrink)
plt.title('Topography after ' + str(total_model_time) + ' years')
# plt.show()
# write_raster_netcdf(
#     f'{model_name}{total_model_time}.nc',
#     mg,
#     format="NETCDF4",
#     names=[
#         'surface_water__discharge',
#         'drainage_area',
#         'bedrock__elevation',
#         'soil__depth',
#         'sediment__flux',
#         'soil_production__rate',
#         'topographic__steepest_slope'
#     ]
# )
# final_topo_img= f'/Users/taranguiz/Research/Lateral_advection/output_model_run/{model_name}/{model_name}{total_model_time}.png'
# plt.savefig(final_topo_img, dpi=300, facecolor='white')
# add_file_to_writer(final_topo_img)
# writer.close()
# mg.save(f'/Users/taranguiz/Research/Lateral_advection/output_model_run/{model_name}/{model_name}_{total_model_time}.asc')

# figsize = [16,4] # size of grid plots
# fig, ax = plt.subplots(figsize=figsize)

# x = mg.node_x[fault_nodes]
# # soil_level = bed + soil
#
# ax.plot(x, mg.at_node['soil__depth'][fault_nodes], 'orange', linewidth=2, markersize=12, label='soil')
# ax.plot(x, mg.at_node['bedrock__elevation'][fault_nodes], linewidth=2, markersize=12, label='bedrock')
# ax.plot(x, mg.at_node['topographic__elevation'][fault_nodes],'red', linewidth=2, markersize=12, label='topo')
#
#
# plt.title('Final cross-Profile topography at fault location')
# #plt.text(480, 9, 'air')
# #plt.text(480, 7, 'soil')
# #plt.text(470, 2, 'bedrock')
#
# # plt.xlim(0,1000)
# # plt.ylim(0, 10)
# ax.set_xlabel('X (m)')
# ax.set_ylabel('Depth (m)')
# ax.legend(loc='upper right')
# plt.show()








