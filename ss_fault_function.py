import numpy as np
import matplotlib.pyplot as plt

from landlab import RasterModelGrid, imshow_grid, imshowhs_grid # landlab grid components
from landlab.io import read_esri_ascii

def ss_fault(grid, fault_loc_y, total_slip, total_time, accumulate):

    #grid: raster model grid
    #fault_loc_y: y value (row) where fault is located
    #total_slip: total slip in meters
    #method: drop or roll
    #accumulate: input that is taking from a loop, is a int and is how many columns are going to be moved.
    #faulted surface faulted_fields=[['topographic__elevation'], ['bedrock__elevation'], ['soil__depth']]

    nrows= grid.number_of_node_rows
    ncols= grid.number_of_node_columns

    z_original= grid.at_node["topographic__elevation"]
    z_original_reshaped = np.reshape(z_original, (nrows, ncols))
    bed_original = grid.at_node["bedrock__elevation"]
    bed_original_reshaped = np.reshape(bed_original, (nrows, ncols))
    soil_original = grid.at_node["soil__depth"]
    soil_original_reshaped = np.reshape(soil_original, (nrows, ncols))


    # INPUT FROM USER ABOUT TECTONICS
    slip_rate = (total_slip / total_time)  # in m/yr
    print(' the slip rate of your fault is ' + str(slip_rate * 1000) + ' in mm/yr')
    number_cols = int(accumulate/grid.dx)   # in m assuming characteristic behavior #number of columns
    print('number of columns dropped per event is: ' + str(number_cols))




    z_top_new = np.roll(z_original_reshaped[fault_loc_y:, 1:-1], number_cols, axis=1)
    z_original_reshaped[fault_loc_y:, 1:-1] = z_top_new
    z_reshaped_after_shift = np.reshape(z_original_reshaped, (nrows * ncols))
    z_original[:] = z_reshaped_after_shift

    soil_top_new = np.roll(soil_original_reshaped[fault_loc_y:, 1:-1], number_cols, axis=1)
    soil_original_reshaped[fault_loc_y:, 1:-1] = soil_top_new
    soil_reshaped_after_shift = np.reshape(soil_original_reshaped, (nrows * ncols))
    soil_original[:] = soil_reshaped_after_shift

    bed_top_new = np.roll(bed_original_reshaped[fault_loc_y:, 1:-1], number_cols, axis=1)
    bed_original_reshaped[fault_loc_y:, 1:-1] = bed_top_new
    bed_reshaped_after_shift = np.reshape(bed_original_reshaped, (nrows * ncols))
    bed_original[:] = bed_reshaped_after_shift

        # imshow_grid(grid, z_original, cmap='viridis')
        # plt.show()

    #displacement += slip_per_event

    # imshow_grid(grid, z_original, cmap='viridis')
    # plt.show()


# This is just an example


# #READING ORIGINAL TOPO
# (grid,z_original)=read_esri_ascii(
#     '/Users/taranguiz/Research/NZ_Tararua_SSF/steady_state_1/finaltopo_topographic__elevation.asc',
#     name="topographic__elevation")
# grid.set_closed_boundaries_at_grid_edges(
#     bottom_is_closed=False,
#     left_is_closed=True,
#     right_is_closed=True,
#     top_is_closed=False)
#
# imshow_grid(grid,z_original, cmap='viridis')
# plt.show()
# #
# grid.add_zeros("node", "soil__depth") #add field to the grid
# grid.at_node["soil__depth"]=grid.at_node["soil__depth"]+1 #2
# grid.at_node["bedrock__elevation"]=grid.at_node["topographic__elevation"] - grid.at_node["soil__depth"]
#
# soil=grid.at_node['soil__depth']
# bed= grid.at_node['bedrock__elevation']
# #
# imshow_grid(grid, soil)
# plt.show()
# displacement = 0
# time_track =0
# model_time=100000
# dt=100
# total_slip_1=1000
# desired_slip_per_event_1=(total_slip_1/model_time)*dt
# print(desired_slip_per_event_1)
# total_slip_2=200
# desired_slip_per_event_2=(total_slip_2/model_time)*dt
#  # effective_slip= grid.dx/desired_slip_per_event
# accumulate_1=0
# accumulate_2=0
# #     # np.copy(desired_slip_per_event)
# #
# while time_track < model_time:
#
#     accumulate_1+= desired_slip_per_event_1
#     accumulate_2+= desired_slip_per_event_2
#     time_track += dt
#     print('is accumulating')
#
#     if accumulate_1 >= grid.dx:
#         ss_fault(grid=grid, fault_loc_y=int(grid.number_of_node_rows/4.),
#                  total_slip=total_slip_1, total_time=model_time, accumulate=accumulate_1)
#         accumulate_1= accumulate_1%grid.dx
#         print('one slip fault 1')
#         imshow_grid(grid, z_original, cmap='viridis', shrink=0.5)
#         plt.title('testing')
#         plt.show()
#         plt.close()
#
#     if accumulate_2 >= grid.dx:
#         ss_fault(grid=grid, fault_loc_y=int(grid.number_of_node_rows /1.5),
#                  total_slip=total_slip_2, total_time=model_time, accumulate=accumulate_2)
#         accumulate_2 = accumulate_2 % grid.dx
#         print('one slip fault 2')
#         imshow_grid(grid, z_original, cmap='viridis', shrink=0.5)
#         plt.title('testing')
#         plt.show()
#         plt.close()
