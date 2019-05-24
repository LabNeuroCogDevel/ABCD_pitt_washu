The scripts and functions in this folder resort matrices based on the Gordon parcellation. Below is a description of each script/function

1. Parcels_Network_Ids.txt
This file contains the ROIs ordered 1-333 and includes coordinates and network affiliations. This path needs to be given as the filename (line 9 in 'read_gordon_coordinates.m' & line 8 in Read_Gordon_Parcel_IDs.m)

2. gordon_matrix_plotting.m
This script will generate the ordered matrices and add network borders. If you pass your matrix as the first argument into the function 'plot_adj_matrix.m' and add the correct path to 'Parcels_Network_IDs.txt' in each of the other 2 scripts, you should get the desired output

3. reorder_gordon_laumann_parcels 
Script to order ROIs by network. Network order is defined by Net_labels. The first column of 'NetworksOrdered' is the proper ordering of ROIs (within networks they are sorted left->right by x-coordinate.

4. Read_Gordon_Parcel_IDs.m
Contained within 'reorder_gordon_laumann_parcels.m'. Just change the path to 'Parcels_Network_IDs.txt' and you should be good to go. This script reads in the parcel text file

5. read_gordon_coordinate.m
Reads the coordinates (centroids) for the parcels to sort them in the x direction within networks. Just change the filename path to 'Parcels_Network_IDs.txt'

6. plot_adj_matrix.m
function to plot matrix according to the roi_order and draw lines at network partitions defined as 'partitionidx'. Should just need to change the first input argument to whatever matrix you are looking to plot. 
