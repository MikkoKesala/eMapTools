
# %%
from osgeo import gdal
import numpy as np
from collections import deque


def d_inf2(dem):
    """
    This calculates D-infinity based on Tarbeton (1997) method
    
    Inputs: dem as raster format
    Outputs: dinf : D-infinity (0..360 where zero points north)

    """
    dem_ds = gdal.Open(dem)
    dem_array = dem_ds.GetRasterBand(1).ReadAsArray()
    # Get the pixel resolution from the geotransform
    geotransform = dem_ds.GetGeoTransform()
    resolution = geotransform[1]  # Assuming square pixels, pixel size in meters
    diag_cell_size = np.sqrt(resolution * resolution + resolution * resolution)
    #adding extra cell for handling boundaries
    dem_array = np.pad(dem_array,pad_width=1, mode='edge')

    dem_ds = None
    
    #print (dem_array)
    e1_col = [-1,0,0,1,1,0,0,-1]
    e1_row = [0,1,1,0,0,-1,-1,0]

    e2_col = [-1,-1,1,1,1,1,-1,-1]
    e2_row = [1,1,1,1,-1,-1,-1,-1]

    ac_vals = {1:0,2:1,3:1,4:2,5:2,6:3,7:3,8:4}
    af_vals = {1:1,2:-1,3:1,4:-1,5:1,6:-1,7:1,8:-1}
    # Get the number of rows and columns in the DEM
    r_list=[]
    s_list=[]

    for i in range(1,9,1):
        e1 = np.roll(dem_array,(e1_row[i-1],e1_col[i-1]),axis=(0,1))
        e2 = np.roll(dem_array,(e2_row[i-1],e2_col[i-1]),axis=(0,1))
        
        ac = ac_vals[i]
        af = af_vals[i]
        
        r = np.where(dem_array>0,0,0)
        s = r
        
        #if e0 > e1 and e0 > e2 cases
        s1 = np.where(((dem_array>e1) & (dem_array>e2)),(dem_array-e1)/resolution,0)
        s2 = np.where(((dem_array>e1) & (dem_array>e2)),(e1-e2)/resolution,0)
        
        r = np.where(((dem_array>e1) & (dem_array>e2) & (s1!=0)),np.arctan(s2/s1),0)
        r = np.where(((dem_array>e1) & (dem_array>e2) & (s1==0)),np.pi,r)
        #s = np.sqrt(abs(s1*s1 + s2*s2))
        s = np.where(((dem_array>e1) & (dem_array>e2)),np.sqrt(s1*s1+s2*s2),0)
        #s = np.where(((s1<0) & (s2<0)),s*-1,s)
        #s = np.where(((s1<0) & (s2==0)),s*-1,s)
        #s = np.where(((s1==0) & (s2<0)),s*-1,s)
        
        s = np.where(r<0,s1,s)
        s = np.where(r>np.arctan(1),(dem_array-e2)/diag_cell_size,s)        
        r = np.where(r<0,0,r)
        r = np.where(r>np.arctan(1),np.arctan(1),r)
        
        #if only e0 > e1 cases
        r = np.where(((dem_array>e1) & (dem_array<=e2)),0,r)
        s = np.where(((dem_array>e1) & (dem_array<=e2)),(dem_array-e1) / resolution,s)

        #if only e0 > e2 cases
        r = np.where(((dem_array>e2) & (dem_array<=e1)),np.arctan(1),r)
        s = np.where(((dem_array>e2) & (dem_array<=e1)),(dem_array-e2) / diag_cell_size,s)

        r = af * r + ac * (np.pi/2)
        r_list.append(r)
        s_list.append(s)


    #direction based on max slope value
    slopes = np.stack(s_list,axis=-1)
    max_slope_index = np.argmax(slopes,axis=-1)
    directions = np.stack(r_list,axis=-1)
    rows, cols = max_slope_index.shape
    dinf = directions[np.arange(rows)[:, None], np.arange(cols), max_slope_index]

    #radians to degrees (0 radian = 90 degree so we rotate 90 cells value)
    #dinf = (180 / np.pi) * dinf
    #dinf = dinf+270
    dinf = 360-np.degrees(dinf)+90
    dinf = np.where(dinf>360,dinf-360,dinf)
    
    #pit cells to -1
    slope = np.max(slopes,axis=-1)
    dinf = np.where(slope>0,dinf,-1)
    dinf = dinf[1:-1, 1:-1]
    slope = slope[1:-1, 1:-1]
    
    #slope to degrees
    angles = [0,90,180,270,360]
    closest_direction = closest_number_2d_optimized(dinf,angles)[0]
    closest_direction = np.where(dinf==-1,-1,abs(closest_direction))
    horizontal_hy = resolution / np.cos(np.radians(abs(closest_direction-dinf)))
    slope_degrees = np.degrees(np.arctan(slope/horizontal_hy))

    #sidedistance = resolution /np.arcsin(90-(closest_direction-dinf))
    #sidedistance = np.where(dinf==-1,-1,sidedistance)
    #slope_degrees = np.tan(slope/sidedistance)


    #cleaning
    closest_direction = None
    horizontal_hy = None
    r_list = None
    s_list = None
    s1 = None
    s2 = None
    e1 = None
    e2 = None
    dem_array = None
    slopes = None
    max_slope_index = None
    directions = None
    slope = None

    return dinf,slope_degrees

def closest_number_2d_optimized(x_2d, y_list):
    """
    Calculation of closest number between 2d array and y-list
    Inputs:
        x_2d = 2d numpy array
        y_list = list type (angles 0,45,...,315,360)
    Outputs:
        closest_array = closest values as 2d numpy array format
        weight = weight coefficient based on closeness to closest direction
    """

    # Convert y_list to a NumPy array for vectorized operations
    y_array = np.array(y_list)
    
    # Calculate the absolute difference between each element in x_2d and every element in y_array
    # This results in a 3D array: (rows, cols, len(y_list))
    differences = np.abs(x_2d[:, :, np.newaxis] - y_array)
    weight = 1  - np.min(differences,axis=2)/45
    # Find the index of the minimum difference along the third axis (which corresponds to y_list)
    min_indices = np.argmin(differences, axis=2)
    
    # Use the indices to retrieve the closest number from y_list
    closest_array = y_array[min_indices]
    
    return closest_array,weight

def second_closest_number_2d_optimized(x_2d, y_list):
    """
    Calculation of second closest number between 2d array and y-list
    Inputs:
        x_2d = 2d numpy array
        y_list = list type (angles 0,45,...,315,360)
    Outputs:
        closest_array = closest values as 2d numpy array format
        weight = weight coefficient based on closeness to second closest direction
    """
    # Convert y_list to a NumPy array for vectorized operations
    y_array = np.array(y_list)
    
    # Calculate the absolute difference between each element in x_2d and every element in y_array
    differences = np.abs(x_2d[:, :, np.newaxis] - y_array)
    
    # Find the index of the closest number (smallest difference)
    closest_indices = np.argmin(differences, axis=2)
    
    # Set the closest number's difference to infinity to exclude it from being selected again
    rows, cols = x_2d.shape
    differences[np.arange(rows)[:, None], np.arange(cols), closest_indices] = np.inf
    
    weight = 1  - np.min(differences,axis=2)/45
    # Now find the index of the second closest number
    second_closest_indices = np.argmin(differences, axis=2)
    
    # Use the indices to retrieve the second closest number from y_list
    second_closest_array = y_array[second_closest_indices]
    
    return second_closest_array,weight


def calculate_flow_accumulation_with_toposort(dem, flow_dir):
    # Initialize accumulation array
    flow_accum = np.ones(dem.shape, dtype=np.float32)
    
    rows, cols = dem.shape
    inflow_count = np.zeros(dem.shape, dtype=int)
    outflow_targets = np.zeros(dem.shape + (2,), dtype=int) - 1  # To store the row, col of outflow

    # Define the direction vectors for 8 primary directions
    drow = [-1, -1, 0, 1, 1, 1, 0, -1, -1]
    dcol = [0, 1, 1, 1, 0, -1, -1, -1, 0]
    angles = [0, 45, 90, 135, 180, 225, 270, 315,360]

    # Determine outflow targets based on flow direction
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            angle = flow_dir[i, j]
            if angle == -1:  # Skip cells without defined direction
                continue
            s_angles = angles
            
            # Find the closest direction based on angle
            closest_direction = min(range(8), key=lambda k: abs(angle - angles[k]))
            s_angles[closest_direction] = 1000
            second_closest = min(range(8), key=lambda k: abs(angle - s_angles[k]))
            out_r, out_c = i + drow[closest_direction], j + dcol[closest_direction]
            outflow_targets[i, j] = (out_r, out_c)

            # Increase inflow count for the outflow target
            inflow_count[out_r, out_c] += 1

    # Initialize a queue for cells with no inflows
    queue = deque([(i, j) for i in range(rows) for j in range(cols) if inflow_count[i, j] == 0])

    # Process cells in topological order
    while queue:
        i, j = queue.popleft()
        
        # Accumulate flow to the outflow target
        out_r, out_c = outflow_targets[i, j]
        if out_r != -1 and out_c != -1:
            flow_accum[out_r, out_c] += flow_accum[i, j]
            
            # Decrement inflow count and add to queue if no more inflows
            inflow_count[out_r, out_c] -= 1
            if inflow_count[out_r, out_c] == 0:
                queue.append((out_r, out_c))

    return flow_accum

def calculate_flow_accumulation_dinf(dem):
    
    dinf,slope = d_inf2(dem)
    dem = gdal.Open(dem).ReadAsArray()
    directions = [i for i in range(0,361,45)]
    (dir1,proportion1) = closest_number_2d_optimized(dinf,directions)
    (dir2,proportion2) = second_closest_number_2d_optimized(dinf,directions)

    rows, cols = dem.shape
    flow_accum = np.ones((rows, cols), dtype=np.float32)
    inflow_count = np.zeros((rows, cols), dtype=int)

    dirs = {0:(-1,0),45:(-1,1),90:(0,1),135:(1,1),
            180:(1,0),225:(1,-1),270:(0,-1),315:(-1,-1),360:(-1,0)}
    # Initialize inflow count based on primary and secondary directions
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            pr, pc = dirs[dir1[i, j]]
            pr,pc = i+pr,j+pc
            sr, sc = dirs[dir2[i, j]]
            sr,sc = i+sr,j+sc
            if pr != -1:
                inflow_count[pr, pc] += 1
            if sr != -1:
                inflow_count[sr, sc] += 1

    # Queue of cells with no inflows
    
    queue = deque([(i, j) for i in range(rows) for j in range(cols) if inflow_count[i, j] == 0])

    # Process cells in topological order
    while queue:
        i, j = queue.popleft()

        # Primary and secondary directions
        pr, pc = dirs[dir1[i, j]]
        sr, sc = dirs[dir2[i, j]]
        
        pr,pc =pr+i,pc+j
        sr,sc=sr+i,sc+j
        # Accumulate flow to primary and secondary directions
        if pr != -1:
            flow_accum[pr, pc] += flow_accum[i, j] * proportion1[i, j]
            inflow_count[pr, pc] -= 1
            if inflow_count[pr, pc] == 0:
                queue.append((pr, pc))

        if sr != -1 and sr<rows and sc<cols:
            flow_accum[sr, sc] += flow_accum[i, j] * proportion2[i, j]
            inflow_count[sr, sc] -= 1
            if inflow_count[sr, sc] == 0:
                queue.append((sr, sc))

    return flow_accum,dinf,slope

def topvalue(array):
    array = np.pad(array,pad_width=1, mode='edge')
    
    dir = {1 : (0,-1),2:(1,-1),3:(1,0),4:(1,1),
           5:(0,1),6:(-1,1),7:(-1,0),8:(-1,-1)}
    
    array_max = array
    for i in range(1,9,1):
        array_max = np.maximum(array_max,np.roll(array,dir[i],(0,1)))
    
    top_array = np.where(array-array_max==0,1,0)
    top_array = top_array[1:-1, 1:-1]
    return top_array

def calculate_massflux_dinf(dem,loading,efficiency,absorption):
    
    dinf = d_inf2(dem)[0]
    dem_ds = gdal.Open(dem)
    dem = dem_ds.GetRasterBand(1).ReadAsArray()
    directions = [i for i in range(0,361,45)]
    (dir1,proportion1) = closest_number_2d_optimized(dinf,directions)
    (dir2,proportion2) = second_closest_number_2d_optimized(dinf,directions)

    rows, cols = dem.shape
    #flow_accum = np.ones((rows, cols), dtype=np.float32)
    loading = gdal.Open(loading).ReadAsArray()
    efficiency  = gdal.Open(efficiency).ReadAsArray()
    if absorption is None:
        absorption =  np.zeros((rows,cols),dtype=float)
    else:
        absorption = gdal.open(absorption).ReadAsArray()

    massflux = loading-absorption
    inflow_count = np.zeros((rows, cols), dtype=int)

    #drow = [-1, -1, 0, 1, 1, 1, 0, -1, -1]
    #dcol = [0, 1, 1, 1, 0, -1, -1, -1, 0]
    #angles = [0, 45, 90, 135, 180, 225, 270, 315,360]
    dirs = {0:(-1,0),45:(-1,1),90:(0,1),135:(1,1),
            180:(1,0),225:(1,-1),270:(0,-1),315:(-1,-1),360:(-1,0)}
    # Initialize inflow count based on primary and secondary directions
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            pr, pc = dirs[dir1[i, j]]
            pr,pc = i+pr,j+pc
            sr, sc = dirs[dir2[i, j]]
            sr,sc = i+sr,j+sc
            if pr != -1:
                inflow_count[pr, pc] += 1
            if sr != -1:
                inflow_count[sr, sc] += 1

    # Queue of cells with no inflows
    
    queue = deque([(i, j) for i in range(rows) for j in range(cols) if inflow_count[i, j] == 0])

    # Process cells in topological order
    while queue:
        i, j = queue.popleft()

        # Primary and secondary directions
        pr, pc = dirs[dir1[i, j]]
        sr, sc = dirs[dir2[i, j]]
        
        pr,pc =pr+i,pc+j
        sr,sc=sr+i,sc+j

        # Accumulate massflux to primary and secondary directions
        if pr != -1:
            massflux[pr, pc] += massflux[i, j] * proportion1[i, j] * efficiency[i,j]
            inflow_count[pr, pc] -= 1
            if inflow_count[pr, pc] == 0:
                queue.append((pr, pc))

        if sr != -1 and sr<rows and sc<cols:
            massflux[sr, sc] += massflux[i, j] * proportion2[i, j] * efficiency[i,j]
            inflow_count[sr, sc] -= 1
            if inflow_count[sr, sc] == 0:
                queue.append((sr, sc))

    return massflux

def ls_factor(dem):
        #get rasters
        #dem = QgsProcessingUtils.mapLayerFromString(parameters['dem'],context)
        #demname = dem.source()

        #lsout = self.parameterAsOutputLayer(parameters,"ls",context)

        fa,aspect,slope = calculate_flow_accumulation_dinf(dem)
        #print (slope)
        #ls_dinf = calculate_flow_accumulation_dinf(demname)
        ds = gdal.Open(dem)
        psize = ds.GetGeoTransform()[1]
        aspect = np.radians(aspect)
        slope_r = np.radians(slope)
        beta = (np.sin(slope_r) / 0.0896) / (3 * (np.sin(slope_r) ** 0.8) + 0.56)
        m = beta / (1 + beta)
        #S = np.where(slope<9,10.8*np.sin(slope_r)+0.03,16.8*np.sin(slope_r)-0.5)
         
          
        S = (65.41 * (np.sin(slope_r) ** 2) +
             4.56 * np.sin(slope_r) + 
             0.065)
        #L = (fa+ pize**2)**m+1
        #x = np.sin(aspect) + np.cos(aspect)
        #L = (fa*2 / 22.13) ** m
        L = ((fa +psize**2)**(m+1) - fa**(m+1)) / ((22.13)**m*psize**m)
        ls  = L*S

        return ls


# %%
dem = r'Z:/Documents/GitHub/eMapTools _dev/processing/algorithms/testfiles_dev/testDEM.tif'
ls = ls_factor(dem)
#dinf,slope = d_inf2(dem)
#print (dinf)
print (ls)

# %%
import numpy as np
slope = np.array([0,2,5,7,10,13,15,20])
slope = np.radians(slope)
fa = 10000
psize = 2
aspect = 30
aspect = np.radians(aspect)
beta = np.sin(slope) / (0.0896 * (3 * (np.sin(slope) ** 0.8) + 0.56))
#beta = (np.sin(slope) / 0.0896) / (3 * (np.sin(slope) ** 0.8) + 0.56)
S = np.where(slope<9,10.8*np.sin(slope)+0.03,16.8*np.sin(slope)-0.5)
m = beta / (1 + beta)
x = np.sin(aspect) + np.cos(aspect)
#L = (slope_length / 22.13) ** m
L = ((fa +psize**2)**(m+1) - fa**(m+1)) / (psize**(m+2)*x*22.13**m)
print (beta)
print (m)
print (S)
print (L)