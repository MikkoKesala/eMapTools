# %%
import geopandas as gpd
import numpy as np
from paras_2 import diameterFromHeight,decay_tree_potential
from geotools2 import pointDelineate
from ecoindices import diversityIndices
from geostats import getisord,spatialWeights
from scipy.spatial import cKDTree
import glob
from osgeo import gdal
import pandas as pd
import os




# %%
#read files to single geodataframe
tree_path = r'C:\Users\mjkesala\pvenv\aanekoski\Tree_circle_shp'
site_path = r'C:\Users\mjkesala\pvenv\aanekoski\Site_boundary_shp'
def combineFiles(folder_path):
    gdfs = []
    for shp_file in glob.glob(os.path.join(folder_path, "*.shp")):
        gdf = gpd.read_file(shp_file)
        
        gdf["site"] = os.path.splitext(os.path.basename(shp_file))[0]
        
        gdfs.append(gdf)

        #print(f"Processed {shp_file}, with {len(gdf)} geometries.")

    # Concatenate all GeoDataFrames into a single GeoDataFrame
    combined_gdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True))
    combined_gdf
    return combined_gdf
# %%
combined_gdf = combineFiles(tree_path)
sites = combineFiles(site_path)


combined_gdf['site'] = combined_gdf['site'].str[:-7]
combined_gdf

sites.set_crs("EPSG:4326", allow_override=True)
sites_utm = sites.to_crs("EPSG:3067")
sites_utm['area_ha'] = np.round(sites_utm.geometry.area / 10000,3)
#sites_utm
combined_gdf = combined_gdf.merge(sites_utm[['site','area_ha']],on='site',how='left')

combined_gdf

# %%
#change geometries

combined_point = combined_gdf.set_crs("EPSG:4326", allow_override=True)
combined_point = combined_point.to_crs("EPSG:3067")
combined_point.geometry = combined_point.geometry.centroid

combined_point



# %%

#unique_trees = combined_point["species"].unique()
#print (unique_trees)

treespecies = {'pine':1,'spruce':2,'birch':3,'aspen':5}
ts4dph = {'pine':1,'spruce':2,'birch':3,'aspen':3}

combined_point['ts_dbh'] = combined_point['species'].apply(lambda x: ts4dph[x])
combined_point['speciesclass'] = combined_point['species'].apply(lambda x: treespecies[x])
combined_point['dbh'] = combined_point.apply(lambda x: diameterFromHeight(x['ts_dbh'],x['chm'],3),axis=1)
combined_point = combined_point.drop(['ts_dbh'],axis=1)

combined_point

# %%
def calcDWP(zone,dbh,fertilityclass,treespecies):

    #calculate ecological indices
    zone = 'zone'+str(zone)
    decay_params = decay_tree_potential('zone3')
    dwp = np.poly1d(decay_params[fertilityclass][treespecies])(dbh)
    if dwp > 2:
        dwp == 2.0

    return dwp

# Function to convert geographic coordinates to pixel coordinates
def geo_to_pixel(geo_x, geo_y,geotransform):
    pixel_x = int((geo_x - geotransform[0]) / geotransform[1])
    pixel_y = int((geo_y - geotransform[3]) / geotransform[5])
    return pixel_x, pixel_y

def sampleValues(gdf,raster,fieldname):
    raster_ds = raster

    # Get the raster band (assuming we want the first band)
    band = raster_ds.GetRasterBand(1)

    # Get the raster geotransform to convert coordinates
    geotransform = raster_ds.GetGeoTransform()
    # Sample values at each point
    sampled_values = []
    for point in gdf.geometry:
        pixel_x, pixel_y = geo_to_pixel(point.x, point.y,geotransform)
    
        # Read the pixel value
        if 0 <= pixel_x < band.XSize and 0 <= pixel_y < band.YSize:
            value = band.ReadAsArray(pixel_x, pixel_y, 1, 1)[0, 0]
            sampled_values.append(value)
        else:
            sampled_values.append(np.nan)  # Handle out-of-bounds points

    # Add sampled values to the GeoDataFrame
    gdf[fieldname] = sampled_values

    # Close the raster dataset
    raster_ds = None

    return gdf

def replaceNull(gdf,target_column,nullvalue):
    # Column with null values (replace 'your_column' with your actual column name)

    # Separate the data into GeoDataFrames with and without null values
    gdf_with_values = gdf[gdf[target_column]!=nullvalue]
    gdf_with_nulls = gdf[gdf[target_column]==nullvalue]

    # Create KDTree from the coordinates of points with values
    tree = cKDTree(gdf_with_values.geometry.apply(lambda geom: (geom.x, geom.y)).tolist())

    # Find nearest neighbor for each null value
    distances, indices = tree.query(gdf_with_nulls.geometry.apply(lambda geom: (geom.x, geom.y)).tolist())

    # Replace nulls with the nearest value
    gdf.loc[gdf[target_column]==nullvalue, target_column] = gdf_with_values.iloc[indices][target_column].values

    return gdf


# %%
# Enrich data raster data

dtw_url = "WCS:https://paituli.csc.fi/geoserver/paituli/wcs?service=WCS&version=1.0.0&request=GetCoverage&coverageId=paituli:luke_dtw_04_2023"
sfclass = "WCS:https://paituli.csc.fi/geoserver/paituli/wcs?service=WCS&version=1.0.0&request=GetCoverage&coverageId=paituli:luke_vmi_kasvupaikka_2019"
dtw = gdal.Open(dtw_url)
sfclass = gdal.Open(sfclass)

combined_point = sampleValues(combined_point,dtw,"dtw")
combined_point = sampleValues(combined_point,sfclass,"sfclass")
combined_point

# %%
combined_point["dtw"] = combined_point["dtw"] / 100
combined_point

combined_point['dtw2'] = np.max(combined_point['dtw']) - combined_point['dtw']
combined_point = replaceNull(combined_point,"sfclass",32767)
print (np.max(combined_point['dtw']),np.min(combined_point['dtw']))
print (combined_point['sfclass'].unique())
# %%
#print (decay_params[3][5])
#value = np.poly1d(decay_params[3][5])(30)
#print (value)
combined_point["dwp"] = combined_point.apply(lambda x:calcDWP(3,x.dbh,x.sfclass,x.speciesclass),axis=1)
combined_point = diversityIndices(combined_point,"speciesclass",20,["shannon"])
combined_point
# %%
combined_point['aspen_prop'] = combined_point['list_speciesclass'].apply(lambda x: len([s for s in x if s == 5]) / len(x))
combined_point['deci_prop'] = combined_point['list_speciesclass'].apply(lambda x: len([s for s in x if s in [3,5]]) / len(x))
combined_point
# %%
#calculating hotspot values
collect = []
indices = ['dtw2','dwp','shannon','aspen_prop','deci_prop']
groupvalues = sorted(combined_point['site'].unique())
for g in groupvalues:
    # selecting rows based on condition 
    #uniques = combined_point.loc[combined_point['site']==g]
    uniques = combined_point[(combined_point['site']==g)]
    for i in indices:
        uniques = getisord(uniques,i,30,"gaussian")
    collect.append(uniques)

combined_point = gpd.GeoDataFrame(pd.concat(collect, ignore_index=True))
combined_point
#combined_point = getisord(combined_point,)

# %%

#combined_point["sum_gi*"] = combined_point.apply(lambda x: np.sum(x['dtw2_gi*'],x['dwp_gi*'],x['shannon_gi*']),axis=1)
combined_point["sum_GiZ"] =combined_point[['dwp_GiZ','shannon_GiZ','aspen_prop_GiZ','deci_prop_GiZ']].sum(axis=1)
combined_point["rtreecount"] = np.ceil(combined_point["area_ha"] * 20).astype("int")
combined_point

#test['test'] = combined_point['dtw2_GiZ'].rank(method='first',ascending=False)
#print (test[['dtw2_GiZ','test']].sort_values(by=["dtw2_GiZ"]))
# %%
#indices = ['dtw2_GiZ','dwp_GiZ','shannon_GiZ','deci_prop_GiZ','sum_GiZ']
# remove small trees
min_dbh = 10
trees_ha = 20
max_tree = 20
min_tree = 5
min_distance = 50
max_size = np.max(combined_point['area_ha'])
#max_iter = max_size*trees_ha / min_tree
#print (max_size,max_iter)
# %%
selection_df = combined_point[combined_point['dbh']>=min_dbh]
selection_df[['site','dtw2_GiZ','dwp_GiZ','shannon_GiZ','aspen_prop_GiZ']]
indices = ['dtw2_GiZ','dwp_GiZ','shannon_GiZ','aspen_prop_GiZ']
test = selection_df.groupby('site')[indices[0]].quantile(0.95)
print (test)
#for i in indices:
# df['90th_percentile'] = df.groupby('group')['value'].transform(lambda x: x.quantile(0.9))
#    selection_df[i+'_p90'] = selection_df.groupby('site')[i].quantile(0.95)
#selection_df
# %%
for i in range(1,10):

    selection_df['max_GiZ'] = selection_df[['dtw2_GiZ','dwp_GiZ','shannon_GiZ','deci_prop_GiZ']].max(axis=1)
    selection_df['max_column'] = selection_df[['dtw2_GiZ','dwp_GiZ','shannon_GiZ','deci_prop_GiZ']].idxmax(axis=1)
    #selection_df['rank']
    #selection_df['selected'] = combined_point.groupby('site')['max_GiZ'].rank(method="first",ascending=False)<=combined_point["rtreecount"]).astype("int")


# %%
agg = {'site':'first',
       'deci_prop':'mean',
       'aspen_prop':'mean',
       'shannon':'mean',
       'dtw':'mean',
       'dbh':'mean',
       'density':'mean',
       'gini':'mean',
       'species':list}
selection_df = combined_point
selection_df['selected'] = (combined_point.groupby('site')['sum_GiZ'].rank(method="first",ascending=False)<=combined_point["rtreecount"]).astype("int")
selection_df = selection_df[selection_df['selected']==1]

groups = sorted(combined_point['site'].unique())
retentiotrees = []
for g in groups:
    group_gdf = selection_df[selection_df['site'] == g]
    ret = pointDelineate(group_gdf,5,agg)
    retentiotrees.append(ret)

retention = gpd.GeoDataFrame(pd.concat(retentiotrees, ignore_index=True))
retention
# %%
#back to wgs84
co = combined_point.to_crs("EPSG:4326")
re = retention.to_crs("EPSG:4326")

# %%
output = r'C:\Users\mjkesala\pvenv\aanekoski\projectdata.gpkg'
name = 'treemap'
rname = 'retentiongroups'
co.to_file(output,layer=name,driver="GPKG")
re.to_file(output,layer=rname,driver="GPKG")

# %%
for i in indices:
    combined_point[i+"_selected"] = (combined_point.groupby('site')[i].rank(method="first",ascending=False)<=combined_point["rtreecount"]).astype("int")

print (combined_point[['dtw2_GiZ','dtw2_GiZ_selected']].sort_values(by='dtw2_GiZ',ascending=False))
#   df['top_n'] = (df['values'].rank(method='first', ascending=False) <= n).astype(int)
#combined_point["treecount"] = combined_point.groupby("site").apply(lambda x: np.ceil(x['area']*10)

# %%
from shapely.ops import unary_union
# %%
#delieanation of retention trees
d = 5
sample = combined_point[(combined_point['sum_GiZ_selected']==1) & (combined_point['site']=='331590520_259')]
print (len(sample))
sample['buffer'] = sample.buffer(d)
sample.set_geometry('buffer')
sample['buffer'].set_crs('EPSG:3067')


sample = sample.dissolve(by=None,aggfunc='first')
#sample = sample.explode(index_parts=True)
print (len(sample))
#area_gdf['buffer'] = area_gdf.buffer(2*d)
#area_gdf.set_geometry('buffer')
#area_gdf['buffer'].set_crs('EPSG:3067')
#area_gdf = area_gdf.dissolve(by="MITKOEALAID",aggfunc='first')
#area_gdf['buffer2'] = area_gdf.buffer((-2*d)-(d/2))

sample.plot()
# %%
import matplotlib.pyplot as plt
import seaborn as sns
# %%
# Create a 3x3 grid for each group
fig, axes = plt.subplots(19, 5, figsize=(20, 80))
indices = ['dwp_GiZ','shannon_GiZ','deci_prop_GiZ','aspen_prop_GiZ','sum_GiZ']
# Define unique groups and values
groups = sorted(combined_point['site'].unique())  # There should be 9 groups
#values = sorted(combined_point['value'].unique())  # There should be 3 values
labels = ['Dead Wood Potential','Tree Diversity','Deciduous Proportion','Aspen Proportion','Sum of indices']
# Define a color palette for different values
#color = sns.color_palette("viridis")
idx=0
# Loop through groups and plot sites within each group
for group in (groups):
    #row, col = divmod(idx, 3)  # Calculate grid position

    # Filter for the current group
    group_gdf = combined_point[combined_point['site'] == group]
    for c,i in enumerate(indices):
        row,col = divmod(idx,5)
        ax=axes[row,col]
        # Filter by site and value for the current plotcol
        #value_gdf = group_gdf[[i]]
        # Plot each point in the filtered GeoDataFrame
        #selection = group_gdf[group_gdf[i+'_selected']==1]
        #notselection = group_gdf[group_gdf[i+'_selected']!=1]
        
        #selection.plot(ax=ax, column=i,legend=False,color='black', markersize=12)
        #notselection.plot(ax=ax, column=i,legend=False,cmap='coolwarm',vmin=-6,vmax=6, marker='x',markersize=3)
        #selection.plot(ax=ax, column=i,legend=False,color='black',marker='o',markersize=8)
        group_gdf.plot(ax=ax, column=i,legend=False,cmap='coolwarm',vmin=-10,vmax=10, markersize=5)
        if col == 0:
            ax.set_ylabel(group)
        if row == 0 or row%3==0:
            ax.set_title(labels[c])
        idx+=1
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        # Set title and legend
        #ax.set_title(f'Group {group}')
        #ax.legend(loc='upper right')

# Adjust layout and show plot
#plt.tight_layout()

fig.suptitle("Hotspot of ecological indices within 30 m radius",loc='upper middle')

plt.show()

# %%

fig, axes = plt.subplots(19, 5, figsize=(20, 80))
indices = ['dtw2_GiZ','dwp_GiZ','shannon_GiZ','deci_prop_GiZ','sum_GiZ']
# Define unique groups and values
groups = sorted(combined_point['site'].unique())  # There should be 9 groups
#values = sorted(combined_point['value'].unique())  # There should be 3 values

# Define a color palette for different values
#color = sns.color_palette("viridis")
idx=0
# Loop through groups and plot sites within each group
for group in (groups):
    #row, col = divmod(idx, 3)  # Calculate grid position

    # Filter for the current group
    group_gdf = combined_point[combined_point['site'] == group]
    for i in indices:
        row,col = divmod(idx,5)
        ax=axes[row,col]
        # Filter by site and value for the current plotcol
        #value_gdf = group_gdf[[i]]
        # Plot each point in the filtered GeoDataFrame
        
        selection = group_gdf[group_gdf[i+'_selected']==1]

        #notselection = group_gdf[group_gdf[i+'_selected']!=1][[i,'dbh']]
        if len(selection)>0:
            sns.histplot(selection,x="dbh",hue='species',multiple="stack",ax=ax)
        #plt.hist(ax=ax,x=selection['dbh'])
        #selection.plot(ax=ax, column=i,legend=False,color='black', markersize=12)
        #sns.lmplot(data=group_gdf,x=i,ax=ax,legend=False,color='blue')
        #selection.plot(ax=ax,x=i,y='dbh',legend=False,color='red',marker='o')
        #notselection.plot(ax=ax,x=i,y='dbh',legend=False,color='blue',marker='o')
        #group_gdf.plot(ax=ax, column=i,legend=False,cmap='coolwarm',vmin=-6,vmax=6, markersize=10)
        if col == 0:
            ax.set_ylabel(group)
        if row == 0:
            ax.set_title(i)
        idx+=1



# %%


distances = [i for i in range(0,31,1)]
#print (distances)

linear = spatialWeights(distances,'linear',30)
guadric = spatialWeights(distances,'quadric',30)
squared = spatialWeights(distances,'squared',30)
gaussian = spatialWeights(distances,'gaussian2',30)
fixed = spatialWeights(distances,'fixed',30)
#print (len(distances),len(linear),len(guadric),len(gaussian),len(fixed))

disttable = {'distances':distances,
             'linear':linear,
             'quadric':guadric,
             'squared':squared,
             'gaussian':gaussian,
             'fixed':fixed}
disttable = pd.DataFrame(disttable)
print(disttable)
#print (weights)
#save output

# %%
plt.figure(figsize=(10,6))

plt.plot(disttable['distances'].values,disttable['linear'].values,label='Linear')
plt.plot(disttable['distances'].values,disttable['quadric'].values,label='Quadric')
plt.plot(disttable['distances'].values,disttable['squared'].values,label='Squared')
plt.plot(disttable['distances'].values,disttable['gaussian'].values,label='Gaussian')
plt.plot(disttable['distances'].values,disttable['fixed'].values,label='Fixed')

plt.xlabel('Distance (m)')
plt.ylabel("Weight value (w)")
plt.title("Distance decay functions of spatial hotspot analysis with 30 m threshold")

plt.legend()
plt.legend(title="Methods")
#disttable.plot(x='distances',y='quadric')

# %%
output = r'C:\Users\mjkesala\pvenv\aanekoski\projectdata.gpkg'
name = 'aanekoski_treemap2'
combined_point.to_file(output,layer=name,driver="GPKG")
# %%
