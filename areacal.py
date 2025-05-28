# -*- coding: utf-8 -*-
# %%
"""
Created on Tue Dec  5 11:28:14 2023

@author: ROG
"""

# %%
import os

def path_convert(inshp,tif_list):
    series = os.path.split(inshp)[1][12:][0:-6]
    intif = os.path.splitext(tif_list[0])[0][0:31]+series+'.tif'
    return intif
    #os.path.exists(intif)
    
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from rasterio.enums import Resampling

def shptoraster(existing_raster_path,shapefile_path,output_raster_path, invert = True, geom = 0):
    # Read the existing raster to use as a template
    with rasterio.open(existing_raster_path) as src:
        # Get raster properties
        crs = src.crs
        transform = src.transform
        width = src.width
        height = src.height
        dtype = src.dtypes[0]  # Assuming a single band raster

    # Read the shapefile
    if geom == 0:
        gdf = gpd.read_file(shapefile_path)
    else:
        gdf = shapefile_path
    if len(gdf) !=0:
        gdf = gdf.to_crs(crs)
        # Create an empty raster with the same properties as the existing raster
        with rasterio.open(
            output_raster_path,
            'w',
            driver='GTiff',
            height=height,
            width=width,
            count=1,  # Number of bands
            dtype=dtype,
            crs=crs,
            transform=transform,
            compress = 'LZW'
        ) as dst:
            # Rasterize the geometry onto the raster
            mask = rasterio.features.geometry_mask(
                gdf.geometry,
                transform=transform,
                invert= invert,
                #all_touched = True,
                out_shape=(height, width)
            )

            # Write the rasterized data to the output raster
            dst.write(mask.astype(dtype), 1)
    else:
        if geom ==0:
            print(shapefile_path)
        else:
            print(existing_raster_path)
        print('no data')

# %%
def grwl(gdbin,dirin,tifpath,txt, mod = 1, extt = '.tif'):
    '''
    gdbin: geopandas read_file
    dirin: dir stored grwl_vector
    txt: txt file to store no exist grwl vector
    mod: control return concated gdb or not, default 1 is return, otherwise, return path 
    
    import numpy as np
    import os
    import geopandas as gpd
    '''
    ext = gdbin.geometry.total_bounds
    zimu="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    zm1 = zimu[int(abs(ext[1])//4)]
    zm2 = zimu[int(abs(ext[3])//4)]
    num1 = int((ext[0]+180)//6)+1
    num2 = int((ext[2]+180)//6)+1
    ns = 'N' if int(ext[1]//4) >= 0 else 'S'
        
    path1 = os.path.join(dirin,ns+zm1+str(num1)+extt)
    path2 = os.path.join(dirin,ns+zm1+str(num2)+extt)
                    
    ns = 'N' if int(ext[3]//4) >= 0 else 'S'
    path3 = os.path.join(dirin,ns+zm2+str(num1)+extt)
    path4 = os.path.join(dirin,ns+zm2+str(num2)+extt)
    
    if num1//10 ==0:
        path1 = os.path.join(dirin,ns+zm1+'0'+str(num1)+extt)
        path3 = os.path.join(dirin,ns+zm2+'0'+str(num1)+extt)
    if num2//10 ==0:
        path2 = os.path.join(dirin,ns+zm1+'0'+str(num2)+extt)
        path4 = os.path.join(dirin,ns+zm2+'0'+str(num2)+extt)    
    
    biao1 = [path1,path2,path3,path4]
    biao = []
    for i in range(len(biao1)):
        if os.path.exists(biao1[i]):
            biao.append(biao1[i])
        else:
            print(biao1[i])
            with open(txt,'a') as f:
                f.write(biao1[i])
    biao = np.unique(biao)
    print('miss')
    print(biao)
    if (mod ==1) & (len(biao) !=0):
        # j = len(biao)-1
        # a = gpd.read_file(biao[j])
        # while j > 0:
        #     b = gpd.read_file(biao[j-1])
        #     a = gpd.GeoDataFrame(pd.concat([a,b], ignore_index=True))
        #     j = j-1
        a = gpdreaddir1(biao,tifpath)
        return a
    else:
        return biao

# %%

def grwlext(ext,dirin,txt, extt = '.tif'):
    '''
    gdbin: geopandas read_file
    dirin: dir stored grwl_vector
    txt: txt file to store no exist grwl vector
    mod: control return concated gdb or not, default 1 is return, otherwise, return path 
    
    import numpy as np
    import os
    import geopandas as gpd
    '''
    #ext = gdbin.geometry.total_bounds
    zimu="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    zm1 = zimu[int(abs(ext[1])//4)]
    zm2 = zimu[int(abs(ext[3])//4)]
    num1 = int((ext[0]+180)//6)+1
    num2 = int((ext[2]+180)//6)+1
    ns = 'N' if int(ext[1]//4) >= 0 else 'S'
        
    path1 = os.path.join(dirin,ns+zm1+str(num1)+extt)
    path2 = os.path.join(dirin,ns+zm1+str(num2)+extt)
                    
    ns = 'N' if int(ext[3]//4) >= 0 else 'S'
    path3 = os.path.join(dirin,ns+zm2+str(num1)+extt)
    path4 = os.path.join(dirin,ns+zm2+str(num2)+extt)
    
    if num1//10 ==0:
        path1 = os.path.join(dirin,ns+zm1+'0'+str(num1)+extt)
        path3 = os.path.join(dirin,ns+zm2+'0'+str(num1)+extt)
    if num2//10 ==0:
        path2 = os.path.join(dirin,ns+zm1+'0'+str(num2)+extt)
        path4 = os.path.join(dirin,ns+zm2+'0'+str(num2)+extt)    
    
    biao1 = [path1,path2,path3,path4]
    biao = []
    for i in range(len(biao1)):
        if os.path.exists(biao1[i]):
            biao.append(biao1[i])
        else:
            print(biao1[i])
            with open(txt,'a') as f:
                f.write(biao1[i])
    biao = np.unique(biao)
    #print('miss')
    print(biao)
    return biao


# %%
import rasterio
import numpy as np
from rasterio.merge import merge
from rasterio.transform import from_origin

# 输入栅格文件路径列表
#raster_paths = ["path/to/raster1.tif", "path/to/raster2.tif", ...]
#output_raster_path = "path/to/output/result.tif"
# 打开所有栅格数据并读取
def raster_read(raster):
    with rasterio.open(raster) as dst:
        data = dst.read(1)
    return data
def mergeraster(raster_paths, output_raster_path):
    #datasets = [rasterio.open(path) for path in raster_paths]
    #raster_data_list = [src.read(1) for src in datasets]
    raster_data_list = [ raster_read(path) for path in raster_paths]
    # 逐个进行逻辑运算
    result = raster_data_list[0]  # 初始化为第一个栅格数据
    for raster_data in raster_data_list[1:]:
        result = np.logical_or(result, raster_data)

    # 获取第一个栅格数据的元数据（空间信息、大小等）
    with rasterio.open(raster_paths[0]) as first_dataset:
        profile = first_dataset.profile

    # 创建输出栅格数据

    with rasterio.open(output_raster_path, 'w', **profile) as dst:
        dst.write(result, 1)

    # 关闭所有栅格数据集
    #for dataset in datasets:
    #    dataset.close()

# %%
import rasterio
import numpy as np
from rasterio.merge import merge
from rasterio.transform import from_origin

# 输入栅格文件路径列表
#raster_paths = ["path/to/raster1.tif", "path/to/raster2.tif", ...]
#output_raster_path = "path/to/output/result.tif"
# 打开所有栅格数据并读取
def raster_read(raster):
    with rasterio.open(raster) as dst:
        data = dst.read(1)
    return data
def mergeraster(raster_paths, output_raster_path):
    #datasets = [rasterio.open(path) for path in raster_paths]
    #raster_data_list = [src.read(1) for src in datasets]
    raster_data_list = [ raster_read(path) for path in raster_paths]
    # 逐个进行逻辑运算
    result = raster_data_list[0]  # 初始化为第一个栅格数据
    for raster_data in raster_data_list[1:]:
        result = np.logical_or(result, raster_data)

    # 获取第一个栅格数据的元数据（空间信息、大小等）
    with rasterio.open(raster_paths[0]) as first_dataset:
        profile = first_dataset.profile

    # 创建输出栅格数据

    with rasterio.open(output_raster_path, 'w', **profile) as dst:
        dst.write(result, 1)

    # 关闭所有栅格数据集
    #for dataset in datasets:
    #    dataset.close()

# %%
import chu
from chu import raster2shp1
#dirgrwl = r"P:\tif\0\temp"

    
def grwltoraster(grwl11, intif,dirgrwl,dirmatch):
    raster_paths = []
    # grwl river shp
    for i in range(len(grwl11)):
        existing_raster_path = grwl11[i]  
        output_raster_path = os.path.join(dirgrwl,os.path.split(existing_raster_path)[1])
        if not os.path.exists(output_raster_path):
            with rasterio.open(existing_raster_path) as dst:
                a = dst.read()
                profile = dst.profile
                transform = dst.transform
                crs = dst.crs
            mask = np.where((a == 255) | (a == 86) | (a == 126), 1, 0)
            with rasterio.open(output_raster_path, 'w', **profile) as dst:
                dst.write(mask)
            raster2shp1(output_raster_path,model=1,dir2 = '')
        
        # grwl river shp match with intif
        shapefile_path = os.path.splitext(output_raster_path)[0]+'.shp'
        num = os.path.split(grwl11[i])[1][0:-4]
        output_raster_path1 = os.path.join(dirmatch, num+'_'+os.path.split(intif)[1])
        
        shptoraster(intif,shapefile_path,output_raster_path1, invert = True)
        
        raster_paths.append(output_raster_path1)
        
        #merge all matched grwl river
    output_raster_path2 =  os.path.join(dirmatch, 'grwl'+'_'+os.path.split(intif)[1])
    mergeraster(raster_paths, output_raster_path2)
        
    return output_raster_path2 

# %%
def countraster(rasterpath, num = 1):
    with rasterio.open(rasterpath) as src:
        # 读取栅格数据
        raster_data = src.read(1)
        cnum = np.sum(raster_data == num)
    return cnum

# %%
def extofraster(tifpath):
    with rasterio.open(tifpath) as src:
        bounds = src.bounds
    ext = [bounds.left,bounds.bottom,bounds.right,bounds.top]
    return ext

# %%
def cdir(name):
    os.makedirs(name, exist_ok=True)
    return name
