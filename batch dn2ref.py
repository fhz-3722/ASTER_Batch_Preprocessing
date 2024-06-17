# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:38:18 2022

@author: Hanze Fu

This script provides two conversions including DN value to radiance-at-sensor and radiance to TOA reflectance. 
"""

from osgeo import gdal, osr
import numpy as np
from datetime import datetime
import os, glob, sys, getopt, argparse, re
import matplotlib.pyplot as plt

from spectral import *
import spectral.io.envi as evni

in_dir = "" # Paste your ASTER dir here, eg. C://Temp or C://Temp//
os.chdir(in_dir)
   
    # Create and set output directory
out_dir = os.path.normpath((os.path.split(in_dir)[0] + os.sep + 'output' ))
rad_out_dir = os.path.normpath((os.path.split(in_dir)[0] + os.sep + 'rad_output'))
ref_out_dir = os.path.normpath((os.path.split(in_dir)[0] + os.sep + 'ref_output'))
rgb_out_dir = os.path.normpath((os.path.split(in_dir)[0] + os.sep + 'rgb_output'))

if not os.path.exists(out_dir): os.makedirs(out_dir)
if not os.path.exists(rad_out_dir): os.makedirs(rad_out_dir)
if not os.path.exists(ref_out_dir): os.makedirs(ref_out_dir)
if not os.path.exists(rgb_out_dir): os.makedirs(rgb_out_dir)
            
    # Create a list of ASTER L1T HDF files in the directory
file_list = glob.glob('AST_L1T_**.hdf')
   
if len(file_list) == 0:
   print('Error: no ASTER L1T hdf files were found in this directory')
   sys.exit(2)
ucc = np.matrix(([[0.676, 1.688, 2.25, 0.0],\
                [0.708, 1.415, 1.89, 0.0],\
                [0.423, 0.862, 1.15, 0.0],\
                [0.1087, 0.2174, 0.2900, 0.2900],\
                [0.0348, 0.0696, 0.0925, 0.4090],\
                [0.0313, 0.0625, 0.0830, 0.3900],\
                [0.0299, 0.0597, 0.0795, 0.3320],\
                [0.0209, 0.0417, 0.0556, 0.2450],\
                [0.0159, 0.0318, 0.0424, 0.2650]]))

    # Thome et al. is used, which uses spectral irradiance values from MODTRAN
    # Ordered b1, b2, b3N, b4, b5...b9
irradiance = [1848, 1549, 1114, 225.4, 86.63, 81.85, 74.85, 66.49, 59.85]
    
def dn2rad (x):
    rad = (x-1.)*ucc1
    return rad
    
def rad2ref (rad):
    ref = (np.pi * rad * (esd * esd)) / (irradiance1 * np.sin(np.pi * sza / 
    180))
    return ref
    # Loop through all ASTER L1T hdf files in the directory


for k in range(len(file_list)):
        
        # Maintains original filename convention    
    file_name = file_list[k] 
        
    print('Processing File: ' + file_name + ' (' + str(k+1) + ' out of ' 
        + str(len(file_list)) + ')')
    
        # Read in the file and metadata
    aster = gdal.Open(file_name)
    aster_sds = aster.GetSubDatasets()
    meta = aster.GetMetadata()
    date = meta['CALENDARDATE']
    dated = datetime.strptime(date, '%Y%m%d')
    day = dated.timetuple()
    doy = day.tm_yday
        
        # Calculate Earth-Sun Distance    
    esd = 1.0 - 0.01672 * np.cos(np.radians(0.9856 * (doy - 4)))    
        
    del date, dated, day, doy
    
    out_filename = '{}\\{}.dat'.format(out_dir,file_name.split('.hdf')[0])
    out_filename_rad = '{}\\{}_rad.dat'.format(rad_out_dir, file_name.split('.hdf')[0])
    out_filename_ref = '{}\\{}_ref.dat'.format(ref_out_dir, file_name.split('.hdf')[0])
    out_filename_rgb = '{}\\{}_rgb.dat'.format(rgb_out_dir, file_name.split('.hdf')[0])
        
        # Need SZA--calculate by grabbing solar elevation info     
    sza = [float(x) for x in meta['SOLARDIRECTION'].split(', ')][1]
      
        # Query gain data for each  band, needed for UCC     
    gain_list = [g for g in meta.keys() if 'GAIN' in g]  ###### AARON HERE
    gain_info = []
    for f in range(len(gain_list)):
        gain_info1 = meta[gain_list[f]].split(', ')#[0] ###### AARON HERE
        gain_info.append(gain_info1)
    gain_dict = dict(gain_info)
    
        # Define UL, LR, UTM zone    
    ul = [float(x) for x in meta['UPPERLEFTM'].split(', ')]
    lr = [float(x) for x in meta['LOWERRIGHTM'].split(', ')]
    utm = int(meta['UTMZONENUMBER'])
    n_s = float(meta['NORTHBOUNDINGCOORDINATE'])
        
        # Create UTM zone code numbers    
    utm_n = [i+32600 for i in range(60)]
    utm_s = [i+32700 for i in range(60)]
        
        # Define UTM zone based on North or South
    if n_s < 0:
        utm_zone = utm_s[utm]
    else:
        utm_zone = utm_n[utm]
            
    del utm_n, utm_s
    #------------------------------------------------------------------------------   
        # Loop through all ASTER L1T SDS (bands)    
    t = aster_sds[0][1]
    t1 = re.findall(r"\d+?\d*", t)
    t2 = [ int(i) for i in t1]
    img_arr = np.empty((t2[0], t2[1], 9))
    img_rad = np.empty((t2[0], t2[1], 9))
    img_ref = np.empty((t2[0], t2[1], 9))
    img_rgb = np.empty((t2[0], t2[1], 3))
        
    for e in range(len(aster_sds)):
        gname = str(aster_sds[e])
            # Maintain original dataset name    
            
        aster_sd = gname.split(',')[0]
            
        vnir = re.search("(VNIR.*)", aster_sd)
        swir = re.search("(SWIR.*)", aster_sd)
        if swir or vnir:
                # Generate output name for tif            
            aster_sd2 = aster_sd.split('(')[1]
            aster_sd3 = aster_sd2[1:-1]
            band = aster_sd3.split(':')[-1]
            #out_filename = '{}/{}_{}.tif'.format(out_dir,file_name.split('.hdf')[0],band)
            #out_filename_rad = '{}_radiance.tif'.format(out_filename.split('.tif')[0])
            #out_filename_ref = '{}_reflectance.tif'.format(out_filename.split('.tif')[0])
                #out_filename = out_dir + file_name.split('.hdf')[0] + '_' + band + '.tif'
                #out_filename_rad = out_filename.split('.tif')[0] + '_radiance.tif'
                #out_filename_ref = out_filename.split('.tif')[0] + '_reflectance.tif'
                # Open SDS and create array            
            band_ds = gdal.Open(aster_sd3, gdal.GA_ReadOnly)
            sds = band_ds.ReadAsArray(buf_xsize = t2[1], buf_ysize = t2[0]).astype(np.uint16)
                
            ncol, nrow = sds.shape
                
            del aster_sd, aster_sd2, aster_sd3
            
            """
        elif vnir:
            aster_sd2 = aster_sd.split('(')[1]
            aster_sd3 = aster_sd2[1:-1]
            band = aster_sd3.split(':')[-1]
            out_filename = '{}/{}_{}.tif'.format(out_dir,file_name.split('.hdf')[0],band)
            out_filename_rad = '{}_radiance.tif'.format(out_filename.split('.tif')[0])
            out_filename_ref = '{}_reflectance.tif'.format(out_filename.split('.tif')[0])
                #out_filename = out_dir + file_name.split('.hdf')[0] + '_' + band + '.tif'
                #out_filename_rad = out_filename.split('.tif')[0] + '_radiance.tif'
                #out_filename_ref = out_filename.split('.tif')[0] + '_reflectance.tif'
                # Open SDS and create array            
            band_ds = gdal.Open(aster_sd3, gdal.GA_ReadOnly)
            sds = band_ds.ReadAsArray(buf_xsize = nrow, buf_ysize = ncol).astype(np.uint16)
            del aster_sd, aster_sd2, aster_sd3
            """


                   
               # Define extent and provide offset for UTM South zones            
            if n_s < 0:
                ul_y = ul[0] + 10000000
                ul_x = ul[1]
                
                lr_y = lr[0] + 10000000
                lr_x = lr[1]
                    
                # Define extent for UTM North zones            
            else:
                ul_y = ul[0] 
                ul_x = ul[1]
                
                lr_y = lr[0] 
                lr_x = lr[1]
                
                # Query raster dimensions and calculate raster x & y resolution
            y_res = -1 * round((max(ul_y, lr_y)-min(ul_y, lr_y))/ncol)
            x_res = round((max(ul_x, lr_x)-min(ul_x, lr_x))/nrow)
                
                # Define UL x and y coordinates based on spatial resolution           
            ul_yy = ul_y - (y_res/2)
            ul_xx = ul_x - (x_res/2)
                
    #------------------------------------------------------------------------------
                # Start conversions by band (1-9)         
            if band == 'ImageData1':
                bn = -1 + 1                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['01'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['01'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                else:
                    ucc1 = ucc[bn, 2] 
                        
            if band == 'ImageData2':
                bn = -1 + 2                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['02'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['02'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                else:
                    ucc1 = ucc[bn, 2] 
                        
            if band == 'ImageData3N':
                bn = -1 + 3                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['3N'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['3N'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                else:
                    ucc1 = ucc[bn, 2] 
                        
            if band == 'ImageData4':
                bn = -1 + 4                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['04'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['04'] == 'NOR':
                     ucc1 = ucc[bn, 1] 
                elif gain_dict['04'] == 'LO1':
                    ucc1 = ucc[bn, 2] 
                else:
                    ucc1 = ucc[bn, 3] 
                        
            if band == 'ImageData5':
                bn = -1 + 5                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['05'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['05'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                elif gain_dict['05'] == 'LO1':
                    ucc1 = ucc[bn, 2] 
                else:
                    ucc1 = ucc[bn, 3] 
                        
            if band == 'ImageData6':
                bn = -1 + 6                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['06'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['06'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                elif gain_dict['06'] == 'LO1':
                    ucc1 = ucc[bn, 2] 
                else:
                    ucc1 = ucc[bn, 3]
                
                        
            if band == 'ImageData7':
                bn = -1 + 7                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['07'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['07'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                elif gain_dict['07'] == 'LO1':
                    ucc1 = ucc[bn, 2] 
                else:
                    ucc1 = ucc[bn, 3] 
                        
            if band == 'ImageData8':
                bn = -1 + 8                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['08'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['08'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                elif gain_dict['08'] == 'LO1':
                     ucc1 = ucc[bn, 2] 
                else:
                    ucc1 = ucc[bn, 3] 
                        
            if band == 'ImageData9':
                bn = -1 + 9                
                    # Query for gain specified in file metadata (by band)            
                if gain_dict['09'] == 'HGH':
                    ucc1 = ucc[bn, 0] 
                elif gain_dict['09'] == 'NOR':
                    ucc1 = ucc[bn, 1] 
                elif gain_dict['09'] == 'LO1':
                    ucc1 = ucc[bn, 2] 
                else:
                    ucc1 = ucc[bn, 3]        
            
    #------------------------------------------------------------------------------
                #Set irradiance value for specific band
            irradiance1 = irradiance[bn]
            band_number = int(re.findall(r"\d", band)[0])
            rad = dn2rad(sds)
            
            rad[rad == dn2rad(0)] = 0
            
            ref = rad2ref(rad)
            
            img_arr[:,:, band_number-1] = sds
            img_rad[:,:,band_number-1] = rad
            img_ref[:,:,band_number-1] = ref
            
            
            del sds, rad, ref
        #img_rgb = img_rad.take([5,2,0], 2)   #Generate rgb color file
    band_wv = [0.556, 0.661, 0.807, 1.656, 2.167, 2.209, 2.262, 2.336, 2.4 ]
            
    #driver = gdal.GetDriverByName("ENVI")
    
    #out_rad = driver.Create(out_filename_rad, img_rad.shape[1], img_rad.shape[0], img_rad.shape[2], 
    #                    gdal.GDT_Float32)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(utm_zone)
    #out_rad.SetProjection(srs.ExportToWkt())
    #out_rad.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))
    #out_rad.SetMetadata({'Wavelength units': 'Micrometers'})
    #for i in range(1, out_rad.RasterCount+1):
    #    outband = out_rad.GetRasterBand(i)
    #    outband.WriteArray(img_rad[:,:, i-1])
    #    outband.SetNoDataValue(0)
        
    out_rad = None
    
    
    md = { "map info": ["UTM", str(1), str(1), str(ul_xx), str(ul_yy), str(30), str(30), str(utm), 'South', 'WGS-84' ],\
          #'coordinate system string': srs.ExportToWkt(),
          'Wavelength units': 'Micrometers',
          'file type': 'ENVI Standard',
          'wavelength': band_wv,
          'bands': 9,
          
        
        }
    envi.save_image('{}\\{}_rad1.hdr'.format(rad_out_dir, file_name.split('.hdf')[0]),\
                    img_rad.astype('float32'), metadata = md, force = True)

    

    
    
    """
    out_rgb = driver.Create(out_filename_rgb, img_rgb.shape[1], img_rgb.shape[0],\
                            img_rgb.shape[2], gdal.GDT_Float32)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(utm_zone)
    out_rgb.SetProjection(srs.ExportToWkt())
    out_rgb.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))
    for i in range(1, out_rgb.RasterCount+1):
        out_band = out_rgb.GetRasterBand(i)
        out_band.WriteArray(img_rgb[:, :, i-1])
        out_band.FlushCache()
        out_band.SetNoDataValue(0)
    out_rgb = None
    
    del out_rgb
    
    
    
    out = driver.Create(out_filename, img_arr.shape[1], img_arr.shape[0], img_arr.shape[2], 
                        gdal.GDT_UInt16)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(utm_zone)
    out.SetProjection(srs.ExportToWkt())
    out.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))
    for i in range(1, out.RasterCount+1):
        outband = out.GetRasterBand(i)
        outband.WriteArray(img_arr[:,:, i-1])
        outband.SetNoDataValue(0)
    out = None
       
    
    out_ref = driver.Create(out_filename_ref, img_ref.shape[1], img_ref.shape[0], img_ref.shape[2], 
                        gdal.GDT_Float32)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(utm_zone)
    out_ref.SetProjection(srs.ExportToWkt())
    out_ref.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))
    for i in range(1, out_ref.RasterCount+1):
        outband = out_ref.GetRasterBand(i)
        outband.WriteArray(img_ref[:,:, i-1])
        outband.SetNoDataValue(0)
    out_ref = None
    
    
    del img_arr, img_rad, img_ref
    """
