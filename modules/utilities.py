# Training: Python and GOES-R Imagery: Script 8 - Functions for download data from AWS
# Required modules
import os  # Miscellaneous operating system interfaces
import numpy as np  # Import the Numpy package
import colorsys  # To make convertion of colormaps
import boto3  # Amazon Web Services (AWS) SDK for Python
from botocore import UNSIGNED  # boto3 config
from botocore.config import Config  # boto3 config
import math  # Mathematical functions
from datetime import datetime, timedelta  # Basic Dates and time types
from osgeo import osr  # Python bindings for GDAL
from osgeo import gdal  # Python bindings for GDAL


def load_cpt(path):
    try:
        f = open(path)
    except:
        print("File ", path, "not found")
        return None

    lines = f.readlines()

    f.close()

    x = np.array([])
    r = np.array([])
    g = np.array([])
    b = np.array([])

    colorModel = 'RGB'

    for l in lines:
        ls = l.split()
        if l[0] == '#':
            if ls[-1] == 'HSV':
                colorModel = 'HSV'
                continue
            else:
                continue
        if ls[0] == 'B' or ls[0] == 'F' or ls[0] == 'N':
            pass
        else:
            x = np.append(x, float(ls[0]))
            r = np.append(r, float(ls[1]))
            g = np.append(g, float(ls[2]))
            b = np.append(b, float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

        x = np.append(x, xtemp)
        r = np.append(r, rtemp)
        g = np.append(g, gtemp)
        b = np.append(b, btemp)

    if colorModel == 'HSV':
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
        r[i] = rr;
        g[i] = gg;
        b[i] = bb

    if colorModel == 'RGB':
        r = r / 255.0
        g = g / 255.0
        b = b / 255.0

    xNorm = (x - x[0]) / (x[-1] - x[0])

    red = []
    blue = []
    green = []

    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])

    colorDict = {'red': red, 'green': green, 'blue': blue}

    return colorDict


def download_cmi(yyyymmddhhmn, band, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'
    product_name = 'ABI-L2-CMIPF'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(band):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")


    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Band-{band}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}'


def download_prod(yyyymmddhhmn, product_name, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")


    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}{file_name}.nc exists')
            else:
                # print(f'Downloading file {path_dest}{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}'


def download_glm(yyyymmddhhmnss, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%H')
    min = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%M')
    seg = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%S')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # File structure
    product_name = "GLM-L2-LCFA"
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}_G16_s{year}{day_of_year}{hour}{min}{seg}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")


    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmnss}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}'


def reproject(file_name, ncfile, array, extent, undef):
    # Read the original file projection and configure the output projection
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4(ncfile.GetProjectionRef())

    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

    # Reproject the data
    GeoT = ncfile.GetGeoTransform()
    driver = gdal.GetDriverByName('MEM')
    raw = driver.Create('raw', array.shape[0], array.shape[1], 1, gdal.GDT_Float32)
    raw.SetGeoTransform(GeoT)
    raw.GetRasterBand(1).WriteArray(array)

    # Define the parameters of the output file  
    kwargs = {'format': 'netCDF',
              'srcSRS': source_prj,
              'dstSRS': target_prj,
              'outputBounds': (extent[0], extent[3], extent[2], extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}

    # Write the reprojected file on disk
    gdal.Warp(file_name, raw, **kwargs)


def download_cmi_joao(yyyymmddhhmn, band, path_dest, logging):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'
    product_name = 'ABI-L2-CMIPF'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(band):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")


    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        ##print(f'No files found for the date: {yyyymmddhhmn}, Band-{band}')
        logging.info(f'No files found for the date: {yyyymmddhhmn}, Band-{band}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            file_name = key.split('/')[-1].split('.')[0]
            #Alteração problema servidor
            new_file_name = file_name.replace('OR','CG')

            # Download the file
            if os.path.exists(f'{path_dest}/{new_file_name}.nc'):
                ##print(f'File {path_dest}/{file_name}.nc exists')
                logging.info(f'File {path_dest}/{file_name}.nc exists')
            else:
                ##print(f'Downloading file {path_dest}/{file_name}.nc')
                logging.info(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')

            # Alteração problema servidor
            old_file = os.path.join(path_dest, f'{file_name}.nc')
            new_file = os.path.join(path_dest, f'{new_file_name}.nc')
            os.rename(old_file, new_file)

    return f'{file_name}'


def reprojectBruno(reproj_file, reproj_var, reproj_extent, reproj_resolution, path_dest):
    
    def get_geot(ex, nlines, ncols):
        # Compute resolution based on data dimension
        resx = (ex[2] - ex[0]) / ncols
        resy = (ex[3] - ex[1]) / nlines
        return [ex[0], resx, 0, ex[3], 0, -resy]

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 '
                               '+lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Abrindo imagem com a biblioteca GDAL
    raw = gdal.Open(f'NETCDF:{reproj_file}:' + reproj_var, gdal.GA_ReadOnly)
    # ds_dqf = gdal.Open(f'NETCDF:{reproj_file}:' + 'DQF', gdal.GA_ReadOnly)

    # Lendo os metadados do cabecalho
    if reproj_var == 'BCM' or 'Mask':  ### O arquivo Clear Sky não possui sacale e offset é um arquivo binário
        metadata = raw.GetMetadata()
        scale = 1
        offset = 0
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    else:  # Alteração João 26/08/2022
        metadata = raw.GetMetadata()
        scale = float(metadata.get(reproj_var + '#scale_factor'))
        offset = float(metadata.get(reproj_var + '#add_offset'))
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    # Setup projection and geo-transformation
    raw.SetProjection(source_prj.ExportToWkt())
    # raw.SetGeoTransform(raw.GetGeoTransform())
    GOES16_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    raw.SetGeoTransform(get_geot(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))

    # Compute grid dimension
    KM_PER_DEGREE = 111.32
    sizex = int(((reproj_extent[2] - reproj_extent[0]) * KM_PER_DEGREE) / reproj_resolution)
    sizey = int(((reproj_extent[3] - reproj_extent[1]) * KM_PER_DEGREE) / reproj_resolution)

    # Get memory driver
    driver = gdal.GetDriverByName('MEM')

    # Create grid
    grid = driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
    # dqf = driver.Create('dqf', sizex, sizey, 1, gdal.GDT_Float32)

    # Setup projection and geo-transformation
    grid.SetProjection(target_prj.ExportToWkt())
    grid.SetGeoTransform(get_geot(reproj_extent, grid.RasterYSize, grid.RasterXSize))


    # Perform the projection/resampling
    gdal.ReprojectImage(raw, grid, source_prj.ExportToWkt(), target_prj.ExportToWkt(), gdal.GRA_NearestNeighbour,
                        options=['NUM_THREADS=ALL_CPUS'])
    
    
    # Close files
    raw = None
    del raw

    # Read grid data
    array = grid.ReadAsArray()

    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)

    # Aplicando scale, offset
    array = array * scale + offset


    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    # Define the parameters of the output file
    kwargs = {'format': 'netCDF',
              'dstSRS': target_prj,
              'outputBounds': (reproj_extent[0], reproj_extent[3], reproj_extent[2], reproj_extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}

    reproj_file = reproj_file.split('\\')
    reproj_file.reverse()
    r_file = reproj_file[0].replace('.nc', f'_reproj_.nc')
    gdal.Warp(f'{path_dest}/{r_file}', grid, **kwargs)

    # os.remove(f'{input}{reproj_file[0]}')
    # return file_dtime, file_satellite, grid
    caminho = F'{path_dest}{r_file}'
    
    return caminho


def download_dmw(yyyymmddhhmn,Ch,path_dest):
    product_name = 'ABI-L2-DMWF'
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(Ch):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")


    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}.nc'