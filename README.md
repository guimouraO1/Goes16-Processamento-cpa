# NOAA GOES 16 on CEPAGRI 

The National Oceanic and Atmospheric Administration (NOAA) operates a constellation of Geostationary Operational Environmental Satellites (GOES) to provide continuous weather imagery and monitoring of meteorological and space environment data for the protection of life e outros fins. GOES satellites provide critical atmospheric, oceanic, climatic and space weather products supporting weather forecasting and warnings, climatologic analysis and prediction, ecosystems management, safe and efficient public and private transportation, and other.

The satellites provide advanced imaging with increased spatial resolution, 16 spectral channels, and up to 1 minute scan frequency for more accurate forecasts and timely warnings.

The real-time feed and full historical archive of original resolution Advanced Baseline Imager (ABI) radiance data (Level 1b) and full resolution Cloud and Moisture Imager (CMI) products (Level 2) are freely available on cpa.unicamp.br.


## About the Data
All data files from GOES-16 are provided in netCDF4 format. The GOES-16 data is hosted in the `noaa-goes16`.

`<Product>/<Year>/<Day of Year>/<Hour>/<Filename>`

where:

- `<Product>` is the product generated from CEPAGRI

  - ABI-L1b-RadF - Advanced Baseline Imager Level 1b Full Disk
  - ABI-L2-CMIPF - Advanced Baseline Imager Level 2 Cloud and Moisture Imagery Full Disk
  - ABI-L2-FDCF  - Advanced Baseline Imager Level 2 Fire (Hot Spot Characterization) Full Disk
  - ABI-L2-NDVI  - Normalized difference vegetation index
  - ABI-L2-TC    - Natural True Color Level 2 (Band 01, Band 02, Band 03)
  - ABI-L2-RRQPEF - Advanced Baseline Imager Level 2 Rainfall Rate (Quantitative Precipitation Estimate) Full Disk
  - GLM-L2-LCFA - Geostationary Lightning Mapper Level 2 Lightning Detection 
  
- `<Year>` is the year the netCDF4 file was created
- `<Day of Year>` is the numerical day of the year (1-365) julian day
- `<Hour>` is the hour the data observation was made
- `<Filename>` is the name of the file containing the data. These are compressed and encapsulated using the netCDF4 standard.

A `<Filename>` is delineated by underscores '_' and looks like this:

`OR_ABI-L1b-RadF-M3C02_G16_s20231671145342_e20231671156109_c20231671156144.nc`

where:

- `OR`: Operational system real-time data
- `ABI-L1b-RadF-M3C02` is delineated by hyphen '-':
  - `ABI`: is ABI Sensor
  - `L1b`: is processing level, L1b data or L2
  - `Rad`: is radiances. Other products include CMIP (Cloud and Moisture Imagery products) and MCMIP (multichannel CMIP).
  - `F`: is full disk (normally every 15 minutes), C is continental U.S. (normally every 5 minutes), M1 and M2 is Mesoscale region 1 and region 2 (usually every minute each)
  - `M3`: is mode 3 (scan operation), M4 is mode 4 (only full disk scans every five minutes – no mesoscale or CONUS)
  - `C02`: is channel or band 02, There will be sixteen bands, 01-16
- `G16`: is satellite id for GOES-16 (future G17)
- `s20231671145342`: is start of scan time
  - 4 digit year
  - 3 digit day of year
  - 2 digit hour
  - 2 digit minute
  - 2 digit second
  - 1 digit tenth of second
- `e20231671156109`: is end of scan time
- `c20231671156144`: is netCDF4 file creation time
- `.nc` is netCDF file extension

---

# GOES-16 NetCDF Image Manipulation

This repository contains a Python script for manipulating images from the GOES-16 satellite in NetCDF format. 
The script processes these NetCDF to generate images, GIFs and other products.

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [About code](#About-code)
- [Processamento.py](#processamento.py)
- [License](#license)

---

## Prerequisites

Before using the script, make sure you have the following packages and tools installed:

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
In linux
```bash
apt install ffmpeg
```
You can create a Conda environment with required packages using the following command:

```bash
conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp pyspectral pyorbital pillow
```
---

To activate this environment, use
```bash
conda activate goes
```
To deactivate an active environment, use
```bash
conda deactivate
```

## Installation

1. Clone this repository to your local machine:

```bash
https://github.com/guimouraO1/processamentopy.git
```

2. Navigate to the project directory:

```bash
cd your-repo
```
---

## Configuration

You can configure the script behavior by modifying the variables in the `dirs.py` script:

```python
dirs.py
```

```python
dir_main = '/home/myDirs/'
dirs = {
    'dir_in': dir_main + 'goes/',
    'dir_main': dir_main,
    'dir_out': dir_main + 'output/',
    'dir_libs': dir_main + 'libs/',
    'dir_shapefiles': dir_main + 'shapefiles/',
    'dir_colortables': dir_main + 'colortables/',
    'dir_logos': dir_main + 'logos/',
    'dir_temp': dir_main + 'temp/',
    'arq_log': dir_main + 'logs/Processamento-GOES_' + str(datetime.date.today()) + '.log'
}
```

Make sure to update the paths in the `dirs` dictionary to match your file system structure and requirements.

---

## Usage

The directory structure of this repository is as follows:

```
.

├── shapefiles
├── output
├── colortables
├── goes
├── logos
├── logs
├── modules/
│   ├── check_new_images.py
│   ├── dirs.py
│   ├── logs.py
│   ├── process_gif.py
│   ├── processamento.py
│   ├── quantity_products.py
│   ├── remap.py
│   ├── remove.py
│   ├── send_products.py
│   ├── utilities.py
├── oldBands.json             
├── download_amazon.py        
├── main.py                   
├── processamento.sh          
└── README.md
```

To use the script, follow these steps:

1. Ensure you have the required NetCDF images in the input directory (`dir_in/band{1-16}`).

2. Run the script:

```bash
/opt/miniconda3/envs/goes/bin/python3 processamento.py
```
3. The script will process the images, generate GIFs, and perform other operations based on the configuration.

---

## Processamento.py

### Code Structure

- The code is highly modularized for ease of maintenance and extensibility. Below, we detail the key functionalities of the code:
---
### Directory Configuration

- The `dirs.py` file defines input, output, and temporary directories where images will be read, processed, and stored.
---
### Band Processing Control

- A dictionary called `bands` is created to represent the processing status of each image band. Initially, all bands are marked as unprocessed (False).
---
### Log Configuration and Time Tracking

- The `conf_log` function is called to configure log generation. Next, the `start` variable is initialized to record the script's start time.
---
### Checking for New Images

- The `checar_imagens` function is called to check if there are new images to be processed. The `bands` dictionary is updated to reflect this information.
---
### Image Processing

- If there is at least one new image to process (any band with a True value in the dictionary), the `processamento_das_imagens` function is called to perform image processing.
---
### Storage of Processed Images

- Processed images are stored in the output directory (`dir_out`).
---
### Removal of Processed Images

- After processing, the `remover_imagens` function is called to remove the `.nc` files that have already been processed from the input directory (`dir_in`).
---
### Product Quantity Control

- The `quantity_products` function is called to control the quantity of products (images) to be retained for producing an animated GIF.
---
### Creation of Animated GIF

- Next, the `process_gif` function is called to create an animated GIF from the processed images.
---
### Sending Processed Products

- The `send_products` function sends the processed images to a specific site (cpa.unicamp.br).


---

## Example: NOAA GOES-16 truecolor Processing on CEPAGRI

![GOES-16 Satellite](shapefiles/truecolor.png)

---

## Acknowledgments

This script is developed as part of an image processing project for the GOES-16 satellite data. It builds upon various open-source libraries and tools.

---

## License

This project is licensed under the []().

---
**Readme.md Author:** [Guilherme de Moura Oliveira]
**Contact:** [guimoura@unicamp.br]
**Last Updated:** [16/08/2023]
