
import numpy as np
import matplotlib.pyplot as plt
import pointCollection as pc
import glob
import re
import os
import shutil
import argparse
import sys
from datetime import datetime, timedelta, timezone
from importlib import resources
import csv
import h5py
import uuid
import subprocess
from ATL11.h5util import create_attribute
from ATL11 import ATL11xo
from ATL11.version import xosoftwareVersion,xosoftwareDate,xosoftwareTitle,xoidentifier,xoseries_version

def make_queue(args):

    for cycle in range(1, args.cycle+1):
        print(f'make_ATL11xo_tiles.py --top_dir {args.top_dir} --dest_dir {args.dest_dir} --release {args.release} --version {args.version} --cycle {cycle} --region {args.region} --ref_cycles {args.ref_cycles[0]} {args.ref_cycles[1]}')

def write_meta_fields(D, h5f, ref_cycles, cycle):
    ''' Write the metadata fields to an hdf5 file handle '''

    maxlon = h5f['/'].attrs['geospatial_lon_max']
    minlon = h5f['/'].attrs['geospatial_lon_min']
    if minlon*maxlon < 0 and maxlon > 150.0:
        eastlon = minlon
        westlon = maxlon
    else:
        eastlon = maxlon
        westlon = minlon
    h5f['METADATA/Extent'].attrs['eastBoundLongitude'] = eastlon
    h5f['METADATA/Extent'].attrs['westBoundLongitude'] = westlon
    h5f['METADATA/Extent'].attrs['northBoundLatitude'] = h5f['/'].attrs['geospatial_lat_max']
    h5f['METADATA/Extent'].attrs['southBoundLatitude'] = h5f['/'].attrs['geospatial_lat_min']
    create_attribute(h5f['/'].id, 'identifier_file_uuid', [], str(uuid.uuid4()))
    create_attribute(h5f['/'].id, 'date_created', [], str(datetime.now().date())+'T'+str(datetime.now().time())+'Z')
    create_attribute(h5f['/'].id, 'history', [], str(datetime.now().date())+'T'+str(datetime.now().time())+'Z')
    create_attribute(h5f['/'].id, 'hdfversion', [], str(h5py.version.hdf5_version))
    create_attribute(h5f.id, 'identifier_product_format_version', [], xosoftwareVersion())

    create_attribute(h5f['METADATA/DatasetIdentification'].id, 'uuid', [], str(uuid.uuid4()))
    create_attribute(h5f['METADATA/DatasetIdentification'].id, 'creationDate', [], str(datetime.now().date())+'T'+str(datetime.now().time())+'Z')
    create_attribute(h5f['METADATA/DatasetIdentification'].id, 'fileName', [], str(os.path.basename(h5f.filename)))
    create_attribute(h5f['METADATA/DatasetIdentification'].id, 'VersionID', [], xosoftwareVersion())
    create_attribute(h5f['METADATA/Lineage/Control'].id, 'control', [], ' '.join(sys.argv))
    create_attribute(h5f['METADATA/Lineage/Control'].id, 'description', [], 'Exact command line execution of ICESat-2/ATL11 algorithm providing all of the conditions required for each individual run of the software.')
    create_attribute(h5f['METADATA/Lineage/Control'].id, 'shortName', [], 'CNTL')
    create_attribute(h5f['METADATA/Lineage/Control'].id, 'version', [], '1')

    create_attribute(h5f['METADATA/Lineage/Control'].id, 'control', [], ' '.join(sys.argv))
    create_attribute(h5f['METADATA/SeriesIdentification'].id, 'VersionID', [], xoseries_version())

    create_attribute(h5f['METADATA/ProcessStep/PGE'].id, 'stepDateTime', [], str(datetime.now().date())+'T'+str(datetime.now().time())+'Z')
    create_attribute(h5f['METADATA/ProcessStep/PGE'].id, 'softwareVersion', [], xosoftwareVersion())
    create_attribute(h5f['METADATA/ProcessStep/PGE'].id, 'runTimeParameters', [], ' '.join(sys.argv))
    create_attribute(h5f['METADATA/ProcessStep/PGE'].id, 'identifier', [], xoidentifier())
    create_attribute(h5f['METADATA/ProcessStep/PGE'].id, 'softwareDate', [], xosoftwareDate())
    create_attribute(h5f['METADATA/ProcessStep/PGE'].id, 'softwareTitle', [], xosoftwareTitle())

    g2 = h5f['ancillary_data']
    g2['atlas_sdp_gps_epoch'][...] = np.array([1198800018.], dtype=np.float64)
    g2['start_delta_time'][...] = np.array([np.nanmin(D.delta_time)])
    g2['end_delta_time'][...] = np.array([np.nanmax(D.delta_time)])
    g2['start_geoseg'][...] = np.array([np.nanmin(D.segment_id).astype(int)])
    g2['end_geoseg'][...] = np.array([np.nanmax(D.segment_id).astype(int)])
    g2['start_rgt'][...] = np.array([np.nanmin(D.rgt).astype(int)])
    g2['end_rgt'][...] = np.array([np.nanmax(D.rgt).astype(int)])
    ctl = h5f['METADATA']['Lineage']['Control'].attrs['control'].decode()
    g2['control'][...] = ctl.encode('ASCII','replace')
    g2['release'][...] = os.path.basename(h5f.filename).split('_')[-2].encode('ASCII','replace')
    g2['version'][...] = os.path.splitext(os.path.basename(h5f.filename))[0].split('_')[-1].encode('ASCII','replace')

    epoch = datetime(2018, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
    start_delta_time = timedelta(seconds=np.array([np.nanmin(D.delta_time)])[0])
    start_utc_time = epoch + start_delta_time
    str_start_utc = (str(start_utc_time.date())+'T'+
                    start_utc_time.strftime("%H:%M:%S.%f")+'Z')
    end_delta_time = timedelta(seconds=np.array([np.nanmax(D.delta_time)])[0])
    end_utc_time = epoch + end_delta_time
    str_end_utc = (str(end_utc_time.date())+'T'+
                    end_utc_time.strftime("%H:%M:%S.%f")+'Z')
    create_attribute(h5f['/'].id, 'time_coverage_start', [], str_start_utc)
    create_attribute(h5f['/'].id, 'time_coverage_end', [], str_end_utc)
    h5f['/'].attrs['time_coverage_duration'] = np.array([np.nanmax(D.delta_time)])[0] - np.array([np.nanmin(D.delta_time)])[0]
    g2['data_start_utc'][...] = str_start_utc.encode('ASCII','replace')
    g2['granule_start_utc'][...] = str_start_utc.encode('ASCII','replace')
    g2['data_end_utc'][...] = str_end_utc.encode('ASCII','replace')
    g2['granule_end_utc'][...] = str_end_utc.encode('ASCII','replace')
    create_attribute(h5f['METADATA/Extent'].id, 'rangeBeginningDateTime', [], str_start_utc)
    create_attribute(h5f['METADATA/Extent'].id, 'rangeEndingDateTime', [], str_end_utc)

    h5f.attrs['ref_surf_cycles'] = ref_cycles
    h5f.attrs['cycle'] = cycle

def parse_attr_file():
    # read in the metadata
    attrfile=os.path.join(str(resources.files('ATL11')), 'package_data', 'ATL11xo_output_attrs.csv')
    all_field_attrs=list(csv.DictReader(open(attrfile, encoding='utf-8-sig')))

    group_attrs={}
    group_descriptions={}
    group_dimensions={}
    # make a dictionary listing attributes for each group
    for field_attrs in all_field_attrs:
        group_name = field_attrs['group']
        if len(group_name)==0:
            # blank line
            continue
        if group_name not in group_attrs:
            group_attrs[group_name] = {}
            group_dimensions[group_name] = {}
        this_group = group_attrs[group_name]
        this_dimensions = group_dimensions[group_name]
        field = field_attrs['field']
        if len(field) == 0:
            # zero-length field indicates a group description
            group_descriptions[group_name] = field_attrs['description']
            continue
        if field not in this_group:
            this_group[field] = {}
            this_dimensions[field] = {}
        for attr, val in field_attrs.items():
            if attr in ['field','group']:
                continue
            if attr == 'dimension':
                group_dimensions[group_name][field]['dimension'] =\
                        val in ['yes','YES','True', 'true','TRUE', 1]
                continue
            if attr == 'dimensions':
                group_dimensions[group_name][field]['dimensions'] = val
                continue
            #'valid_min' and 'valid_max' fields are numeric
            if 'valid_m' in attr:
                # make sure valid_max and valid_min match the variable's datatype
                if val == '':
                    continue
                this_group[field][attr] = np.dtype(field_attrs['datatype'].lower()).type(val)
            else:
                # otherwise write a string
                this_group[field][attr] = str(val)
        if field_attrs['datatype'].startswith('int'):
            this_group[field]['_FillValue'] = np.iinfo(np.dtype(field_attrs['datatype'])).max
        elif field_attrs['datatype'].startswith('float'):
            this_group[field]['_FillValue'] = np.finfo(np.dtype(field_attrs['datatype'])).max
        if 'standard_name' not in field_attrs or field_attrs['standard_name' ] is None:
            this_group[field]['standard_name'] = field
    return group_attrs, group_descriptions, group_dimensions

def make_dimensions(h5f, group_dimensions):
    # first make the dimension(s)
    for group_key in ['ROOT','datum_track','crossing_track']:
        if group_key=='ROOT':
            dst=h5f
        else:
            dst=h5f[group_key]
        for field, dim_dict in group_dimensions[group_key].items():
            if not dim_dict['dimension']:
                continue
            dst[field].make_scale()
    # next attach the dimensions to the variables:
    for group_key in ['ROOT','datum_track','crossing_track']:
        if group_key=='ROOT':
            dst=h5f
            group=''
        else:
            dst=h5f[group_key]
            group=group_key
        for field, dim_dict in group_dimensions[group_key].items():
            # don't add a field to the dimension dictionary if it's a dimension itself
            if dim_dict['dimension']:
                continue
            dims = dim_dict['dimensions']
            if isinstance(dims, str):
                dims = dims.split(',')
            dset = dst[field]
            for ind, dim in enumerate(dims):
                if '../' in dim:
                    group_path = group.split('/')
                    while '../' in dim:
                        dim=dim.lstrip('../')
                        group_path=group_path[:-1]
                    try:
                        dset.dims[ind].attach_scale(h5f["/"+'/'.join(group_path+[dim])])
                    except Exception as e:
                        print("-----")
                        print([group,field])
                        print('/'.join(group_path+[dim]))
                        print("------")
                        raise e
                else:
                    dset.dims[ind].attach_scale(h5f[dim])
                dset.dims[ind].label=dim

def write_data(out_file, xyT, D_cache, args, group_attrs, group_descriptions, group_dimensions):
    atl11xo_template = str(resources.files('ATL11').joinpath("package_data/atl11xo_template.h5"))

    if os.path.isfile(out_file):
        os.remove(out_file)
    try:
        shutil.copyfile(atl11xo_template, out_file)
    except PermissionError:
        print("Error: Permission denied. Cannot write to {out_file}.")
    except FileNotFoundError:
        print(f"Error: Template file not found at {atl11xo_template}.")
    except shutil.Error as e:
        print(f"Error during template copy operation to {out_file}: {e}")
    except Exception as e:
        print(f'Error: {e} Failed to copy template to {out_file}')
    # write the data
    with h5py.File(out_file,'a') as fh:
    # Clean up group not needed
        del fh['orbit_info']
        for group in ['ROOT','datum_track','crossing_track']:
            out_group = '/' if group=='ROOT' else group
            Dsub = D_cache[group][xyT]
            Dsub.to_h5(out_file, 
#                       h5f_out=fh,
                       group=out_group,
                       replace = False,
                       meta_dict = group_attrs[group])
            if group in group_descriptions:
                fh[out_group].attrs['description'.encode('ascii')] =\
                    group_descriptions[group].encode('ascii')
            if group=='crossing_track':
                # this group contains delta_time, segment, and rgt
                write_meta_fields(Dsub, fh, args.ref_cycles, args.cycle)
            if 'delta_time' in fh[out_group]:
                fh[out_group]['delta_time'].attrs['standard_name'] = 'time'
                fh[out_group]['delta_time'].attrs['calendar'] = 'standard'
            if group=='ROOT':
                # do this for the last group written
                fh.attrs['geospatial_lon_min'] = np.nanmin(Dsub.longitude)
                fh.attrs['geospatial_lon_max'] = np.nanmax(Dsub.longitude)
                fh.attrs['geospatial_lat_min'] = np.nanmin(Dsub.latitude)
                fh.attrs['geospatial_lat_max'] = np.nanmax(Dsub.latitude)
                create_attribute(fh['/'].id, 'geospatial_ellipsoid', [], 'WGS84')

        make_dimensions(fh, group_dimensions)

def post_process(xover_output_dir):
    project_bin="/discover/nobackup/bjelley/bin"
    atlas_meta = project_bin+"/atlas_meta"
    atl11xo_qa_util = project_bin+"/atl11xo_qa_util"

    if not os.path.isdir(xover_output_dir):
        print('ERROR: xover_output_dir does not exist')
        exit(1)


    ATL11xo_glob=xover_output_dir+'/cycle_??/ATL11xo_*_??.h5'

    ATL11xo_list = []
    ATL11xo_list = glob.glob(ATL11xo_glob)

# Read controlfile template
    with open(project_bin+"/ATL11xo_AX_EXXXX_NXXXX_cXX_0XX_XX.ctl", "r") as ctltemplate:
        ctllines = ctltemplate.readlines()

    for file in ATL11xo_list:
        control = file.replace('.h5', '.ctl')
        with open(control, "w") as ctlfile:
            for line in ctllines:
                ctlfile.write(re.sub(r'_atl11file_', file, line))
# Run atlas_meta
        try:
            result = subprocess.run([atlas_meta, control], capture_output=True, text= True)
        except subprocess.CalledProcessError as e:
            print(f"{atlas_meta} failed with error: {e}")
# Run qa utility
        try:
            result = subprocess.run([atl11xo_qa_util, control], capture_output=True, text= True)
        except subprocess.CalledProcessError as e:
            print(f"{atl11xo_qa_util} failed with error: {e}")

def main():
    parser=argparse.ArgumentParser(description='Generate a set of ATL11xo tiles from a directory of ATL11_atxo along-track crossover files')
    parser.add_argument('--top_dir', type=str, required=True, help='top directory containing ATL11_atxo files')
    parser.add_argument('--dest_dir', type=str, required=False, help='output directory for ATL11xo files')
    parser.add_argument('--EPSG', type=int, required=False, help='EPSG string for output')
    parser.add_argument('--release', type=int, required=True, help='ATL11 release number')
    parser.add_argument('--version', type=int, required=True, help='ATL11xo version number')
    parser.add_argument('--cycle', type=int, required=True, help='cycle number')
    parser.add_argument('--ref_cycles', type=int, nargs=2, required=True, help='first and last reference-track cycles included in the fit')
    parser.add_argument('--queue','-q', action="store_true", help='if set, a queue of commands will be ouput that make tiles for cycles 1...args.cycle')
    parser.add_argument('--region', type=str, required=True, help='region for output, AA=Antarctic, AR=Arctic')
    parser.add_argument('--tile_spacing', type=float, default=200000, help='tile spacing,  m')
    parser.add_argument('--max_files', type=int, default=-1, help='if specified, limit number of ATL11 files to this, default is -1 (all_files)')
    parser.add_argument('--post_process', type=bool, default=False, help='if true, run atlas_meta and QA utility')
    args=parser.parse_args()

    if args.dest_dir is None:
        args.dest_dir = args.top_dir

    if args.queue:
        make_queue(args)
        sys.exit(0)

    if args.EPSG is None:
        if args.region=='AA':
            args.EPSG=3031
        else:
            args.EPSG=3413

    tile_out_dir = os.path.join(args.dest_dir, f'cycle_{args.cycle:02d}')
    try:
        os.mkdir(tile_out_dir)
    except FileExistsError:
        pass

    tS = pc.tilingSchema(mapping_function_name='floor', tile_spacing=args.tile_spacing, EPSG=args.EPSG,
                        format_str = f'ATL11xo_{args.region}_E%d_N%d_c{args.cycle:02d}_{args.release:03d}_{args.version:02d}')
    schema_file = os.path.join(tile_out_dir, f'{int(args.tile_spacing/1000)}km_tiling_{args.region}.json')
    if not os.path.isfile(schema_file):
        tS.to_json(schema_file)
    tS.directory=tile_out_dir

    group_attrs, group_descriptions, group_dimensions = parse_attr_file()

    index_for_xyT = {}
    files_list=[]
    # store the outputs (N.B: this will mean keeping a second full copy of the data in memory}
    D_cache={'crossing_track':{},'datum_track':{},'ROOT':{}}

    for group in ['crossing_track', 'datum_track','ROOT']:
        D=[]
        for file in glob.glob(args.top_dir+f'/cycle_{args.cycle:02d}/ATL11_atxo*_{args.cycle:02d}_*_*.h5')[:args.max_files]:
            for pair in [1, 2, 3]:
                try:
                    D += [pc.data().from_h5(file, group=f'pt{pair}/{group}')]
                except Exception:
                    pass
        # D is all the data from the input group
        D = ATL11xo.data().from_list(D)
        if group=='crossing_track':
            # Dxy is the location data (stored in crossing_track)
            Dxy = ATL11xo.data().from_dict({'latitude':D.latitude.copy(),
                                       'longitude':D.longitude.copy(),
                                       'delta_time':D.delta_time.copy()}).get_xy(args.EPSG)
            # bin_dict is the spatial index for the data
            bin_dict = tS.tile_xy(data=Dxy, return_dict=True)
        # loop over the spatial index:
        for xyT, ii in bin_dict.items():
            for this_group in D_cache.keys():
                if xyT not in D_cache[this_group]:
                    D_cache[this_group][xyT]={}

            if  group=='ROOT':
                # D_cache['ROOT'] will have already been populated with everything except the spatial information
                Dxy_sub = Dxy[ii]
                Dxy_sub = Dxy_sub[index_for_xyT[xyT]]
                Dxy_sub.assign(xo_index=np.arange(0, Dxy_sub.size, dtype=int))
                for field in ['latitude','longitude','x','y', 'xo_index']:
                    if field not in D_cache['ROOT'][xyT]:
                        D_cache['ROOT'][xyT][field] = getattr(Dxy_sub, field)
                D_cache[group][xyT] = ATL11xo.data().from_dict(D_cache['ROOT'][xyT])
            else:
                # subset the data to the bin
                Dsub=D[ii]
                # make an index that sorts the data by floor(y/10k), then floor(x/10k), then delta_time
                if xyT not in index_for_xyT:
                    Dxy_sub = Dxy[ii]
                    index_for_xyT[xyT] = np.lexsort((Dxy_sub.delta_time,
                                                     np.floor(Dxy_sub.x/1.e4),
                                                     np.floor(Dxy_sub.y/1.e4)))
                # sort the data by the index
                Dsub = Dsub[index_for_xyT[xyT]]
                # copy fields from this group to D_cache['ROOT'] as needed:
                for field in group_attrs['ROOT']:
                    if field in Dsub.fields and field not in Dxy.fields:
                        D_cache['ROOT'][xyT][field] = getattr(Dsub, field)
                        Dsub.fields.remove(field)
                D_cache[group][xyT] = Dsub
    # now loop over output files:
    for xyT in D_cache['ROOT'].keys():
        out_file = tS.tile_filename(xyT)
        write_data(out_file, xyT, D_cache, args, group_attrs, group_descriptions, group_dimensions)

    # Run post processing steps
    if args.post_process:
        post_process(args.dest_dir)

if __name__=="__main__":
    main()
