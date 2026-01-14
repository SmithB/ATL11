
import numpy as np
import matplotlib.pyplot as plt
import pointCollection as pc
import glob
import re
import os
import shutil
import argparse
import sys
from importlib import resources
import csv
import h5py
import uuid
from ATL11.h5util import create_attribute

atl11xo_template = str(resources.files('ATL11').joinpath("package_data/atl11xo_template.h5"))

def make_queue(args):
    for cycle in range(1, args.cycle+1):
        print(f'make_ATL11xo_tiles.py --top_dir {args.top_dir} --dest_dir {args.dest_dir} --release {args.release} --version {args.version} --cycle {cycle} --region {args.region} --ref_cycles {args.ref_cycles[0]} {args.ref_cycles[1]}')

def write_meta_fields(D, h5f, ref_cycles, cycle):
    ''' Write the metadata fields to an hdf5 file handle '''

#    g0=h5f.create_group('METADATA')
#    g1=h5f.create_group('METADATA/DatasetIdentification')
#    create_attribute(g1.id, 'uuid', [], str(uuid.uuid4()))
    create_attribute(h5f['METADATA/DatasetIdentification'].id, 'uuid', [], str(uuid.uuid4()))
#    g2 = h5f.create_group('ancillary_data')
    g2 = h5f['ancillary_data']
    g2['start_delta_time'][...] = np.array([np.nanmin(D.delta_time)])
    g2['end_delta_time'][...] = np.array([np.nanmax(D.delta_time)])
    g2['start_geoseg'][...] = np.array([np.nanmin(D.segment_id).astype(int)])
    g2['end_geoseg'][...] = np.array([np.nanmax(D.segment_id).astype(int)])
    g2['start_rgt'][...] = np.array([np.nanmin(D.rgt).astype(int)])
    g2['end_rgt'][...] = np.array([np.nanmax(D.rgt).astype(int)])
#    g2.create_dataset('start_delta_time',data=np.array([np.nanmin(D.delta_time)]))
#    g2.create_dataset('end_delta_time',data=np.array([np.nanmax(D.delta_time)]))
#    g2.create_dataset('start_geoseg',data=np.array([np.nanmin(D.segment_id).astype(int)]))
#    g2.create_dataset('end_geoseg',data=np.array([np.nanmax(D.segment_id).astype(int)]))
#    g2.create_dataset('start_rgt',data=np.array([np.nanmin(D.rgt).astype(int)]))
#    g2.create_dataset('end_rgt',data=np.array([np.nanmax(D.rgt).astype(int)]))
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
                # boolean flag to tell us whether a field is a dimension
                val = val in ['yes','YES','True', 'true','TRUE', 1]
                group_dimensions[group_name][field]['dimension'] = val
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
            field_attrs['_FillValue'] = np.iinfo(np.dtype(field_attrs['datatype'])).max
        elif field_attrs['datatype'].startswith('float'):
            field_attrs['_FillValue'] = np.finfo(np.dtype(field_attrs['datatype'])).max
    return group_attrs, group_descriptions, group_dimensions

def make_dimensions(out_file, group_dimensions):
    with h5py.File(out_file,'a') as h5f:
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

    replace=True
    index_for_xyT = {}
    D_root = {}
    files_list=[]
    for group in ['crossing_track', 'datum_track','ROOT']:
        D=[]
        for file in glob.glob(args.top_dir+f'/ATL11_atxo*_*_{args.cycle:02d}_*.h5')[:args.max_files]:
            for pair in [1, 2, 3]:
                try:
                    D += [pc.data().from_h5(file, group=f'pt{pair}/{group}')]
                except Exception:
                    pass
        D=pc.data().from_list(D)
        if group=='crossing_track':
            Dxy = pc.data().from_dict({'latitude':D.latitude.copy(),
                                       'longitude':D.longitude.copy(),
                                       'delta_time':D.delta_time.copy()}).get_xy(args.EPSG)
            bin_dict = tS.tile_xy(data=Dxy, return_dict=True)
        for xyT, ii in bin_dict.items():
            # choose the out file
            out_file = tS.tile_filename(xyT)
            if out_file not in files_list:
                files_list.append(out_file)
            # If starting new file, copy template (replace=False to avoid overwrite)
            if os.path.isfile(out_file) and group=='crossing_track' :
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
                replace=False
            else:
                replace=False

            if xyT not in D_root:
                D_root[xyT]={}

            if  group=='ROOT':
                # D_root will have already been populated
                out_group='/'
                Dxy_sub = Dxy[ii]
                Dxy_sub = Dxy_sub[index_for_xyT[xyT]]
                Dxy_sub.assign(xo_index=np.arange(0, Dxy_sub.size, dtype=int))
                for field in ['latitude','longitude','x','y', 'xo_index']:
                    if field not in D_root[xyT]:
                        D_root[xyT][field] = getattr(Dxy_sub, field)
                Dsub=pc.data().from_dict(D_root[xyT])
            else:
                out_group=group
                # subset the data to the bin
                Dsub=D[ii]

                # make an index that sorts the data by floor(y/10k), then floor(x/10k), then delta_time
                if xyT not in index_for_xyT:
                    Dxy_sub = Dxy[ii]
                    index_for_xyT[xyT] = np.lexsort((Dxy_sub.delta_time, np.floor(Dxy_sub.x/1.e4), np.floor(Dxy_sub.y/1.e4)))
                # sort the data by the index
                Dsub = Dsub[index_for_xyT[xyT]]
                # copy fields from this group to D_root as needed:
                for field in group_attrs['ROOT']:
                    if field in Dsub.fields and field not in Dxy.fields:
                        D_root[xyT][field] = getattr(Dsub, field)
                        Dsub.fields.remove(field)

            # write the data
            Dsub.to_h5(out_file, group=out_group,
                       replace = replace,
                       meta_dict = group_attrs[group])
            with h5py.File(out_file,'a') as fh:
                if group in group_descriptions:
                    fh[out_group].attrs['description'.encode('ascii')] =\
                        group_descriptions[group].encode('ascii')
                if group=='crossing_track':
                    # this group contains delta_time, segment, and rgt
                    write_meta_fields(Dsub, fh, args.ref_cycles, args.cycle)
    # setup dimensions:
    for out_file in files_list:
        make_dimensions(out_file, group_dimensions)
if __name__=="__main__":
    main()
