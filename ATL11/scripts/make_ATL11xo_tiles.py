
import numpy as np
import matplotlib.pyplot as plt
import pointCollection as pc
import glob
import re
import os
import argparse

# Temporary: setup output names for reference-track and crossing-track groups (should not be necessary in production
out_names = {'reference_track_data':'datum_track',
             'crossing_track_data':'crossing_track'}

parser=argparse.ArgumentParser(description='Generate a set of ATL11xo tiles from a directory of ATL11_atxo along-track crossover files')
parser.add_argument('--top_dir', type=str, required=True, help='top directory containing ATL11_atxo files')
parser.add_argument('--dest_dir', type=str, required=False, help='output directory for ATL11xo files')
parser.add_argument('--EPSG', type=int, required=False, help='EPSG string for output')
parser.add_argument('--release', type=int, required=True, help='ATL11 release number')
parser.add_argument('--version', type=int, required=True, help='ATL11xo version number')
parser.add_argument('--cycle', type=int, required=True, help='cycle number')
parser.add_argument('--region', type=str, required=True, help='region for output, AA=Antarctic, AR=Arctic')
parser.add_argument('--tile_spacing', type=float, default=200000, help='tile spacing,  m')
args=parser.parse_args()

if args.EPSG is None:
    if args.region=='AA':
        args.EPSG=3031
    else:
        args.EPSG=3413

if args.dest_dir is None:
    args.dest_dir = args.top_dir
tile_out_dir = os.path.join(args.dest_dir, f'cycle_{args.cycle:02d}')
if not os.path.isdir(tile_out_dir):
    os.mkdir(tile_out_dir)
tS = pc.tilingSchema(mapping_function_name='floor', tile_spacing=args.tile_spacing, EPSG=args.EPSG,
                    format_str = f'ATL11xo_{args.region}_E%d_N%d_c{args.cycle:02d}_{args.release:03d}_{args.version:02d}')
schema_file = os.path.join(tile_out_dir, f'{int(args.tile_spacing/1000)}km_tiling.json')
if not os.path.isfile(schema_file):
    tS.to_json(schema_file)
tS.directory=tile_out_dir

replace=True
for group in ['reference_track_data','crossing_track_data']:
    D=[]
    for file in glob.glob(args.top_dir+f'/ATL11_atxo*_*_{args.cycle:02d}_*.h5'):
        for pair in [1, 2, 3]:
            try:
                D += [pc.data().from_h5(file, group=f'pt{pair}/{group}')]
            except Exception:
                pass
    D=pc.data().from_list(D).get_xy(args.EPSG)
    bins, bin_dict = tS.tile_xy(data=D, return_dict=True)
    for xyT, ii in bin_dict.items():
        out_file = tS.tile_filename(xyT)
        if os.path.isfile(out_file) and group=='reference_track_data' :
            os.remove(out_file)
        D[ii].to_h5(out_file, group=out_names[group], replace=replace)
        replace=False
