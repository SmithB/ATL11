
import os
os.environ['MKL_NUM_THREADS']="1"
os.environ['NUMEXPR_NUM_THREADS']="1"
os.environ['OMP_NUM_THREADS']="1"
os.environ['OPENBLAS_NUM_THREADS']="1"

import numpy as np
import ATL11
import glob
import sys, h5py
import matplotlib.pyplot as plt
import matplotlib
#from mpl_toolkits.basemap import Basemap
from PointDatabase import geo_index

#591 10 -F /Volumes/ice2/ben/scf/AA_06/001/cycle_02/ATL06_20190205041106_05910210_001_01.h5 -b -101. -76. -90. -74.5 -o test.h5 -G "/Volumes/ice2/ben/scf/AA_06/001/cycle*/index/GeoIndex.h5" 
#591 10 -F /Volumes/ice2/ben/scf/AA_06/001/cycle_02/ATL06_20190205041106_05910210_001_01.h5 -o test.h5 -G "/Volumes/ice2/ben/scf/AA_06/001/cycle*/index/GeoIndex.h5" 

def get_proj4(hemisphere):
    if hemisphere==-1:
        return'+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' 
    if hemisphere==1:
        return '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '

def main(argv):
    # account for a bug in argparse that misinterprets negative agruents
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    # command-line interface: run ATL06_to_ATL11 on a list of ATL06 files
    import argparse
    parser=argparse.ArgumentParser(description='generate an ATL11 file from a collection of ATL06 files.')
    parser.add_argument('rgt', type=int, help="reference ground track number")
    parser.add_argument('subproduct', type=int, help="ICESat-2 subproduct number (latltude band)")
    parser.add_argument('--directory','-d', default=os.getcwd(), help="directory in which to search for ATL06 files")
    parser.add_argument('--pair','-p', type=int, default=None, help="pair number to process (default is all three)")
    parser.add_argument('--Release','-R', type=int, default=2, help="Release number")
    parser.add_argument('--Version','-V', type=str, default='001')
    parser.add_argument('--cycles', '-c', type=int, nargs=2, default=[3, 4], help="first and last cycles")
    parser.add_argument('--GI_file_glob','-G', type=str, default=None, help="Glob (wildcard) string used to match geoindex files for crossing tracks")
    parser.add_argument('--out_dir','-o', default=None, required=True, help="Output directory")
    parser.add_argument('--first_point','-f', type=int, default=None, help="First reference point")
    parser.add_argument('--last_point','-l', type=int, default=None, help="Last reference point")
    parser.add_argument('--num_points','-N', type=int, default=None, help="Number of reference points to process")
    parser.add_argument('--Hemisphere','-H', type=int, default=-1)
    parser.add_argument('--bounds', '-b', type=float, nargs=4, default=None, help="latlon bounds: west, south, east, north")
    parser.add_argument('--test_plot', action='store_true', help="plots locations, elevations, and elevation differences between cycles")
    parser.add_argument('--verbose','-v', action='store_true')
    args=parser.parse_args()

    # output file format is ATL11_RgtSubprod_c1c2_rel_vVer.h5
    out_file="%s/ATL11_%04d%02d_%02d%02d_%02d_v%s.h5" %( \
            args.out_dir,args.rgt, args.subproduct, args.cycles[0], \
            args.cycles[1], args.Release, args.Version)
    if os.path.isfile(out_file):
        os.remove(out_file)

    if args.verbose:
        print(out_file)
    glob_str='%s/*ATL06*_*_%04d??%02d_*.h5' % (args.directory, args.rgt, args.subproduct)
    files=glob.glob(glob_str)

    print("found ATL06 files:" + str(files))

    if args.pair is None:
        pairs=[1, 2, 3]
    else:
        pairs=[args.pair]
    
    if args.GI_file_glob is not None:
        GI_files=glob.glob(args.GI_file_glob)
    else:
        GI_files=None   
    print("found GI files:"+str(GI_files))
    
    for pair in pairs:
        print('files in ',files)
        D6 = ATL11.read_ATL06_data(files, beam_pair=pair, cycles=args.cycles)
        if D6 is None:
            continue
        D6, ref_pt_numbers, ref_pt_x = ATL11.select_ATL06_data(D6, first_ref_pt=args.first_point, last_ref_pt=args.last_point, lonlat_bounds=args.bounds, num_ref_pts=args.num_points)
#        D6.get_xy(EPSG=3413)
#        plt.plot(D6.x, D6.y,'r.')
#        plt.show()
#        exit(-1)

        if D6 is None or len(ref_pt_numbers)==0: 
            continue
        D11=ATL11.data().from_ATL06(D6, ref_pt_numbers=ref_pt_numbers, ref_pt_x=ref_pt_x,\
                      cycles=args.cycles, beam_pair=pair, verbose=args.verbose, \
                      GI_files=GI_files, hemisphere=args.Hemisphere) # defined in ATL06_to_ATL11
        # fill cycle_number list in cycle_stats and corrected_h
        setattr(D11.cycle_stats,'cycle_number',list(range(args.cycles[0],args.cycles[1]+1)))
        setattr(D11.corrected_h,'cycle_number',list(range(args.cycles[0],args.cycles[1]+1)))
        
        if D11 is not None:
            D11.write_to_file(out_file)

    # create a geo index for the current file.  This gets saved in the '/index' group
    if os.path.isfile(out_file):
        GI=geo_index(SRS_proj4=get_proj4(args.Hemisphere), delta=[1.e4, 1.e4]).for_file(out_file, 'ATL11', dir_root=args.out_dir)
        GI.attrs['bin_root']=None

        # the 'file' attributes of the geo_index are of the form :pair1, :pair2, :pair3, which means that the 
        # data for each bin are to be read from the current file
        for file in ['file_0','file_1','file_2']:
            if file in GI.attrs:
                temp = ':'+GI.attrs[file].split(':')[1]
                GI.attrs[file] = temp
        GI.to_file(out_file)
    
    # copy METADATA group from ATL06. Make lineage/ group for each ATL06 file, where the ATL06 filenames and their unique metadata are saved. 
    if os.path.isfile(out_file):        
        g = h5py.File(out_file,'r+')
        for ii,infile in enumerate(sorted(files)):
            if os.path.isfile(infile):
                f = h5py.File(infile,'r')         
                if ii==0:
                    f.copy('METADATA',g) # get all METADATA groups except Lineage
                    if 'Lineage' in list(g['METADATA'].keys()):
                        del g['METADATA']['Lineage']
                    g['METADATA'].create_group('Lineage')
                # make ATL06 file group
                gf = g['METADATA']['Lineage'].create_group('ATL06-{:02d}'.format(ii+1))
                gf.attrs['fileName'] = os.path.basename(infile)
                # fill ATL06 file group with unique file metadata
                for fgrp in list(f['METADATA']['Lineage']):
                    f.copy('METADATA/Lineage/{}'.format(fgrp), g['METADATA']['Lineage']['ATL06-{:02d}'.format(ii+1)])

                f.close()
        g.close()
    print("ATL06_to_ATL11: done with "+out_file)
        
    if args.test_plot:
        print(out_file)
        D = ATL11.data().from_file(out_file, field_dict=None)
        cm = matplotlib.cm.get_cmap('jet')
        colorslist = ['blue','green','red','orange','purple','brown','pink','gray','olive','cyan','black','yellow']
        
        fig = plt.figure(1)
        im = plt.scatter(D.corrected_h.longitude, D.corrected_h.latitude, c=D.corrected_h.delta_time[:,0], s=35, marker='.', cmap=cm)
        im = plt.scatter(D.crossing_track_data.longitude, D.crossing_track_data.latitude, c=D.crossing_track_data.delta_time, s=35, marker='.', cmap=cm)
        plt.title('Delta times of ATL11 data')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        fig.colorbar(im)
        #m = Basemap(projection='stere',resolution='i',lat_ts=70, lat_0=90, lon_0=45,
        #            llcrnrlon=-56, llcrnrlat=58, urcrnrlon=11, urcrnrlat=80, rsphere=(6378137.0, 6356752.3142))
        #m.drawcoastlines()
        
        ref, xo, delta = D.get_xovers()
        
        fig = plt.figure(2)
        for ii in range(len(ref.h_corr[:])):
            im = plt.errorbar(ref.x_atc,ref.h_corr[ii],ref.h_corr_sigma[ii],fmt='.',capsize=4,color=colorslist[ii]) 
        im = plt.errorbar(xo.x_atc,xo.h_corr[:],xo.h_corr_sigma[:],fmt='.',capsize=4,color=colorslist[ii+1])
        plt.title('Corrected Heights: cyc3(b), cyc4(g), crossing(r)')
        plt.xlabel('Along Track Distance [m]')
        plt.ylabel('Heights [m]')
        
        fig = plt.figure(3)
        for ii in range(len(ref.h_corr[:])-1):
            im = plt.scatter(ref.x_atc,ref.h_corr[ii+1]-ref.h_corr[ii],c=colorslist[ii],marker='.') 
        plt.title('Difference in Corrected Heights btn sequential cycles: later minus earlier')
        plt.xlabel('Along Track Distance [m]')
        plt.ylabel('Heights [m]')
        plt.grid()
        
        fig = plt.figure(4)
        print
        for ii in range(len(ref.h_corr[:])):
            im = plt.scatter(ref.x_atc,ref.h_corr[ii].ravel()-xo.h_corr[:],c=colorslist[ii],marker='.') 
        plt.title('Diff Corrected Heights: cyc3-xo(b), cyc4-xo(g)')
        plt.xlabel('Along Track Distance [m]')
        plt.ylabel('Heights [m]')
        plt.grid()
        print('type Control-C after viewing figures, to continue')
        plt.show()

if __name__=="__main__":
    main(sys.argv)
