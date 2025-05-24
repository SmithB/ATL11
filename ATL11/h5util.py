#! /usr/bin/env python3

from h5py import h5a, h5t, h5s
import h5py
import numpy as np


def str_to_h5(_s):
  """ Converts a unicode string (the Python3 default) to ASCII (bytes) """
  _bytes = _s.encode('ascii', 'ignore')
  return _bytes
# enddef -----------------------------------------------------------------------
 
def h5_to_str(_bytes):
  """ Converts ASCII bytes to a unicode string (the Python3 default) """
  _string = _bytes.decode('utf-8')
  return _string
# enddef -----------------------------------------------------------------------
 
def create_attribute(_id, _name, _dims, _value):
  """
  Writes a HDF5 string attribute, ASCII, NULLTERM
 
  _id should be something like dset.id
 
  _dims should be a list.  For a scalar, use an empty list []
 
  """
 
# Make sure we don't have a unicode name
  _name=str_to_h5(_name)

# This routine for string attributes
  _dtype = h5t.FORTRAN_S1
# Create a scalar space (if dims len=0); otherwise a simple space
  if len(_dims) == 0:
    _sid=h5s.create(h5s.SCALAR)
  elif len(_dims) == 1 and _dims[0] == 0 :
    _sid=h5s.create(h5s.SCALAR)
  else:
    _sid=h5s.create_simple(tuple(_dims))
# endif
 
# Create the memory & file datatypes. Adjust if datatype is string.
  _mdtype = _dtype.copy()
  _fdtype = _dtype.copy()
  _classtype = _dtype.get_class()
  if _classtype == h5t.STRING:
    if isinstance(_value, list):
      _strlen=0
      for _part in _value: _strlen=max(_strlen, len(_part))
    else:
      _strlen = len(_value)
#   endif
    if _strlen < 1: return None
    _mdtype.set_size(_strlen)
    _mdtype.set_strpad(h5t.STR_SPACEPAD)
    _fdtype.set_size(_strlen+1)
    _fdtype.set_strpad(h5t.STR_NULLTERM)
# endif
 
## Either add or replace the attribute
#  if h5a.exists(_id, _name):
#    _aid = h5a.open(_id, name=_name)
#  else:
#    _aid=h5a.create(_id, _name, _fdtype, _sid)
# endif
# Either add or replace the attribute
  if h5a.exists(_id, _name):
    _aid = h5a.delete(_id, name=_name)
# endif
  _aid=h5a.create(_id, _name, _fdtype, _sid)

  if _classtype == h5t.STRING:
    if isinstance(_value, list):
      _value = np.array(_value, dtype=np.bytes_)
    else:
      _value = np.array(str_to_h5(_value))
#   endif
  else:
    _pytype = _fdtype.dtype
    _value = np.array(_value, dtype=_pytype)
# endif
  _aid.write(_value)
  return _aid
# enddef -----------------------------------------------------------------------

def duplicate_group(f_in, f_out, ingroup):
  """
  Duplicates group from one h5 object to another
 
  f_in should the result of something like f = h5py.File(infile,'r')

  f_out should the result of something like f = h5py.File(infile,'r+')
 
  ingroup should be a string.
 
  """
# Check that ingroup exists in f_in
  if ingroup not in f_in:
      print("[duplicate_group] Error, requested input group not present")
      return
# Remove outgroup if it already exists
  if ingroup in f_out:
      del f_out[ingroup]
  f_out.require_group(ingroup)
  for dset in f_in[ingroup].values():
      f_out[ingroup].create_dataset_like(dset.name, dset, data=np.array(dset))
      for att in dset.attrs:
        if att == 'DIMENSION_SCALE': continue
        if att == 'CLASS': continue
        if att == 'NAME': continue
        if att == 'DIMENSION_LIST': continue
        if att == 'REFERENCE_LIST': continue
        f_out[ingroup][dset.name].attrs[att] = dset.attrs[att]
  return
