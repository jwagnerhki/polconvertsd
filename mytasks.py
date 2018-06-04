#
# User defined tasks setup.
# Generated from buildmytask.
#

import sys
from casa_stack_manip import stack_frame_find

if sys.path[1] != '/home/marti/WORKAREA/ARC_TOOLS/PolConvert/LaunchPad/trunk':
  sys.path.insert(1, '/home/marti/WORKAREA/ARC_TOOLS/PolConvert/LaunchPad/trunk')
from odict import odict
if not globals().has_key('mytasks') :
  mytasks = odict()

mytasks['polconvert'] = '\n\nVersion 1.7.3\n\nConverts VLBI visibilities from mixed-polarization basis (i.e., linear-to-circular) into circular basis. Works with single VLBI stations as well as with phased arrays (i.e., phased ALMA).\n\n'

if not globals().has_key('task_location') :
  task_location = odict()

task_location['polconvert'] = '/home/marti/WORKAREA/ARC_TOOLS/PolConvert/LaunchPad/trunk'
myglobals = stack_frame_find( )
tasksum = myglobals['tasksum'] 
for key in mytasks.keys() :
  tasksum[key] = mytasks[key]

from polconvert_cli import  polconvert_cli as polconvert
