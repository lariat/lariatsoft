# stash_spilltrees.py
#
# Given a file like MergedAtStartLinesim_LAriaT_13degProdxn_10degAna_SurveyedGeom_5000jobsof60k_64GeV_pos100Amps.root:
# - Add a unique (timestamp) string to the name, and 
# - Create the JSON metadata file for it, with some human input.
# - Copy both to the dropbox for spilltree files.
###################################################################################################################

import optparse
import json
import datetime as dt
import os
import subprocess


parser = optparse.OptionParser("usage: %prog [options] <input file.ROOT>\n")
parser.add_option ('-v', dest='debug', action="store_true", default=False,
                   help="Turn on verbose debugging.")
options, args = parser.parse_args()
debug = options.debug
if not len(args) >= 1: exit('Please supply a filename.')
infilename = args[0]
outdir = '/pnfs/lariat/scratch/fts/dropbox/'

## create_date
## data_tier generated-beam-only
## file_format: 'root'
## file_name [[Must be unique ==> include datetime]]
## secondary.intensity 300,000
## secondary.momentum
## tertiary.magnet_current 100
## tertiary.magnet_polarity Negative

metadump = {} # Dictionary to hold all the values before we write JSON
# Hard-coded for merged trees
metadump['data_tier'] = 'generated-beam-only'
metadump['file_format'] = 'root'

# Time we run this script. Needed for create_date and filename
## e.g. 2016-06-04T20:21:27+00:00. 
now_str = dt.datetime.now().strftime('%Y-%m-%dT%H:%M:%S+00:00') # Hardcoding UTC. Ugly but simple.
metadump['create_date'] = now_str

now_clean = now_str.replace('+00:00','').replace(':','').replace('-','')
# Making a unique-for-all-time file name:
newfilename = infilename.replace('.root','__'+now_clean+'.root')
metadump['file_name'] = newfilename

# Pions per target simulated intensity:
intensity = raw_input('Pions on target per spill? (Default: 300e3)')
if intensity == '': intensity = int(300e3)
else: intensity = int(float(intensity))
metadump['secondary.intensity'] = intensity

# Extract the bit about the secondary beam momentum on the target
parts = infilename.split('_')
momentum_str = '' # Establish scope
for part in parts:
    if part.count('GeV') == 1:
        momentum_str = part.rsplit('GeV')[0]
        break
metadump['secondary.momentum'] = int(momentum_str)

# Extract the bit about the magnet current
current_str  = '' # Establish scope
for part in parts:
    if part.count('Amps') == 1:
        current_str = part.rsplit('Amps')[0]
        break
magnet_polarity = '' # Establish scope
if current_str.count('neg') == 1: magnet_polarity = 'Negative'
elif current_str.count('pos') == 1: magnet_polarity = 'Positive'
else: exit('Unable to find a single neg/pos in '+current_str)
metadump['tertiary.magnet_polarity'] = magnet_polarity

magnet_current = current_str.replace('neg','').replace('pos','')
metadump['tertiary.magnet_current'] = int(magnet_current)

jsonfilename = newfilename+'.json'
cmd_json = 'ifdh cp '+jsonfilename+' '+outdir+jsonfilename+'\n'
cmd_root = 'ifdh cp '+infilename+' '+outdir+newfilename+'\n' # Change file name
print '\n  Move commands will be these:\n',cmd_json,cmd_root

print json.dumps(metadump, sort_keys=True, indent=4)
response = raw_input('Everything look ok? Enter y or script exits without doing anything.')
if response != 'y': exit('Exiting.')

jsonfile = open(newfilename+'.json', 'w')
jsonfile.write(json.dumps(metadump, sort_keys=True, indent=4) + '\n')
jsonfile.close()

tscriptname = 'tinymvscript.sh'
# Clear it away if file exists already
if os.path.isfile(tscriptname): subprocess.check_output('rm '+tscriptname, shell=True)
tscript = open(tscriptname,'w')
tscript.write('source /grid/fermiapp/lariat/setup_lariat.sh\n')
tscript.write('setup lariatsoft v06_05_00 -q e10:prof\n')
tscript.write(cmd_json)
tscript.write(cmd_root)
tscript.write('\n')
tscript.close()

if debug: subprocess.call('cat '+tscriptname, shell=True)
sourceme = 'source '+tscriptname
# Move file and metadata file to dropbox
out = subprocess.check_output(sourceme, shell=True)
