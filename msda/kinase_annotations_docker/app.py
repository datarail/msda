import subprocess


out_folder = '/output'
fasfile = '/input/substrate.fas'
resfile = '/input/psite.res'
outfile = '%s/Network_predictions.txt' % out_folder

nk = 'resources/NetworKIN_release3.0/NetworKIN.py'
np = 'resources/NetPhorest/netphorest'
bl = 'resources/blast-2.2.17/bin/blastall'

f = open(outfile, 'wb')
command = [nk, '-n', np, '-b', bl, '9606', fasfile, resfile]
subprocess.call(command, stdout=f)
