import os
import dadi
import pylab

NomPop1 = "India"
NomPop2 = "PR"
N_ind1 = 24
N_ind2 = 30


subfiles = []
for dirpath, subdirs, files in os.walk("/work/user/kdurand/dgimi/DADI/boot2024/"):
    for x in files:
        if x.endswith(".txt"):
        	dd = dadi.Misc.make_data_dict(x)
		fs = dadi.Spectrum.from_data_dict(dd, [NomPop1,NomPop2], [N_ind1, N_ind2], polarized = False)
		fs.to_file(x + ".boot2024.fs")
