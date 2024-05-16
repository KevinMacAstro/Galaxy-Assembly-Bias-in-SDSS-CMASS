import os
TYPE1=raw_input('SFH or sSFR: ')
bins=int(input('12-bin=0 19-bin=1:  '))
jack_samp=int(input('Number of jackknife samples:  '))
data=raw_input('early or late? ')
choice=int(input('Data-data=0, data-random=1, random-random=2: '))

if data=='early':
	TYPE='vespaBC_1_gals_SFH6HSsig_dr7_sp2'
else:
	TYPE='LINksmhigh22_{}zcut_parent'.format(TYPE1)

TYPE1='dr7'

if bins==0:
	if choice==0:
		cmd ='./xi_rppi_tree_jack {0}_jack_rppiCALC_12tree_{2}.dat {0}_jack_rppiCALC_12tree_{2}.dat {0}_rppiGALind_12tree_{2}.dat {0}_rppiGALind_12tree_{2}.dat rppiCELLnei_{2}_12tree_1Rp.dat -1.0 1.4 0 40 12 20 250 250 250 {0}DD_rppi_12bin_jack_tree.dat {1} &'.format(TYPE,jack_samp,TYPE1)
		os.system(cmd) 
	if choice==1:
        	cmd ='./xi_rppi_tree_jack  {0}_jack_rppiCALC_12tree_{2}.dat random_{0}_jack_rppiCALC_12tree_{2}.dat {0}_rppiGALind_12tree_{2}.dat random_{0}_rppiGALind_12tree_{2}.dat rppiCELLnei_{2}_12tree_1Rp.dat  -1.0 1.4 0 40 12 20 250 250 250 {0}DR_rppi_12bin_jack_tree.dat {1} &'.format(TYPE,jack_samp,TYPE1)
        	os.system(cmd)
	if choice==2:
        	cmd ='./xi_rppi_tree_jack random_{0}_jack_rppiCALC_12tree_{2}.dat random_{0}_jack_rppiCALC_12tree_{2}.dat random_{0}_rppiGALind_12tree_{2}.dat random_{0}_rppiGALind_12tree_{2}.dat rppiCELLnei_{2}_12tree_1Rp.dat  -1.0 1.4 0 40 12 20 250 250 250 {0}RR_rppi_12bin_jack_tree.dat {1} &'.format(TYPE,jack_samp,TYPE1)
        	os.system(cmd)
if bins==1:
	if choice==0:
                #cmd ='./xi_rppi_pairs_jack {0}_jack_rppiCALC_19tree_SFH.dat {0}_jack_rppiCALC_19tree_SFH.dat -2. 1.8 0 60 19 30 250 250 250 {0}DD_rppi_19bin_jack_SFH1.dat {1} &'.format(data,jack_samp)

        	cmd ='./xi_rppi_tree_jack {0}_jack_rppiCALC_19tree_{2}.dat {0}_jack_rppiCALC_19tree_{2}.dat {0}_rppiGALind_19tree_{2}.dat {0}_rppiGALind_19tree_{2}.dat rppiCELLnei_{2}_19tree_1Rp.dat -2. 1.8 0 60 19 30 250 250 250 {0}DD_rppi_19bin_jack_tree.dat {1} &'.format(TYPE,jack_samp,TYPE1)
        	os.system(cmd)
	if choice==1:
        	cmd ='./xi_rppi_tree_jack {0}_jack_rppiCALC_19tree_{2}.dat random_{0}_jack_rppiCALC_19tree_{2}.dat {0}_rppiGALind_19tree_{2}.dat random_{0}_rppiGALind_19tree_{2}.dat rppiCELLnei_{2}_19tree_1Rp.dat -2. 1.8 0 60 19 30 250 250 250 {0}DR_rppi_19bin_jack_tree.dat {1} &'.format(TYPE,jack_samp,TYPE1)
        	os.system(cmd)
	if choice==2:
        	cmd ='./xi_rppi_tree_jack  random_{0}_jack_rppiCALC_19tree_{2}.dat random_{0}_jack_rppiCALC_19tree_{2}.dat random_{0}_rppiGALind_19tree_{2}.dat random_{0}_rppiGALind_19tree_{2}.dat rppiCELLnei_{2}_19tree_1Rp.dat -2. 1.8 0 60 19 30 250 250 250 {0}RR_rppi_19bin_jack_tree.dat {1} &'.format(TYPE,jack_samp,TYPE1)
        	os.system(cmd)
                    
