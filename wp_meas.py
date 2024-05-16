import os
TYPE1=raw_input('SFH or sSFR: ')
bins=int(input('12-bin=0 19-bin=1:  '))
data=raw_input('early or late? ')
choice=int(input('Data-data=0, data-random=1, random-random=2: '))

if data=='early':
	TYPE='vespaBC_1_gals_SFH6HSsig_dr7_sp2'
else:
	TYPE='LINksmhigh22_{}zcut_parent'.format(TYPE1)

TYPE1='dr7'

if bins==0:
	if choice==0:
		cmd ='./xi_rppi_tree {0}_rppiCALC_12tree_{1}.dat {0}_rppiCALC_12tree_{1}.dat {0}_rppiGALind_12tree_{1}.dat {0}_rppiGALind_12tree_{1}.dat rppiCELLnei_{1}_12tree_1Rp.dat -1.0 1.4 0.0 40 12 20 250 250 250 {0}DD_rppi_12bin_tree.dat &'.format(TYPE,TYPE1)
		os.system(cmd) 
	if choice==1:
        	cmd ='./xi_rppi_tree {0}_rppiCALC_12tree_{1}.dat random_{0}_rppiCALC_12tree_{1}.dat {0}_rppiGALind_12tree_{1}.dat random_{0}_rppiGALind_12tree_{1}.dat rppiCELLnei_{1}_12tree_1Rp.dat -1.0 1.4 0.0 40.0 12 20 250 250 250 {0}DR_rppi_12bin_tree.dat &'.format(TYPE,TYPE1)
        	os.system(cmd)
	if choice==2:
        	cmd ='./xi_rppi_tree random_{0}_rppiCALC_12tree_{1}.dat random_{0}_rppiCALC_12tree_{1}.dat random_{0}_rppiGALind_12tree_{1}.dat random_{0}_rppiGALind_12tree_{1}.dat rppiCELLnei_{1}_12tree_1Rp.dat -1.0 1.4 0 40 12 20 250 250 250 {0}RR_rppi_12bin_tree.dat &'.format(TYPE,TYPE1)
        	os.system(cmd)
if bins==1:
	if choice==0:
        	cmd ='./xi_rppi_tree {0}_rppiCALC_30tree_{1}.dat {0}_rppiCALC_30tree_{1}.dat {0}_rppiGALind_30tree_{1}.dat {0}_rppiGALind_30tree_{1}.dat rppiCELLnei_{1}_30tree_1Rp.dat 0 60 0 60 30 30 250 250 250 {0}DD_rppi_30bin_tree_{1}.dat &'.format(TYPE,TYPE1)
        	os.system(cmd)
	if choice==1:
        	cmd ='./xi_rppi_tree {0}_rppiCALC_30tree_{1}.dat random_{0}_rppiCALC_30tree_{1}.dat {0}_rppiGALind_30tree_{1}.dat random_{0}_rppiGALind_30tree_{1}.dat rppiCELLnei_{1}_30tree_1Rp.dat 0 60 0 60 30 30 250 250 250 {0}DR_rppi_30bin_tree_{1}.dat &'.format(TYPE,TYPE1)
        	os.system(cmd)
	if choice==2:
        	cmd ='./xi_rppi_tree random_{0}_rppiCALC_30tree_{1}.dat random_{0}_rppiCALC_30tree_{1}.dat random_{0}_rppiGALind_30tree_{1}.dat random_{0}_rppiGALind_30tree_{1}.dat rppiCELLnei_{1}_30tree_1Rp.dat 0 60 0 60 30 30 250 250 250 {0}RR_rppi_30bin_tree_{1}.dat &'.format(TYPE,TYPE1)
        	os.system(cmd)
                    
