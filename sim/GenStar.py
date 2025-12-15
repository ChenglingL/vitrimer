import sys
sys.path.append('/home/sciarella/H-B/buildCompass')
from hoomd import *
from hoomd import md
from hoomd import group
from hoomd import deprecated

context.initialize('--mode=cpu')

import math

length_box=float(sys.argv[2])
NAstar=600
NBstar=300
star_size=25
# compute the box volume
varSysVol = length_box * length_box * length_box
totrun=1000



#define the star structure
varBondStruct = []
varBondStruct.append((0, 1, 'polymerCha'))      # 0 is the central node
varBondStruct.append((0, 4, 'polymerCha'))
varBondStruct.append((0, 7, 'polymerCha'))
varBondStruct.append((0, 10, 'polymerCha'))
varBondStruct.append((0, 13, 'polymerCha'))
varBondStruct.append((0, 16, 'polymerCha'))
varBondStruct.append((0, 19, 'polymerCha'))
varBondStruct.append((0, 22, 'polymerCha'))

varBondStruct.append((1, 2, 'polymerCha'))
varBondStruct.append((2, 3, 'polymerCha'))

varBondStruct.append((4, 5, 'polymerCha'))
varBondStruct.append((5, 6, 'polymerCha'))

varBondStruct.append((7, 8, 'polymerCha'))
varBondStruct.append((8, 9, 'polymerCha'))

varBondStruct.append((11, 12, 'polymerCha'))
varBondStruct.append((10, 11, 'polymerCha'))

varBondStruct.append((13, 14, 'polymerCha'))
varBondStruct.append((14, 15, 'polymerCha'))

varBondStruct.append((16, 17, 'polymerCha'))
varBondStruct.append((17, 18, 'polymerCha'))

varBondStruct.append((19, 20, 'polymerCha'))
varBondStruct.append((20, 21, 'polymerCha'))

varBondStruct.append((22, 23, 'polymerCha'))
varBondStruct.append((23, 24, 'polymerCha'))


polymerStarA = dict(bond_len=0.5, type=(['C'] + 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['A']+ 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['B']), bond=varBondStruct, count=NAstar)
polymerStarB = dict(bond_len=0.5, type=(['C'] + 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['B']+ 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['A']),  bond=varBondStruct, count=NBstar)




# generate the polymer system
system = deprecated.init.create_random_polymers(box=data.boxdim(volume=varSysVol), polymers=[polymerStarA,polymerStarB], separation=dict(C=0.2,B=0.2, N=0.2, A=0.2),seed=int(sys.argv[3]))


nl = md.nlist.cell()
all = group.all()

#3 body potential between A and B
potC = md.pair.revcross(r_cut=1.8,nlist=nl)
potC.pair_coeff.set(['N','A','B','C'],['N','A','B','C'],sigma=0,n=0,epsilon=0,lambda3=0)
potC.pair_coeff.set('A','B',sigma=0.5,n=10,epsilon=100,lambda3=1)

#the respulsion is  WAC  (lj with cut and shift in the minimum)   (supposed diam=0.5)      
sigmalj=0.9
my_rcut=math.pow(2,1/6)*sigmalj    
potRep = md.pair.lj(r_cut=my_rcut,nlist=nl)
potRep.pair_coeff.set(['N','A','B','C'],['N','A','B','C'],epsilon=1.,sigma=sigmalj)
potRep.pair_coeff.set('A','B',epsilon=0.0,sigma=0.01)
potRep.set_params(mode="shift")



harm=md.bond.harmonic()
harm.bond_coeff.set('polymerCha', k=1000.0, r0=1)   


NAME='energyT%s_L%s_%s.log'%(sys.argv[1],sys.argv[2],sys.argv[3])
og1 = analyze.log(filename=NAME, quantities=['time', 'potential_energy', 'temperature','bond_harmonic_energy','pair_lj_energy'], period=500, overwrite=True)

T=float(sys.argv[1])

#### DUMP FOR VMD
#deprecated.dump.xml(group=all, filename="dump.xml", position=True, velocity=True, type=True, period=None, restart=False)
#dump.dcd(filename="movie.dcd",period=1000)

##FIRST RUN WITH LIMIT TO EQUILIBRATE
md.integrate.mode_standard(dt=0.0001)
a=md.integrate.nve(group=all, limit=0.001)
run(int(totrun))

md.integrate.nve.disable(a)

# output the snapshot
#deprecated.dump.xml(group=group.all(), filename="_Init.xml", vis=True)
#NAME='%scl/restart%s.gsd'%(sys.argv[2],sys.argv[1])
NAME='restartT%s_L%s_%s.gsd'%(sys.argv[1],sys.argv[2],sys.argv[3])
dump.gsd(filename=NAME, group=group.all(), truncate=True, period=5000, phase=0)


##temperature quench
md.integrate.mode_standard(dt=0.0001)
b=md.integrate.nvt(group=all, kT=T, tau=0.1) 
run(int(totrun)*1000)
md.integrate.nvt.disable(b)

##SECOND RUN NVE WITHOUT LIMIT
md.integrate.mode_standard(dt=0.001)
md.integrate.nve(group=all)

run(int(totrun)*1000)

# output the snapshot
#deprecated.dump.xml(group=group.all(), filename="_Init.xml", vis=True)
NAME='restartT%s_L%s_%s.gsd'%(sys.argv[1],sys.argv[2],sys.argv[3])
dump.gsd(filename=NAME, group=group.all(), truncate=True, period=None, phase=0)



