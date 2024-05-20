import MDAnalysis as mda
import sys
import numpy as np

rsw=1.25    # switching function parameters for heavy atoms (currently in A)
dsw=0.04    # switching function parameters for heavy atoms (currently in A)
sigma=0.025 # parameter for switching functionfor CV (currently in A)

# if no arguments passed, print help
if len(sys.argv) == 1:
    print('usage: python pt_setup.py ref.pdb --initial (resid) --target (resid)')


    print('\n\t--initial\tinitial site residue ID') 
    print('\t--target\ttarget site residue ID')    
    print('\t--trajectory\ttrajectory file .xtc or .pdb (optional)')
    print('\t--debug\t\tprint out debug information (optional)')
    print('\t--write\t\twrite pdb and ndx files for frames where CV > argument (optional)\n')
    exit()

# raise error if reference structure is not a pdb file or a gro file
if sys.argv[1][-4:] != '.pdb':
    print('reference structure must be a pdb file, i havent tested this with gro files')
    exit()
# check that user has specified initial and target sites
if '--initial' not in sys.argv or '--target' not in sys.argv:
    print('please specify initial and target sites using flags --initial and --target')
    exit()
# flag for debugging
if '--debug' in sys.argv:
    debug = True
else:
    debug = False
if '--write' in sys.argv:
    write = True
else:
    write = False

# get arguments (reference structure, initial site, target site)
reference_structure = sys.argv[1]
initial_site_resid = sys.argv[sys.argv.index('--initial')+1]
target_site_resid = sys.argv[sys.argv.index('--target')+1]
# otherwise, set trajectory to reference structure
if '--trajectory' not in sys.argv:
    trajectory = reference_structure
else :
    trajectory = sys.argv[sys.argv.index('--trajectory')+1]

# md analysis universe object
u = mda.Universe(reference_structure, trajectory)

# run below for each frame
for ts in u.trajectory:

    QMheavy = u.select_atoms('((resid %s or resid %s) and protein and (name O*) and not backbone)' % (initial_site_resid, target_site_resid)) \
    + u.select_atoms('(name O* and (resname TIP3 or resname SOL) and cyzone 8 8 -8 ((resid %s or resid %s) and protein))' % (initial_site_resid, target_site_resid)) 
    #+ u.select_atoms('((resid 281 or resid 280 or resid 134 or resid 335 or resid 173) and protein and (name O* or name N*) and not backbone) ', updating=True)

    QMhydrogen = u.select_atoms('(resid %s or resid %s) and protein and (name HD* or name HE*) and not backbone' % (initial_site_resid, target_site_resid)) \
    + u.select_atoms('(byres (name O* and (resname TIP3 or resname SOL) and cyzone 8 8 -8 ((resid %s or resid %s) and protein))) and name H*' % (initial_site_resid, target_site_resid)) 
    #+ u.select_atoms('((resid 281 or resid 280 or resid 134 or resid 335 or resid 173) and protein and (name HZ* or name HH) and not backbone)', updating=True)

    # QMall
    QMall = QMheavy + QMhydrogen

    # select relevant (protonable) heavy atoms in initial and target sites
    for site in [initial_site_resid, target_site_resid]:
        if len(u.select_atoms('resid %s' % site)) == 0:
            print('residue %s is not present in reference structure' % site)
            exit()
            # can also be ASPH or GLUH
        if u.select_atoms('resid %s' % site).residues[0].resname in ('ASP', 'ASPH', 'ASH'):
            cv_atoms = u.select_atoms('resid %s and name OD*' % site)
        elif u.select_atoms('resid %s' % site).residues[0].resname in ('GLU', 'GLUH'):
            cv_atoms = u.select_atoms('resid %s and name OE*' % site)
        elif u.select_atoms('resid %s' % site).residues[0].resname in ('HIS', 'HSD', 'HSE', 'HSP'):
            cv_atoms = u.select_atoms('resid %s and name NE* or name ND*' % site)
        elif u.select_atoms('resid %s' % site).residues[0].resname in ('LYS', 'LYSH'):
            cv_atoms = u.select_atoms('resid %s and name NZ' % site)
        else:
            print('residue %s named' % site, u.select_atoms('resid %s' % site).residues[0].resname, 'is not a valid selection (e.g. ASP, GLU, HIS, LYS)')
            exit()
        if site == initial_site_resid:
            initial_cv_atoms = cv_atoms
        elif site == target_site_resid:
            target_cv_atoms = cv_atoms
        else:
            print('problem with initial and target site selection')
            exit()

    # raise error if one of the sites is not an ASP or GLU because i haven't implemented the other residues yet
    if initial_cv_atoms.residues[0].resname not in ('ASP', 'ASPH', 'GLU', 'GLUH', 'ASH'):
        print('initial site must be an ASP or GLU residue (currently only ASP and GLU are supported)')
        exit()
    if target_cv_atoms.residues[0].resname not in ('ASP', 'ASPH', 'GLU', 'GLUH', 'ASH'):
        print('target site must be an ASP or GLU residue (currently only ASP and GLU are supported)')
        exit()

    # heavy atom weight factor is equal to the coordination number of hydrogens in the deprotonated state - maybe the weight should be one for hte protonated oxygen in asp
    weight_factors = {}
    for atom in QMheavy:
        if atom.name == 'OW':                           # water
            weight_factors[atom.index] = 2.0
        elif atom.name == 'OD2' or atom.name == 'OE2':  # aspartate oxygens
            weight_factors[atom.index] = 0.0
        elif atom.name == 'OD1' or atom.name == 'OE1':  # glutamate oxygens
            weight_factors[atom.index] = 0.0
        elif atom.name == 'NZ':                         # lysine
            weight_factors[atom.index] = 3.0               
        elif atom.name == 'OH':                         # tyrosine
            weight_factors[atom.index] = 1.0
        elif atom.name == 'OH2':                        # water (alternate)
            weight_factors[atom.index] = 2.0
        # ligand specific
        elif atom.name == 'N' and atom.index == 42583:     # workaround for a ligand amide N
            weight_factors[atom.index] = 1.0
        elif atom.name == 'N' and atom.index == 42571:     # workaround for a ligand N-terminus
            weight_factors[atom.index] = 3.0
        elif atom.name == 'O' and atom.index == 42582:     # workaround for a ligand
            weight_factors[atom.index] = 0.0
        elif atom.name == 'OC1' or atom.name == 'OC2':  # workaround for a ligand
            weight_factors[atom.index] = 0.0

        else:
            print('atom name not recognized')
            print(atom.name)
            exit()
            
    # get all hydrogen positions
    QMhydrogen_positions = QMhydrogen.positions
    QMheavy_positions = QMheavy.positions

    # get heavy atom COG
    QMheavy_COG = np.mean(QMheavy_positions, axis=0)
    QMhydrogen_COG = np.mean(QMhydrogen_positions, axis=0)

    # do the above using nested loops
    QMhydrogen_to_QMheavy_vectors = np.zeros((len(QMhydrogen), len(QMheavy), 3))
    QMhydrogen_to_QMheavy_distances = np.zeros((len(QMhydrogen), len(QMheavy)))
    QMhydrogen_to_QMheavy_vectors_scaled = np.zeros((len(QMhydrogen), len(QMheavy), 3))
    for i in range(len(QMhydrogen)):
        for j in range(len(QMheavy)):
            QMhydrogen_to_QMheavy_vectors[i,j,:] = QMhydrogen_positions[i,:] - QMheavy_positions[j,:]
            QMhydrogen_to_QMheavy_distances[i,j] = np.linalg.norm(QMhydrogen_to_QMheavy_vectors[i,j,:])
            QMhydrogen_to_QMheavy_vectors_scaled[i,j,:] = QMhydrogen_to_QMheavy_vectors[i,j,:] * (1.0 / (1.0 + np.exp((QMhydrogen_to_QMheavy_distances[i,j] - rsw) / dsw)))

    # sum of all hydrogen vectors (A)
    QMhydrogen_sum_vector = np.sum(QMhydrogen.positions, axis=0)
    # do the above using a loop
    QMhydrogen_sum_vector = np.zeros(3)
    for atom in QMhydrogen:
        QMhydrogen_sum_vector += atom.position
    
    # weighted sum of all heavy atom vectors (B)
    QMheavy_weighted_sum_vector = np.zeros(3)
    for atom in QMheavy:
        QMheavy_weighted_sum_vector += weight_factors[atom.index] * atom.position

    # sum of all scaled pairwise vectors (C)
    QMhydrogen_to_QMheavy_vectors_scaled_sum = np.sum(QMhydrogen_to_QMheavy_vectors_scaled, axis=(0,1))

    # center of excess charge = A - B - C
    CEC = (QMhydrogen_sum_vector - QMheavy_weighted_sum_vector - QMhydrogen_to_QMheavy_vectors_scaled_sum) #+ QMheavy_COG

    # project the CEC onto the vector between the initial and target sites
    # vector between initial and target sites
    initial_to_target_vector = np.mean(target_cv_atoms.positions, axis=0) - np.mean(initial_cv_atoms.positions, axis=0)
    # normalise this vector
    initial_to_target_vector = initial_to_target_vector / np.linalg.norm(initial_to_target_vector)

    # distances between all atoms in initial site and the CEC
    initial_site_to_CEC_distances = np.linalg.norm(initial_cv_atoms.positions - CEC, axis=1)

    # distances between all atoms in target site and the CEC
    target_site_to_CEC_distances = np.linalg.norm(target_cv_atoms.positions - CEC, axis=1)

    # sort both these arrays from smallest to largest
    initial_site_to_CEC_distances_sorted = np.sort(initial_site_to_CEC_distances)
    target_site_to_CEC_distances_sorted = np.sort(target_site_to_CEC_distances)

    # update the identities of the atoms in the initial and target sites to be the closest two to the CEC
    initial_cv_atoms = initial_cv_atoms[np.argsort(initial_site_to_CEC_distances)]
    target_cv_atoms = target_cv_atoms[np.argsort(target_site_to_CEC_distances)]
    

    # weight for each distance is 1 / (exp((RD1 - RD2)/sigma) + 1))
    w1 = 1.0 / (1.0 + np.exp((initial_site_to_CEC_distances_sorted[0] - initial_site_to_CEC_distances_sorted[1]) / sigma))
    w2 = 1.0 / (1.0 + np.exp((initial_site_to_CEC_distances_sorted[1] - initial_site_to_CEC_distances_sorted[0]) / sigma))
    w3 = 1.0 / (1.0 + np.exp((target_site_to_CEC_distances_sorted[0] - target_site_to_CEC_distances_sorted[1]) / sigma))
    w4 = 1.0 / (1.0 + np.exp((target_site_to_CEC_distances_sorted[1] - target_site_to_CEC_distances_sorted[0]) / sigma))

    # first two colvar terms 
    RD = w1*initial_site_to_CEC_distances_sorted[0] + w2*initial_site_to_CEC_distances_sorted[1]
    RE = w3*target_site_to_CEC_distances_sorted[0] + w4*target_site_to_CEC_distances_sorted[1]

    RD1E1 = np.linalg.norm(initial_cv_atoms[0].position - target_cv_atoms[0].position)
    RD1E2 = np.linalg.norm(initial_cv_atoms[0].position - target_cv_atoms[1].position)
    RD2E1 = np.linalg.norm(initial_cv_atoms[1].position - target_cv_atoms[0].position)
    RD2E2 = np.linalg.norm(initial_cv_atoms[1].position - target_cv_atoms[1].position)

    # cross term
    RDE = (w1*w3*RD1E1) + (w1*w4*RD1E2) + (w2*w3*RD2E1) + (w2*w4*RD2E2)

    # final CV (projection)
    CV = ((RD * RD) + (RDE * RDE) - (RE * RE)) / (2.0 * (RDE * RDE))
 
    # if xtc was passed as argument, print out frame number, CV, time and number of atoms, else simpler printout
    if trajectory[-4:] == '.xtc':
        print('frame:', ts.frame, '\tCV: %.6f' % CV, '\ttime:', ts.time, '\tnumber of atoms:', len(QMheavy) + len(QMhydrogen), end='\r')
        # write CV to a xvg file
        with open('cv_value.xvg', 'a') as f:
            f.write('%s %.6f\n' % (ts.time, CV))
    else:
        print('\n\tCV: %.6f' % CV, '\tnumber of atoms:', len(QMheavy) + len(QMhydrogen), end='\n\r') 

    # if flag -w is passed, write pdb and ndx file for frames where CV > argument
    if write == True:
        if CV > float(sys.argv[sys.argv.index('--write')+1]):
           # write pdb file with all the QM CEC atoms
            with mda.Writer('frame-%s.pdb' % ts.frame, n_atoms=len(u.atoms)) as W:
                W.write(u.atoms)

            # write a gromacs index file with all the QM CEC atoms
            with open('frame-%s.ndx' % ts.frame, 'w') as f:
                f.write('[QM-CEC-atoms]\n')
                # must sort so they are ascendiung
                QMlist = np.sort(np.array(QMall.atoms.indices))
                for atom in QMlist:
                    f.write('%s ' % int(atom + 1))
                f.write('\n\n')

    ############################################################################################################################

    # if debug flag is set to true, print out everything below
    if debug == True :

        # check atoms
        print('\n')
        for atom in QMall:
            print('\tIndex:', atom.index, 'ID:', atom.index, ' resID:', atom.resid, 'resname', atom.resname)
        print('\n')

        # DEBUG PRINTOUTS
        #set print preferences to 3 dp
        np.set_printoptions(precision=3)
        print('\n\tinitial site COG:', np.mean(initial_cv_atoms.positions, axis=0))
        print('\ttarget site COG:', np.mean(target_cv_atoms.positions, axis=0))
        print('\n\tX COG:', QMheavy_COG)
        print('\tH COG:', QMhydrogen_COG)
        print('\n\tQM hydrogens:', len(QMhydrogen), 'QM heavy:', len(QMheavy))    
        print('\tmatrix shape:', QMhydrogen_to_QMheavy_vectors.shape)
        print('\tN cross terms:', len(QMhydrogen) * len(QMheavy))
        print('\n         H vector sum:', QMhydrogen_sum_vector)
        print('        wX vector sum:', QMheavy_weighted_sum_vector)
        print('           cross term:', QMhydrogen_to_QMheavy_vectors_scaled_sum)
        print('            CEC coord:', CEC)
        print('\n\tRD CALCULATION:', '(', '%.2f' % w1, '*', '%.2f' % initial_site_to_CEC_distances_sorted[0], ') + (', '%.2f' % w2,'*', '%.2f' % initial_site_to_CEC_distances_sorted[1], ')')
        print('\tRE CALCULATION:',   '(', '%.2f' % w3, '*', '%.2f' % target_site_to_CEC_distances_sorted[0],  ') + (', '%.2f' % w4,'*', '%.2f' % target_site_to_CEC_distances_sorted[1],  ')')
        print('\n\tRDE = (w1*w3*RD1E1) + (w1*w4*RD1E2) + (w2*w3*RD2E1) + (w2*w4*RD2E2)')
        print('\n\tw1: %.2f' % w1 ,'\n\tw2: %.2f' % w2 ,'\n\tw3: %.2f' % w3 ,'\n\tw4: %.2f' % w4)
        print('\n\tRD1E1: %.2f' % RD1E1 ,'\n\tRD1E2: %.2f' % RD1E2 ,'\n\tRD2E1: %.2f' % RD2E1 ,'\n\tRD2E2: %.2f' % RD2E2)
        print('\n\tRDE CALCULATION: (', '%.2f' % w1,'%.2f' % w3,'%.2f' % RD1E1, ') + (', '%.2f' % w1,'%.2f' % w4,'%.2f' % RD1E2, ') + (', '%.2f' % w2,'%.2f' % w3,'%.2f' % RD2E1, ') + (', '%.2f' % w2,'%.2f' % w4,'%.2f' % RD2E2, ')')
        print('\n\tRD:', '%.6f' % RD)
        print('\tRE:', '%.6f' % RE)
        print('\tRDE:', '%.6f' % RDE)
        print('\n')

        # DEBUGGING PLOT
        import matplotlib.pyplot as plt
        wf=[]
        for atom in QMheavy:
            wf.append(weight_factors[atom.index]+1)
        wf=np.array(wf)
        wf=wf/np.max(wf)
        ax = plt.axes(projection='3d')
        ax.scatter3D(QMheavy_positions[:,0], QMheavy_positions[:,1], QMheavy_positions[:,2], color='red', marker='o', s=wf*100 )
        # now plot the hydrogens
        ax.scatter3D(QMhydrogen_positions[:,0], QMhydrogen_positions[:,1], QMhydrogen_positions[:,2], color='grey', marker='o')
        # now plot the initial and target sites
        ax.scatter3D(initial_cv_atoms.positions[:,0], initial_cv_atoms.positions[:,1], initial_cv_atoms.positions[:,2], color='yellow', marker='o')
        ax.scatter3D(target_cv_atoms.positions[:,0], target_cv_atoms.positions[:,1], target_cv_atoms.positions[:,2], color='yellow', marker='o')
        # label the initial and target sites at the average position of the atoms
        ax.text3D(np.mean(initial_cv_atoms.positions[:,0]), np.mean(initial_cv_atoms.positions[:,1]), np.mean(initial_cv_atoms.positions[:,2]), 'initial site', color='black')
        ax.text3D(np.mean(target_cv_atoms.positions[:,0]), np.mean(target_cv_atoms.positions[:,1]), np.mean(target_cv_atoms.positions[:,2]), 'target site', color='black')
        # connect the hydrogens to the heavy atoms with lines, colour their alpha based on the switching function
        for i in range(len(QMhydrogen)):
            for j in range(len(QMheavy)):
                ax.plot3D([QMhydrogen_positions[i,0], QMheavy_positions[j,0]], [QMhydrogen_positions[i,1], QMheavy_positions[j,1]], [QMhydrogen_positions[i,2], QMheavy_positions[j,2]], color='grey', alpha=(1.0 / (1.0 + np.exp((QMhydrogen_to_QMheavy_distances[i,j] - rsw) / dsw))))
        ## a notate the heavy atoms with their weight factors
        for i in range(len(QMheavy)):
            ax.text3D(QMheavy_positions[i,0], QMheavy_positions[i,1]+0.2, QMheavy_positions[i,2]+0.2, '%.0f' % weight_factors[QMheavy[i].index], color='black')
        ## now plot the CEC
        ax.scatter3D(CEC[0], CEC[1], CEC[2], color='green', marker='o')
        # label the CEC
        ax.text3D(CEC[0], CEC[1], CEC[2], 'CEC', color='black')
        # draw lines twixt the CEC and the initial and target sites
        ax.plot3D([CEC[0], initial_cv_atoms.positions[0,0]], [CEC[1], initial_cv_atoms.positions[0,1]], [CEC[2], initial_cv_atoms.positions[0,2]], color='black', linestyle='--')
        # label this line "RD"
        ax.text3D((CEC[0] + initial_cv_atoms.positions[0,0]) / 2.0, (CEC[1] + initial_cv_atoms.positions[0,1]) / 2.0, (CEC[2] + initial_cv_atoms.positions[0,2]) / 2.0, 'RD', color='black')
        ax.plot3D([CEC[0], target_cv_atoms.positions[0,0]], [CEC[1], target_cv_atoms.positions[0,1]], [CEC[2], target_cv_atoms.positions[0,2]], color='black', linestyle='--')
        # label this line "RE"
        ax.text3D((CEC[0] + target_cv_atoms.positions[0,0]) / 2.0, (CEC[1] + target_cv_atoms.positions[0,1]) / 2.0, (CEC[2] + target_cv_atoms.positions[0,2]) / 2.0, 'RE', color='black')
        # draw a line between the initial and target sites
        ax.plot3D([initial_cv_atoms.positions[0,0], target_cv_atoms.positions[0,0]], [initial_cv_atoms.positions[0,1], target_cv_atoms.positions[0,1]], [initial_cv_atoms.positions[0,2], target_cv_atoms.positions[0,2]], color='black', linestyle='--')
        # label this line "RDE"
        ax.text3D((initial_cv_atoms.positions[0,0] + target_cv_atoms.positions[0,0]) / 2.0, (initial_cv_atoms.positions[0,1] + target_cv_atoms.positions[0,1]) / 2.0, (initial_cv_atoms.positions[0,2] + target_cv_atoms.positions[0,2]) / 2.0, 'RDE', color='black')
        # plot cecp
        #ax.scatter3D(CECp[0], CECp[1], CECp[2], color='blue', marker='o')
        plt.show()

