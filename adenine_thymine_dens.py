from pyscf import gto, scf, dft
import sys
from pyscf.tools import cubegen
import density_functional_approximation_dm21 as dm21
from pylibnxc.pyscf import UKS

mol = gto.Mole()
mol.atom = """
C            0.79798578248978        0.27519057470284       -0.01890509695066
C            0.70702747312535       -1.06575607739564       -0.01327835807200
N           -0.50035980543377       -1.70730936525323        0.06704320769620
C           -1.69568922338892       -1.03103345715763        0.14618410761476
N           -1.59915072204173        0.31979902406806        0.14076456708071
C           -0.43384804095291        1.06625807754303        0.06401971185273
O           -2.75362515578970       -1.63531483072713        0.21529615221125
C            2.08908976790964        1.01960035707866       -0.10541973075319
O           -0.47822793311336        2.27762741012063        0.06829630572069
H            1.57627045408281       -1.70496020473060       -0.07152939155485
H           -0.55050831848261       -2.71559440368892        0.06999153163645
H           -2.50968320821782        0.87257126191031        0.20082453496904
H            2.93282644091015        0.33737998664090       -0.15972749421850
H            2.19829952583287        1.66643057215990        0.76341252118527
H            2.08323641963996        1.66285560893889       -0.98372747616555
N           -7.17714822826316        3.98156165640549        0.48384323632164
C           -5.95643680340869        3.38069640978421        0.41225779243010
C           -6.21836587989208        2.00091690901628        0.42425064620424
N           -7.57164658925112        1.78664179542793        0.50161649555365
C           -8.09733331733828        2.97045848850909        0.53445411326732
C           -5.10170291295142        1.14164068588014        0.35787665821786
N           -3.89336892284171        1.72204544690860        0.28915771001715
C           -3.77953075533591        3.04671488063477        0.28551402597380
N           -4.74411706466496        3.93731710283537        0.34350216370909
H           -2.76736953251811        3.42383151702245        0.22726326863966
N           -5.19140570895963       -0.18743467682696        0.36007719484137
H           -7.34835812456523        4.97266827085858        0.49578953243705
H           -9.15126094504415        3.15822050439344        0.59538233129986
H           -6.10345688775085       -0.60425355030555        0.41244310182460
H           -4.35452178378429       -0.76653997475585        0.31127663700995
"""
mol.basis = 'def2-tzvp'
mol.charge = 1
mol.spin = 1
#mol.verbose = 4
mol.build()

#func = sys.argv[1]
func = 'GGA_HM'  #or the functional which you want to have

if func in ['DM21', 'DM21m', 'DM21mc', 'DM21mu']:

    mf = dft.UKS(mol, xc='B3LYP').density_fit()
    mf.run()
    dm0 = mf.make_rdm1()
    if func == 'DM21':
        mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
    elif func == 'DM21m':
        mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21m)
    elif func == 'DM21mc':
        mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
    elif func == 'DM21mu':
        mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mu)
    mf.grids.level = 3
    # mf.damp = 0.5  #uncomment if you need damping for SCF convergence
    # mf.init_guess = 'huckel'
    mf.conv_tol = 1E-5
    mf.conv_tol_grad = 1E-3
    # mf.verbose = 4
    mf.max_cycle = 250
    mf.kernel(dm0=dm0)
    rho_up, rho_down = mf.make_rdm1()

    # calculate total spindensity
    total_spin_density = rho_up - rho_down
    cubegen.density(mol, f"./resultsDFT/density_cubes/def2tzvp/{func}-ionized_SPIN_def2tzvp.cube", total_spin_density)


if func in ['GGA_XC_PBE', 'MGGA_XC_SCAN', 'GGA_HM', 'MGGA_HM']:
    mf = UKS(mol, nxc=func, nxc_kind='grid')
    mf.grids.level = 3
    mf.damp = 0.5  # uncomment if you need damping for SCF convergence
    mf.level_shift = 0.5
    # mf.init_guess = 'huckel'
    mf.conv_tol = 1E-5
    mf.conv_tol_grad = 1E-3
    mf.max_cycle = 250
    mf.kernel()
    rho_up, rho_down = mf.make_rdm1()

    # calculate total spindensity
    total_spin_density = rho_up - rho_down

    cubegen.density(mol, f"./resultsDFT/density_cubes/def2tzvp/{func}-ionized_SPIN_def2tzvp.cube", total_spin_density)

else:
    mf = dft.UKS(mol).density_fit()
    mf.xc = func
    mf.grids.level = 3
    # mf.damp = 0.5  #uncomment if you need damping for SCF convergence
    mf.conv_tol = 1E-5
    mf.conv_tol_grad = 1E-3
    # mf.level_shift = 0.5  #tested levelshift if HOMO-LUMO gap is small
    mf.verbose = 4
    mf.max_cycle = 250
    mf.kernel()
    rho_up, rho_down = mf.make_rdm1()

    # calculate total spindensity
    total_spin_density = rho_up - rho_down
    cubegen.density(mol, f"./resultsDFT/density_cubes/def2tzvp/{func}-ionized_SPIN_def2tzvp.cube", total_spin_density)

















































