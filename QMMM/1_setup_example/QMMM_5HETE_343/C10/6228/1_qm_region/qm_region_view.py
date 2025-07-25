import chimera
from chimera import runCommand

runCommand('sel solvent & @EPW')
runCommand('~disp sel')
runCommand('sel #1 & protein & element.C')
runCommand('color pink sel')
runCommand('sel #1:LIG & element.C')
runCommand('color pink sel')
runCommand('sel #1 ')
runCommand('transp 85,a sel')
runCommand('transp 96,r sel')
runCommand('sel #0 & element.C')
runCommand('color #bfffbf sel')
runCommand('sel solvent')
runCommand('transp 70,a sel')
runCommand('disp #0:LIG')
runCommand('sel #1:LIG')
runCommand('transp 100,a sel')
runCommand('sel #1:LIG@c1,c2,c3,c4,c5,c15,c16,c17,c18,c19,c20,o1,o2,o3,o4,h31')
runCommand('transp 0,a sel')
runCommand('background solid white')
runCommand('save /Volumes/white_HDD/qmmm/QMMM_5HETE_343/2_calculations/6228/1_qm_region/view.py')