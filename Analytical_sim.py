from Cit_par import *
from EigFuncVer import *



print("Dimension-less for short period are:" , short_period(muc , KY2 , CZa , Cmadot, Cmq, Cma, V0 , c ))
print()
print("Dimension-less for phugoid are:" , phugoid(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0 ,V0 , c))
print()
print("Dimension-less for Dutch roll are:" , dutch_roll(mub , KZ2 , Cnr , CYb , Cnb, V0 , b))
print ("-----------------------------------------------")
print("Dimension-having eigenvalues for short period are:" , short_period_d(muc , KY2 , CZa , Cmadot, Cmq, Cma , V0 , c))
print()
print("Dimension-having eigenvalues for phugoid are:" , phugoid_d(muc, CZa, Cmq, Cma,CXu , Cmu, CXa, CZu , CZ0, V0 , c))
print()
print("Dimension-having eigenvalues for Dutch roll are:" , dutch_roll_d(mub , KZ2 , Cnr , CYb , Cnb , V0 , b))




