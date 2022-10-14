from function import *

# This is a test for the value of Rs
# The expected result for i=55 deg, e=0, w=0 deg is Rs = 6371.009 km.

r = Rs(i=55, e=0, w=0)
error_rel = (r - 6371.009)/6371.009
print(error_rel)

# This is a test for the value of a
# The expected result for Ap = 500 km, e = 0.1, Rs = 6371.009 km is a = 7634.454444 km
a = semimayor_axis(Ap = 500, e = 0.1, Rs = 6371.009)
a_err_rel = (a - 7634.454444)/7634.454444*100
print("The porcentual relative error  obtained for the value of the semi-major axis is" , a_err_rel, "%")

# This is a test for the value of MM
# The expected result for a = 7634.4544 km is MM =  13.0147405 rev/day
MM = MeanMotion(a = 7634.4544)
MM_err_rel = (MM - 13.0147405)/13.0147405*100
print("The porcentual relative error  obtained for the value of the mean motion is" , MM_err_rel, "%")