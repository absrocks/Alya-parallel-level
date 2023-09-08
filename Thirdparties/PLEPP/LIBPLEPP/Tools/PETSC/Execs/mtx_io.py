from scipy.io     import *
from scipy.sparse import *
from numpy        import * 
from sys          import argv 

A = mmread(argv[1])
print ".        A:", A.shape

A  = A.tocsr() 
b = ones(A.shape[0])
 

from scipy.sparse.linalg import spsolve, gmres, bicgstab 
from scipy.sparse.linalg import aslinearoperator, LinearOperator   
from scipy import * 

x2, info2 = bicgstab(A, b)
print ". bicgstab:", average(x2), x2.std(), where(~info2, "ok!", "error")
print ".        x:", x2 
print "ok!" 
