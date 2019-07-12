import numpy as np
import l1_procs as l1

x=np.array([1,2,3,0],dtype="double")
UT=np.array([[1,3,0,0.1],[-1,2,-1,0],[-1,-1,0,0]],dtype="double")
v0=np.zeros((3))
l1.BR_solve(x,UT,v0)
print 'Solving linear robust regression'
print 'The correct solution is [ 5. -3. 7.]'
print 'Calculated value (transposed case):' 
print v0
UN=np.array([[1,-1,-1],[3,2,-1],[0,-1,0],[0.1,0,0]],dtype="double")
v1=np.zeros((3))
l1.BR_solve(x,UN,v1)
print 'Calculated value (normal case):' 
print v1
v2=np.zeros((3))
l1.BR_solve(x,UN.T,v2)
print 'Calculated value (normal case that is transposed directly before calling):' 
print v2
U4N=np.zeros((4,2,3,2));
U4N[:,1,:,1]=UN[:,:]
v3=np.zeros((3))
l1.BR_solve(x,U4N[:,1,:,1],v3)
print 'Calculated value (normal case, not continuous):'
print v3 
xf=np.array([1,2,3,0],dtype="float32")
UNf=np.array([[1,-1,-1],[3,2,-1],[0,-1,0],[0.1,0,0]],dtype="float32")
vf=np.zeros((3),dtype="float32")
l1.BR_solve(xf,UNf,vf)
print 'Calculated value (normal case float32):'
print vf 
xi=np.array([1,2,3,0],dtype="int")
UTi=np.array([[1,3,0,0.1],[-1,2,-1,0],[-1,-1,0,0]],dtype="int")
vi=np.zeros((3),dtype="int")
print 'Processing Integer arrays. Not allowed. Exception should be raised:'
try:
  l1.BR_solve(xi,UTi,vi)
except Exception, e:
    print 'Exception catched: '+ str(e)  
else:
    print 'Failed to catch an exception'
print 'Passing different types. Not allowed. Exception should be raised:'
try:
  l1.BR_solve(x,UN,vf)
except Exception, e:
    print 'Exception catched: '+ str(e)  
else:
    print 'Failed to catch an exception'
print 'Passing matrix instead of vector. Not allowed. Exception should be raised:'
m=np.empty((2,2))
try:
  l1.BR_solve(x,UN,m)
except Exception, e:
    print 'Exception catched: '+ str(e)  
else:
    print 'Failed to catch an exception'
print 'Passing too long output vector. Not allowed. Exception should be raised:'
v5=np.empty(5)
try:
  l1.BR_solve(x,UN,v5)
except Exception, e:
    print 'Exception catched: '+ str(e)  
else:
    print 'Failed to catch an exception'
print 'Passing arrays of incompatible sizes. Not allowed. Exception should be raised:'
v2=np.empty(2)
try:
  l1.BR_solve(x,UN,v2)
except Exception, e:
    print 'Exception catched: '+ str(e)  
else:
    print 'Failed to catch an exception'

