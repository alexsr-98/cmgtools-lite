import sys, math
print ('argument list', sys.argv)
n1 = float(sys.argv[1])
n2 = float(sys.argv[2])
print ("n1  = {}".format(n1))
print ("n2  = {}".format(n2))

err = math.sqrt(n1**2 + n2**2)
print ("err = {:.2f}".format(err))
