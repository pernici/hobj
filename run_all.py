
import os

pyt = 'python'
#pyt = 'python3'
os.system('cd tests; %s test_hobj.py' % pyt)
os.system('cd tests; %s test_active_nodes.py' % pyt)
os.system('cd bench; %s bench.py' % pyt)
os.system('cd examples; %s d40nano.py 4' % pyt)
os.system('cd examples; %s ms_sqnp.py 12 12' % pyt)
os.system('cd examples; %s rand_reg6.py 6' % pyt)
os.system('cd examples; %s rnp_gen.py 4' % pyt)

