from gae import gae
import argparse
import numpy as np
import scipy.sparse as sp
parser = argparse.ArgumentParser()
parser.add_argument('--avg', type=int, default=5, help='averge times')
parser.add_argument('--nodes', type=int, default=0, help='Number of nodes.')
parser.add_argument('--fin', type=str, metavar='FILE',help='input edge file')
parser.add_argument('--fout', type=str, metavar='FILE',help='output reconstructed adjacency file')
parser.add_argument('--epochs', type=int, default=400, help='Number of epochs to train.')
parser.add_argument('--hidden1', type=int, default=128, help='Number of units in hidden layer 1.')
parser.add_argument('--lr', type=float, default=0.001, help='Initial learning rate.')
parser.add_argument('--dropout', type=float, default=0., help='Dropout rate (1 - keep probability).')
args = parser.parse_args()
nod = args.nodes
t = args.avg
print(t)
new_data = np.zeros([nod,nod],dtype = float, order = 'F')
for x in range(t):
	data1 = gae(args)
	new_data += data1
new_data = new_data/t
fout = args.fout
np.savetxt(fout, new_data, fmt='%1.5f')
#print(new_data)

