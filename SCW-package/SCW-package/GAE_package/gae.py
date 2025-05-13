from __future__ import division
from __future__ import print_function
import os
import argparse
import numpy as np
import scipy.sparse as sp
import torch
from torch import optim
import torch.nn.modules.loss
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau

from model import GAE
from utils import split_edges, preprocess_graph, get_roc_score

parser = argparse.ArgumentParser()

parser.add_argument('--nodes', type=int, default=0, help='Number of nodes.')
parser.add_argument('--fin', type=str, metavar='FILE',help='input edge file')
parser.add_argument('--fout', type=str, metavar='FILE',help='output reconstructed adjacency file')
parser.add_argument('--epochs', type=int, default=800, help='Number of epochs to train.')
parser.add_argument('--hidden1', type=int, default=128, help='Number of units in hidden layer 1.')
parser.add_argument('--lr', type=float, default=0.001, help='Initial learning rate.')
parser.add_argument('--dropout', type=float, default=0., help='Dropout rate (1 - keep probability).')
args = parser.parse_args()


def loss_function_GAE(preds, labels, pos_weight, norm):
	cost = norm * F.binary_cross_entropy_with_logits(preds, labels, pos_weight=pos_weight)
	return cost

def loss_function_VGAE(preds, labels, mu, logvar, n_nodes, norm, pos_weight):
	cost = norm * F.binary_cross_entropy_with_logits(preds, labels, pos_weight=pos_weight)
	# see Appendix B from VAE paper:
	# Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
	# https://arxiv.org/abs/1312.6114
	# 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
	KLD = -0.5 / n_nodes * torch.mean(torch.sum(1 + 2 * logvar - mu.pow(2) - logvar.exp().pow(2), 1))
	return cost + KLD

def gae(args):

	n_nodes = args.nodes
	edges = np.loadtxt(args.fin, dtype='int')
	fout = args.fout
	num_epochs = args.epochs
	lr = args.lr
	ndim = args.hidden1
	droprate = args.dropout

	train_edges, val_edges, test_edges, val_edges_false, test_edges_false = split_edges(n_nodes, edges)
	#print("Done. split edges.")
	# rebuild adj matrix
	adj = sp.csr_matrix((np.ones(edges.shape[0]), (edges[:, 0], edges[:, 1])), shape=(n_nodes, n_nodes))
	# make it symmetric and remove dupliates of diagonal
	adj = adj + adj.T- sp.diags(adj.diagonal())
	# store original adj matrix for evaluation
	adj_orig = adj
	adj_orig.eliminate_zeros()
	# rebuild adj matrix for training
	# TODO: change weights for training edges based on Hi-C contacts, not just set all to 1s
	weight_train_edges = np.ones(train_edges.shape[0])
	adj_train = sp.csr_matrix((weight_train_edges, (train_edges[:, 0], train_edges[:, 1])), shape=(n_nodes, n_nodes))
	adj_train = adj_train + adj_train.T- sp.diags(adj_train.diagonal())
	adj = adj_train
	# normalize training network
	adj_norm = preprocess_graph(adj)
	#adj_label = adj_train + sp.eye(adj_train.shape[0])
	adj_label = adj_train
	adj_label = torch.FloatTensor(adj_label.toarray())
	# some weights for loss function
	pos_weight = float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()
	pos_weight = torch.from_numpy(np.array(pos_weight)).float()
	norm = adj.shape[0] * adj.shape[0] / float((adj.shape[0] * adj.shape[0] - adj.sum()) * 2)

	# Graph auto-encoder model
	model = GAE(n_nodes, ndim, droprate)
	optimizer = optim.Adam(model.parameters(), lr=lr)

	best_hidden_emb, best_model = None, None
	best_val_ap, best_epoch = -1, -1
	#print("Start training.")
	for epoch in range(num_epochs):
		# training
		model.train()
		optimizer.zero_grad()
		recovered = model(adj_norm)
		loss = loss_function_GAE(preds=recovered, labels=adj_label, pos_weight=pos_weight, norm=norm)
		loss.backward()
		optimizer.step()
		# validation
		model.eval()
		hidden_emb_tmp = recovered.data.numpy()
		hidden_emb_tmp = 1. / (1. + np.exp(-hidden_emb_tmp))
		roc_curr, ap_curr,val_acc = get_roc_score(hidden_emb_tmp, adj_orig, val_edges, val_edges_false)
		if ap_curr > best_val_ap:
			best_val_ap = ap_curr
			best_epoch = epoch
			best_hidden_emb = hidden_emb_tmp
			#best_model = model.state_dict()
	# testing
	#roc, ap, acc, pre, recall, f1 = get_roc_score(best_hidden_emb, adj_orig, test_edges, test_edges_false)
	roc, ap, acc = get_roc_score(best_hidden_emb, adj_orig, test_edges, test_edges_false)
	#return best_hidden_emb
	#return best_hidden_emb, roc, ap, acc
	print('AUC: ', roc, ' AP: ', ap)
	np.savetxt(fout, best_hidden_emb, fmt='%1.5f')
	#return best_hidden_emd
if __name__ == '__main__':
	gae(args)

