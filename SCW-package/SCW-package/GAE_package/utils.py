import pickle as pkl
import networkx as nx
import numpy as np
import scipy.sparse as sp
import torch
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import accuracy_score, precision_score
from sklearn.metrics import f1_score, recall_score

def split_edges(num_nodes, ip0):
	ip = ip0[:,0:2]
	npos = ip.shape[0]
	# training : valiation: testing = 7:2:1
	np_val, np_test = int(np.floor(npos / 20)), int(np.floor(npos / 10))
	# positive edges
	ipi_test = np.random.choice(range(npos), np_test, replace=False)
	ip2 = list(set(range(npos)) - set(ipi_test))
	ipi_val = np.random.choice(ip2, np_val, replace=False)
	ipi_train = list(set(range(npos)) - set(ipi_test) - set(ipi_val))
	ip_train = ip[ipi_train, :]
	ip_val, ip_test = ip[ipi_val, :], ip[ipi_test, :]

	nall = num_nodes*(num_nodes-1)/2
	nper = int(nall - npos)

	def ismember(a, b, tol=5):
		rows_close = np.all(np.round(a - b[:, None], tol) == 0, axis=-1)
		return np.any(rows_close)
	# randomly select negative edges
	n_neg = np_test + np_val
	if n_neg > nper:
		n_neg = nper
		np_test = int(n_neg/3)
		np_val = n_neg - np_test
	test_edges_false = []
	while len(test_edges_false) < n_neg:
		idx_i = np.random.randint(0, num_nodes)
		idx_j = np.random.randint(0, num_nodes)
		if idx_i >= idx_j:
			continue
		if ismember([idx_i, idx_j], ip):
			continue
		if test_edges_false:
			if ismember([idx_i, idx_j], np.array(test_edges_false)):
				continue
		test_edges_false.append([idx_i, idx_j])

	ineg = np.array(test_edges_false)
	ini_test = np.random.choice(range(n_neg), np_test, replace=False)
	ini_val = list(set(range(n_neg)) - set(ini_test))
	in_val, in_test = ineg[ini_val, :], ineg[ini_test, :]
	return ip_train, ip_val, ip_test, in_val, in_test

def split_edges_fast(num_nodes, ip0):
	ip = ip0[:,0:2]
	npos = ip.shape[0]
	# training : valiation: testing = 7:2:1
	np_val, np_test = int(np.floor(npos / 20)), int(np.floor(npos / 10))
	# positive edges
	ipi_test = np.random.choice(range(npos), np_test, replace=False)
	ip2 = list(set(range(npos)) - set(ipi_test))
	ipi_val = np.random.choice(ip2, np_val, replace=False)
	ipi_train = list(set(range(npos)) - set(ipi_test) - set(ipi_val))
	ip_train = ip[ipi_train, :]
	ip_val, ip_test = ip[ipi_val, :], ip[ipi_test, :]

	iall = np.asarray(np.triu_indices(num_nodes)).T
	idiag = iall[iall[:,0] == iall[:,1]]
	ineg = np.array(list(set(map(tuple, iall)) - set(map(tuple, ip)) - set(map(tuple, idiag))))
	
	nper = ineg.shape[0]
	# randomly select negative edges
	n_neg = np_test + np_val
	if n_neg > nper:
		n_neg = nper
		np_test = int(n_neg/3)
		np_val = n_neg - np_test

	ini_test = np.random.choice(range(n_neg), np_test, replace=False)
	ini_val = list(set(range(n_neg)) - set(ini_test))
	in_val, in_test = ineg[ini_val, :], ineg[ini_test, :]
	return ip_train, ip_val, ip_test, in_val, in_test

def parse_index_file(filename):
	index = []
	for line in open(filename):
		index.append(int(line.strip()))
	return index

def sparse_to_tuple(sparse_mx):
	if not sp.isspmatrix_coo(sparse_mx):
		sparse_mx = sparse_mx.tocoo()
	coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
	values = sparse_mx.data
	shape = sparse_mx.shape
	return coords, values, shape

def sparse_mx_to_torch_sparse_tensor(sparse_mx):
	"""Convert a scipy sparse matrix to a torch sparse tensor."""
	sparse_mx = sparse_mx.tocoo().astype(np.float32)
	indices = torch.from_numpy(np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
	values = torch.from_numpy(sparse_mx.data)
	shape = torch.Size(sparse_mx.shape)
	return torch.sparse.FloatTensor(indices, values, shape)

def preprocess_graph(adj):
	adj = sp.coo_matrix(adj)
	adj_ = adj + sp.eye(adj.shape[0])
	rowsum = np.array(adj_.sum(1))
	degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
	adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
	# return sparse_to_tuple(adj_normalized)
	return sparse_mx_to_torch_sparse_tensor(adj_normalized)

def get_roc_score(adj_rec, adj_orig, edges_pos, edges_neg, acc_cutoff=0.7):
	def sigmoid(x):
		return 1 / (1 + np.exp(-x))
	#adj_rec = np.dot(emb, emb.T)
	preds, pos, preds_neg, neg = [], [], [], []
	for e in edges_pos:
		#preds.append(sigmoid(adj_rec[e[0], e[1]]))
		preds.append(adj_rec[e[0], e[1]])
		pos.append(adj_orig[e[0], e[1]])
	for e in edges_neg:
		#preds_neg.append(sigmoid(adj_rec[e[0], e[1]]))
		preds_neg.append(adj_rec[e[0], e[1]])
		neg.append(adj_orig[e[0], e[1]])
	preds_all = np.hstack([preds, preds_neg])
	labels_all = np.hstack([np.ones(len(preds)), np.zeros(len(preds_neg))]) 
	# ROC-AUC
	roc_score = roc_auc_score(labels_all, preds_all)
	ap_score = average_precision_score(labels_all, preds_all)
	# traditional
	y_true = labels_all
	y_pred = np.where(preds_all > acc_cutoff, 1, 0)
	acc = accuracy_score(y_true, y_pred)
	#precision, recall, f1 = precision_score(y_true, y_pred), recall_score(y_true, y_pred), f1_score(y_true, y_pred)
	#return roc_score, ap_score, acc, precision, recall, f1
	return roc_score, ap_score, acc


