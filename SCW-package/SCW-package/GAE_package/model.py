import torch
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module

class ConvGraph(Module):
	""" Simple GCN layer: A*W """
	def __init__(self, in_features, out_features, dropout=0., act=F.relu):
		super(ConvGraph, self).__init__()
		self.in_features = in_features
		self.out_features = out_features
		self.dropout = dropout
		self.act = act
		self.weight = Parameter(torch.FloatTensor(in_features, out_features))
		self.reset_parameters()

	def reset_parameters(self):
		torch.nn.init.xavier_uniform_(self.weight)
	
	def forward(self, adj):
		#x = F.dropout(input, self.dropout, training=self.training)
		#x = torch.mm(x, self.weight)
		#x = torch.spmm(adj, x)
		x = torch.spmm(adj, self.weight)
		#x = self.act(x)
		return x

class InnerProductDecoder(Module):
	""" Decoder = Z*(Z^T) """
	def __init__(self, dropout=0., act=torch.sigmoid):
		super(InnerProductDecoder, self).__init__()
		self.dropout = dropout
		self.act = act

	def forward(self, z):
		z = F.dropout(z, self.dropout, training=self.training)
		adj_predicted = self.act(torch.mm(z, z.t()))
		return adj_predicted

class DistMultDecoder(Module):
	""" Decoder = Z*R*(Z^T) R is a diagonal matrix """
	def __init__(self, size, dropout=0., act=torch.sigmoid):
		super(DistMultDecoder, self).__init__()
		self.dropout = dropout
		self.act = act
		self.weight = Parameter(torch.randn(size))

	def forward(self, z):
		z = F.dropout(z, self.dropout, training=self.training)
		tmp = torch.diag(self.weight)
		z = torch.mm(z, tmp)
		adj_predicted = self.act(torch.mm(z, z.t()))
		return adj_predicted

class MatDecoder(Module):
	""" Decoder = Z*M*(Z^T) M is a n-by-n matrix """
	def __init__(self, size, dropout=0., act=torch.sigmoid):
		super(MatDecoder, self).__init__()
		self.dropout = dropout
		self.act = act
		self.weight = Parameter(torch.diag(torch.randn(size)))

	def forward(self, z):
		z = F.dropout(z, self.dropout, training=self.training)
		z = torch.mm(z, self.weight)
		adj_predicted = self.act(torch.mm(z, z.t()))
		return adj_predicted

class DEDICOMDecoder(Module):
	""" Decoder = Z*R*M*R*(Z^T) M is a n-by-n matrix, R is a diagonal matrix """
	def __init__(self, size, dropout=0., act=torch.sigmoid):
		super(DEDICOMDecoder, self).__init__()
		self.dropout = dropout
		self.act = act
		self.weightR = Parameter(torch.randn(size))
		self.weightM = Parameter(torch.diag(torch.randn(size)))

	def forward(self, z):
		z = F.dropout(z, self.dropout, training=self.training)
		tmp = torch.diag(self.weightR)
		z = torch.mm(z, tmp)
		z = torch.mm(z, self.weightM)
		z = torch.mm(z, tmp)
		adj_predicted = self.act(torch.mm(z, z.t()))
		return adj_predicted

class Encoder(Module):
	"""   """
	def __init__(self, in_features, out1_features, dropout):
		super(Encoder, self).__init__()
		self.layer1 = ConvGraph(in_features, out1_features, dropout, act=F.relu)
		#self.layer2 = ConvGraph(out1_features, out2_features, dropout, act=lambda x: x)
	
	def forward(self, adj):
		#x = input
		x = self.layer1(adj)
		#x = self.layer2(x, adj)
		return x

class GAE(Module):
	""" Graph auto-encoder """
	def __init__(self, in_features, out1_features, dropout):
		super(GAE, self).__init__()
		self.layer1 = ConvGraph(in_features, out1_features, dropout, act=F.relu)
		self.decoder = InnerProductDecoder(dropout, act=lambda x: x)
		#self.decoder_DistMult = DistMultDecoder(out2_features, dropout, act=lambda x: x)
		#self.decoder_Mat = MatDecoder(out2_features, dropout, act=lambda x: x)
		#self.decoder_DEDICOM = DEDICOMDecoder(out2_features, dropout, act=lambda x: x)
	
	def forward(self, adj):
		#x = input
		# encoder
		x = self.layer1(adj)
		# decoder
		x = self.decoder(x)
		#x = self.decoder_DistMult(x)
		#x = self.decoder_Mat(x)
		#x = self.decoder_DEDICOM(x)
		return x
		

class VGAE(Module):
	""" Variational graph auto-encoder """
	def __init__(self, in_features, out1_features, out2_features, dropout):
		super(VGAE, self).__init__()
		self.layer1 = ConvGraph(in_features, out1_features, dropout, act=F.relu)
		self.layer2 = ConvGraph(out1_features, out2_features, dropout, act=lambda x: x)
		self.layer3 = ConvGraph(out1_features, out2_features, dropout, act=lambda x: x)
		self.decoder = InnerProductDecoder(dropout, act=lambda x: x)
		self.decoder_DistMult = DistMultDecoder(out2_features, dropout, act=lambda x: x)
		self.decoder_Mat = MatDecoder(out2_features, dropout, act=lambda x: x)
		self.decoder_DEDICOM = DEDICOMDecoder(out2_features, dropout, act=lambda x: x)

	def reparameterize(self, zmean, zlogstd):
		if self.training:
			zstd = torch.exp(zlogstd)
			return torch.randn_like(zstd).mul(zstd).add_(zmean)
		else:
			return zmean

	def forward(self, input, adj):
		x = input
		# encoder
		x = self.layer1(x, adj)
		z_mean = self.layer2(x, adj)
		z_log_std = self.layer3(x, adj)
		z = self.reparameterize(z_mean, z_log_std)
		# decoder
		#x = self.decoder(z)
		#x = self.decoder_DistMult(z)
		#x = self.decoder_Mat(z)
		x = self.decoder_DEDICOM(z)
		return x, z_mean, z_log_std


