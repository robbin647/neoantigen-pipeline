import numpy as np


def Gaussian_ll(subdata):
	mu = np.mean(subdata)
	sigma = np.std(subdata)

	if sigma == 0:
		ll = -1e7
	else:
		ll = 0
		for i in subdata:
			ll += -1 / 2 * np.log(2 * np.pi) - np.log(sigma) - (i - mu) ** 2 / (2 * (sigma ** 2))
	print(ll)
	return ll


class GaussianMixtureDP:
	def __init__(self, n_components):
		self.K = n_components

	def Gaussian_ll(self, subdata):
		mu = np.mean(subdata)
		sigma = np.std(subdata)

		if sigma == 0:
			ll = -1e7
		else:
			ll = 0
			for i in subdata:
				ll += -1 / 2 * np.log(2 * np.pi) - np.log(sigma) - (i - mu) ** 2 / (2 * (sigma ** 2))

		return ll

	def fit(self, data):
		self.data = np.sort(data)

		dp_table = np.zeros((len(self.data), self.K)) - 1e7
		dp_index_table = np.zeros((len(self.data), self.K)) - 1

		for i in range(len(self.data)):
			dp_table[i, 0] = self.Gaussian_ll(self.data[:i+1])

		for j in range(1, self.K):
			for i in range(j, len(self.data)):
				ll_list = []
				for n in range(0, i):
					ll_list.append(dp_table[n, j-1] + self.Gaussian_ll(self.data[n+1:i+1]))

				dp_table[i, j] = max(ll_list)
				dp_index_table[i, j] = ll_list.index(dp_table[i, j])

		idx = len(self.data) - 1
		breakpoints_list = []
		breakpoints_list.append(len(self.data)-1)
		for j in range(self.K-1, 0, -1):
			breakpoints_list.append(int(dp_index_table[int(idx)][int(j)]))
			idx = breakpoints_list[-1]

		breakpoints_list.reverse()

		self.weights_ = []
		self.means_ = []
		self.variances_ = []
		left_bp = 0
		for right_bp in breakpoints_list:
			data_c = self.data[left_bp: right_bp+1]
			self.weights_.append(len(data_c)/len(self.data))
			self.means_.append(np.mean(data_c))
			self.variances_.append(np.var(data_c))
			left_bp = right_bp + 1

		self.weights_ = np.array(self.weights_)
		self.means_ = np.array(self.means_)
		self.variances_ = np.array(self.variances_)

		# print(self.data)
		# print(breakpoints_list)
		# print(dp_table)
		# print(dp_index_table)


