import os
import sys
import copy
import time
import argparse
import numpy as np
import pandas as pd


class Simplex():
	def __init__(self, folder, verbose, show_dual):
		self.folder = folder
		self.verbose = verbose
		self.show_dual = show_dual
		self.pivot_cnt = 0
		self.B = None					# Base updated in real time
		self.dual_list = None			# Indices of dual, set at the start of phase two
		
		# Check if all necessary input files exist under the specified path
		mandatory_files = ["A.csv", "b.csv", "c.csv"]
		mandatory_file_paths = list(map(lambda f: os.path.join(folder, f), mandatory_files))
		if any(list(map(lambda x: not os.path.exists(x), mandatory_file_paths))):
			raise Exception("Lack of input files! Please check if 'A.csv', 'b.csv' and 'c.csv' exist under the specified path!")

		A_df = pd.read_csv(mandatory_file_paths[0], header=None)
		self.A = np.array(A_df, dtype=np.float64)

		b_df = pd.read_csv(mandatory_file_paths[1], header=None)
		self.b = np.array(b_df, dtype=np.float64).reshape(-1)
		self.cons_cnt = self.b.size

		c_df = pd.read_csv(mandatory_file_paths[2], header=None)
		self.c = np.array(c_df, dtype=np.float64).reshape(-1)
		self.var_cnt = self.c.size

	def solve(self):
		phase_one_out = self.firstPhaseSimplex(self.A, self.b)
		if phase_one_out == None:
			return
		
		phase_two_out = self.secondPhaseSimplex(*phase_one_out, self.c)
		if phase_two_out == None:
			return
		
		self.printFinal(*phase_two_out)

	# Build AUXILIARY problem for initial solution
	def firstPhaseSimplex(self, A, b):
		# Build the AUXILIARY matrices
		A_aux = np.hstack((A, np.eye(self.cons_cnt)))
		A_aux[b < 0] *= -1														# Modify A_aux so that the respective b is positive
		b_aux = np.where(b < 0, -b, b)											# b_aux should be positive
		self.B = B_aux = np.array([(i + self.var_cnt) for i in range(self.cons_cnt)])	# List of basis for the auxiliary problem
		c_aux = np.hstack((np.zeros(self.var_cnt), np.ones(self.cons_cnt)))		# Cost for the auxiliary problem

		# Calculate the (negative) objective value for the AUXILIARY problem
		obj_aux = -np.sum(b_aux)

		# Calculate reduced cost for the AUXILIARY problem
		r_aux = c_aux - np.dot(c_aux[B_aux].reshape(1, -1), A_aux)[0]

		if self.verbose:
			print("---------------------------------------------")
			print("- Starting phase ONE")
			self.printTableau(self.buildTableau(A_aux, b_aux, r_aux, obj_aux))
		
		# Choose the entering and leaving basis for the AUXILIARY problem, then PIVOT
		# Ignore small errors
		while (r_aux < -1e-12).any():
    		# Choose the smallest index of negative reduce cost to avoid degeneration
			enter_var_idx = np.where(r_aux < 0)[0][0]

			# Choose the smallest index of the minimal b/y ratio to avoid degeneration
			b_y_ratio = np.divide(b_aux, A_aux[:, enter_var_idx], out=np.full_like(b_aux, np.inf), where=A_aux[:, enter_var_idx]>0)
			b_y_ratio_min = b_y_ratio[np.argmin(b_y_ratio[np.where(b_y_ratio >= 0)])]
			leave_var_idx = np.where(b_y_ratio == b_y_ratio_min)[0][0]

			A_aux, B_aux, b_aux, r_aux, obj_aux = self.doPivot(A_aux, B_aux, b_aux, r_aux, obj_aux, enter_var_idx, leave_var_idx)
			
		# If the final objective value for the AUXILIARY problem is not zero,
		# then the primal problem is INFEASIBLE. Ignore small errors.
		if abs(obj_aux) >= 1e-12:
			print("Infeasible!")
			print("---------------------------------------------")
			return

		# Continue pivoting if the base contains AUXILIARY variable
		for i, basis_aux in enumerate(B_aux):
			if basis_aux >= self.var_cnt:
				enter_var_idx = np.where(A_aux[i] != 0)[0][0]
				if (A_aux[i, :self.var_cnt] == 0).all():
					B_aux = np.delete(B_aux, i, axis=0)
					A_aux = np.delete(A_aux, i, axis=0)
					b_aux = np.delete(b_aux, i, axis=0)
					self.cons_cnt -= 1
					continue
				A_aux, B_aux, b_aux, r_aux, obj_aux = self.doPivot(A_aux, B_aux, b_aux, r_aux, obj_aux, enter_var_idx, i)

		return (A_aux[:, 0:self.var_cnt], B_aux, b_aux)

	# Solve the PRIMAL problem with initial solution
	def secondPhaseSimplex(self, A, B, b, c):
		# Calculate the (negative) objective value for the PRIMAL problem
		obj = -np.dot(c[B], b)

		# Calculate reduced cost for the PRIMAL problem
		r = c - np.dot(c[B].reshape(1, -1), A)[0]

		# Set base for dual
		self.dual_list = copy.deepcopy(B)

		if self.verbose:
			print("---------------------------------------------")
			print("- Starting phase TWO")
			self.printTableau(self.buildTableau(A, b, r, obj))

		# Choose the entering basis and leaving basis for the PRIMAL problem, then PIVOT
		while (r < 0).any():
    		# Choose the smallest index of the smallest reduce cost to avoid degeneration
			enter_var_idx = np.where(r < 0)[0][0]
			
			if (A[:, enter_var_idx] <= 0).all():
				print("Unbounded!")
				print("---------------------------------------------")
				return
			
			# Choose the smallest index of the minimal b/y ratio to avoid degeneration
			b_y_ratio = np.divide(b, A[:, enter_var_idx], out=np.full_like(b, np.inf), where=A[:, enter_var_idx]>0)
			b_y_ratio_min = b_y_ratio[np.argmin(b_y_ratio[np.where(b_y_ratio >= 0)])]
			leave_var_idx = np.where(b_y_ratio == b_y_ratio_min)[0][0]

			A, B, b, r, obj = self.doPivot(A, B, b, r, obj, enter_var_idx, leave_var_idx)

		sol = np.zeros(self.var_cnt, dtype=np.float64)
		sol[B] = b
		return (sol, obj, r)

	# Pivoting two basis
	def doPivot(self, A, B, b, r, obj, enter_var_idx, leave_var_idx):
		self.pivot_cnt += 1
		pivot_elem = A[leave_var_idx, enter_var_idx]
		if self.verbose:
			print("-- %d pivot: x_%d enter, x_%d leave" % (self.pivot_cnt, enter_var_idx+1, B[leave_var_idx]+1))

		# Update the line for leaving variable
		A[leave_var_idx], b[leave_var_idx] = A[leave_var_idx] / pivot_elem, b[leave_var_idx] / pivot_elem

		# Update other lines
		for i in range(self.cons_cnt):
			if i != leave_var_idx:
				b[i] -= A[i, enter_var_idx] * b[leave_var_idx]
				A[i] -= A[i, enter_var_idx] * A[leave_var_idx]

		obj -= b[leave_var_idx] * r[enter_var_idx]	# Update objective value
		r -= r[enter_var_idx] * A[leave_var_idx]	# Update reduced cost
		B[leave_var_idx] = enter_var_idx			# Update the base
		self.B = B
		
		if self.verbose:
			self.printTableau(self.buildTableau(A, b, r, obj))
		return (A, B, b, r, obj)

	def buildTableau(self, A, b, r, obj):
		tableau = np.zeros((A.shape[0]+1, A.shape[1]+1), dtype=np.float64)
		tableau[:-1, :-1] = A
		tableau[-1, :-1] = r
		tableau[:-1, -1] = b
		tableau[-1, -1] = obj
		return tableau

	def printTableau(self, tableau):
		print("---------------------------------------------")
		
		# Print variable names and b
		print("{:<6}".format(""), end='')
		for j in range(tableau.shape[1]-1):
			print("x_{:<8}".format(j+1), end='')
		print("b\t")

		# Print variable values
		for i in range(tableau.shape[0]-1):
			print("x_{:<4}".format(self.B[i]+1), end='')
			for j in range(tableau.shape[1]):
				print("{:<10.4f}".format(tableau[i, j]), end='')
			print("")
		
		# Print reduced costs
		print("{:<6}".format("r"), end='')
		for j in range(tableau.shape[1]):
			print("{:<10.4f}".format(tableau[-1, j]), end='')
		print("")

		print("---------------------------------------------")
	
	def printFinal(self, x, obj, r):
		print("---------------------------------------------")
		print("Objective value: %.4f" % -obj)
		print("---------------------------------------------")
		print("Value of each variable:")
		for i, sol in enumerate(x):
			print("x_%d = %.4f" % (i+1, sol))
		if self.show_dual:
			print("---------------------------------------------")
			dual_vals = r[self.dual_list]
			print("Value of dual variable:")
			for i, d_val in enumerate(dual_vals):
				print("d_%d = %.4f" % (i+1, -d_val))
		print("---------------------------------------------")
		print("Total number of pivots: %d" % self.pivot_cnt)
		print("---------------------------------------------")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Commands for the Simplex solver")
	parser.add_argument("folder", type=str, default="data", help="The path of the folder that contains inputs")
	parser.add_argument("-s", "--step", help="Turn this on to see step-by-step solution. Turn this off to test running speed.", action="store_true")
	parser.add_argument("-d", "--dual", help="Turn this on to see dual values.", action="store_true")
	args = parser.parse_args()
	s = Simplex(folder=args.folder, verbose=args.step, show_dual=args.dual)

	time_start = time.time()
	s.solve()
	time_end = time.time()
	print("Run time:", time_end - time_start, "second")
