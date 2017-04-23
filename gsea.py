import numpy as np
from scipy.stats import ks_2samp
from collections import deque
import sys, os



class gsea(object):
	# will need the class object to handle and store all that rickety data
	def __init__(self):
		self.seed_val = 20170422
		self.data_expression_profile = {}
		# store if the patient is AML or ALL. This can operate later a key to the lists
		# upon construction, it will look like: "data_expression_profile['patient'] = ["AML", "ALL", "ALL", "ALL", "AML"]
		# upon construction, it will look like: "data_expression_profile['STAT1'] = [-732, 34, 55, 123, -5]"
		self.data_expression_profile['patient_type'] = []
		self.gene_count = 0

		# raw strings from the pathways.txt
		self.gene_set_pathways = []
		self.pathways_count = 0

		self.examined_genes = []
		# an iterable list of the expression profile keys
		self.gene_samples_found = 0

		# store all of the test results here
		# how will you actually store your findings? 
		# test_result['GENE'] = [ results_int, results_int, results_int, ...]
		# test_result['GENE'][0] will ALWAYS contains the actual score. [1...n] is null hypothesis testing 
		self.initial_test_results = []
		self.test_results = {}

	def error_handler(self, content):
		print ("Error: {}\n".format(content))
		exit()

	def init_expressions(self, source_file):
		print ("reading expressions from file {}... ".format(source_file))

		with open(source_file) as fd:
			for line in fd:
				#print ("Here's a line in the expressions")
				#print ("{}".format(line))
				self._format_init_expression(line)

		self.gene_count = len(self.data_expression_profile) - 1

	def _format_init_expression(self, line):
		expression_data = line.split()

		if 'gene/patient' in expression_data[0]:
			for patient_type in expression_data[1:]:
				self.data_expression_profile['patient_type'].append(patient_type)
		else:
			# this is a new gene marker
			self.data_expression_profile[expression_data[0]] = []
			for expression_value in expression_data[1:]:
				self.data_expression_profile[expression_data[0]].append(float(expression_value))

	# store each pathway profile as a string. This should handle the weird formatting, especially with entries like "Manually Curated" 
	def init_pathways_profile(self, source_file, search_term):
		print ("reading gene_expression_profile from file {}... ".format(source_file))
		with open(source_file) as fd:
			for line in fd:
				if search_term.lower() in line.lower():
					#print line
					#self.gene_set_pathways.append(line)
					self._refine_pathways_profile(line)

		self.pathways_count = len(self.gene_set_pathways)

	# modify the contents found through the search terms to get the list of genes
	def _refine_pathways_profile(self, gene_input_pathway):
		split_contents = deque(gene_input_pathway.split())

		while "Curated" not in split_contents[0]:
			split_contents.popleft()
		split_contents.popleft()

		self.gene_set_pathways = self.gene_set_pathways + list(split_contents)






	'''
	# locate union of data_expression_profile from gene_set_pathways // narrow the list
	def narrow_set(self, search_term):
		desired_gene_set = self.show_profile(search_term)
		for k in self.data_expression_profile.keys():
			if k in desired_gene_set:
				#print k
				#self.narrowed_gene_expressions.append(k)
				self.gene_samples_found += 1

		dummy = desired_gene_set.split("Manually Curated")
		print "HERES THE DUMMY INFO"
		print dummy[1].split()

		#print desired_gene_set
		#print self.narrowed_gene_expressions
	'''


	# run diference calculation
	def test_gene_set(self):
		# test set A and set B

		# iterate through each dict, using the narrowed_gene expressions shit
		# that AML/ ALL key array
		# the unique typs
		# the values within a given gene array

		# return a tuple {NAME, differential expression}
		for gene in self.data_expression_profile[1:]:
			self.initial_test_results = self.calculate_mean_differential_expression(self.data_expression_profile.get('patient_type'), self.data_expression_profile.get(gene))
			#self.test_results[gene] = (self._perform_gene_set_test(self.data_expression_profile.get(gene), self.data_expression_profile.get('patient_type')))



		for d in self.initial_test_results:
			print d

		# self.data_expression_profile['patient_type']
		# self.narrowed_gene_expressions


	'''
	# we'll need the patient types, and their values	
	def _perform_gene_set_test(self, gene_set, patient_key):
		# build the sets

		unique_patient_types = list(set(patient_key))

		gene_dict = {}
		for pos in range(len(patient_key)):
			if patient_key[pos] in gene_dict:
				gene_dict[patient_key[pos]].append(gene_set[pos])
			else:
				gene_dict[patient_key[pos]] = []
				gene_dict[patient_key[pos]].append(gene_set[pos])


		np.random.seed(self.seed_val)
		n1 = len(gene_dict[unique_patient_types[0]])
		n2 = len(gene_dict[unique_patient_types[1]])


		return ks_2samp(gene_dict.get(unique_patient_types[0]), gene_dict.get(unique_patient_types[1]))
	'''

	#Calculate the mean differential expression for each gene between the two groups.

	def calculate_mean_differential_expression(self, patient_key, gene_set):
	#The first step is to rank the genes on the basis of their level of 
	#differential expression. Here, we simply averaged the expression of the 
	#two samples and then ranked them on the basis of the difference in expression across samples. 
		unique_patient_types = list(set(patient_key))
		group_a = []
		group_b = []

		# sort into two groups based on the patient_key
		for pos in range(len(patient_key)):
			if patient_key[pos] is unique_patient_types[0]:
				group_a.append(gene_set[pos])
			else:
				group_b.append(gene_set[pos])

		de_a = _calculate_individual_expression(group_a)
		de_b = _calculate_individual_expression(group_b)

		return (de_a - de_b)

	def _calculate_individual_expression(self, gene_set):
		gene_set.sort()
		difference = (numpy.diff(gene_set, n=len(gene_set)))
		
		return (difference / len(gene_set))


	#
	#
	#
	#
	#
	#
	#
	#



	# testing functions to read the output of the dictionaries
	def show_genes(self):
		genes = []
		for g in self.data_expression_profile.keys():
			genes.append(g)

		print ("Genes: {}".format(genes))

	def show_profile(self, search_term):
		for p in self.gene_set_pathways:
			if search_term.lower() in p.lower():
				print p
				return p

	def show_test_results(self):
		for (k, v) in self.test_results.items():
			print (k, v)




def main():
	gs = gsea()

	if len(sys.argv) != 3:
		gs.error_handler("Incorrect input given. \n $python gsea.py { gene_expression } { gene_profile }")
	else:
		gs.init_expressions(sys.argv[1])
		gs.init_pathways_profile(sys.argv[2], "leukemia")

		#print("Patient_expression_profile: {}".format(gs.data_expression_profile['patient_type']))

		#gs.narrow_set("leukemia")
		gs.test_gene_set()
		#gs.show_test_results()



main()











































































































































