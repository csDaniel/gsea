import numpy as np
import math
import random
from collections import deque
import sys, os

class gsea(object):
	# will need the class object to handle and store all that rickety data
	def __init__(self):
		self.data_expression_profile = {}
		# store if the patient is AML or ALL. This can operate later a key to the lists
		# upon construction, it will look like: "data_expression_profile['patient'] = ["AML", "ALL", "ALL", "ALL", "AML"]
		# upon construction, it will look like: "data_expression_profile['STAT1'] = [-732, 34, 55, 123, -5]"
		self.data_expression_profile['patient_type'] = []
		self.gene_count = 0

		# raw strings from the pathways.txt
		self.gene_set_pathways = []
		self.pathways_count = 0
		self.initial_test_results = []
		self.test_results = {}
		self.enrichment_scores = []

	def error_handler(self, content):
		print ("Error: {}\n".format(content))
		exit()

	def init_expressions(self, source_file):
		print ("reading expressions from file {}... ".format(source_file))

		with open(source_file) as fd:
			for line in fd:
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
					self._refine_pathways_profile(line)

		self.pathways_count = len(self.gene_set_pathways)

	# modify the contents found through the search terms to get the list of genes
	def _refine_pathways_profile(self, gene_input_pathway):
		split_contents = deque(gene_input_pathway.split())

		# I was unsure how to 'clean up' the pathways.txt file. This was a quick work-around.
		while "Curated" not in split_contents[0]:
			split_contents.popleft()
		split_contents.popleft()

		self.gene_set_pathways = self.gene_set_pathways + list(split_contents)

	def test_gene_set_controller(self):
		base_type = self.data_expression_profile.get('patient_type')
		self.permutation_test_gene_set(base_type, 0)

		print "Testing Null Hypothesis. Generating Permuation Tests Now."
		random_patient_types = np.array(base_type)

		for n in range(1,3):
			np.random.seed(n)
			self.permutation_test_gene_set(np.random.permutation(random_patient_types), n)
			if n % 100 is 0:		
				print ("Testing {}% complete.".format(n/10))

		for h in self.enrichment_scores:
			print h

	# randomize the 'patient_type' key-set 
	def permutation_test_gene_set(self, patient_types, test_sequence):

		gene_results = []

		for gene in self.data_expression_profile.keys():
			if 'patient_type' not in gene:
				result =  self.calculate_mean_differential_expression(patient_types, self.data_expression_profile.get(gene))
				gene_results.append((gene,result))
		
		#self.initial_test_results.append(gene_results)

		#self.initial_test_results[test_sequence].sort(key= lambda x: x[1])
		#self.significance_of_difference(self.initial_test_results[test_sequence])

		gene_results.sort(key= lambda x: x[1])
		self.significance_of_difference(gene_results)

	# sub-sections of the test itself

	def significance_of_difference(self, test_sequence):
		score = 0
		max_score = 0
		max_location = 0


		up_step = math.sqrt ((self.gene_count - self.pathways_count) / self.pathways_count)
		down_step = math.sqrt (self.pathways_count / (self.gene_count - self.pathways_count))

		print ("upstep = {}".format(up_step))
		print ("downstep = {}".format(down_step))
		print self.pathways_count
		print self.gene_count

		for loc in range(self.gene_count):
			#print score
			if test_sequence[loc][0] in self.gene_set_pathways:
				score += up_step
			else:
				score -= down_step
			if score > max_score:
				max_score = score
				max_location = loc

		results = {}
		results['score'] = score
		results['es_static'] = max_location
		self.enrichment_scores.append(results)

	#Using the ordered list generated in Problem 2, calculate the running sum, adding the appropriate up- 
	#and down-steps as you move down the ranked list of genes. What is the supremum (enrichment score) of this running sum?


	def calculate_mean_differential_expression(self, patient_key, gene_set):
		#The first step is to rank the genes on the basis of their level of 
		#differential expression. Here, we simply averaged the expression of the 
		#two samples and then ranked them on the basis of the difference in expression across samples. 
		unique_patient_types = list(set(patient_key))
		group_a = []
		group_b = []

		# sort into two groups based on the patient_key
		for pos in range(len(patient_key)):
			if patient_key[pos] == unique_patient_types[0]:
				group_a.append(gene_set[pos])
			else:
				group_b.append(gene_set[pos])

		de_a = self._calculate_individual_expression(group_a)
		de_b = self._calculate_individual_expression(group_b)

		return (de_a - de_b)

	def _calculate_individual_expression(self, gene_set):
		gene_set.sort()

		difference = (np.diff(gene_set, n=(len(gene_set)-1)))
		result = float(difference) / len(gene_set)
		return result

def main():
	gs = gsea()

	if len(sys.argv) != 3:
		gs.error_handler("Incorrect input given. \n $python gsea.py { gene_expression } { gene_profile }")
	else:
		gs.init_expressions(sys.argv[1])
		gs.init_pathways_profile(sys.argv[2], "leukemia")

		#print("Patient_expression_profile: {}".format(gs.data_expression_profile['patient_type']))

		#gs.narrow_set("leukemia")
		gs.test_gene_set_controller()
		#gs.show_test_results()

main()











































































































































