import numpy as np
import sys, os



# https://www.lynda.com/NumPy-tutorials/NumPy-data-science-IMQAV/508873/543332-4.html

def testing():
	print "Hello World"
	a = np.arange(1,100000001)
	print a.sum()
#testing()


def sales_ex():
	sales_a = np.array([
		[5.0,2000.0],
		[10.0,500.0],
		[20.0,200.0]])
	
	print sales_a[:,0]
#sales_ex()




class gsea(object):
	# will need the class object to handle and store all that rickety data
	def __init__(self):
		self.data_expression_profile = {}
		# store if the patient is AML or ALL. This can operate later a key to the lists
		# upon construction, it will look like: "data_expression_profile['patient'] = ["AML", "ALL", "ALL", "ALL", "AML"]
		# upon construction, it will look like: "data_expression_profile['STAT1'] = [-732, 34, 55, 123, -5]"
		self.data_expression_profile['patient_type'] = []

		# raw strings from the pathways.txt
		self.gene_set_pathways = []

		# an iterable list of the expression profile keys
		self.narrowed_gene_expressions = []
		self.gene_samples_found = 0

		# store all of the test results here
		# how will you actually store your findings? 
		# test_result['GENE'] = [ results_int, results_int, results_int, ...]
		# test_result['GENE'][0] will ALWAYS contains the actual score. [1...n] is null hypothesis testing 
		self.test_result = {}

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
	def init_pathways_profile(self, source_file):
		print ("reading gene_expression_profile from file {}... ".format(source_file))
		with open(source_file) as fd:
			for line in fd:
				self.gene_set_pathways.append(line)


	# locate union of data_expression_profile from gene_set_pathways // narrow the list
	def narrow_set(self, search_term):
		desired_gene_set = self.show_profile(search_term)
		for k in self.data_expression_profile.keys():
			if k in desired_gene_set:
				#print k
				self.narrowed_gene_expressions.append(k)
				self.gene_samples_found += 1

		# testing print statement for myself~~~
		print ("Narrowed Set: {}".format(len(self.narrowed_gene_expressions)))
		#print ("Narrowed Set: {}".format(self.gene_samples_found))



	# run diference calculation
	def test_gene_set(self):
		# test set A and set B

		# test for null hypothesis, so test it a bunch
	# shuffle data set
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





def main():
	gs = gsea()

	if len(sys.argv) != 3:
		gs.error_handler("Incorrect input given. \n $python gsea.py { gene_expression } { gene_profile }")
	else:
		gs.init_expressions(sys.argv[1])
		gs.init_pathways_profile(sys.argv[2])

		#print("Patient_expression_profile: {}".format(gs.data_expression_profile['patient_type']))

		gs.narrow_set("leukemia")




main()











































































































































