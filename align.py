import sys
import copy


#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """
    def __init__(self):
        self.match_matrix = {}

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it to
        """
        self.match_matrix[(a, b)] = score

    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        if (a, b) in self.match_matrix.keys():
            return self.match_matrix[(a, b)]
        else:
            print('Your key is not in the match matrix!')

    def print_matches(self):
        """
        Returns a nicely formatted string containing the scores in the match matrix

        Example: ('A', 'A') 1

        """
        for key in self.match_matrix.keys():
            print(str(key) + '\t' + str(self.match_matrix[key]))

class ScoreEntry(object):
    """
    Object to represent a score in the score matrix.
    """
    def __init__(self):
        self.score = 0
        self.pointers = []

    def get_score(self):
        """
        Returns the score of the score entry
        """
        return self.score

    def set_score(self, score):
        """
        Updates or adds a score for the score entry

        Inputs:
            score = the score to set it to
        """
        self.score = score

    def get_pointers(self):
        """
        Returns the pointers of the score entry
        """
        return self.pointers

    def set_pointers(self, pointers):
        """
        Updates or adds pointers for the score entry

        Inputs:
            pointers = the pointers to set it to
        """
        self.pointers = pointers

class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix = [[ScoreEntry() for b in range(ncol)] for a in range(nrow)]

    def get_score(self, row, col):
        """
        Returns the score (a float) of the score entry object at the given row and column

        Inputs:
            row = the row of the score entry
            col = the col of the score entry
        """
        if row >= 0 & col >= 0:
            return self.score_matrix[row][col].get_score()
        else:
            print("Your indices are negative!")
        
    def set_score(self, row, col, score):
        """
        Updates or adds the score (a float) to the score entry object at the given row and column

        Inputs:
            row = the row of the score entry
            col = the col of the score entry
        """
        if row >=0 & col >= 0:  
            self.score_matrix[row][col].set_score(score)
        else:
            print("Your indices are negative!")

    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to (formatted as a list of tuples)
        
        Example:
        [(1,1), (1,0)]
        """
        if row >= 0 & col >= 0:
            return self.score_matrix[row][col].get_pointers()
        else:
            print("Your indices are negative!")

    def set_pointers(self, row, col, pointers):
        """
        Sets the pointers to the score entry at the given row and column
        Inputs:
            row = the row of the score entry being set
            col = the column of the score entry being set
        """
        if row >= 0 & col >= 0:
            self.score_matrix[row][col].set_pointers(pointers)
        else:
            print("Your indices are negative!")

    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """
        print(self.name + '=')
        print('\t')
        for i in range(self.nrow):
            for j in range(self.ncol):
                print(str(round(self.get_score(i, j), 1)) + ', ', end = '')
            print('')

    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """
        print(self.name + ' pointers =')
        for i in range(self.nrow):
            for j in range(self.ncol):
                print('(' + str(i) + ', ' + str(j) + '): ', end = '')
                print(self.get_pointers(i, j), end = '')
            print('')

class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """

    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """

        with open(input_file, 'r') as file: #
        	split_file = file.read().splitlines()

        self.seq_a = split_file[0]
        self.seq_b = split_file[1]

        if split_file[2] == '0':
        	self.global_alignment = True

        gap_penalties = split_file[3].split()
        self.dx = float(gap_penalties[0])
        self.ex = float(gap_penalties[1])
        self.dy = float(gap_penalties[2])
        self.ey = float(gap_penalties[3])
        self.len_alphabet_a = int(split_file[4])
        self.alphabet_a = split_file[5]
        self.len_alphabet_b = int(split_file[6])
        self.alphabet_b = split_file[7]
        
        #going through each entry in the input match matrix
        #key is the pair of letters being matched, (a tuple)
        #value is the score of the pair, (a float)
        end = self.len_alphabet_a * self.len_alphabet_b
        for i in range(8, 8 + end):
            match = split_file[i].split()
            a = str(match[2])
            b = str(match[3])
            score = float(match[4])
            self.match_matrix.set_score(a, b, score)

class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters() 

        nrow = 50
        ncol = 50

        self.m_matrix = ScoreMatrix('M', nrow, ncol)
        self.ix_matrix = ScoreMatrix('Ix', nrow, ncol)
        self.iy_matrix = ScoreMatrix('Iy', nrow, ncol)

    def align(self):
        """
        Main method for running alignment.
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        nrow = len(self.align_params.seq_a) + 1
        ncol = len(self.align_params.seq_b) + 1

        self.m_matrix = ScoreMatrix('M', nrow, ncol)
        self.ix_matrix = ScoreMatrix('Ix', nrow, ncol)
        self.iy_matrix = ScoreMatrix('Iy', nrow, ncol)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # perform a traceback and write the output to an output file
        maxes = self.find_traceback_start()
        paths = self.traceback()
        self.write_output(maxes, paths)

    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """

        nrow = len(self.align_params.seq_a) + 1
        ncol = len(self.align_params.seq_b) + 1

        for i in range(ncol):
            self.m_matrix.set_score(0, i, float(0))
            self.ix_matrix.set_score(0, i, float(0))
            self.iy_matrix.set_score(0, i, float(0))

        for j in range(nrow):
            self.m_matrix.set_score(j, 0, float(0))
            self.ix_matrix.set_score(j, 0, float(0))
            self.iy_matrix.set_score(j, 0, float(0))

        for i in range(1, nrow):
            for j in range(1, ncol):
                self.update(i, j)

    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row, col):
        """
        Updates the m matrix based on the function defined in Durbin and lecture:
        Take the score from the previous diagonal entry in the m, ix, and iy matrices
        Add the offset of matching the two current letters in sequence a and sequence b

        Set the score entry to the max score of the above method and its associated pointers 

        Input:
            row = the row index to update
            col = the column index to update
        """
        seq_a = self.align_params.seq_a
        seq_b = self.align_params.seq_b
        letter_in_a = seq_a[(row - 1):row]
        letter_in_b = seq_b[(col - 1):col]
        match_matrix = self.align_params.match_matrix

        score_m = self.m_matrix.get_score(row - 1, col - 1) + match_matrix.get_score(letter_in_a, letter_in_b)
        score_ix = self.ix_matrix.get_score(row - 1, col - 1) + match_matrix.get_score(letter_in_a, letter_in_b)
        score_iy = self.iy_matrix.get_score(row - 1, col - 1) + match_matrix.get_score(letter_in_a, letter_in_b)
        
        scores_and_pointers = [(score_m, ('M', row - 1, col - 1)), (score_ix, ('Ix', row - 1, col - 1)), (score_iy, ('Iy', row - 1, col - 1))]

        #find the max score
        max_score_and_pointer = max(scores_and_pointers)
        max_score = max_score_and_pointer[0]
        max_pointers = []

        #add pointers if there several scores that are the max
        for i in range(3):
            if fuzzy_equals(scores_and_pointers[i][0], max_score):
                max_pointers.append(scores_and_pointers[i][1])

        #if aligning locally, discard all the negative scores
        if self.align_params.global_alignment == False:
            if max_score < 0:
                max_score = 0
                max_pointers = []

        #update corresponding score entry
        self.m_matrix.set_score(row, col, max_score)
        self.m_matrix.set_pointers(row, col, max_pointers)

    def update_ix(self, row, col):
        """
        Updates the ix matrix based on the function defined in Durbin and lecture:
        Take the score from the previous vertical entry in the m and ix matrices
        Subtract the penalty of adding a gap in sequence b

        Set the score entry to the max score of the above method and its associated pointers 

        Input:
            row = the row index to update
            col = the column index to update
        """
        score_m = self.m_matrix.get_score(row - 1, col) - self.align_params.dy
        score_ix = self.ix_matrix.get_score(row - 1, col) - self.align_params.ey
        
        scores_and_pointers = [(score_m, ('M', row - 1, col)), (score_ix, ('Ix', row - 1, col))]

        #find the max score
        max_score_and_pointer = max(scores_and_pointers)
        max_score = max_score_and_pointer[0]
        max_pointers = []

        #add pointers if there several scores that are the max
        for i in range(2):
            if fuzzy_equals(scores_and_pointers[i][0], max_score):
                max_pointers.append(scores_and_pointers[i][1])

        #if aligning locally, discard all the negative scores
        if self.align_params.global_alignment == False:
            if max_score < 0:
                max_score = 0
                max_pointers = []

        #update corresponding score entry
        self.ix_matrix.set_score(row, col, max_score)
        self.ix_matrix.set_pointers(row, col, max_pointers)

    def update_iy(self, row, col):
        """
        Updates the iy matrix based on the function defined in Durbin and lecture:
        Take the score from the previous horizontal entry in the m and iy matrices
        Subtract the penalty of introducing a gap in sequence a

        Set the score entry to the max score of the above method and its associated pointers 

        Input:
            row = the row index to update
            col = the column index to update
        """
        score_m = self.m_matrix.get_score(row, col - 1) - self.align_params.dx
        score_iy = self.iy_matrix.get_score(row, col - 1) - self.align_params.ex
        
        scores_and_pointers = [(score_m, ('M', row, col - 1)), (score_iy, ('Iy', row, col - 1))]
        
        #find the max score
        max_score_and_pointer = max(scores_and_pointers)
        max_score = max_score_and_pointer[0]
        max_pointers = []

        #add pointers if there several scores that are the max
        for i in range(2):
            if fuzzy_equals(scores_and_pointers[i][0], max_score):
                max_pointers.append(scores_and_pointers[i][1])

        #if aligning locally, discard all the negative scores
        if self.align_params.global_alignment == False:
            if max_score < 0:
                max_score = 0
                max_pointers = []

        #update corresponding score entry
        self.iy_matrix.set_score(row, col, max_score)
        self.iy_matrix.set_pointers(row, col, max_pointers)

    def find_traceback_start(self):
        """
        Finds the location to start the traceback.

        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a list [] containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """

        #global alignment
        if self.align_params.global_alignment == True:
            nrow = len(self.align_params.seq_a) + 1
            ncol = len(self.align_params.seq_b) + 1

            max_entries = []

            #only check entries in the last row
            for i in range(ncol):
                score_m = self.m_matrix.get_score(nrow - 1, i)
                score_ix = self.ix_matrix.get_score(nrow - 1, i)
                score_iy = self.iy_matrix.get_score(nrow - 1, i)

                max_entries.append(score_m)
                max_entries.append(score_ix)
                max_entries.append(score_iy)

            #only check entries in the last column
            for j in range(nrow - 1):
                score_m = self.m_matrix.get_score(j, ncol - 1)
                score_ix = self.ix_matrix.get_score(j, ncol - 1)
                score_iy = self.iy_matrix.get_score(j, ncol - 1)

                max_entries.append(score_m)
                max_entries.append(score_ix)
                max_entries.append(score_iy)

            max_val = max(max_entries)
            max_loc = []

            #find all locations that correspond to the max score in the last row and column
            for i in range(ncol):
                if fuzzy_equals(max_val, self.m_matrix.get_score(nrow - 1, i)):
                    max_loc.append(('M', nrow - 1, i))

                if fuzzy_equals(max_val, self.ix_matrix.get_score(nrow - 1, i)):
                    max_loc.append(('Ix', nrow - 1, i))

                if fuzzy_equals(max_val, self.iy_matrix.get_score(nrow - 1, i)):
                    max_loc.append(('Iy', nrow - 1, i))

            for j in range(nrow - 1):
                if fuzzy_equals(max_val, self.m_matrix.get_score(j, ncol - 1)):
                    max_loc.append(('M', j, ncol - 1))

                if fuzzy_equals(max_val, self.ix_matrix.get_score(j, ncol - 1)):
                    max_loc.append(('Ix', j, ncol - 1))

                if fuzzy_equals(max_val, self.iy_matrix.get_score(j, ncol - 1)):
                    max_loc.append(('Iy', j, ncol - 1))

        #local alignment
        else:
            max_entries = []

            #check all entries the matrices
            for i in range(len(self.align_params.seq_a) + 1):
                for j in range(len(self.align_params.seq_b) + 1):
                    score_m = self.m_matrix.get_score(i, j)
                    score_ix = self.ix_matrix.get_score(i, j)
                    score_iy = self.iy_matrix.get_score(i, j)

                    max_entries.append(score_m)
                    max_entries.append(score_ix)
                    max_entries.append(score_iy)

            max_val = max(max_entries)
            max_loc = []

            #find all locations that correspond to the max score
            for i in range(len(self.align_params.seq_a) + 1):
                for j in range(len(self.align_params.seq_b) + 1):
                    if fuzzy_equals(max_val, self.m_matrix.get_score(i, j)):
                        max_loc.append(('M', i, j))

                    if fuzzy_equals(max_val, self.ix_matrix.get_score(i, j)):
                        max_loc.append(('Ix', i, j))

                    if fuzzy_equals(max_val, self.iy_matrix.get_score(i, j)):
                        max_loc.append(('Iy', i, j))
        
        return (round(max_val, 1), max_loc)

    def traceback_global_rec(self, curr_entry, curr_path, best_paths):
    	"""
    	Recursive helper function for performing a global traceback

    	Inputs:
    		curr_entry is the current ScoreEntry being visited
			curr_path is the current path being built out
			best_paths is a list of all best paths
    	"""
    	matrix_name = curr_entry[0]
    	entry_row = curr_entry[1]
    	entry_col = curr_entry[2]
    	copy_path = copy.deepcopy(curr_path) #deep copy since python does not pass by value

        #stop in the m matrix if either row or column are 1
    	if fuzzy_equals(entry_row, 1) or fuzzy_equals(entry_col, 1) and matrix_name == 'M':
    		copy_path.append(curr_entry)
    		best_paths.append(copy_path)
    		return

        #stop in the ix matrix only if row is 1: traceback can still move vertically
    	if fuzzy_equals(entry_row, 1) and matrix_name == 'Ix':
    		copy_path.append(curr_entry)
    		best_paths.append(copy_path)
    		return

        #stop in the iy matrix only if column is 1: traceback can still move horizontally
    	if fuzzy_equals(entry_col, 1) and matrix_name == 'Iy':
    		copy_path.append(curr_entry)
    		best_paths.append(copy_path)
    		return

    	else:
    		copy_path.append(curr_entry)

    		#get score and pointers from m matrix based on pointer
    		if matrix_name == "M":
    			pointers_m = self.m_matrix.get_pointers(entry_row, entry_col)
    			for pointer in pointers_m:
    				self.traceback_global_rec(pointer, copy_path, best_paths)

    		#get score and pointers from ix matrix based on pointer
    		if matrix_name == 'Ix':
    			pointers_ix = self.ix_matrix.get_pointers(entry_row, entry_col)
    			for pointer in pointers_ix:
    				self.traceback_global_rec(pointer, copy_path, best_paths)

    		#get score and pointers from iy matrix based on pointer
    		if matrix_name == 'Iy':
    			pointers_iy = self.iy_matrix.get_pointers(entry_row, entry_col)
    			for pointer in pointers_iy:
    				self.traceback_global_rec(pointer, copy_path, best_paths)

    def traceback_local_rec(self, curr_entry, curr_path, best_paths):
    	"""
    	Recursive helper function for performing a local traceback

    	Inputs:
    		curr_entry is the current ScoreEntry
			curr_path is the current path being built out
			best_paths is a list of all best paths
    	"""
    	matrix_name = curr_entry[0]
    	entry_row = curr_entry[1]
    	entry_col = curr_entry[2]
    	copy_path = copy.deepcopy(curr_path) #deep copy since python does not pass by value

        #get score and pointers from m matrix based on pointer
    	if matrix_name == 'M':
        	score_m = self.m_matrix.get_score(entry_row, entry_col)

        	if fuzzy_equals(score_m, 0):
        		best_paths.append(copy_path)
        		return

        	else:
        		copy_path.append(curr_entry)
        		pointers_m = self.m_matrix.get_pointers(entry_row, entry_col)
        		for pointer in pointers_m:
        			self.traceback_local_rec(pointer, copy_path, best_paths)    		
    	
        #get score and pointers from ix matrix based on pointer
    	if matrix_name == 'Ix':

        	score_ix = self.ix_matrix.get_score(entry_row, entry_col)

        	if fuzzy_equals(score_ix, 0):
        		best_paths.append(copy_path)
        		return

        	else:
        		copy_path.append(curr_entry)
        		pointers_ix = self.ix_matrix.get_pointers(entry_row, entry_col)
        		for pointer in pointers_ix:
        			self.traceback_local_rec(pointer, copy_path, best_paths)

        #get score and pointers from iy matrix based on pointer
    	if matrix_name == 'Iy':
        	score_iy = self.iy_matrix.get_score(entry_row, entry_col)

        	if fuzzy_equals(score_iy, 0):
        		best_paths.append(copy_path)
        		return
        	
        	else: 
        		copy_path.append(curr_entry)
        		pointers_iy = self.iy_matrix.get_pointers(entry_row, entry_col)
        		for pointer in pointers_iy:
        			self.traceback_local_rec(pointer, copy_path, best_paths)

    def traceback(self):
        """
        Performs a traceback.
        Hint: include a way to printing the traceback path. This will be helpful for debugging!
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)

        """
        best_paths = []
        maxes = self.find_traceback_start()
        max_loc = maxes[1]

        #global alignment
        if self.align_params.global_alignment == True:
        	for loc in max_loc:
        		self.traceback_global_rec(loc, [], best_paths)
        
        #local alignment
        else:
        	for loc in max_loc:
        		self.traceback_local_rec(loc, [], best_paths)

        #building up alignments (derived from the given sequences) based on traceback paths
        best_alignments = []
        
        for path in best_paths:
            
            alignment_a = ''
            alignment_b = ''
            
            for entry in path:

                if entry[0] == 'M':
                    alignment_a += self.align_params.seq_a[entry[1] - 1]
                    alignment_b += self.align_params.seq_b[entry[2] - 1]
                
                if entry[0] == 'Ix':
                    alignment_a += self.align_params.seq_a[entry[1] - 1]
                    alignment_b += '_'

                if entry[0] == 'Iy':
                    alignment_a += '_'
                    alignment_b += self.align_params.seq_b[entry[2] - 1]
            
            #reversing strings since they are built backwards
            alignment_a = alignment_a[::-1]
            alignment_b = alignment_b[::-1] 
            best_alignments.append((alignment_a, alignment_b))

        #discarding any alignments that are repeated
        best_alignments_no_reps = []

        for alignment in best_alignments:
        	if alignment not in best_alignments_no_reps:
        		best_alignments_no_reps.append(alignment)

        best_alignments = best_alignments_no_reps
        return best_alignments

    def write_output(self, maxes, alignments):
    	"""
        Writes to the given output file the max score of the traceback and all the alignments with that score

        Inputs:
            maxes = a tuple containing the max score (a float) as its first element
            and the location(s) of these max scores (a list of tuples) as its second element
            
            alignments = a list of alignments (tuples) that were computed after traceback
            the first element of the tuple is the alignment of the first sequence
            the second element of the tuple is the alignment of the second sequence
    	"""
    	with open(self.output_file, 'w') as file:
            file.write(str(maxes[0]))
            file.write('\n')
            file.write('\n')
            for alignment in alignments:
            	file.write(str(alignment[0]))
            	file.write('\n')
            	file.write(str(alignment[1]))
            	file.write('\n')
            	file.write('\n')

def main():

    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()

if __name__=="__main__":
    main()
