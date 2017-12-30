import sys
import string
import copy as copy
import types
from collections import defaultdict
col_comb = []
row_comb = []
board_set = []
def get_permutations(counts, length):
    if len(counts) == 0:
        row = []
        for x in xrange(length):
            row.append(False)
        return [row]
    if counts[0] == '':
    	row = []
        for x in xrange(length):
            row.append(False)
    	return [row]
    permutations = []
    for start in xrange(length - int(float(counts[0])) + 1):
        permutation = []
        for x in xrange(start):
            permutation.append(False)
        for x in xrange(start, start + int(float(counts[0]))):
            permutation.append(True)
        x = start + int(float(counts[0]))
        if x < length:
            permutation.append(False)
            x += 1
        if x == length and len(counts) == 0:
            permutations.append(permutation)
            break
        sub_start = x
        sub_rows = get_permutations(counts[1:len(counts)], length - sub_start)
        for sub_row in sub_rows:
            sub_permutation = copy.deepcopy(permutation)
            for x in xrange(sub_start, length):
                sub_permutation.append(sub_row[x-sub_start])
            permutations.append(sub_permutation)
    return permutations
def generate_col(col, index, n):
	all_possible = get_permutations(col[index], n)
	col_comb.append(all_possible)
def generate_row(row, index, n):
	all_possible = get_permutations(row[index], n)
	row_comb.append(all_possible)
def check(board, col_counts, row_counts, n):
	final_col = []
	for j in range(n):
		tmp = []
		con_black_count = 0
		last_cell = False
		for i in range(n):
			if board[i][j]:
				if last_cell:
					con_black_count += 1
				last_cell = True
			elif not board[i][j]:
				if last_cell:
					tmp.append(con_black_count + 1)
					con_black_count = 0
				last_cell = False
		if last_cell:
			tmp.append(con_black_count + 1)
		final_col.append(tmp)
	# print final_col
	correct = True
	for index, col in enumerate(final_col):
		if col != col_counts[index]:
			correct = False
	if correct:
		for i in range(n):
			print board[i]
	return correct

def BFS_find_solution(col_counts, row_counts, n_col_comb, n_row_comb, board, current, n):
	tmp_board = []
	correct = False
	if current >= n:
		tmp_board = copy.deepcopy(board)
		board_set.append(tmp_board)
		correct = check(tmp_board, col_counts, row_counts, n)
		if correct:
			return True
		else:
			return False
	for i in range(len(n_row_comb[current])):
		print current, i
		board[current] = n_row_comb[current][i]
		stop = BFS_find_solution(col_counts, row_counts, n_col_comb, n_row_comb, board, current+1, n)
		if stop:
			break
	return stop
if __name__ == "__main__":
	col = []
	row = []
	board = []
	count = 0
	n = 0
	DATA_FILE = ""
	try:
		if sys.argv[1] and sys.argv[2]:
			n = int(float(sys.argv[1]))
			DATA_FILE = sys.argv[2]
	except: 
		print "need two arguments"
	# read nonogram data
	with open(DATA_FILE, "r") as f:
		for line in f:
			line = line.replace('\r','').replace('\n','')
			line = line.split('\t')
			# print line
			if count == 0:
				count+=1
				continue
			if count <= n: 
				col.append(line)
				count+=1
			elif count > n: 
				row.append(line)
				count+=1
	# initialize the board
	for i in range(n):
		tmp = []
		for j in range(n):
			tmp.append(None)
		board.append(tmp)
	for i in range(n):
		generate_col(col, i, n)
		generate_row(row, i, n)
	stop = BFS_find_solution(col, row, col_comb, row_comb,board, 0, n)
	print "-----------------final result------------------"
	for i in range(n):
		print ''.join(['1' if cell else '0' for cell in board[i]])