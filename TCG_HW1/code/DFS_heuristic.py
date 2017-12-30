import sys
import string
import copy as copy
import types
import timeit

# DATA_FILE = sys.argv[2]
# n = sys.argv[1]
col_comb = []
row_comb = []
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
        # print permutation
        for x in xrange(start, start + int(float(counts[0]))):
            permutation.append(True)
        # print permutation
        x = start + int(float(counts[0]))
        if x < length:
            permutation.append(False)
            x += 1
        # print permutation, x
        if x == length and len(counts) == 0:
            print permutation, "??"
            permutations.append(permutation)
            break
        sub_start = x
        sub_rows = get_permutations(counts[1:len(counts)], length - sub_start)
       
        # print "sub rows", sub_rows
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
def delete_some_cols(board, n_col_comb, n):
	scanning = False
	for i in range(n):
		index, indexB = 0, 0
		while indexB < len(n_col_comb[i]):
			same = True
			for j in range(n):
				if board[j][i] != None:
					if n_col_comb[i][indexB][j] != board[j][i]:
						same = False
						scanning = True
			if same:
				n_col_comb[i][index] = n_col_comb[i][indexB]
				index+=1
			indexB+=1
		while len(n_col_comb[i]) > index:
			n_col_comb[i].pop()
	return scanning, n_col_comb
def delete_some_rows(board, n_row_comb, n):
	scanning = False
	for i in range(n):
		index, indexB = 0, 0
		while indexB < len(n_row_comb[i]):
			same = True
			for j in range(n):
				if board[i][j] != None:
					if n_row_comb[i][indexB][j] != board[i][j]:
						same = False
						scanning = True
			if same:
				n_row_comb[i][index] = n_row_comb[i][indexB]
				index+=1
			indexB+=1
		while len(n_row_comb[i]) > index:
			n_row_comb[i].pop()
	return scanning, n_row_comb
def fill_the_board_by_rows(n_row_comb, board, n):
	for i in range(n):
		sames = [True for k in range(n)]
		first_row = n_row_comb[i][0]
		for row_possible in n_row_comb[i]:
			for j in range(n):
				if row_possible[j] != first_row[j]:
					sames[j] = False
		for index, same in enumerate(sames):
			if same:
				if board[i][index] == None:
					board[i][index] = first_row[index]
	return board
def fill_the_board_by_cols(n_col_comb, board, n):
	for i in range(n):
		sames = [True for k in range(n)]
		first_col = n_col_comb[i][0]
		for col_possible in n_col_comb[i]:
			for j in range(n):
				if col_possible[j] != first_col[j]:
					sames[j] = False
		for index, same in enumerate(sames):
			if same and board[index][i] == None:
				board[index][i] = first_col[index]
	return board
def check_board(row_counts, col_counts, n_col_comb, n_row_comb, board, n):
	scanning = True
	no_solution = False
	temp_board = []
	temp_board = copy.deepcopy(board)
	# for i in range(n):
	# 	print board[i]
	while scanning:
		for i in range(n):
			if len(n_col_comb[i]) == 0 or len(n_row_comb[i]) == 0: 
				no_solution = True
				return n_row_comb, no_solution
		temp_board = fill_the_board_by_rows(n_row_comb, temp_board, n)
		scanning, n_col_comb = delete_some_cols(temp_board, n_col_comb, n)
		for i in range(n):
			if len(n_col_comb[i]) == 0 or len(n_row_comb[i]) == 0: 
				no_solution = True
				return n_row_comb, no_solution
		if scanning == False:
			break
		temp_board = fill_the_board_by_cols(n_col_comb, temp_board, n)
		scanning, new_row_comb = delete_some_rows(temp_board, n_row_comb, n)
		if scanning == False:
			break
	return n_row_comb, no_solution
def DFS_find_solution(col_counts, row_counts, n_col_comb, n_row_comb, board, current, n):
	stop = False
	no_solution = False
	new_col_comb = copy.deepcopy(n_col_comb)
	new_row_comb = copy.deepcopy(n_row_comb)
	# print "----------------DFS visiting-------------------"
	new_row_comb, no_solution = check_board(col_counts, row_counts, new_col_comb, new_row_comb, board, n)
	if no_solution:
		return board, not no_solution
	if current >= n:
		# print "================================end======================================="
		if no_solution:
			return board, False
		else:
			return board, True
	for i in range(len(new_row_comb[current])):
		board[current] = new_row_comb[current][i]
		tmp_board = copy.deepcopy(board)
		tmp_board, stop = DFS_find_solution(col_counts, row_counts, new_col_comb, new_row_comb, tmp_board, current+1, n)
		if stop:
			break
	if stop:
		return tmp_board, True
	else: 
		return board, False
if __name__ == "__main__":
	n = 0
	DATA_FILE = ""
	try:
		if sys.argv[1] and sys.argv[2]:
			n = int(float(sys.argv[1]))
			DATA_FILE = sys.argv[2]
	except: 
		print "need two arguments"
	print n, DATA_FILE
	start = timeit.default_timer()
	col = []
	row = []
	board = []
	count = 0
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
	board, stop = DFS_find_solution(col, row, col_comb, row_comb,board, 0, n)
	print "-----------------final result------------------"
	final_col = []
	final_row = []
	for i in range(n):
		tmp = []
		con_black_count = 0
		last_cell = False
		for j in range(n):
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
		final_row.append(tmp)
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
	print final_col
	print final_row
	for i in range(n):
		print ''.join([' 1' if cell else ' 0' for cell in board[i]])
	stop = timeit.default_timer()
	print stop - start 
