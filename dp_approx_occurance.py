# Implemented dynamic programming to find the approximate
# occurrence of a given DNA sequence in a reference genome.

def editDistance(x, y):
    # Create distance matrix
	D = []
	for i in range(len(x)+1):
		D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
	for i in range(len(x)+1):
		D[i][0] = i
    # Fill in the rest of the matrix
	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			distHor = D[i][j-1] + 1
			distVer = D[i-1][j] + 1
			if x[i-1] == y[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + 1
			D[i][j] = min(distHor, distVer, distDiag)
	last_row= D[-1]
    # Edit distance is the value in the bottom right corner of the matrix
	return min(last_row)