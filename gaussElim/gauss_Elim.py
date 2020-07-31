# -*- coding: utf-8 -*-
"""
Spyder Editor

Gaussian Elimination: A Toy Example (for square systems of dimension 3 or lower)
"""
import numpy as np
import string
import sys

def gaussElim(dim, tol=1e-6):
          """ 
          Inputs:
          - dim (integer): Dimension of the matrix (if you enter 3, it means a 3x3 matrix)
          - A (real): Matrix of coefficients (has to be a square matrix of dimension 3 or lower)
                         Enter A as a series of numbers separated by single spaces 
                         starting from A[1,1] and proceeding row-wise
                         (e.g. The 2x2 matrix [1 2 
                                               4 5] would be entered
                          as 1 2 4 5)
          - b (real): Coefficients on the right hand side in Ax = b
                         Enter b as a series of numbers starting from b[1,1] and 
                         proceeding down the column 
                         (e.g. The 3x1 column [1 
                                               2 
                                               3] 
                          would be entered as 1 2 3)
          
          Output: 
          - The final form of the augmented matrix, [U|b]
          - Solution for the system, if it exists
          
          Description:
          This function performs Gaussian elimination on square matrices of 
          dimension 3 or lower and calculates the soluion for the system Ax = b provided 
          the matrix is non-singular. It has the capability to do row exchanges 
          whenever they are necessary. The default tolerance used to check if a quanitity 
          is zero is 1x10^-6 i.e. if the absolute value of an entry is less than 1x10^-6, 
          it is considered to be equal to zero. Also, row exchanges will only 
          be done if the element below is greater than the tolerance. To change 
          the tolerance, pass it as the second argument to the function 
          (e.g. gaussElim(3,1e-4)).
          """
          ################## Inputs ##################
          
          # Check if dimension is okay
          try:
               int(dim)
          except ValueError:
               sys.exit("Please enter a valid integer for dimension.")
          else:
               if dim < 2 or dim > 3:
                    print("Dimension too low or too high. Check function help and try again.")
                    return
          
          # Messages for A and b inputs
          msg_A = "Enter matrix A (square with dimension 3 or lower) as a series of numbers separated by a space. Start from the first element and proceed row-wise.\n"
          msg_b = "Enter column vector b as a series of numbers separated by a space. Start from the first element proceed row-wise.\n"
          
          ### Input error handling ###
          while True:
               
               # Prompt for input. We want A and b
               A = input(msg_A)
               b = input(msg_b)
               
               # Check if A is okay
               try:
                    A = list(map(float, \
                        A.translate(str.maketrans('', '', \
                                        string.punctuation.replace("-",""))).split(" ")))
               except:
                    print("Invalid input for A. Check help for accepted input types.")
                    inp = input("Do you want to try again? y or n: \n")
                    if inp == "y":
                         continue
                    if inp == "n":
                         print("Goodbye!")
                         return 
                    else:
                         print("I didn't understand that. Goodbye!")
                         return
               
               # Check if b is okay
               try:
                    b = list(map(float, \
                        b.translate(str.maketrans('', '', \
                                        string.punctuation.replace("-",""))).split(" ")))
                    break
               except:
                    print("Invalid input for b. Check help for accepted input types.")
                    inp = input("Do you want to try again? y or n: \n")
                    if inp == "y":
                         continue
                    if inp == "n":
                         print("Goodbye!")
                         return
                    else:
                         print("I did not understand that. Goodbye!")
                         return
                    
          ################## Calculations ##################
          
          ### Convert matrix and column vector to numpy arrays ###
          
          # Matrix of coefficients
          try:
               A = np.array(A).reshape((dim, dim))
          except ValueError:
               sys.exit("Dimension entered and input for A are incompatible. Try again.")
          
          # Right hand side coefficients (form: Ax = b)
          try:
               b = np.array(b).reshape((dim,1))
          except ValueError:
               sys.exit("Dimension entered and input for b are incompatible. Try again.")
          
          # Augmented matrix (form: [AUG|b])
          AUG = np.concatenate((A,b),axis=1)
          print("\nAugmented Matrix is:\n", AUG)
          print("\n")
          
          ### Start looping row-wise ###
          
          # Remember the current pivot and have a temp row for exchanges
          curr_pivot = 0
          temp_row = np.zeros(AUG.shape[1])
          
          # Start loop
          for i in range(AUG.shape[0] - 1):
               
                    # Check if pivot is non-zero
                    if abs(AUG[i,i]) < tol: 
                         
                         # Try exchanging rows
                         if i == 0:
                              if abs(AUG[i+1,i]) > tol:
                                   temp_row = AUG[i,:].copy()
                                   AUG[i,:] = AUG[i+1,:]
                                   AUG[i+1,:] = temp_row
                              elif abs(AUG[i+2,i]) > tol: 
                                   temp_row = AUG[i,:].copy()                        
                                   AUG[i,:] = AUG[i+2,:]
                                   AUG[i+2,:] = temp_row
                              elif dim > 3 and abs(AUG[i+3,i]) > tol:
                                   temp_row = AUG[i,:].copy()                        
                                   AUG[i,:] = AUG[i+3,:]
                                   AUG[i+3,:] = temp_row
                              else:
                                   print("Reduced matrix:\n",AUG)
                                   print("\nI couldn't find a pivot. Possible singularity.\n")
                                   return
                         elif i == 1:
                              if abs(AUG[i+1,i]) > tol:
                                   temp_row = AUG[i,:].copy()
                                   AUG[i,:] = AUG[i+1,:]
                                   AUG[i+1,:] = temp_row
                              if abs(AUG[i+1,i]) > tol:
                                   temp_row = AUG[i,:].copy()
                                   AUG[i,:] = AUG[i+1,:]
                                   AUG[i+1,:] = temp_row
                              else:
                                   print("Reduced matrix:\n",AUG)
                                   print("\nI couldn't find a pivot. Possible singularity.\n")
                                   return
                         elif i == 2:
                              if abs(AUG[i+1,i]) > tol:
                                   temp_row = AUG[i,:].copy()
                                   AUG[i,:] = AUG[i+1,:]
                                   AUG[i+1,:] = temp_row
                              else:
                                   print("Reduced matrix:\n",AUG)
                                   print("\nI couldn't find a pivot. Possible singularity.\n")
                                   return
          
                    
                    # The current pivot
                    curr_pivot = AUG[i,i]
                    print("Pivot for row",i+1,"is:",curr_pivot,"\n")
                    
                    # Do row operations
                    zero_pivot = (abs(curr_pivot) < tol)
                    if zero_pivot:
                         print("At least one pivot is zero. The matrix is singular.\n")
                         return
                    elif i < 1:
                         if abs(AUG[i+1,i]) > tol: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]
                         if dim > 2 and abs(AUG[i+2,i]) > tol: 
                              AUG[i+2,:] = AUG[i+2,:] - (AUG[i+2,i]/curr_pivot)*AUG[i,:]
                    elif i < 2:
                         if abs(AUG[i+1,i]) > tol: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]
          
          ################## Output ##################
          
          # Matrix output               
          print("Upper triangular (augmented) matrix, [U|b]:\n", AUG)
          
          ### Calculate solution by back substitution ###
          
          x = np.zeros(dim)
          # Iterate in reverse for back substitution
          for i in range(dim-1,-1,-1):
               if i == dim-1:
                    x[i] = AUG[i,-1]/AUG[i,i]
               else:
                    x[i] = (AUG[i,-1] - np.dot(AUG[i,i+1:-1],x[i+1:]))/AUG[i,i]
               
          # Print solution
          if dim == 2:
               print("\nSolution is: x = {:.2f}, y = {:.2f}".format(x[0],x[1]))
          else: 
               print("\nSolution is: x = {:.2f}, y = {:.2f}, z = {:.2f}".format(x[0],x[1],x[2])) 
               
               
          return
          