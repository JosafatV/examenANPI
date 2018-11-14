function lookupExample ()
  A = [1, 2, 3, 4, 5, 6, 7.1, 8, 9, 9, 9, 10]
  
  display ("length")
    length = length(A)
 
  display ("expected")
    lookup(A, 5) # returns index
    lookup(A, 9) # returns largest possible index
    
  display ("out of bounds")
    lookup(A,-10) # returns index 0: NE
    lookup(A, 0)  # returns index 0: NE
    lookup(A,100) # returns index n: length
    
  display ("in between")
    lookup(A, 7.2)    # return index 7 truncates
    lookup(A, 7.5)    # return index 7 truncates
    lookup(A, 7.9999) # return index 7 truncates
  
  display ("m")
    lookup(A, 5, "m") # 5 is contained in A -> 5 (self)
    lookup(A, 7, "m") # 5 is contained in A -> 7 (self)
    lookup(A, 0, "m") # 0 is not contained in A -> 0
  
  display ("b")
    lookup(A, 5, "b") # 5 is contained in A -> true
    lookup(A, 0, "b") # 0 is not contained in A -> false
  
  display ("l")
    lookup(A,-2, "l") # return 1 -> leftmost index
    lookup(A, 9, "l") # return largest index
    lookup(A,30, "l") # return n -> length
  
  display ("r")
    lookup(A,-2, "r") # returns index 0 -> NE
    lookup(A, 9, "r") # return largest index
    lookup(A,30, "r") # return 11 -> last index

endfunction