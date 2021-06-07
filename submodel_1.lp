\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  zeta_sub
Subject To
 R0: zeta_sub - 9 y_0 - 12 y_1 - 3 y_2 >= 0
 R1: zeta_sub - 9 y_0 - y_1 - 10 y_2 >= 0
 R2: zeta_sub - 9 y_0 - 12 y_1 - 10 y_2 >= 0
 R3: y_0 + y_1 + y_2 = 1
 R4: - y_0 - y_1 - y_2 = -1
Bounds
End
