\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  zeta_sub
Subject To
 R0: zeta_sub - 9 y_0_0 - 12 y_0_1 - 3 y_0_2 >= 0
 R1: zeta_sub - 9 y_1_0 - y_1_1 - 10 y_1_2 >= 0
 R2: zeta_sub - 9 y_2_0 - 12 y_2_1 - 10 y_2_2 >= 0
 R3: y_0_0 + y_0_1 + y_0_2 = 1
 R4: - y_0_0 - y_0_1 - y_0_2 = -1
 R5: y_1_0 + y_1_1 + y_1_2 = 1
 R6: - y_1_0 - y_1_1 - y_1_2 = -1
 R7: y_2_0 + y_2_1 + y_2_2 = 1
 R8: - y_2_0 - y_2_1 - y_2_2 = -1
Bounds
End
