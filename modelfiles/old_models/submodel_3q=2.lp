\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  zeta_sub_2
Subject To
 R0: zeta_sub_2 - 100009 y_2_0 - 12 y_2_1 - 10 y_2_2 >= 0
 R1: y_2_0 + y_2_1 + y_2_2 = 1
 R2: - y_2_0 - y_2_1 - y_2_2 = -1
Bounds
 y_2_0 <= 1
 y_2_1 <= 1
 y_2_2 <= 1
End
