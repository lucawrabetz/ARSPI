\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  zeta_sub_1
Subject To
 R0: zeta_sub_1 - 9 y_1_0 - 100001 y_1_1 - 100010 y_1_2 >= 0
 R1: y_1_0 + y_1_1 + y_1_2 = 1
 R2: - y_1_0 - y_1_1 - y_1_2 = -1
Bounds
 y_1_0 <= 1
 y_1_1 <= 1
 y_1_2 <= 1
End
