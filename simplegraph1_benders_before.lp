\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - zeta
Subject To
 R0: x_0 + x_1 + x_2 <= 2
 R1: - 100000 x_2 + zeta <= 3
 R2: - 100000 x_1 + zeta <= 1
 R3: - 100000 x_0 + zeta <= 9
Bounds
Binaries
 x_0 x_1 x_2
End
