\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - z
Subject To
 R0: x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 <= 2
 R1: z + pi_0_0 - pi_0_5 <= 0
 R2: z + pi_1_0 - pi_1_5 <= 0
 R3: z + pi_2_0 - pi_2_5 <= 0
 R4: - 100 x_0 - pi_0_0 + pi_0_1 <= 3
 R5: - 100 x_0 - pi_1_0 + pi_1_1 <= 3
 R6: - 100 x_0 - pi_2_0 + pi_2_1 <= 3
 R7: - 100 x_1 - pi_0_0 + pi_0_2 <= 4
 R8: - 100 x_1 - pi_1_0 + pi_1_2 <= 1
 R9: - 100 x_1 - pi_2_0 + pi_2_2 <= 4
 R10: - 100 x_2 - pi_0_1 + pi_0_3 <= 3
 R11: - 100 x_2 - pi_1_1 + pi_1_3 <= 3
 R12: - 100 x_2 - pi_2_1 + pi_2_3 <= 3
 R13: - 100 x_3 - pi_0_2 + pi_0_4 <= 4
 R14: - 100 x_3 - pi_1_2 + pi_1_4 <= 1
 R15: - 100 x_3 - pi_2_2 + pi_2_4 <= 4
 R16: - 100 x_4 - pi_0_3 + pi_0_5 <= 3
 R17: - 100 x_4 - pi_1_3 + pi_1_5 <= 3
 R18: - 100 x_4 - pi_2_3 + pi_2_5 <= 3
 R19: - 100 x_5 - pi_0_4 + pi_0_5 <= 4
 R20: - 100 x_5 - pi_1_4 + pi_1_5 <= 1
 R21: - 100 x_5 - pi_2_4 + pi_2_5 <= 4
 R22: - 100 x_6 - pi_0_0 + pi_0_5 <= 3
 R23: - 100 x_6 - pi_1_0 + pi_1_5 <= 10
 R24: - 100 x_6 - pi_2_0 + pi_2_5 <= 10
 R25: pi_0_0 = 0
 R26: pi_1_0 = 0
 R27: pi_2_0 = 0
Bounds
 pi_0_0 free
 pi_0_1 free
 pi_0_2 free
 pi_0_3 free
 pi_0_4 free
 pi_0_5 free
 pi_1_0 free
 pi_1_1 free
 pi_1_2 free
 pi_1_3 free
 pi_1_4 free
 pi_1_5 free
 pi_2_0 free
 pi_2_1 free
 pi_2_2 free
 pi_2_3 free
 pi_2_4 free
 pi_2_5 free
Binaries
 x_0 x_1 x_2 x_3 x_4 x_5 x_6
End
