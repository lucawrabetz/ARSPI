\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - z
Subject To
 R0: x_0 + x_1 + x_2 <= 2
 R1: z + pi_0_0 - pi_0_1 + lambda_0_0 + lambda_0_1 + lambda_0_2 <= 0
 R2: z + pi_1_0 - pi_1_1 + lambda_1_0 + lambda_1_1 + lambda_1_2 <= 0
 R3: z + pi_2_0 - pi_2_1 + lambda_2_0 + lambda_2_1 + lambda_2_2 <= 0
 R4: - 100000 x_0 - pi_0_0 + pi_0_1 - lambda_0_0 <= 9
 R5: - 100000 x_1 - pi_0_0 + pi_0_1 - lambda_0_1 <= 12
 R6: - 100000 x_2 - pi_0_0 + pi_0_1 - lambda_0_2 <= 3
 R7: - 100000 x_0 - pi_1_0 + pi_1_1 - lambda_1_0 <= 9
 R8: - 100000 x_1 - pi_1_0 + pi_1_1 - lambda_1_1 <= 1
 R9: - 100000 x_2 - pi_1_0 + pi_1_1 - lambda_1_2 <= 10
 R10: - 100000 x_0 - pi_2_0 + pi_2_1 - lambda_2_0 <= 9
 R11: - 100000 x_1 - pi_2_0 + pi_2_1 - lambda_2_1 <= 12
 R12: - 100000 x_2 - pi_2_0 + pi_2_1 - lambda_2_2 <= 10
 R13: pi_0_0 = 0
 R14: pi_1_0 = 0
 R15: pi_2_0 = 0
Bounds
 pi_0_0 free
 pi_0_1 free
 pi_1_0 free
 pi_1_1 free
 pi_2_0 free
 pi_2_1 free
Binaries
 x_0 x_1 x_2
End
