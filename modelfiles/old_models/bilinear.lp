\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  pi_5
Subject To
 R0: lambda_0 + lambda_1 + lambda_2 = 1
 R1: x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 <= 2
 R2: pi_0 = 0
 qc0: - 3 lambda_0 - 3 lambda_1 - 3 lambda_2 - pi_0 + pi_1 + [
   - 100 x_0 * lambda_0 - 100 x_0 * lambda_1 - 100 x_0 * lambda_2 ] <= -0
 qc1: - 4 lambda_0 - lambda_1 - 4 lambda_2 - pi_0 + pi_2 + [
   - 100 x_1 * lambda_0 - 100 x_1 * lambda_1 - 100 x_1 * lambda_2 ] <= -0
 qc2: - 3 lambda_0 - 3 lambda_1 - 3 lambda_2 - pi_1 + pi_3 + [
   - 100 x_2 * lambda_0 - 100 x_2 * lambda_1 - 100 x_2 * lambda_2 ] <= -0
 qc3: - 4 lambda_0 - lambda_1 - 4 lambda_2 - pi_2 + pi_4 + [
   - 100 x_3 * lambda_0 - 100 x_3 * lambda_1 - 100 x_3 * lambda_2 ] <= -0
 qc4: - 3 lambda_0 - 3 lambda_1 - 3 lambda_2 - pi_3 + pi_5 + [
   - 100 x_4 * lambda_0 - 100 x_4 * lambda_1 - 100 x_4 * lambda_2 ] <= -0
 qc5: - 4 lambda_0 - lambda_1 - 4 lambda_2 - pi_4 + pi_5 + [
   - 100 x_5 * lambda_0 - 100 x_5 * lambda_1 - 100 x_5 * lambda_2 ] <= -0
 qc6: - 3 lambda_0 - 10 lambda_1 - 10 lambda_2 - pi_0 + pi_5 + [
   - 100 x_6 * lambda_0 - 100 x_6 * lambda_1 - 100 x_6 * lambda_2 ] <= -0
Bounds
 lambda_0 <= 1
 lambda_1 <= 1
 lambda_2 <= 1
 pi_0 free
 pi_1 free
 pi_2 free
 pi_3 free
 pi_4 free
 pi_5 free
Binaries
 x_0 x_1 x_2 x_3 x_4 x_5 x_6
End
