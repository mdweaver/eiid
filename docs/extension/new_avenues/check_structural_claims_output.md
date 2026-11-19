# Structural checks

C1: 300 random unit tables (R in 2..4, C1 in 2..4, C2 in 2..3): max |DD(one-way) - DD(joint)| = 2.22e-16  -> PASS
C2: R=2, C=2: dim K = 1, rank of K -> B* = 1, predicted min(d,(R-1)(C-1)) = 1; free directions left by the true table = 0
C2: R=2, C=3: dim K = 2, rank of K -> B* = 2, predicted min(d,(R-1)(C-1)) = 2; free directions left by the true table = 0
C2: R=3, C=2: dim K = 3, rank of K -> B* = 2, predicted min(d,(R-1)(C-1)) = 2; free directions left by the true table = 1
C2: R=3, C=3: dim K = 6, rank of K -> B* = 4, predicted min(d,(R-1)(C-1)) = 4; free directions left by the true table = 2
C2: R=4, C=3: dim K = 12, rank of K -> B* = 6, predicted min(d,(R-1)(C-1)) = 6; free directions left by the true table = 6
C2: R=4, C=4: dim K = 18, rank of K -> B* = 9, predicted min(d,(R-1)(C-1)) = 9; free directions left by the true table = 9
C2: R=5, C=3: dim K = 20, rank of K -> B* = 8, predicted min(d,(R-1)(C-1)) = 8; free directions left by the true table = 12
C2: PASS
