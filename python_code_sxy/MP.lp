\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  400 y[0] + 414 y[1] + 326 y[2] + 18 z[0] + 25 z[1] + 20 z[2] + η
Subject To
 column1[0]: - 800 y[0] + z[0] <= 0
 column1[1]: - 800 y[1] + z[1] <= 0
 column1[2]: - 800 y[2] + z[2] <= 0
 z: z[0] + z[1] + z[2] >= 772
 column5: g[0] + g[1] <= 1.2
 column6: g[0] + g[1] + g[2] <= 1.8
Bounds
 y[0] <= 1
 y[1] <= 1
 y[2] <= 1
 g[0] <= 1
 g[1] <= 1
 g[2] <= 1
Generals
 y[0] y[1] y[2]
End
