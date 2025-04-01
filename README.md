# Problem Description
使用有限體積法離散SteadyNS，並使用SIMPLE algorithm求得流體速度與壓力

# Geometry and Fluid properties
+ `2D channel: Length L = 0.1 m, Width 2h = 0.01 m `
+ `Origin located at the inlet along the centerline `
+ `Uniform velocity of 0.01 m/s at the inlet `
+ `Uniform pressure of 0 Pa at the outlet`
+ `Density of ρ = 1.205 kg/m3 `
+ `Dynamic viscosity of μ = 1.82 x 10-5 Pa-s`

# Anylytical Solutuin
$$
\begin{aligned}
u(y) &= \frac{1}{2\mu} \frac{dP}{dx} (y^2 - h^2) \\
u_{avg} &= -\frac{h^2}{3\mu} \frac{dP}{dx} \\
u_{max} &= -\frac{h^2}{2\mu} \frac{dP}{dx}
\end{aligned}
$$

#Mesh format
|numRow|numColume|
|-|-|
|x(0)|y(0)|
|x(1)|y(2)|
|.|.|
|.|.|
|x(N)|y(N)|
##  
網格點順序為下到上、左至右。如下圖所示  
  
<img src="https://github.com/KWGHG/FVM_SIMPLE_2D_SteadyLaminarNS/blob/main/mesh%20format.jpg" width="400" />

# Dependences
+ `Eigen-3.4.0`
  安裝至Dependences，記得要更改VS庫設定
