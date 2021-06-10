# 三阶静倒立摆的LQR控制器设计

## 1、windows端：

## 用visual studio 2017 x64 

## 依赖于外部库Eigen、boost-odeint 和 gnuplot

## 在windows端最方便的方法是利用visual studio（2015及以上的版本） 和 vcpkg 来编译原代码

## 2、linux 端（recommanded）：

## 先安装必要的工具g++和cmake

## 安装g++：

```shell
sudo apt-get install gcc
```

## 验证g++：

```shell
gcc --version
```

## 安装cmake: 

```shell
sudo apt-get install cmake
```

## 验证cmake:

```bash
cmake -version
```

## 将单独的matrixoperation.cpp文件和CMakeLists.txt文件复制到linux根目录下的新创建空间cppSpace下。

## 依次执行以下命令：

```shell
cd cppSpace
cmake .
make
./matrixoperation
```

下载gnuplot

## 即可得到图像



## 结果解释：

![](C:\Users\guixiangyu\Desktop\tuxiang.png)

第一行左是第一个关节的角度；

第一行右是第二个关节的角度；

第二行左是第三个关节的角度；

第二行右是第一个关节的角速度；

第三行左是第二个关节的角速度；

第三行右是第三个关节的角速度；

直接运行可以进入Desktop\lqr_cpp\x64\Debug 运行exe文件

MATLAB版本可以参考另一个文件（完全与C++无关）
