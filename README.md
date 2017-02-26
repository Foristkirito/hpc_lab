# Labs of High Performance Computing

The official introduction of the class is below:

>This course consists of three parts, parallel programming technology, programming optimization and grid computing technology. Parallel programming technology includes parallel programming over OpenMP and MPI respectively. Programming optimization covers usage of VTune and popular compilers, optimization method for high performance computing users. Grid computing is composed grid environment construction and usage, data grid deployment, grid application development and so on. This course will last 48 hours, 24 hours for parallel programming, 12 hours for programming optimization, and 12 hours for grid computing. Each parts of this course will be practiced through small course project at least.

## Table of Content

- [Lab1: Optimize Image Rotating](#lab1)
- [Lab2: Distributed Stencil Computing](#lab2)
- [lab3: Solve large - scale linear equations in parallel](#lab3)
- [lab4: K-means on Spark](#lab4)
- [lab5: Implement two Graph Algoritm on GridGraph](#lab4)
- [Contact](#contact)

## <a name="lab1"></a>Lab1: Optimize Image Rotating
`cpp`  

In this lab, you are given a program showing moving and rotating image. Your task is to optimize the program to make it execute faster. [More Detail](https://github.com/Foristkirito/hpc_lab/tree/master/lab_1)

---

## <a name="lab2"></a>Lab2: Distributed Stencil Computing
`cpp`  

In this lab, you have to implement the [stencil computing](https://en.wikipedia.org/wiki/Stencil_code) using multi-machines. To get the fastes speed, you need to consider the communication and local in-momery computation. [More detail](https://github.com/Foristkirito/hpc_lab/tree/distributed/lab_2)

## <a name="lab3"></a>Lab3: Solve large - scale linear equations in parallel
`cpp`  

In this lab, you are ask to solve a linear equations. The equations are abstraction of three dimension features of the Earth about climate. The theory of problem solving is [here](http://www.netlib.org/templates/templates.pdf). Here uses the RGC algorithm and the preconditioner is ILU(0). [More detail](https://github.com/Foristkirito/hpc_lab/tree/distributed/lab_3)

## <a name="lab4"></a>Lab4: K-means on Spark
`Scala`

Using datasets from [kddcup99](http://kdd.ics.uci.edu/databases/kddcup99/kddcup99.html) to do k-means. The aim of this lab is to get famililar with Saprk. [More detail](https://github.com/Foristkirito/hpc_lab/tree/distributed/lab_4)

## <a name="lab5"></a>Lab5: Implement two Graph Algoritm on GridGraph
`cpp`

In this lab, you have to implement two graph algorithm using [GridGraph ATC '15](https://www.usenix.org/system/files/conference/atc15/atc15-paper-zhu.pdf) engine. The two is:

- PageRank-Delta
- [HyperANF](http://dl.acm.org/citation.cfm?id=1963493)

[More detail]()


## <a name="contact"></a> Contact

You can contact me using following methods:
- Email: xiaxin0202@foxmail.com
- WechatID: xiaxin0202
