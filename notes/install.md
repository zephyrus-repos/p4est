# p4est 编译安装指南

## 1. 概述

本文介绍在 Linux 环境下使用 GNU Autotools 编译安装 p4est 的完整流程，包括依赖准备、编译测试、安装验证。

p4est 是一个用于并行管理自适应四叉树和八叉树森林的 C 语言库，常用于自适应网格细化（AMR）和并行数值计算。其核心依赖 sc 库，启用并行功能时需要 MPI。

本文使用以下占位符，实际操作时应替换为对应目录：

| 占位符      | 含义             |
| -------- | -------------- |
| `<源码目录>` | p4est 源代码所在目录  |
| `<构建目录>` | 独立编译目录         |
| `<安装前缀>` | p4est 最终安装目录   |
| `<版本号>`  | 实际使用的 p4est 版本 |

## 2. 安装编译依赖

以 Ubuntu / Debian 为例：

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    autoconf automake libtool pkg-config \
    zlib1g-dev libjansson-dev
```

同时本机上要有MPI编译器。

主要依赖说明：

| 依赖                            | 用途                   |
| ----------------------------- | -------------------- |
| GCC / Make                    | C 语言编译与构建            |
| Autoconf / Automake / Libtool | GNU Autotools 构建系统   |
| MPI                           | 并行计算支持               |
| zlib                          | 压缩、校验和及 VTK 压缩输出相关功能 |
| Jansson                       | JSON 配置选项读取          |

检查 MPI 是否可用：

```bash
which mpicc
mpicc --version
which mpirun
```

## 3. 获取源码

从 [p4est 官方网站](https://www.p4est.org/) 下载正式发布的 tarball。

注意：GitHub 自动生成的 Source code 压缩包可能缺少 `sc` 子目录、`configure` 脚本和其他生成文件。建议使用官方网站提供的正式发布包。

```bash
tar -xzf p4est-<版本号>.tar.gz
cd p4est-<版本号>
```

## 4. 创建独立构建目录

建议采用源码目录外构建，避免生成文件与源码混杂：

```bash
cd <源码目录>
mkdir -p build
cd build
```

此时当前目录即为 `<构建目录>`，下文使用 `../configure` 调用上一级源码目录中的配置脚本。

## 5. 配置 p4est

推荐的 Release 构建配置如下：

```bash
../configure --prefix=<安装前缀> \
    --enable-mpi \
    --disable-shared \
    CFLAGS="-O2 -Wall -Wno-unused-parameter" \
    LIBS="-lm -lz -ljansson"
```

参数说明：

| 参数                         | 含义                     |
| -------------------------- | ---------------------- |
| `--prefix=<安装前缀>`          | 指定最终安装目录               |
| `--enable-mpi`             | 启用 MPI 并行功能            |
| `--disable-shared`         | 不构建共享库，使用静态库           |
| `CFLAGS="-O2 ..."`         | 启用优化及常用编译警告            |
| `LIBS="-lm -lz -ljansson"` | 显式链接数学库、zlib 和 Jansson |

如果未安装 zlib 或 Jansson，则不要加入对应的 `-lz` 或 `-ljansson`。其中 `-lm` 是解决本文所记录数学函数链接错误的关键选项。

### 5.1 调试构建

开发或调试时可以使用：

```bash
../configure --prefix=<安装前缀> \
    --enable-mpi \
    --enable-debug \
    CFLAGS="-O0 -g -Wall -Wextra" \
    LIBS="-lm -lz -ljansson"
```

调试构建启用断言和额外检查，性能通常低于优化构建。

### 5.2 配置结果检查

配置日志中应重点确认：

```text
checking whether we are using MPI... yes
configure: CC set to mpicc
checking compile/link for MPI C program... successful
checking whether to build shared libraries... no
checking whether to build static libraries... yes
configure: Building with source of package sc
```

如果 MPI 检测失败，可显式指定编译器：

```bash
../configure CC=mpicc \
    --prefix=<安装前缀> \
    --enable-mpi
```

## 6. 编译与测试

### 6.1 编译

```bash
make -j$(nproc)
```

`$(nproc)` 表示使用当前系统可用的处理器数量进行并行编译。

### 6.2 运行测试

```bash
make check
```

测试过程中可能调用 MPI，例如：

```bash
mpirun -np 2
```

应检查测试总结，确认是否存在失败项目。若 MPI 运行环境存在限制，应先排查运行环境问题，再判断是否为 p4est 本身的测试失败。

## 7. 安装

编译和测试通过后执行：

```bash
make install
```

如果安装前缀属于系统目录且当前用户没有写入权限，可使用：

```bash
sudo make install
```

建议仅在安装阶段使用必要的管理员权限，不要使用 root 权限进行整个编译过程。

## 8. 验证安装

检查安装目录下的头文件和库文件：

```bash
find <安装前缀> -name 'p4est.h'
find <安装前缀> -name 'p8est.h'
find <安装前缀> -name 'libp4est*'
find <安装前缀> -name 'libsc*'
```

安装成功后，通常可以找到 p4est、p8est 的头文件以及 p4est、sc 对应的静态库或共享库。

这些检查仅用于确认文件已安装。完整可用性仍应以测试结果及实际应用程序的编译运行结果为准。
