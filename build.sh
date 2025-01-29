cwd=$(pwd)
dir_temp="$(dirname -- $(readlink -fn -- "$0"; echo x))"
dir="${dir_temp%x}"
cd $dir
cmake -B ./build
cd build
make
cd $cwd