#!/bin/bash

sed -i s/a2fac/afac2/g module_bl_mynn_1d_cc.cpp
max=10
for ((i = 0 ; i < max ; i++ )); do
    sed -i 's/\([0-9]f\)/\1_rt/' module_bl_mynn_1d_cc.cpp
    sed -i s/f_rt/_rt/g module_bl_mynn_1d_cc.cpp
    sed -i 's/\(\.f\)/\1_rt/' module_bl_mynn_1d_cc.cpp
    sed -i s/f_rt/_rt/g module_bl_mynn_1d_cc.cpp
    sed -i s/'float\*'/'Real\*'/g  module_bl_mynn_1d_cc.cpp
    sed -i s/'float('/'Real('/g  module_bl_mynn_1d_cc.cpp
    sed -i s/'float '/'Real '/g  module_bl_mynn_1d_cc.cpp
    sed -i s/'float\&'/'Real\&'/g  module_bl_mynn_1d_cc.cpp
    sed -i s/'<float>'/'<Real>'/g  module_bl_mynn_1d_cc.cpp
    sed -i s/' abs'/' std::abs'/g  module_bl_mynn_1d_cc.cpp
    sed -i s/'=abs'/'=std::abs'/g  module_bl_mynn_1d_cc.cpp
    sed -i s/'(abs'/'(std::abs'/g  module_bl_mynn_1d_cc.cpp
done
sed -i s/afac2/a2fac/g module_bl_mynn_1d_cc.cpp
#grep [0-9]f ~/codes/WRF_mynn_changes/phys/module_bl_mynn_1d_cc.cpp | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g  | sed 's/\([0-9]f\)/\1_rt/' | sed s/f_rt/_rt/g

